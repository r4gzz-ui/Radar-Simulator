#include "core/tracker.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace radar {

// ---------------------------------------------------------------------------
// Kalman
// ---------------------------------------------------------------------------

void Kalman::predict(double dt, double q) {
    if (dt <= 0)
        return;
    // x <- F x avec F = [[1,0,dt,0],[0,1,0,dt],[0,0,1,0],[0,0,0,1]]
    x[0] += dt * x[2];
    x[1] += dt * x[3];

    // P <- F P F^T : calcul direct, F etant tres creuse.
    std::array<double, 16> A = P;   // A = F P
    for (int j = 0; j < 4; ++j) {
        A[0 * 4 + j] = P[0 * 4 + j] + dt * P[2 * 4 + j];
        A[1 * 4 + j] = P[1 * 4 + j] + dt * P[3 * 4 + j];
    }
    std::array<double, 16> B = A;   // B = A F^T
    for (int i = 0; i < 4; ++i) {
        B[i * 4 + 0] = A[i * 4 + 0] + dt * A[i * 4 + 2];
        B[i * 4 + 1] = A[i * 4 + 1] + dt * A[i * 4 + 3];
    }

    // Q : bruit d'acceleration blanc, modele "vitesse constante" discretise.
    const double q2 = q * q;
    const double dt2 = dt * dt, dt3 = dt2 * dt, dt4 = dt3 * dt;
    B[0 * 4 + 0] += q2 * dt4 / 4;
    B[1 * 4 + 1] += q2 * dt4 / 4;
    B[0 * 4 + 2] += q2 * dt3 / 2;
    B[2 * 4 + 0] += q2 * dt3 / 2;
    B[1 * 4 + 3] += q2 * dt3 / 2;
    B[3 * 4 + 1] += q2 * dt3 / 2;
    B[2 * 4 + 2] += q2 * dt2;
    B[3 * 4 + 3] += q2 * dt2;
    P = B;
}

double Kalman::mahalanobis2(const std::array<double, 2>& z, const std::array<double, 4>& R) const {
    // Innovation nu = z - H x, covariance S = H P H^T + R (2x2).
    const double nx = z[0] - x[0], ny = z[1] - x[1];
    const double s00 = P[0] + R[0], s01 = P[1] + R[1];
    const double s10 = P[4] + R[2], s11 = P[5] + R[3];
    const double det = s00 * s11 - s01 * s10;
    if (det <= 0)
        return std::numeric_limits<double>::infinity();
    // nu^T S^-1 nu
    return (s11 * nx * nx - (s01 + s10) * nx * ny + s00 * ny * ny) / det;
}

void Kalman::update(const std::array<double, 2>& z, const std::array<double, 4>& R) {
    const double nx = z[0] - x[0], ny = z[1] - x[1];
    const double s00 = P[0] + R[0], s01 = P[1] + R[1];
    const double s10 = P[4] + R[2], s11 = P[5] + R[3];
    const double det = s00 * s11 - s01 * s10;
    if (det <= 0)
        return;
    const double i00 = s11 / det, i01 = -s01 / det, i10 = -s10 / det, i11 = s00 / det;

    // K = P H^T S^-1 : H selectionne les deux premieres colonnes de P.
    double K[4][2];
    for (int i = 0; i < 4; ++i) {
        const double p0 = P[i * 4 + 0], p1 = P[i * 4 + 1];
        K[i][0] = p0 * i00 + p1 * i10;
        K[i][1] = p0 * i01 + p1 * i11;
    }
    for (int i = 0; i < 4; ++i)
        x[i] += K[i][0] * nx + K[i][1] * ny;

    // P <- P - K H P  (H P = les deux premieres lignes de P)
    std::array<double, 16> newP = P;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            newP[i * 4 + j] -= K[i][0] * P[0 * 4 + j] + K[i][1] * P[1 * 4 + j];
    // Symetrisation pour eviter la derive numerique.
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j) {
            const double m = 0.5 * (newP[i * 4 + j] + newP[j * 4 + i]);
            newP[i * 4 + j] = newP[j * 4 + i] = m;
        }
    P = newP;
}

// ---------------------------------------------------------------------------
// Tracker
// ---------------------------------------------------------------------------

Tracker::Tracker(const Params& params, TrackerConfig cfg) : p_(params), cfg_(cfg) {}

void Tracker::clear() {
    tracks_.clear();
    nextId_ = 1;
}

std::array<double, 4> Tracker::measurementNoise(double range, double azimuth) const {
    // Incertitudes en polaire : ~1/2 case en distance, ~1/4 de lobe en azimut.
    const double sr = 0.5 * p_.rangeResolution();
    const double sa = 0.25 * p_.beamWidth();
    // x = r sin(az), y = r cos(az) -> J = d(x,y)/d(r,az)
    const double s = std::sin(azimuth), c = std::cos(azimuth);
    const double j00 = s, j01 = range * c;
    const double j10 = c, j11 = -range * s;
    // R = J diag(sr^2, sa^2) J^T
    const double a = sr * sr, b = sa * sa;
    const double r00 = j00 * j00 * a + j01 * j01 * b;
    const double r01 = j00 * j10 * a + j01 * j11 * b;
    const double r11 = j10 * j10 * a + j11 * j11 * b;
    return {r00, r01, r01, r11};
}

void Tracker::process(const std::vector<Detection>& plots) {
    for (const Detection& d : plots) {
        const std::array<double, 2> z{d.x, d.y};
        const std::array<double, 4> R = measurementNoise(d.range, d.azimuth);

        // Association "plus proche voisin" dans la porte de validation.
        Track* best = nullptr;
        double bestD2 = cfg_.gateChi2;
        for (Track& t : tracks_) {
            Kalman k = t.kf;
            k.predict(d.time - t.lastUpdate, cfg_.accelNoise);
            const double d2 = k.mahalanobis2(z, R);
            if (d2 < bestD2) {
                bestD2 = d2;
                best = &t;
            }
        }

        if (best) {
            best->kf.predict(d.time - best->lastUpdate, cfg_.accelNoise);
            best->kf.update(z, R);
            best->lastUpdate = d.time;
            best->hits++;
            if (best->hits >= cfg_.confirmHits)
                best->confirmed = true;
        } else {
            // Nouvelle piste "tentative". La vitesse radiale mesuree par
            // Doppler sert d'estimation initiale ; la vitesse tangentielle
            // est inconnue, d'ou une grande incertitude.
            Track t;
            t.id = nextId_++;
            const double ux = std::sin(d.azimuth), uy = std::cos(d.azimuth);
            t.kf.x = {d.x, d.y, d.radialVelocity * ux, d.radialVelocity * uy};
            t.kf.P.fill(0.0);
            t.kf.P[0] = R[0];
            t.kf.P[1] = R[1];
            t.kf.P[4] = R[2];
            t.kf.P[5] = R[3];
            const double vVar = 20.0 * 20.0;
            t.kf.P[10] = vVar;
            t.kf.P[15] = vVar;
            t.lastUpdate = d.time;
            t.hits = 1;
            tracks_.push_back(t);
            best = &tracks_.back();
        }

        best->lastSnrDb = d.snrDb;
        best->lastRadialVelocity = d.radialVelocity;
        best->history.push_back({best->kf.x[0], best->kf.x[1], d.time});
        best->rawPlots.push_back({d.x, d.y, d.time});
        while (best->history.size() > cfg_.historyLength)
            best->history.pop_front();
        while (best->rawPlots.size() > cfg_.historyLength)
            best->rawPlots.pop_front();
    }
}

void Tracker::prune(double now) {
    const double period = p_.rotationPeriod;
    tracks_.erase(std::remove_if(tracks_.begin(), tracks_.end(),
                                 [&](const Track& t) {
                                     const double limit = t.confirmed ? cfg_.deleteConfirmedScans
                                                                      : cfg_.deleteTentativeScans;
                                     return now - t.lastUpdate > limit * period;
                                 }),
                  tracks_.end());
}

std::array<double, 4> Tracker::predicted(const Track& t, double now) const {
    Kalman k = t.kf;
    k.predict(now - t.lastUpdate, cfg_.accelNoise);
    return k.x;
}

} // namespace radar
