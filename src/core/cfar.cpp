#include "core/cfar.hpp"

#include <algorithm>
#include <cmath>

#include "core/target.hpp"

namespace radar {

double cfarAlpha(int nTrain, double pfa) {
    // Pour un bruit exponentiel (loi de |bruit gaussien complexe|^2), la
    // probabilite de fausse alarme du CA-CFAR vaut (1 + alpha / N)^-N.
    return nTrain * (std::pow(pfa, -1.0 / nTrain) - 1.0);
}

namespace {

// Image integrale : somme d'un rectangle quelconque en 4 lectures.
struct Integral {
    int w = 0, h = 0;
    std::vector<double> s;   // (h + 1) x (w + 1)

    void build(const std::vector<float>& data, int width, int height) {
        w = width;
        h = height;
        s.assign(static_cast<size_t>(w + 1) * (h + 1), 0.0);
        for (int y = 0; y < h; ++y) {
            double row = 0;
            for (int x = 0; x < w; ++x) {
                row += data[static_cast<size_t>(y) * w + x];
                s[static_cast<size_t>(y + 1) * (w + 1) + x + 1] = s[static_cast<size_t>(y) * (w + 1) + x + 1] + row;
            }
        }
    }
    // Somme sur [x0, x1] x [y0, y1] (bornes incluses, deja bornees a l'image).
    double sum(int x0, int y0, int x1, int y1) const {
        const auto at = [&](int x, int y) { return s[static_cast<size_t>(y) * (w + 1) + x]; };
        return at(x1 + 1, y1 + 1) - at(x0, y1 + 1) - at(x1 + 1, y0) + at(x0, y0);
    }
};

// Interpolation parabolique du sommet a partir de 3 valeurs (en dB) :
// renvoie un decalage dans [-0.5, 0.5] autour de la case centrale.
double parabolicOffset(double left, double center, double right) {
    const double denom = left - 2.0 * center + right;
    if (std::abs(denom) < 1e-12)
        return 0.0;
    return std::clamp(0.5 * (left - right) / denom, -0.5, 0.5);
}

} // namespace

std::vector<Detection> cfar2d(const RangeDopplerMap& map, const Params& p, const CfarConfig& cfg,
                              double azimuth, double time, std::vector<uint8_t>* hitMask) {
    const int nr = map.nRange;
    const int nd = map.nDoppler;
    Integral integ;
    integ.build(map.power, nr, nd);

    std::vector<uint8_t> hits(static_cast<size_t>(nr) * nd, 0);
    std::vector<float> noise(static_cast<size_t>(nr) * nd, 0.0f);

    const int outerR = cfg.guardRange + cfg.trainRange;
    const int outerD = cfg.guardDoppler + cfg.trainDoppler;

    for (int d = 0; d < nd; ++d) {
        for (int r = cfg.minRangeBin; r < nr; ++r) {
            const int r0 = std::max(0, r - outerR), r1 = std::min(nr - 1, r + outerR);
            const int d0 = std::max(0, d - outerD), d1 = std::min(nd - 1, d + outerD);
            const int gr0 = std::max(0, r - cfg.guardRange), gr1 = std::min(nr - 1, r + cfg.guardRange);
            const int gd0 = std::max(0, d - cfg.guardDoppler), gd1 = std::min(nd - 1, d + cfg.guardDoppler);

            const int nOuter = (r1 - r0 + 1) * (d1 - d0 + 1);
            const int nInner = (gr1 - gr0 + 1) * (gd1 - gd0 + 1);
            const int nTrain = nOuter - nInner;
            if (nTrain < 8)
                continue;

            const double trainSum = integ.sum(r0, d0, r1, d1) - integ.sum(gr0, gd0, gr1, gd1);
            const double localNoise = trainSum / nTrain;
            const double threshold = cfarAlpha(nTrain, cfg.pfa) * localNoise;
            const size_t idx = static_cast<size_t>(d) * nr + r;
            noise[idx] = static_cast<float>(localNoise);
            if (map.power[idx] > threshold)
                hits[idx] = 1;
        }
    }

    // Une cible occupe plusieurs cases (fenetrage) : on ne garde que les
    // maxima locaux parmi les cases detectees.
    std::vector<Detection> out;
    for (int d = 0; d < nd; ++d) {
        for (int r = cfg.minRangeBin; r < nr; ++r) {
            const size_t idx = static_cast<size_t>(d) * nr + r;
            if (!hits[idx])
                continue;
            const float v = map.power[idx];
            bool isMax = true;
            for (int dd = -1; dd <= 1 && isMax; ++dd) {
                for (int dr = -1; dr <= 1; ++dr) {
                    if (dd == 0 && dr == 0)
                        continue;
                    const int rr = r + dr;
                    const int di = (d + dd + nd) % nd;   // l'axe Doppler est circulaire
                    if (rr < 0 || rr >= nr)
                        continue;
                    if (map.at(di, rr) > v) {
                        isMax = false;
                        break;
                    }
                }
            }
            if (!isMax)
                continue;

            const auto db = [](float x) { return 10.0 * std::log10(std::max(x, 1e-20f)); };
            double rOff = 0, dOff = 0;
            if (r > 0 && r < nr - 1)
                rOff = parabolicOffset(db(map.at(d, r - 1)), db(v), db(map.at(d, r + 1)));
            dOff = parabolicOffset(db(map.at((d - 1 + nd) % nd, r)), db(v), db(map.at((d + 1) % nd, r)));

            Detection det;
            det.rangeBin = r;
            det.dopplerBin = d;
            det.range = map.rangeOfBin(r + rOff, p);
            det.radialVelocity = map.velocityOfBin(d + dOff, p);
            det.azimuth = azimuth;
            det.time = time;
            det.snrDb = db(v) - db(std::max(noise[idx], 1e-20f));
            det.x = det.range * std::sin(azimuth);
            det.y = det.range * std::cos(azimuth);
            out.push_back(det);
        }
    }

    // Suppression des lobes secondaires, du pic le plus fort au plus faible.
    std::sort(out.begin(), out.end(), [](const Detection& a, const Detection& b) { return a.snrDb > b.snrDb; });
    std::vector<Detection> kept;
    for (const Detection& d : out) {
        const float pw = map.at(d.dopplerBin, d.rangeBin);
        bool sidelobe = false;
        for (const Detection& k : kept) {
            int dd = std::abs(d.dopplerBin - k.dopplerBin);
            dd = std::min(dd, nd - dd);
            if (std::abs(d.rangeBin - k.rangeBin) <= cfg.sidelobeRange && dd <= cfg.sidelobeDoppler &&
                10.0 * std::log10(map.at(k.dopplerBin, k.rangeBin) / std::max(pw, 1e-20f)) > cfg.sidelobeDb) {
                sidelobe = true;
                break;
            }
        }
        if (!sidelobe)
            kept.push_back(d);
    }

    if (hitMask)
        *hitMask = std::move(hits);
    return kept;
}

// ---------------------------------------------------------------------------

PlotExtractor::PlotExtractor(const Params& params, int minDwells) : p_(params), minDwells_(minDwells) {}

void PlotExtractor::clear() { pending_.clear(); }

std::vector<Detection> PlotExtractor::push(const std::vector<Detection>& dets, double beamAzimuth) {
    const double gateRange = 2.5 * p_.rangeResolution();
    const double gateVel = 2.5 * p_.velocityResolution();
    const double gateAz = 2.0 * p_.beamWidth();

    for (const Detection& d : dets) {
        const double w = std::pow(10.0, d.snrDb / 10.0);
        Pending* match = nullptr;
        for (Pending& pl : pending_) {
            const double r = pl.range / pl.wSum;
            const double v = pl.velocity / pl.wSum;
            if (std::abs(r - d.range) < gateRange && std::abs(v - d.radialVelocity) < gateVel &&
                std::abs(angleDiff(d.azimuth, pl.lastAzimuth)) < gateAz) {
                match = &pl;
                break;
            }
        }
        if (!match) {
            pending_.emplace_back();
            match = &pending_.back();
            match->rangeBin = d.rangeBin;
            match->dopplerBin = d.dopplerBin;
        }
        match->count++;
        match->wSum += w;
        match->range += w * d.range;
        match->velocity += w * d.radialVelocity;
        match->sinSum += w * std::sin(d.azimuth);
        match->cosSum += w * std::cos(d.azimuth);
        match->time = d.time;
        match->lastAzimuth = d.azimuth;
        if (d.snrDb > match->bestSnrDb) {
            match->bestSnrDb = d.snrDb;
            match->rangeBin = d.rangeBin;
            match->dopplerBin = d.dopplerBin;
        }
    }

    // Un plot est termine quand le lobe s'en est suffisamment eloigne.
    std::vector<Detection> done;
    std::vector<Pending> still;
    for (const Pending& pl : pending_) {
        if (std::abs(angleDiff(beamAzimuth, pl.lastAzimuth)) > gateAz) {
            if (pl.count < minDwells_)
                continue;   // detection isolee : tres probablement du bruit
            Detection d;
            d.range = pl.range / pl.wSum;
            d.radialVelocity = pl.velocity / pl.wSum;
            d.azimuth = wrapAngle(std::atan2(pl.sinSum, pl.cosSum));
            d.snrDb = pl.bestSnrDb;
            d.time = pl.time;
            d.x = d.range * std::sin(d.azimuth);
            d.y = d.range * std::cos(d.azimuth);
            d.rangeBin = pl.rangeBin;
            d.dopplerBin = pl.dopplerBin;
            done.push_back(d);
        } else {
            still.push_back(pl);
        }
    }
    pending_ = std::move(still);
    return done;
}

} // namespace radar
