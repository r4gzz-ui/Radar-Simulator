#pragma once

#include <array>
#include <deque>
#include <vector>

#include "core/cfar.hpp"
#include "core/params.hpp"

namespace radar {

// Filtre de Kalman lineaire a vitesse constante.
// Etat x = [px, py, vx, vy], mesure z = [px, py] (le plot converti en cartesien).
struct Kalman {
    std::array<double, 4> x{};
    std::array<double, 16> P{};   // covariance 4x4, ligne par ligne

    // Prediction : x <- F x, P <- F P F^T + Q
    // q = ecart-type de l'acceleration non modelisee (m/s^2).
    void predict(double dt, double q);

    // Distance de Mahalanobis au carre de la mesure z (sans modifier l'etat).
    double mahalanobis2(const std::array<double, 2>& z, const std::array<double, 4>& R) const;

    // Correction : x <- x + K (z - H x), P <- (I - K H) P
    void update(const std::array<double, 2>& z, const std::array<double, 4>& R);
};

struct TrackPoint {
    double x, y, time;
};

struct Track {
    int id = 0;
    Kalman kf;
    double lastUpdate = 0;
    int hits = 0;
    bool confirmed = false;
    double lastSnrDb = 0;
    double lastRadialVelocity = 0;
    std::deque<TrackPoint> history;   // positions filtrees successives
    std::deque<TrackPoint> rawPlots;  // plots bruts associes (pour comparer)
};

struct TrackerConfig {
    double accelNoise = 2.0;      // m/s^2 : manoeuvrabilite supposee des cibles
                                  // (petit = trajectoires lissees, grand = reactif aux virages)
    double gateChi2 = 13.8;       // seuil de Mahalanobis (chi2 2 ddl, 99.9 %)
    int confirmHits = 3;          // nb de plots pour confirmer une piste
    double deleteConfirmedScans = 2.5;   // tours sans plot avant suppression
    double deleteTentativeScans = 1.5;
    size_t historyLength = 40;
};

class Tracker {
public:
    explicit Tracker(const Params& params, TrackerConfig cfg = {});

    // Associe chaque plot a la piste la plus proche (au sens de Mahalanobis),
    // met a jour cette piste, ou cree une nouvelle piste si aucune ne convient.
    void process(const std::vector<Detection>& plots);

    // Supprime les pistes qui n'ont plus ete vues depuis trop longtemps.
    void prune(double now);

    // Etat predit d'une piste a l'instant now (pour l'affichage).
    std::array<double, 4> predicted(const Track& t, double now) const;

    const std::vector<Track>& tracks() const { return tracks_; }
    void clear();

    // Bruit de mesure d'un plot, converti de (distance, azimut) en (x, y).
    std::array<double, 4> measurementNoise(double range, double azimuth) const;

    TrackerConfig& config() { return cfg_; }

private:
    const Params& p_;
    TrackerConfig cfg_;
    std::vector<Track> tracks_;
    int nextId_ = 1;
};

} // namespace radar
