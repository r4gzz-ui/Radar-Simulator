#pragma once

#include <vector>

#include "core/cfar.hpp"
#include "core/params.hpp"
#include "core/range_doppler.hpp"
#include "core/target.hpp"
#include "core/tracker.hpp"
#include "raylib.h"

namespace app {

// Ecran PPI (Plan Position Indicator) : vue de dessus, radar au centre,
// nord en haut, azimut dans le sens horaire.
class Ppi {
public:
    struct Options {
        bool showTruth = false;
        bool showTracks = true;
        bool showRaw = true;
        bool spotlight = true;
        int selectedTarget = -1;
    };

    explicit Ppi(const radar::Params& params);

    void setLayout(Vector2 center, float radiusPx);
    float radius() const { return radius_; }
    Vector2 center() const { return center_; }

    // "Video brute" : l'energie recue le long du lobe, case par case.
    void addRawVideo(const radar::RangeDopplerMap& rd, double azimuth, double time);
    // Plots (detections regroupees) : les echos francs.
    void addPlots(const std::vector<radar::Detection>& plots, double time);
    void prune(double now);
    void clear();

    void draw(double now, double beamAzimuth, const std::vector<radar::Target>& targets,
              const radar::Tracker& tracker, const Options& opt) const;

    Vector2 toScreen(double x, double y) const;
    // Renvoie false si le point ecran est hors du cercle.
    bool toWorld(Vector2 screen, double& x, double& y) const;
    double metersPerPixel() const;

private:
    struct Blip {
        double x, y;
        float strength;   // 0..1
        double time;
        bool plot;
    };

    void drawGrid() const;
    void drawSweep(double beamAzimuth) const;
    void drawBlips(double now, bool showRaw) const;
    void drawTracks(double now, const radar::Tracker& tracker) const;
    void drawTruth(const std::vector<radar::Target>& targets, const Options& opt) const;

    const radar::Params& p_;
    Vector2 center_{0, 0};
    float radius_ = 100;
    std::vector<Blip> blips_;
};

} // namespace app
