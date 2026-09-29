#pragma once

#include <cstdint>
#include <vector>

#include "core/cfar.hpp"
#include "core/params.hpp"
#include "core/range_doppler.hpp"
#include "core/simulator.hpp"
#include "raylib.h"

namespace app {

// Tout ce qu'il faut pour dessiner les panneaux "ondes" d'un dwell.
struct DwellView {
    const radar::Params* params = nullptr;
    const radar::Frame* frame = nullptr;
    const radar::RangeDopplerMap* rd = nullptr;
    const std::vector<radar::EchoInfo>* echoes = nullptr;
    const std::vector<radar::Detection>* detections = nullptr;
    int focusTarget = -1;          // echo a mettre en avant
    bool showTruth = false;
    const char* source = "";       // "projecteur #2", "balayage"...
};

// Panneaux de droite :
//  1. formes d'onde emises / recues et frequence instantanee des chirps
//  2. signal de battement (IF) d'un chirp
//  3. profil distance (FFT n1)
//  4. carte distance / vitesse (FFT n2) et detections CFAR
class Plots {
public:
    void init(const radar::Params& params);
    void unload();
    void draw(Rectangle area, const DwellView& view);

private:
    void drawWaves(Rectangle r, const DwellView& v);
    void drawBeat(Rectangle r, const DwellView& v);
    void drawRangeProfile(Rectangle r, const DwellView& v);
    void drawRangeDoppler(Rectangle r, const DwellView& v);

    Texture2D rdTexture_{};
    std::vector<Color> pixels_;
    int texW_ = 0, texH_ = 0;
    bool ready_ = false;
};

} // namespace app
