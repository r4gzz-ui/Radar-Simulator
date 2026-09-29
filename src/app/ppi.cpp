#include "app/ppi.hpp"

#include <algorithm>
#include <cmath>

#include "app/theme.hpp"

namespace app {

using radar::degToRad;
using radar::kPi;

Ppi::Ppi(const radar::Params& params) : p_(params) {}

void Ppi::setLayout(Vector2 center, float radiusPx) {
    center_ = center;
    radius_ = radiusPx;
}

double Ppi::metersPerPixel() const { return p_.maxRange() / radius_; }

Vector2 Ppi::toScreen(double x, double y) const {
    const double s = radius_ / p_.maxRange();
    return {static_cast<float>(center_.x + x * s), static_cast<float>(center_.y - y * s)};
}

bool Ppi::toWorld(Vector2 screen, double& x, double& y) const {
    const double s = p_.maxRange() / radius_;
    x = (screen.x - center_.x) * s;
    y = (center_.y - screen.y) * s;
    return std::hypot(x, y) <= p_.maxRange();
}

void Ppi::clear() { blips_.clear(); }

void Ppi::addRawVideo(const radar::RangeDopplerMap& rd, double azimuth, double time) {
    // On affiche les cases dont l'energie depasse le bruit de 3 dB.
    const double s = std::sin(azimuth), c = std::cos(azimuth);
    for (int r = 3; r < rd.nRange; ++r) {
        const double db = 10.0 * std::log10(std::max(rd.rangeProfile[r] / rd.noiseLevel, 1e-12f));
        if (db < 3.0)
            continue;
        const double range = rd.rangeOfBin(r, p_);
        const float strength = static_cast<float>(std::clamp((db - 3.0) / 30.0, 0.08, 1.0));
        blips_.push_back({range * s, range * c, strength, time, false});
    }
}

void Ppi::addPlots(const std::vector<radar::Detection>& plots, double time) {
    for (const auto& d : plots) {
        const float strength = static_cast<float>(std::clamp(d.snrDb / 30.0, 0.4, 1.0));
        blips_.push_back({d.x, d.y, strength, time, true});
    }
}

void Ppi::prune(double now) {
    const double keep = 1.3 * p_.rotationPeriod;
    blips_.erase(std::remove_if(blips_.begin(), blips_.end(),
                                [&](const Blip& b) { return now - b.time > keep; }),
                 blips_.end());
}

// ---------------------------------------------------------------------------

namespace {

Vector2 polar(Vector2 c, float r, double az) {
    return {static_cast<float>(c.x + r * std::sin(az)), static_cast<float>(c.y - r * std::cos(az))};
}

} // namespace

void Ppi::drawGrid() const {
    DrawCircleV(center_, radius_, Color{3, 16, 8, 255});

    // Anneaux de distance tous les 100 m.
    for (double r = 100; r < p_.maxRange(); r += 100) {
        const float px = static_cast<float>(r / p_.maxRange() * radius_);
        DrawCircleLines(static_cast<int>(center_.x), static_cast<int>(center_.y), px, theme::grid);
        const char* label = TextFormat("%.0f m", r);
        const Vector2 at = polar(center_, px, degToRad(45));
        DrawText(label, static_cast<int>(at.x) + 3, static_cast<int>(at.y) - 12, 10, theme::gridText);
    }
    DrawCircleLines(static_cast<int>(center_.x), static_cast<int>(center_.y), radius_, theme::border);

    // Graduations d'azimut.
    for (int deg = 0; deg < 360; deg += 5) {
        const double az = degToRad(deg);
        const bool major = deg % 30 == 0;
        const float len = major ? 12.0f : (deg % 10 == 0 ? 7.0f : 3.0f);
        DrawLineV(polar(center_, radius_ - len, az), polar(center_, radius_, az), theme::gridText);
        if (major) {
            DrawLineV(center_, polar(center_, radius_, az), Fade(theme::grid, 0.6f));
            const char* label = deg == 0 ? "N" : TextFormat("%03d", deg);
            const int w = MeasureText(label, 12);
            const Vector2 at = polar(center_, radius_ + 16, az);
            DrawText(label, static_cast<int>(at.x) - w / 2, static_cast<int>(at.y) - 6, 12, theme::gridText);
        }
    }
}

void Ppi::drawSweep(double beamAzimuth) const {
    // Trainee lumineuse derriere la ligne de balayage.
    const int steps = 40;
    const double span = degToRad(30);
    for (int i = 0; i < steps; ++i) {
        const double a0 = beamAzimuth - span * i / steps;
        const double a1 = beamAzimuth - span * (i + 1) / steps;
        const float alpha = 0.22f * (1.0f - static_cast<float>(i) / steps);
        // Ordre des sommets adapte a raylib (sens anti-horaire a l'ecran).
        DrawTriangle(center_, polar(center_, radius_, a0), polar(center_, radius_, a1), Fade(theme::phosphor, alpha));
    }
    DrawLineEx(center_, polar(center_, radius_, beamAzimuth), 2.0f, theme::phosphor);
}

void Ppi::drawBlips(double now, bool showRaw) const {
    const double fade = p_.rotationPeriod;
    const float binPx = static_cast<float>(p_.rangeResolution() / metersPerPixel());
    for (const Blip& b : blips_) {
        const double age = now - b.time;
        const float life = static_cast<float>(1.0 - age / (b.plot ? 1.2 * fade : fade));
        if (life <= 0 || (!b.plot && !showRaw))
            continue;
        const Vector2 s = toScreen(b.x, b.y);
        if (b.plot) {
            DrawCircleV(s, 4.0f, Fade(theme::phosphor, 0.35f * life));
            DrawCircleV(s, 2.2f, Fade(WHITE, 0.9f * life * b.strength));
        } else {
            const float size = std::max(2.0f, binPx);
            DrawRectangleV({s.x - size / 2, s.y - size / 2}, {size, size},
                           Fade(theme::phosphor, 0.75f * life * b.strength));
        }
    }
}

void Ppi::drawTracks(double now, const radar::Tracker& tracker) const {
    for (const radar::Track& t : tracker.tracks()) {
        const auto x = tracker.predicted(t, now);
        const Vector2 pos = toScreen(x[0], x[1]);
        if (!t.confirmed) {
            DrawCircleLines(static_cast<int>(pos.x), static_cast<int>(pos.y), 5, theme::tentative);
            continue;
        }

        // Historique filtre.
        const size_t n = t.history.size();
        for (size_t i = 0; i + 1 < n; ++i) {
            const float a = 0.15f + 0.6f * static_cast<float>(i) / n;
            DrawLineV(toScreen(t.history[i].x, t.history[i].y), toScreen(t.history[i + 1].x, t.history[i + 1].y),
                      Fade(theme::track, a));
            DrawCircleV(toScreen(t.history[i].x, t.history[i].y), 1.5f, Fade(theme::track, a));
        }

        // Vecteur vitesse : ou sera la cible dans 15 s.
        const Vector2 lead = toScreen(x[0] + 15 * x[2], x[1] + 15 * x[3]);
        DrawLineEx(pos, lead, 1.5f, theme::track);

        // Symbole de piste.
        DrawRectangleLinesEx({pos.x - 5, pos.y - 5, 10, 10}, 1.5f, theme::track);
        const double speed = std::hypot(x[2], x[3]);
        double heading = std::atan2(x[2], x[3]) * 180.0 / kPi;
        if (heading < 0)
            heading += 360;
        DrawText(TextFormat("T%d", t.id), static_cast<int>(pos.x) + 9, static_cast<int>(pos.y) - 14, 14, theme::track);
        DrawText(TextFormat("%.0f m/s %03.0f", speed, heading), static_cast<int>(pos.x) + 9,
                 static_cast<int>(pos.y) + 1, 10, Fade(theme::track, 0.8f));
    }
}

void Ppi::drawTruth(const std::vector<radar::Target>& targets, const Options& opt) const {
    for (const radar::Target& t : targets) {
        const Vector2 s = toScreen(t.x, t.y);
        const bool sel = t.id == opt.selectedTarget;
        if (opt.showTruth) {
            DrawCircleLines(static_cast<int>(s.x), static_cast<int>(s.y), 3, theme::truth);
            DrawText(TextFormat("#%d", t.id), static_cast<int>(s.x) - 18, static_cast<int>(s.y) + 6, 10, theme::truth);
        }
        if (sel) {
            DrawRing(s, 10, 12, 0, 360, 32, Fade(theme::selected, 0.9f));
            if (opt.spotlight)
                DrawLineEx(center_, s, 1.0f, Fade(theme::selected, 0.35f));
        }
    }
}

void Ppi::draw(double now, double beamAzimuth, const std::vector<radar::Target>& targets,
               const radar::Tracker& tracker, const Options& opt) const {
    drawGrid();
    drawSweep(beamAzimuth);
    drawBlips(now, opt.showRaw);
    if (opt.showTracks)
        drawTracks(now, tracker);
    drawTruth(targets, opt);
    DrawCircleV(center_, 3, theme::phosphor);
}

} // namespace app
