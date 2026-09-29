// Simulateur de radar FMCW : generation des echos, double FFT, CFAR,
// pistage par filtre de Kalman et affichage PPI.
//
// Usage : ./radar_sim [--warmup SECONDES] [--screenshot fichier.png [--after SECONDES]]

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "app/plots.hpp"
#include "app/ppi.hpp"
#include "app/theme.hpp"
#include "core/cfar.hpp"
#include "core/params.hpp"
#include "core/range_doppler.hpp"
#include "core/simulator.hpp"
#include "core/target.hpp"
#include "core/tracker.hpp"
#include "raylib.h"
#include "rlgl.h"

using namespace radar;
namespace theme = app::theme;

namespace {

struct App {
    Params params;
    Simulator sweepSim{params, 1};
    Simulator spotSim{params, 2};
    RangeDopplerProcessor processor{params};
    CfarConfig cfar;
    PlotExtractor extractor{params};
    Tracker tracker{params};
    app::Ppi ppi{params};
    app::Plots plots;

    std::vector<Target> targets;
    int nextTargetId = 1;
    std::mt19937 rng{1234};

    double time = 0;
    double beamAzimuth = 0;
    bool paused = false;
    bool showHelp = true;
    app::Ppi::Options ppiOpt;

    // Dernier dwell du balayage, et dernier dwell qui contenait un echo.
    Frame sweepFrame;
    RangeDopplerMap sweepRd;
    std::vector<EchoInfo> sweepEchoes;
    std::vector<Detection> sweepDets;
    Frame heldFrame;
    RangeDopplerMap heldRd;
    std::vector<EchoInfo> heldEchoes;
    std::vector<Detection> heldDets;
    double heldTime = -1;

    // Dwell "projecteur" : un faisceau pointe en permanence sur la cible
    // selectionnee, pour voir ses ondes en continu.
    Frame spotFrame;
    RangeDopplerMap spotRd;
    std::vector<EchoInfo> spotEchoes;
    std::vector<Detection> spotDets;

    int plotsThisScan = 0;
    int plotsLastScan = 0;
    double scanStart = 0;

    // Creation de cible a la souris.
    bool dragging = false;
    Vector2 dragStart{};
};

const Target* findTarget(const App& a, int id) {
    for (const auto& t : a.targets)
        if (t.id == id)
            return &t;
    return nullptr;
}

void addTarget(App& a, double x, double y, double vx, double vy, double rcs, double turnRate = 0) {
    Target t;
    t.id = a.nextTargetId++;
    t.x = x;
    t.y = y;
    t.vx = vx;
    t.vy = vy;
    t.rcs = rcs;
    t.turnRate = turnRate;
    a.targets.push_back(t);
}

void resetScenario(App& a) {
    a.targets.clear();
    a.nextTargetId = 1;
    a.tracker.clear();
    a.extractor.clear();
    a.ppi.clear();
    a.heldTime = -1;
    // Quelques cibles de depart, avec des trajectoires variees.
    addTarget(a, -250, 250, 12, -4, 10);                 // traverse vers l'est
    addTarget(a, 180, 330, -3, -14, 30);                 // s'approche
    addTarget(a, 300, -150, 0, 18, 5, degToRad(3.0));    // vire doucement
    addTarget(a, -120, -380, 6, 9, 50);                  // gros et lointain
    addTarget(a, 60, 120, -9, 2, 1);                     // petit et proche
    a.ppiOpt.selectedTarget = 1;
}

void processDwell(App& a, double az, double t) {
    a.sweepSim.generate(a.targets, az, t, a.sweepFrame, &a.sweepEchoes);
    a.processor.process(a.sweepFrame, a.sweepRd);
    a.sweepDets = cfar2d(a.sweepRd, a.params, a.cfar, az, t);
    const std::vector<Detection> plots = a.extractor.push(a.sweepDets, az);
    a.tracker.process(plots);
    a.ppi.addRawVideo(a.sweepRd, az, t);
    a.ppi.addPlots(plots, t);
    a.plotsThisScan += static_cast<int>(plots.size());

    if (!a.sweepDets.empty() || !a.sweepEchoes.empty()) {
        a.heldFrame = a.sweepFrame;
        a.heldRd = a.sweepRd;
        a.heldEchoes = a.sweepEchoes;
        a.heldDets = a.sweepDets;
        a.heldTime = t;
    }
}

void step(App& a, double dt) {
    const double t0 = a.time;
    a.time += dt;
    for (auto& t : a.targets)
        moveTarget(t, dt, 40.0, 0.75 * a.params.maxRange());

    // Plusieurs dwells par image si besoin, pour ne jamais "sauter" une cible.
    const double azStep = kTwoPi * dt / a.params.rotationPeriod;
    const int n = std::max(1, static_cast<int>(std::ceil(azStep / degToRad(a.params.maxAzStepDeg))));
    for (int i = 0; i < n; ++i) {
        const double before = a.beamAzimuth;
        a.beamAzimuth = wrapAngle(a.beamAzimuth + azStep / n);
        if (a.beamAzimuth < before) {   // passage au nord : nouveau tour
            a.plotsLastScan = a.plotsThisScan;
            a.plotsThisScan = 0;
            a.scanStart = a.time;
        }
        processDwell(a, a.beamAzimuth, t0 + dt * (i + 1) / n);
    }
    a.tracker.prune(a.time);
    a.ppi.prune(a.time);

    if (const Target* sel = findTarget(a, a.ppiOpt.selectedTarget)) {
        const double az = sel->azimuth();
        a.spotSim.generate(a.targets, az, a.time, a.spotFrame, &a.spotEchoes);
        a.processor.process(a.spotFrame, a.spotRd);
        a.spotDets = cfar2d(a.spotRd, a.params, a.cfar, az, a.time);
    }
}

// --- Entrees ---------------------------------------------------------------

void selectNext(App& a) {
    if (a.targets.empty()) {
        a.ppiOpt.selectedTarget = -1;
        return;
    }
    auto it = std::find_if(a.targets.begin(), a.targets.end(),
                           [&](const Target& t) { return t.id == a.ppiOpt.selectedTarget; });
    it = (it == a.targets.end() || std::next(it) == a.targets.end()) ? a.targets.begin() : std::next(it);
    a.ppiOpt.selectedTarget = it->id;
}

int targetNear(const App& a, Vector2 screen, float maxPx) {
    int best = -1;
    float bestD = maxPx;
    for (const auto& t : a.targets) {
        const Vector2 s = a.ppi.toScreen(t.x, t.y);
        const float d = std::hypot(s.x - screen.x, s.y - screen.y);
        if (d < bestD) {
            bestD = d;
            best = t.id;
        }
    }
    return best;
}

Vector2 dragVelocity(const App& a, Vector2 from, Vector2 to) {
    // 1 pixel de glisse = deplacement parcouru en 5 s.
    const double mpp = a.ppi.metersPerPixel();
    double vx = (to.x - from.x) * mpp / 5.0;
    double vy = -(to.y - from.y) * mpp / 5.0;
    // On autorise un peu plus que la vitesse max non ambigue, pour pouvoir
    // observer le repliement Doppler.
    const double vmax = 70.0;
    const double v = std::hypot(vx, vy);
    if (v > vmax) {
        vx *= vmax / v;
        vy *= vmax / v;
    }
    return {static_cast<float>(vx), static_cast<float>(vy)};
}

void handleInput(App& a) {
    const Vector2 mouse = GetMousePosition();
    double wx, wy;
    const bool inPpi = a.ppi.toWorld(mouse, wx, wy);

    if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && inPpi) {
        a.dragging = true;
        a.dragStart = mouse;
    }
    if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT) && a.dragging) {
        a.dragging = false;
        const float dist = std::hypot(mouse.x - a.dragStart.x, mouse.y - a.dragStart.y);
        const int near = targetNear(a, a.dragStart, 14);
        if (dist < 5 && near >= 0) {
            a.ppiOpt.selectedTarget = near;   // simple clic sur une cible : selection
        } else {
            double sx, sy;
            a.ppi.toWorld(a.dragStart, sx, sy);
            const Vector2 v = dragVelocity(a, a.dragStart, mouse);
            static const double rcsChoices[] = {1, 5, 10, 30, 100};
            std::uniform_int_distribution<int> pick(0, 4);
            addTarget(a, sx, sy, v.x, v.y, rcsChoices[pick(a.rng)]);
            a.ppiOpt.selectedTarget = a.targets.back().id;
        }
    }
    if (IsMouseButtonPressed(MOUSE_BUTTON_RIGHT)) {
        const int near = targetNear(a, mouse, 20);
        if (near >= 0) {
            a.targets.erase(std::remove_if(a.targets.begin(), a.targets.end(),
                                           [&](const Target& t) { return t.id == near; }),
                            a.targets.end());
            if (a.ppiOpt.selectedTarget == near)
                selectNext(a);
        }
    }

    if (IsKeyPressed(KEY_SPACE))
        a.paused = !a.paused;
    if (IsKeyPressed(KEY_TAB))
        selectNext(a);
    if (IsKeyPressed(KEY_KP_ADD))
        a.params.rotationPeriod = std::max(1.0, a.params.rotationPeriod / 1.25);
    if (IsKeyPressed(KEY_KP_SUBTRACT))
        a.params.rotationPeriod = std::min(16.0, a.params.rotationPeriod * 1.25);

    // Touches "caracteres" : independantes de la disposition du clavier (AZERTY/QWERTY).
    for (int c = GetCharPressed(); c > 0; c = GetCharPressed()) {
        switch (c) {
        case '+': a.params.rotationPeriod = std::max(1.0, a.params.rotationPeriod / 1.25); break;
        case '-': a.params.rotationPeriod = std::min(16.0, a.params.rotationPeriod * 1.25); break;
        case 'n': case 'N': a.params.noiseSigma = std::min(20.0, a.params.noiseSigma * 1.41); break;
        case 'b': case 'B': a.params.noiseSigma = std::max(0.1, a.params.noiseSigma / 1.41); break;
        case 'g': case 'G': a.ppiOpt.showTruth = !a.ppiOpt.showTruth; break;
        case 't': case 'T': a.ppiOpt.showTracks = !a.ppiOpt.showTracks; break;
        case 'v': case 'V': a.ppiOpt.showRaw = !a.ppiOpt.showRaw; break;
        case 'p': case 'P': a.ppiOpt.spotlight = !a.ppiOpt.spotlight; break;
        case 'c': case 'C': a.tracker.clear(); break;
        case 'r': case 'R': resetScenario(a); break;
        case 'h': case 'H': a.showHelp = !a.showHelp; break;
        default: break;
        }
    }
}

// --- Affichage -------------------------------------------------------------

void drawHud(const App& a, Rectangle ppiArea) {
    const Params& p = a.params;
    int x = static_cast<int>(ppiArea.x) + 12, y = static_cast<int>(ppiArea.y) + 10;
    DrawText("RADAR FMCW 77 GHz", x, y, 20, theme::phosphor);
    y += 26;
    // TextFormat() reutilise un petit tampon interne : on dessine chaque ligne
    // tout de suite au lieu de les preparer a l'avance.
    const auto line = [&](const char* txt) {
        DrawText(txt, x, y, 12, theme::text);
        y += 16;
    };
    line(TextFormat("B %.0f MHz   Tc %.0f us", p.bandwidth / 1e6, p.chirpTime * 1e6));
    line(TextFormat("Ns %d   Nc %d", p.samplesPerChirp, p.chirpsPerFrame));
    line(TextFormat("res. %.1f m  portee %.0f m", p.rangeResolution(), p.maxRange()));
    line(TextFormat("res. v %.2f m/s", p.velocityResolution()));
    line(TextFormat("vmax +/-%.1f m/s", p.maxVelocity()));
    line(TextFormat("lobe %.1f deg  tour %.1f s", p.beamWidthDeg, p.rotationPeriod));
    line(TextFormat("bruit x%.2f", p.noiseSigma));

    int confirmed = 0;
    for (const auto& t : a.tracker.tracks())
        confirmed += t.confirmed;
    const int rx = static_cast<int>(ppiArea.x + ppiArea.width) - 200;
    int ry = static_cast<int>(ppiArea.y) + 10;
    DrawText(TextFormat("t = %.1f s   %d FPS", a.time, GetFPS()), rx, ry, 12, theme::text);
    ry += 16;
    DrawText(TextFormat("cibles : %zu", a.targets.size()), rx, ry, 12, theme::text);
    ry += 16;
    DrawText(TextFormat("pistes : %d confirmees / %zu", confirmed, a.tracker.tracks().size()), rx, ry, 12,
             theme::track);
    ry += 16;
    DrawText(TextFormat("plots au dernier tour : %d", a.plotsLastScan), rx, ry, 12, theme::text);
    ry += 16;
    if (a.paused)
        DrawText("PAUSE", rx, ry + 4, 20, theme::selected);

    const int by = static_cast<int>(ppiArea.y + ppiArea.height);
    if (a.showHelp) {
        const char* help[] = {
            "clic + glisser : creer une cible (le trait = sa vitesse)",
            "clic sur une cible : la selectionner   TAB : suivante",
            "clic droit : supprimer   R : scenario initial   C : effacer les pistes",
            "ESPACE pause   +/- vitesse de rotation   N/B bruit +/-",
            "G verite terrain   T pistes   V video brute   P projecteur   H aide",
        };
        int hy = by - 12 - 15 * 5;
        for (const char* h : help) {
            DrawText(h, x, hy, 11, theme::dimText);
            hy += 15;
        }
    } else {
        DrawText("H : aide", x, by - 22, 11, theme::dimText);
    }

    // Legende.
    const int lx = static_cast<int>(ppiArea.x + ppiArea.width) - 150;
    int ly = by - 12 - 15 * 4;
    DrawCircle(lx + 4, ly + 5, 2.2f, WHITE);
    DrawText("plot (detection)", lx + 14, ly, 11, theme::dimText);
    ly += 15;
    DrawRectangleLines(lx, ly + 1, 9, 9, theme::track);
    DrawText("piste confirmee", lx + 14, ly, 11, theme::dimText);
    ly += 15;
    DrawCircleLines(lx + 4, ly + 5, 3, theme::truth);
    DrawText("verite terrain (G)", lx + 14, ly, 11, theme::dimText);
    ly += 15;
    DrawRing({static_cast<float>(lx + 4), static_cast<float>(ly + 5)}, 4, 5.5f, 0, 360, 16, theme::selected);
    DrawText("cible selectionnee", lx + 14, ly, 11, theme::dimText);
}

void drawDragPreview(const App& a) {
    if (!a.dragging)
        return;
    const Vector2 m = GetMousePosition();
    DrawCircleLines(static_cast<int>(a.dragStart.x), static_cast<int>(a.dragStart.y), 5, theme::selected);
    DrawLineEx(a.dragStart, m, 2, theme::selected);
    const Vector2 v = dragVelocity(a, a.dragStart, m);
    DrawText(TextFormat("%.1f m/s", std::hypot(v.x, v.y)), static_cast<int>(m.x) + 8, static_cast<int>(m.y) - 8,
             12, theme::selected);
}

void render(App& a) {
    const float w = static_cast<float>(GetScreenWidth());
    const float h = static_cast<float>(GetScreenHeight());
    const float ppiSize = std::min(h, w * 0.55f);
    const Rectangle ppiArea{0, 0, ppiSize, h};
    a.ppi.setLayout({ppiSize / 2, h / 2}, ppiSize / 2 - 42);
    const Rectangle plotArea{ppiSize + 6, 8, w - ppiSize - 14, h - 16};

    BeginDrawing();
    ClearBackground(theme::background);

    a.ppi.draw(a.time, a.beamAzimuth, a.targets, a.tracker, a.ppiOpt);
    drawDragPreview(a);
    drawHud(a, ppiArea);

    app::DwellView view;
    view.params = &a.params;
    view.focusTarget = a.ppiOpt.selectedTarget;
    view.showTruth = a.ppiOpt.showTruth;
    static std::string source;
    if (a.ppiOpt.spotlight && findTarget(a, a.ppiOpt.selectedTarget)) {
        view.frame = &a.spotFrame;
        view.rd = &a.spotRd;
        view.echoes = &a.spotEchoes;
        view.detections = &a.spotDets;
        source = TextFormat("projecteur sur la cible #%d (P : balayage)", a.ppiOpt.selectedTarget);
    } else if (a.heldTime >= 0) {
        view.frame = &a.heldFrame;
        view.rd = &a.heldRd;
        view.echoes = &a.heldEchoes;
        view.detections = &a.heldDets;
        source = TextFormat("balayage : dernier echo il y a %.1f s (P : projecteur)", a.time - a.heldTime);
    } else {
        view.frame = &a.sweepFrame;
        view.rd = &a.sweepRd;
        view.echoes = &a.sweepEchoes;
        view.detections = &a.sweepDets;
        source = "balayage";
    }
    view.source = source.c_str();
    a.plots.draw(plotArea, view);

    EndDrawing();
}

} // namespace

int main(int argc, char** argv) {
    double warmup = 0;
    const char* screenshot = nullptr;
    double screenshotAfter = 0;   // delai (temps reel) avant la capture
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--warmup") && i + 1 < argc)
            warmup = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--screenshot") && i + 1 < argc)
            screenshot = argv[++i];
        else if (!std::strcmp(argv[i], "--after") && i + 1 < argc)
            screenshotAfter = std::atof(argv[++i]);
    }

    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE | FLAG_VSYNC_HINT);
    InitWindow(1600, 900, "Radar FMCW - simulateur");
    SetTargetFPS(60);
    SetExitKey(KEY_ESCAPE);
    rlDisableBackfaceCulling();   // la trainee du balayage est faite de triangles

    App a;
    a.plots.init(a.params);
    resetScenario(a);

    // Avance rapide optionnelle (utile pour arriver directement a des pistes confirmees).
    for (double t = 0; t < warmup; t += 1.0 / 60)
        step(a, 1.0 / 60);

    int frame = 0;
    while (!WindowShouldClose()) {
        handleInput(a);
        if (!a.paused)
            step(a, std::min(static_cast<double>(GetFrameTime()), 1.0 / 30));
        render(a);

        if (screenshot && ++frame >= 3 && GetTime() >= screenshotAfter) {
            Image img = LoadImageFromScreen();
            ExportImage(img, screenshot);
            UnloadImage(img);
            break;
        }
    }

    a.plots.unload();
    CloseWindow();
    return 0;
}
