// Tests du coeur de calcul (sans affichage).
// Lancement : make test

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "core/cfar.hpp"
#include "core/fft.hpp"
#include "core/params.hpp"
#include "core/range_doppler.hpp"
#include "core/simulator.hpp"
#include "core/tracker.hpp"

using namespace radar;

static int failures = 0;

#define CHECK(cond, ...)                                                   \
    do {                                                                   \
        if (cond) {                                                        \
            std::printf("  [OK]   ");                                      \
        } else {                                                           \
            std::printf("  [FAIL] ");                                      \
            ++failures;                                                    \
        }                                                                  \
        std::printf(__VA_ARGS__);                                          \
        std::printf("\n");                                                 \
    } while (0)

static void testFft() {
    std::printf("FFT (Cooley-Tukey) contre DFT naive\n");
    std::mt19937 rng(1);
    std::normal_distribution<float> g;
    const int n = 64;
    std::vector<cfloat> x(n);
    for (auto& v : x)
        v = cfloat(g(rng), g(rng));

    std::vector<cfloat> ref(n);
    for (int k = 0; k < n; ++k) {
        std::complex<double> s = 0;
        for (int i = 0; i < n; ++i)
            s += std::complex<double>(x[i]) * std::polar(1.0, -kTwoPi * k * i / n);
        ref[k] = cfloat(s);
    }
    std::vector<cfloat> y = x;
    fft(y);
    double err = 0;
    for (int k = 0; k < n; ++k)
        err = std::max(err, static_cast<double>(std::abs(y[k] - ref[k])));
    CHECK(err < 1e-3, "ecart max FFT/DFT = %.2e", err);

    fft(y, true);
    double errInv = 0;
    for (int i = 0; i < n; ++i)
        errInv = std::max(errInv, static_cast<double>(std::abs(y[i] - x[i])));
    CHECK(errInv < 1e-4, "FFT inverse retrouve le signal (ecart %.2e)", errInv);
}

static void testRangeVelocity() {
    std::printf("Simulation d'un echo puis double FFT + CFAR\n");
    Params p;
    Simulator sim(p, 7);
    RangeDopplerProcessor rdp(p);
    CfarConfig cfg;

    struct Case {
        double range, vr;
    } cases[] = {{200.0, 10.0}, {87.0, -23.0}, {420.0, 3.0}};

    for (const Case& c : cases) {
        Target t;
        t.id = 1;
        t.rcs = 10;
        t.x = 0;
        t.y = c.range;          // plein nord
        t.vy = c.vr;            // vitesse purement radiale
        Frame f;
        std::vector<EchoInfo> echoes;
        sim.generate({t}, 0.0, 0.0, f, &echoes);
        RangeDopplerMap map;
        rdp.process(f, map);
        auto dets = cfar2d(map, p, cfg, 0.0, 0.0);

        const Detection* best = nullptr;
        for (const auto& d : dets)
            if (!best || d.snrDb > best->snrDb)
                best = &d;
        CHECK(best != nullptr, "cible R=%.0f m v=%+.0f m/s detectee (%zu detection(s))", c.range, c.vr,
              dets.size());
        if (!best)
            continue;
        CHECK(std::abs(best->range - c.range) < p.rangeResolution(), "distance estimee %.1f m (vraie %.1f, resolution %.1f)",
              best->range, c.range, p.rangeResolution());
        CHECK(std::abs(best->radialVelocity - c.vr) < p.velocityResolution(),
              "vitesse estimee %+.2f m/s (vraie %+.2f, resolution %.2f)", best->radialVelocity, c.vr,
              p.velocityResolution());
    }
}

static void testCfarFalseAlarms() {
    std::printf("CFAR sur du bruit seul\n");
    Params p;
    Simulator sim(p, 99);
    RangeDopplerProcessor rdp(p);
    CfarConfig cfg;
    size_t total = 0;
    const int frames = 20;
    for (int i = 0; i < frames; ++i) {
        Frame f;
        sim.generate({}, 0.0, 0.0, f);
        RangeDopplerMap map;
        rdp.process(f, map);
        total += cfar2d(map, p, cfg, 0.0, 0.0).size();
    }
    const double cells = static_cast<double>(frames) * p.samplesPerChirp * p.chirpsPerFrame;
    CHECK(total <= 10, "%zu fausse(s) alarme(s) sur %.0f cases (pfa visee %.0e)", total, cells, cfg.pfa);
}

static void testKalman() {
    std::printf("Pistage : Kalman sur des plots bruites\n");
    Params p;
    Tracker tracker(p);
    std::mt19937 rng(3);
    std::normal_distribution<double> g;

    // Cible en ligne droite, vue une fois par tour (4 s).
    double x = -150, y = 300;
    const double vx = 4, vy = -3;
    double measSq = 0, estSq = 0;
    int nStats = 0;
    for (int scan = 0; scan < 25; ++scan) {
        const double t = scan * p.rotationPeriod;
        const double px = x + vx * t, py = y + vy * t;
        const double r = std::hypot(px, py);
        const double az = std::atan2(px, py);
        // Bruit de mesure realiste : distance +/- 1 m, azimut +/- 0.75 deg
        const double rm = r + 1.0 * g(rng);
        const double am = az + degToRad(0.75) * g(rng);
        Detection d;
        d.range = rm;
        d.azimuth = am;
        d.time = t;
        d.x = rm * std::sin(am);
        d.y = rm * std::cos(am);
        d.radialVelocity = (px * vx + py * vy) / r;
        d.snrDb = 20;
        tracker.process({d});
        tracker.prune(t);
        // Sur les 15 derniers tours (filtre converge), on compare l'erreur
        // du plot brut et celle de la position filtree.
        if (scan >= 10 && !tracker.tracks().empty()) {
            const auto& k = tracker.tracks().front().kf;
            measSq += (d.x - px) * (d.x - px) + (d.y - py) * (d.y - py);
            estSq += (k.x[0] - px) * (k.x[0] - px) + (k.x[1] - py) * (k.x[1] - py);
            ++nStats;
        }
    }
    CHECK(tracker.tracks().size() == 1, "une seule piste creee (%zu)", tracker.tracks().size());
    if (tracker.tracks().empty())
        return;
    const Track& tr = tracker.tracks().front();
    const double tEnd = 24 * p.rotationPeriod;
    const double ex = tr.kf.x[0] - (x + vx * tEnd), ey = tr.kf.x[1] - (y + vy * tEnd);
    const double evx = tr.kf.x[2] - vx, evy = tr.kf.x[3] - vy;
    CHECK(tr.confirmed, "piste confirmee apres %d plots", tr.hits);
    const double rmsMeas = std::sqrt(measSq / nStats), rmsEst = std::sqrt(estSq / nStats);
    CHECK(std::hypot(ex, ey) < 8.0, "erreur de position finale %.2f m", std::hypot(ex, ey));
    CHECK(std::hypot(evx, evy) < 1.5, "erreur de vitesse finale %.2f m/s", std::hypot(evx, evy));
    CHECK(rmsEst < rmsMeas, "le filtre lisse : erreur RMS %.2f m (plots bruts : %.2f m)", rmsEst, rmsMeas);
}

// Chaine complete, comme dans l'application : 5 cibles, antenne tournante,
// simulation -> double FFT -> CFAR -> extracteur de plots -> pistage.
static void testScenario() {
    std::printf("Scenario complet : 5 cibles pendant 40 s\n");
    Params p;
    Simulator sim(p, 11);
    RangeDopplerProcessor rdp(p);
    CfarConfig cfg;
    PlotExtractor extractor(p);
    Tracker tracker(p);

    std::vector<Target> targets;
    const double init[5][6] = {
        {-250, 250, 12, -4, 10, 0}, {180, 330, -3, -14, 30, 0}, {300, -150, 0, 18, 5, degToRad(3.0)},
        {-120, -380, 6, 9, 50, 0},  {60, 120, -9, 2, 1, 0},
    };
    for (int i = 0; i < 5; ++i) {
        Target t;
        t.id = i + 1;
        t.x = init[i][0];
        t.y = init[i][1];
        t.vx = init[i][2];
        t.vy = init[i][3];
        t.rcs = init[i][4];
        t.turnRate = init[i][5];
        targets.push_back(t);
    }

    // (Au-dela, la plus petite cible (1 m2) s'eloigne a plus de 400 m, ou son
    // echo passe sous le seuil de detection : sa piste finit par etre perdue.)
    const double dt = 1.0 / 60, duration = 40;
    double az = 0;
    size_t plots = 0;
    Frame f;
    RangeDopplerMap map;
    for (double t = 0; t < duration; t += dt) {
        for (auto& tg : targets)
            moveTarget(tg, dt, 40.0, 0.75 * p.maxRange());
        az = wrapAngle(az + kTwoPi * dt / p.rotationPeriod);
        sim.generate(targets, az, t, f);
        rdp.process(f, map);
        const auto dets = cfar2d(map, p, cfg, az, t);
        const auto pl = extractor.push(dets, az);
        plots += pl.size();
        tracker.process(pl);
        tracker.prune(t);
    }
    const double scans = duration / p.rotationPeriod;
    int confirmed = 0;
    for (const auto& tr : tracker.tracks())
        confirmed += tr.confirmed;
    CHECK(plots / scans < 6.0, "%.1f plots par tour pour 5 cibles", plots / scans);
    CHECK(confirmed == 5, "%d pistes confirmees (5 attendues)", confirmed);
    CHECK(tracker.tracks().size() <= 7, "%zu pistes au total (confirmees + tentatives)", tracker.tracks().size());

    // Chaque cible doit avoir une piste confirmee proche.
    int matched = 0;
    for (const auto& tg : targets) {
        for (const auto& tr : tracker.tracks()) {
            const auto x = tracker.predicted(tr, duration);
            if (tr.confirmed && std::hypot(x[0] - tg.x, x[1] - tg.y) < 15.0) {
                ++matched;
                break;
            }
        }
    }
    CHECK(matched == 5, "%d cible(s) sur 5 suivie(s) a moins de 15 m", matched);
}

int main() {
    Params p;
    std::printf("Radar FMCW : resolution %.2f m, portee max %.0f m, vitesse max +/-%.1f m/s, resolution vitesse %.2f m/s\n\n",
                p.rangeResolution(), p.maxRange(), p.maxVelocity(), p.velocityResolution());
    testFft();
    testRangeVelocity();
    testCfarFalseAlarms();
    testKalman();
    testScenario();
    std::printf("\n%s (%d echec(s))\n", failures ? "ECHEC" : "TOUT EST OK", failures);
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
