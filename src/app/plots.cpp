#include "app/plots.hpp"

#include <algorithm>
#include <cmath>

#include "app/theme.hpp"

namespace app {

using radar::kTwoPi;

namespace {

struct Axis {
    double min, max;
    double norm(double v) const { return (v - min) / (max - min); }
};

// Cadre d'un panneau ; renvoie la zone de trace interieure.
Rectangle panelFrame(Rectangle r, const char* title, const char* subtitle) {
    DrawRectangleRec(r, theme::panel);
    DrawRectangleLinesEx(r, 1.0f, theme::border);
    DrawText(title, static_cast<int>(r.x) + 8, static_cast<int>(r.y) + 5, 14, theme::text);
    if (subtitle && *subtitle) {
        const int w = MeasureText(subtitle, 10);
        DrawText(subtitle, static_cast<int>(r.x + r.width) - w - 8, static_cast<int>(r.y) + 8, 10, theme::dimText);
    }
    return {r.x + 48, r.y + 24, r.width - 58, r.height - 42};
}

float mapX(const Rectangle& p, const Axis& a, double v) {
    return p.x + static_cast<float>(a.norm(v)) * p.width;
}
float mapY(const Rectangle& p, const Axis& a, double v) {
    return p.y + p.height - static_cast<float>(a.norm(v)) * p.height;
}

void dashedV(float x, float y0, float y1, Color c) {
    for (float y = y0; y < y1; y += 8)
        DrawLineV({x, y}, {x, std::min(y + 4, y1)}, c);
}

void dashedH(float y, float x0, float x1, Color c) {
    for (float x = x0; x < x1; x += 8)
        DrawLineV({x, y}, {std::min(x + 4, x1), y}, c);
}

// Palette de la carte distance/vitesse : noir -> bleu -> vert -> jaune -> rouge.
Color heat(float v) {
    struct Stop {
        float at;
        float r, g, b;
    };
    static const Stop stops[] = {
        {0.00f, 0, 0, 0},     {0.25f, 0, 35, 85},     {0.50f, 0, 150, 95},
        {0.75f, 210, 220, 40}, {1.00f, 255, 70, 30},
    };
    v = std::clamp(v, 0.0f, 1.0f);
    for (int i = 0; i < 4; ++i) {
        if (v <= stops[i + 1].at) {
            const float t = (v - stops[i].at) / (stops[i + 1].at - stops[i].at);
            return Color{static_cast<unsigned char>(stops[i].r + t * (stops[i + 1].r - stops[i].r)),
                         static_cast<unsigned char>(stops[i].g + t * (stops[i + 1].g - stops[i].g)),
                         static_cast<unsigned char>(stops[i].b + t * (stops[i + 1].b - stops[i].b)), 255};
        }
    }
    return Color{255, 70, 30, 255};
}

const radar::EchoInfo* mainEcho(const DwellView& v) {
    if (!v.echoes || v.echoes->empty())
        return nullptr;
    const radar::EchoInfo* best = nullptr;
    for (const auto& e : *v.echoes) {
        if (e.targetId == v.focusTarget)
            return &e;
        if (!best || e.amplitude > best->amplitude)
            best = &e;
    }
    return best;
}

// Texte sur un fond sombre, lisible par-dessus les courbes.
void label(const char* txt, float x, float y, int size, Color c) {
    const int w = MeasureText(txt, size);
    DrawRectangle(static_cast<int>(x) - 2, static_cast<int>(y) - 1, w + 4, size + 2, Fade(theme::panel, 0.85f));
    DrawText(txt, static_cast<int>(x), static_cast<int>(y), size, c);
}

void noData(Rectangle p, const char* msg) {
    const int w = MeasureText(msg, 14);
    DrawText(msg, static_cast<int>(p.x + p.width / 2) - w / 2, static_cast<int>(p.y + p.height / 2) - 7, 14,
             theme::dimText);
}

} // namespace

// ---------------------------------------------------------------------------

void Plots::init(const radar::Params& params) {
    texW_ = params.samplesPerChirp;
    texH_ = params.chirpsPerFrame;
    pixels_.assign(static_cast<size_t>(texW_) * texH_, BLACK);
    Image img = GenImageColor(texW_, texH_, BLACK);
    rdTexture_ = LoadTextureFromImage(img);
    UnloadImage(img);
    SetTextureFilter(rdTexture_, TEXTURE_FILTER_POINT);
    ready_ = true;
}

void Plots::unload() {
    if (ready_)
        UnloadTexture(rdTexture_);
    ready_ = false;
}

void Plots::draw(Rectangle area, const DwellView& view) {
    const float gap = 8;
    const float h = area.height - 3 * gap;
    const float fr[] = {0.27f, 0.19f, 0.20f, 0.34f};
    float y = area.y;
    Rectangle rects[4];
    for (int i = 0; i < 4; ++i) {
        rects[i] = {area.x, y, area.width, h * fr[i]};
        y += h * fr[i] + gap;
    }
    drawWaves(rects[0], view);
    drawBeat(rects[1], view);
    drawRangeProfile(rects[2], view);
    drawRangeDoppler(rects[3], view);
}

// 1. Ondes emises et recues -------------------------------------------------

void Plots::drawWaves(Rectangle r, const DwellView& v) {
    const radar::Params& p = *v.params;
    const radar::EchoInfo* echo = mainEcho(v);
    Rectangle in = panelFrame(r, "1. Emission / reception", v.source);

    const float split = in.height * 0.42f;
    Rectangle top{in.x, in.y, in.width, split - 6};
    Rectangle bot{in.x, in.y + split + 4, in.width, in.height - split - 4};

    // --- Haut : les ondes elles-memes, sur un chirp. Les vraies frequences
    // (77 GHz) sont invisibles a cette echelle : on dessine un chirp "ralenti"
    // qui monte de 3 a 17 oscillations, avec le vrai rapport retard / Tc.
    const double tauRatio = echo ? echo->delay / p.chirpTime : 0.0;
    const auto chirpWave = [](double u) { return std::cos(kTwoPi * (3.0 * u + 0.5 * 14.0 * u * u)); };
    const int n = 500;
    Vector2 prevTx{}, prevRx{};
    for (int i = 0; i <= n; ++i) {
        const double u = static_cast<double>(i) / n;
        const float x = top.x + static_cast<float>(u) * top.width;
        const float yTx = top.y + top.height * 0.27f - static_cast<float>(chirpWave(u)) * top.height * 0.22f;
        const bool received = echo && u >= tauRatio;
        const double rx = received ? 0.8 * chirpWave(u - tauRatio) : 0.0;
        const float yRx = top.y + top.height * 0.76f - static_cast<float>(rx) * top.height * 0.22f;
        const Vector2 ptx{x, yTx}, prx{x, yRx};
        if (i > 0) {
            DrawLineV(prevTx, ptx, theme::tx);
            DrawLineV(prevRx, prx, received ? theme::rx : theme::dimText);
        }
        prevTx = ptx;
        prevRx = prx;
    }
    DrawText("Tx", static_cast<int>(in.x) - 26, static_cast<int>(top.y + top.height * 0.27f) - 6, 12, theme::tx);
    DrawText("Rx", static_cast<int>(in.x) - 26, static_cast<int>(top.y + top.height * 0.76f) - 6, 12, theme::rx);
    if (echo) {
        const float xTau = top.x + static_cast<float>(tauRatio) * top.width;
        dashedV(xTau, top.y, top.y + top.height, Fade(theme::rx, 0.6f));
    } else {
        label("aucun echo dans le lobe", top.x + top.width - 140, top.y + top.height - 10, 10, theme::dimText);
    }

    // --- Bas : frequence instantanee sur deux chirps (dents de scie).
    const Axis ax{0, 2 * p.chirpTime};
    const Axis ay{0, p.bandwidth};
    DrawRectangleLinesEx(bot, 1, theme::grid);
    DrawText("fc+B", static_cast<int>(bot.x) - 40, static_cast<int>(bot.y) - 2, 10, theme::gridText);
    DrawText("fc", static_cast<int>(bot.x) - 20, static_cast<int>(bot.y + bot.height) - 8, 10, theme::gridText);
    for (int k = 0; k <= 2; ++k) {
        const float x = mapX(bot, ax, k * p.chirpTime);
        DrawText(TextFormat("%.0f us", k * p.chirpTime * 1e6), static_cast<int>(x) - 10,
                 static_cast<int>(bot.y + bot.height) + 3, 10, theme::gridText);
    }

    const auto ramp = [&](double t0, Color c, float thick) {
        // Rampe emise/recue qui commence a t0 ; decoupee a la fenetre [0, 2Tc].
        const double a = std::max(t0, 0.0), b = std::min(t0 + p.chirpTime, 2 * p.chirpTime);
        if (b <= a)
            return;
        const double fa = p.bandwidth * (a - t0) / p.chirpTime;
        const double fb = p.bandwidth * (b - t0) / p.chirpTime;
        DrawLineEx({mapX(bot, ax, a), mapY(bot, ay, fa)}, {mapX(bot, ax, b), mapY(bot, ay, fb)}, thick, c);
    };
    for (int m = 0; m < 2; ++m)
        ramp(m * p.chirpTime, theme::tx, 2.0f);
    if (v.echoes) {
        for (const auto& e : *v.echoes) {
            const bool focus = echo && e.targetId == echo->targetId;
            for (int m = -1; m < 2; ++m)
                ramp(m * p.chirpTime + e.delay, focus ? theme::rx : Fade(theme::otherRx, 0.7f), focus ? 2.0f : 1.0f);
        }
    }
    if (echo) {
        // Le battement = ecart vertical constant entre Tx et Rx.
        const double t = 0.6 * p.chirpTime;
        const float x = mapX(bot, ax, t);
        const float y0 = mapY(bot, ay, p.bandwidth * 0.6);
        const float y1 = mapY(bot, ay, p.bandwidth * (0.6 - echo->delay / p.chirpTime));
        DrawLineEx({x, y0}, {x, y1}, 1.5f, WHITE);
        label(TextFormat("fb = S.tau = %.3f MHz", echo->beatFreq / 1e6), x + 6, (y0 + y1) / 2 - 5, 10, WHITE);
        // Coin haut-gauche : zone laissee libre par les rampes.
        label(TextFormat("tau = 2R/c = %.2f us  (cible #%d a %.1f m)", echo->delay * 1e6, echo->targetId, echo->range),
              bot.x + 6, bot.y + 4, 10, theme::rx);
        label(TextFormat("Doppler fd = 2v/lambda = %+.1f kHz (invisible a cette echelle)", echo->dopplerFreq / 1e3),
              bot.x + 6, bot.y + 17, 10, theme::dimText);
    }
}

// 2. Signal de battement ----------------------------------------------------

void Plots::drawBeat(Rectangle r, const DwellView& v) {
    const radar::Params& p = *v.params;
    Rectangle in = panelFrame(r, "2. Signal de battement (IF) - chirp n0", "I (vert)  Q (sombre)");
    if (!v.frame || v.frame->data.empty()) {
        noData(in, "pas de donnees");
        return;
    }
    const radar::cfloat* s = v.frame->chirp(0);
    const int ns = v.frame->ns;
    float peak = 1e-6f;
    for (int i = 0; i < ns; ++i)
        peak = std::max({peak, std::abs(s[i].real()), std::abs(s[i].imag())});

    dashedH(in.y + in.height / 2, in.x, in.x + in.width, theme::grid);
    const auto px = [&](int i) { return in.x + in.width * i / (ns - 1); };
    const auto py = [&](float val) { return in.y + in.height / 2 - val / peak * in.height * 0.48f; };
    for (int i = 1; i < ns; ++i) {
        DrawLineV({px(i - 1), py(s[i - 1].imag())}, {px(i), py(s[i].imag())}, Fade(theme::phosphor, 0.3f));
        DrawLineV({px(i - 1), py(s[i - 1].real())}, {px(i), py(s[i].real())}, theme::phosphor);
    }
    DrawText("0", static_cast<int>(in.x) - 4, static_cast<int>(in.y + in.height) + 3, 10, theme::gridText);
    DrawText(TextFormat("%d ech. (%.0f us)", ns, p.chirpTime * 1e6), static_cast<int>(in.x + in.width) - 90,
             static_cast<int>(in.y + in.height) + 3, 10, theme::gridText);

    if (const radar::EchoInfo* e = mainEcho(v)) {
        const double cycles = e->beatFreq * p.chirpTime;
        const double snr = 20.0 * std::log10(e->amplitude / (p.noiseSigma * std::sqrt(2.0)));
        DrawText(TextFormat("cible #%d : %.1f oscillations par chirp -> R = fb.c/(2S) = %.1f m   |   SNR par echantillon %.1f dB",
                            e->targetId, cycles, e->range, snr),
                 static_cast<int>(in.x) + 4, static_cast<int>(in.y) + 2, 10, theme::rx);
    } else {
        DrawText("bruit seul", static_cast<int>(in.x) + 4, static_cast<int>(in.y) + 2, 10, theme::dimText);
    }
}

// 3. Profil distance --------------------------------------------------------

void Plots::drawRangeProfile(Rectangle r, const DwellView& v) {
    const radar::Params& p = *v.params;
    Rectangle in = panelFrame(r, "3. FFT n1 : profil distance", "puissance moyenne sur les chirps (dB / bruit)");
    if (!v.rd || v.rd->rangeProfile.empty()) {
        noData(in, "pas de donnees");
        return;
    }
    const Axis ax{0, p.maxRange()};
    const Axis ay{-5, 50};
    for (double m = 0; m <= p.maxRange(); m += 100) {
        const float x = mapX(in, ax, m);
        DrawLineV({x, in.y}, {x, in.y + in.height}, theme::grid);
        DrawText(TextFormat("%.0f m", m), static_cast<int>(x) - 12, static_cast<int>(in.y + in.height) + 3, 10,
                 theme::gridText);
    }
    for (int db = 0; db <= 50; db += 20) {
        const float y = mapY(in, ay, db);
        dashedH(y, in.x, in.x + in.width, theme::grid);
        DrawText(TextFormat("%d dB", db), static_cast<int>(in.x) - 40, static_cast<int>(y) - 5, 10, theme::gridText);
    }

    const auto& prof = v.rd->rangeProfile;
    const int nr = static_cast<int>(prof.size());
    const float noise = v.rd->noiseLevel;
    const auto dbAt = [&](int i) {
        return std::clamp(10.0 * std::log10(std::max(prof[i] / noise, 1e-12f)), ay.min, ay.max);
    };
    for (int i = 1; i < nr; ++i) {
        DrawLineV({mapX(in, ax, v.rd->rangeOfBin(i - 1, p)), mapY(in, ay, dbAt(i - 1))},
                  {mapX(in, ax, v.rd->rangeOfBin(i, p)), mapY(in, ay, dbAt(i))}, theme::phosphor);
    }

    if (v.echoes) {
        for (const auto& e : *v.echoes) {
            if (!v.showTruth && e.targetId != v.focusTarget)
                continue;
            dashedV(mapX(in, ax, e.range), in.y, in.y + in.height, Fade(theme::rx, 0.7f));
        }
    }
    if (v.detections) {
        for (const auto& d : *v.detections) {
            const float x = mapX(in, ax, d.range);
            const float y = mapY(in, ay, dbAt(d.rangeBin));
            DrawTriangle({x, y - 4}, {x - 5, y - 12}, {x + 5, y - 12}, theme::detection);
            DrawText(TextFormat("%.1f m", d.range), static_cast<int>(x) + 6, static_cast<int>(y) - 14, 10,
                     theme::detection);
        }
    }
}

// 4. Carte distance / vitesse ------------------------------------------------

void Plots::drawRangeDoppler(Rectangle r, const DwellView& v) {
    const radar::Params& p = *v.params;
    Rectangle in = panelFrame(r, "4. FFT n2 : carte distance / vitesse + CFAR",
                              "cercles : detections CFAR   x : position vraie");
    if (!v.rd || v.rd->power.empty() || !ready_) {
        noData(in, "pas de donnees");
        return;
    }
    const auto& rd = *v.rd;
    if (rd.nRange == texW_ && rd.nDoppler == texH_) {
        for (int d = 0; d < rd.nDoppler; ++d) {
            for (int rr = 0; rr < rd.nRange; ++rr) {
                const float db = 10.0f * std::log10(std::max(rd.at(d, rr) / rd.noiseLevel, 1e-12f));
                // Vitesses positives en haut : on retourne l'axe vertical.
                pixels_[static_cast<size_t>(rd.nDoppler - 1 - d) * texW_ + rr] = heat((db + 3.0f) / 40.0f);
            }
        }
        UpdateTexture(rdTexture_, pixels_.data());
    }
    DrawTexturePro(rdTexture_, {0, 0, static_cast<float>(texW_), static_cast<float>(texH_)}, in, {0, 0}, 0, WHITE);
    DrawRectangleLinesEx(in, 1, theme::border);

    const double vmax = p.maxVelocity();
    // Une case = un pixel de texture ; la case d est centree sur (d + 0.5).
    const auto xOf = [&](double range) {
        return in.x + static_cast<float>((range / p.rangeResolution() + 0.5) / rd.nRange) * in.width;
    };
    const auto yOf = [&](double vel) {
        const double bin = vel / p.velocityResolution() + rd.nDoppler / 2;
        return in.y + in.height - static_cast<float>((bin + 0.5) / rd.nDoppler) * in.height;
    };

    for (double m = 0; m <= p.maxRange(); m += 100)
        DrawText(TextFormat("%.0f m", m), static_cast<int>(xOf(m)) - 12, static_cast<int>(in.y + in.height) + 3, 10,
                 theme::gridText);
    for (int vel = -40; vel <= 40; vel += 20) {
        const float y = yOf(vel);
        DrawText(TextFormat("%+d", vel), static_cast<int>(in.x) - 30, static_cast<int>(y) - 5, 10, theme::gridText);
        DrawLineV({in.x - 4, y}, {in.x, y}, theme::gridText);
    }
    DrawText("m/s", static_cast<int>(in.x) - 30, static_cast<int>(yOf(-40)) + 10, 10, theme::gridText);
    dashedH(yOf(0), in.x, in.x + in.width, Fade(WHITE, 0.25f));

    if (v.echoes) {
        for (const auto& e : *v.echoes) {
            if (!v.showTruth && e.targetId != v.focusTarget)
                continue;
            if (std::abs(e.radialVelocity) > vmax)
                continue;   // repliee : on n'affiche pas la vraie position
            const Vector2 c{xOf(e.range), yOf(e.radialVelocity)};
            DrawLineEx({c.x - 5, c.y - 5}, {c.x + 5, c.y + 5}, 1.5f, theme::rx);
            DrawLineEx({c.x - 5, c.y + 5}, {c.x + 5, c.y - 5}, 1.5f, theme::rx);
        }
    }
    if (v.detections) {
        for (const auto& d : *v.detections) {
            const Vector2 c{xOf(d.range), yOf(d.radialVelocity)};
            DrawCircleLines(static_cast<int>(c.x), static_cast<int>(c.y), 9, theme::detection);
            const char* txt = TextFormat("%.1f m  %+.1f m/s  SNR %.0f dB", d.range, d.radialVelocity, d.snrDb);
            // Etiquette a gauche du cercle si elle depasserait du panneau.
            const float tx = c.x + 12 + MeasureText(txt, 10) > in.x + in.width ? c.x - 12 - MeasureText(txt, 10) : c.x + 12;
            label(txt, tx, c.y - 5, 10, theme::detection);
        }
    }
}

} // namespace app
