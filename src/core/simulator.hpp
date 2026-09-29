#pragma once

#include <random>
#include <vector>

#include "core/fft.hpp"
#include "core/params.hpp"
#include "core/target.hpp"

namespace radar {

// Un "dwell" : Nc chirps emis dans une meme direction d'antenne.
// data est rangee chirp par chirp : data[m * Ns + n] = echantillon n du chirp m.
struct Frame {
    int ns = 0;
    int nc = 0;
    double azimuth = 0;   // direction du lobe (rad)
    double time = 0;      // instant de la mesure (s)
    std::vector<cfloat> data;

    const cfloat* chirp(int m) const { return data.data() + m * ns; }
};

// Ce que "sait" le simulateur sur chaque echo (verite terrain, pour l'affichage).
struct EchoInfo {
    int targetId = 0;
    double range = 0;           // m
    double radialVelocity = 0;  // m/s
    double delay = 0;           // retard aller-retour tau = 2R/c (s)
    double beatFreq = 0;        // fb = S * tau (Hz)
    double dopplerFreq = 0;     // fd = 2 v / lambda (Hz)
    double amplitude = 0;       // amplitude de l'echo dans le signal IF
    double beamGain = 0;        // gain du lobe dans la direction de la cible (0..1)
};

// Gain aller-retour (en amplitude) du lobe gaussien pour un ecart d'azimut
// delta : vaut 1 dans l'axe et 1/sqrt(2) (-3 dB) a +/- beamWidth / 2.
double beamGain(double delta, double beamWidth);

class Simulator {
public:
    explicit Simulator(const Params& params, unsigned seed = 42);

    // Genere le signal IF recu pendant un dwell dans la direction azimuth :
    // somme des echos des cibles eclairees par le lobe + bruit blanc gaussien.
    void generate(const std::vector<Target>& targets, double azimuth, double time,
                  Frame& out, std::vector<EchoInfo>* echoes = nullptr);

private:
    const Params& p_;
    std::mt19937 rng_;
    std::normal_distribution<float> gauss_{0.0f, 1.0f};
};

} // namespace radar
