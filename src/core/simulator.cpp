#include "core/simulator.hpp"

#include <algorithm>
#include <cmath>

namespace radar {

double beamGain(double delta, double beamWidth) {
    const double u = delta / beamWidth;
    return std::exp(-2.0 * std::log(2.0) * u * u);
}

Simulator::Simulator(const Params& params, unsigned seed) : p_(params), rng_(seed) {}

void Simulator::generate(const std::vector<Target>& targets, double azimuth, double time,
                         Frame& out, std::vector<EchoInfo>* echoes) {
    const int ns = p_.samplesPerChirp;
    const int nc = p_.chirpsPerFrame;
    out.ns = ns;
    out.nc = nc;
    out.azimuth = azimuth;
    out.time = time;
    out.data.resize(static_cast<size_t>(ns) * nc);
    if (echoes)
        echoes->clear();

    // 1. Bruit thermique : gaussien, independant sur I et Q.
    const float sigma = static_cast<float>(p_.noiseSigma);
    for (cfloat& s : out.data)
        s = cfloat(sigma * gauss_(rng_), sigma * gauss_(rng_));

    const double fs = p_.sampleRate();
    const double tc = p_.chirpTime;
    const double lambda = p_.lambda();
    const double bw = p_.beamWidth();

    // 2. Un echo par cible dans le lobe.
    for (const Target& t : targets) {
        const double delta = angleDiff(t.azimuth(), azimuth);
        if (std::abs(delta) > 2.5 * bw)
            continue;
        const double r = t.range();
        if (r < 1.0 || r > p_.maxRange())
            continue;

        const double gain = beamGain(delta, bw);
        // Equation radar (en amplitude) : la puissance recue decroit en 1/R^4,
        // donc l'amplitude en 1/R^2. Elle croit avec la racine de la RCS.
        const double amp = p_.echoGain * std::sqrt(t.rcs) / (std::max(r, 5.0) * std::max(r, 5.0)) * gain;
        const double vr = t.radialVelocity();

        const double tau = 2.0 * r / kC;          // retard aller-retour
        const double fb = p_.slope() * tau;        // frequence de battement
        const double fd = 2.0 * vr / lambda;       // decalage Doppler

        if (echoes) {
            EchoInfo e;
            e.targetId = t.id;
            e.range = r;
            e.radialVelocity = vr;
            e.delay = tau;
            e.beatFreq = fb;
            e.dopplerFreq = fd;
            e.amplitude = amp;
            e.beamGain = gain;
            echoes->push_back(e);
        }

        // Phase du signal IF (emission x conjugue de la reception) :
        //   phi(n, m) = 2pi (fb + fd) n / fs      <- dans le chirp : donne la distance
        //             + 2pi fd m Tc              <- de chirp en chirp : donne la vitesse
        //             + 4pi R / lambda           <- phase initiale
        // La phase initiale est calculee en double : R / lambda vaut ~1e5.
        const double phase0 = std::fmod(4.0 * kPi * r / lambda, kTwoPi);
        const double fastStep = kTwoPi * (fb + fd) / fs;
        const std::complex<double> fastRot(std::cos(fastStep), std::sin(fastStep));

        for (int m = 0; m < nc; ++m) {
            const double chirpPhase = phase0 + kTwoPi * fd * m * tc;
            std::complex<double> ph = std::polar(amp, chirpPhase);
            cfloat* row = out.data.data() + static_cast<size_t>(m) * ns;
            for (int n = 0; n < ns; ++n) {
                row[n] += cfloat(ph);
                ph *= fastRot;
            }
        }
    }
}

} // namespace radar
