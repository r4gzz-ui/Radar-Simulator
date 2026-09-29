#pragma once

#include <vector>

#include "core/fft.hpp"
#include "core/params.hpp"
#include "core/simulator.hpp"

namespace radar {

// Carte distance / vitesse d'un dwell.
// power[d * nRange + r] : puissance (lineaire) de la case distance r, vitesse d.
// L'axe Doppler est recentre (fftshift) : d = nDoppler / 2 correspond a v = 0,
// d < nDoppler / 2 aux cibles qui s'approchent, d > nDoppler / 2 a celles qui s'eloignent.
struct RangeDopplerMap {
    int nRange = 0;
    int nDoppler = 0;
    std::vector<float> power;
    std::vector<float> rangeProfile;  // moyenne sur les chirps de |FFT distance|^2
    float noiseLevel = 1.0f;          // puissance moyenne estimee du bruit

    float at(int d, int r) const { return power[static_cast<size_t>(d) * nRange + r]; }

    double rangeOfBin(double r, const Params& p) const { return r * p.rangeResolution(); }
    double velocityOfBin(double d, const Params& p) const {
        return (d - nDoppler / 2) * p.velocityResolution();
    }
};

class RangeDopplerProcessor {
public:
    explicit RangeDopplerProcessor(const Params& params);

    // FFT n1 : sur chaque chirp (temps rapide) -> distance.
    // FFT n2 : sur chaque case distance, a travers les chirps (temps lent) -> vitesse.
    void process(const Frame& frame, RangeDopplerMap& out);

private:
    const Params& p_;
    std::vector<float> winRange;
    std::vector<float> winDoppler;
    std::vector<cfloat> cube;    // resultat de la FFT distance, [m][r]
    std::vector<cfloat> column;  // tampon pour la FFT Doppler
    std::vector<float> median;   // tampon pour l'estimation du bruit
};

} // namespace radar
