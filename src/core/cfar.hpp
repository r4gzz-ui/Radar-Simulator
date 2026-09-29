#pragma once

#include <cstdint>
#include <vector>

#include "core/params.hpp"
#include "core/range_doppler.hpp"

namespace radar {

// CFAR "Cell Averaging" en 2D.
// Pour chaque case testee, on estime le bruit local en moyennant les cases
// d'entrainement qui l'entourent (en excluant les cases de garde, trop proches,
// qui contiennent peut-etre encore un bout de la cible). La case est declaree
// "detection" si sa puissance depasse alpha * bruit_local, ou alpha est choisi
// pour garantir une probabilite de fausse alarme pfa constante.
struct CfarConfig {
    int guardRange = 2;
    int guardDoppler = 2;
    int trainRange = 8;
    int trainDoppler = 4;
    double pfa = 1e-6;
    int minRangeBin = 3;   // ignore les toutes premieres cases (fuite de l'emission)

    // Suppression des lobes secondaires : une cible tres forte "bave" sur les
    // cases voisines (lobes secondaires de la fenetre de Hann), et le CFAR
    // y verrait de fausses cibles. On ecarte donc tout pic situe pres d'un pic
    // beaucoup plus fort (plus de sidelobeDb dB d'ecart).
    int sidelobeRange = 12;
    int sidelobeDoppler = 8;
    double sidelobeDb = 20.0;
};

// Une detection = un pic de la carte distance/vitesse qui passe le seuil CFAR.
struct Detection {
    double range = 0;           // m (interpole entre les cases)
    double radialVelocity = 0;  // m/s
    double azimuth = 0;         // rad (direction du lobe)
    double snrDb = 0;           // rapport signal / bruit local
    double time = 0;            // s
    double x = 0, y = 0;        // position cartesienne (m)
    int rangeBin = 0;
    int dopplerBin = 0;
};

// hitMask (optionnel) : 1 pour chaque case au-dessus du seuil, [d][r].
std::vector<Detection> cfar2d(const RangeDopplerMap& map, const Params& p, const CfarConfig& cfg,
                              double azimuth, double time, std::vector<uint8_t>* hitMask = nullptr);

// Facteur de seuil du CA-CFAR pour n cases d'entrainement.
double cfarAlpha(int nTrain, double pfa);

// Regroupe les detections d'une meme cible vue sur plusieurs dwells
// consecutifs (le lobe balaye la cible) en un seul "plot", avec un azimut
// estime par barycentre pondere par le SNR. C'est ce que fait l'extracteur
// de plots d'un vrai radar.
//
// Regle "M sur N" : un vrai echo est vu sur plusieurs dwells consecutifs
// (le lobe met du temps a passer sur la cible), alors qu'une fausse alarme
// du bruit est isolee. On exige donc au moins minDwells detections.
class PlotExtractor {
public:
    explicit PlotExtractor(const Params& params, int minDwells = 2);

    // A appeler a chaque dwell. Renvoie les plots termines (le lobe est passe).
    std::vector<Detection> push(const std::vector<Detection>& dets, double beamAzimuth);
    void clear();

private:
    struct Pending {
        double wSum = 0;
        double range = 0, velocity = 0, snrLin = 0, bestSnrDb = 0;
        double sinSum = 0, cosSum = 0;   // barycentre d'angle (evite le saut 360 -> 0)
        double time = 0;
        double lastAzimuth = 0;
        int rangeBin = 0, dopplerBin = 0;
        int count = 0;
    };
    const Params& p_;
    int minDwells_;
    std::vector<Pending> pending_;
};

} // namespace radar
