#pragma once

// Parametres du radar FMCW simule.
//
// Un radar FMCW emet des "chirps" : une rampe de frequence de fc a fc + B
// pendant Tc. L'echo d'une cible a la distance R revient avec un retard
// tau = 2R/c. En melangeant l'emission et la reception, on obtient un signal
// de battement (IF) de frequence fb = S * tau = 2 * S * R / c, ou S = B / Tc
// est la pente du chirp. Mesurer fb (par FFT) donne donc la distance.
//
// Sur plusieurs chirps successifs, une cible qui bouge fait tourner la phase
// de l'echo de 4 * pi * v * Tc / lambda par chirp : une seconde FFT sur les
// chirps donne la frequence Doppler fd = 2v / lambda, donc la vitesse radiale.

namespace radar {

constexpr double kC = 299792458.0;
constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;

inline double degToRad(double d) { return d * kPi / 180.0; }
inline double radToDeg(double r) { return r * 180.0 / kPi; }

struct Params {
    // --- Forme d'onde ---
    double fc = 77e9;            // frequence porteuse (Hz)
    double bandwidth = 75e6;     // excursion du chirp B (Hz)
    double chirpTime = 20e-6;    // duree d'un chirp Tc (s)
    int samplesPerChirp = 256;   // Ns : echantillons par chirp (puissance de 2)
    int chirpsPerFrame = 64;     // Nc : chirps par "dwell" (puissance de 2)

    // --- Antenne / balayage ---
    double beamWidthDeg = 3.0;   // ouverture du lobe a -3 dB (aller-retour)
    double rotationPeriod = 4.0; // duree d'un tour d'antenne (s)
    double maxAzStepDeg = 1.0;   // pas d'azimut max entre deux dwells

    // --- Bilan de liaison simplifie ---
    // Amplitude de l'echo = echoGain * sqrt(RCS) / R^2 * gain_de_lobe.
    // Avec noiseSigma = 1, une cible de 1 m2 a 500 m sort ~10 dB au-dessus
    // du bruit apres les deux FFT.
    double echoGain = 12500.0;
    double noiseSigma = 1.0;     // ecart-type du bruit (par voie I et Q)

    // --- Grandeurs derivees ---
    double lambda() const { return kC / fc; }
    double slope() const { return bandwidth / chirpTime; }
    double sampleRate() const { return samplesPerChirp / chirpTime; }

    // Une case de FFT distance = fs / Ns = 1 / Tc en frequence de battement,
    // soit c / (2B) en distance.
    double rangeResolution() const { return kC / (2.0 * bandwidth); }

    // Signal IF complexe (I/Q) : toutes les cases 0..Ns-1 sont exploitables.
    double maxRange() const { return rangeResolution() * (samplesPerChirp - 1); }

    // Une case Doppler = 1 / (Nc * Tc) Hz, soit lambda / (2 * Nc * Tc) m/s.
    double velocityResolution() const {
        return lambda() / (2.0 * chirpsPerFrame * chirpTime);
    }

    // Au-dela de +/- lambda / (4 Tc), la phase tourne de plus de pi par
    // chirp : la vitesse se replie (ambiguite Doppler).
    double maxVelocity() const { return lambda() / (4.0 * chirpTime); }

    double beamWidth() const { return degToRad(beamWidthDeg); }
};

} // namespace radar
