#pragma once

#include <vector>

namespace radar {

// Repere : le radar est a l'origine, x vers l'est, y vers le nord (metres).
// L'azimut se compte depuis le nord, dans le sens horaire, comme sur un PPI.
struct Target {
    int id = 0;
    double x = 0, y = 0;     // position (m)
    double vx = 0, vy = 0;   // vitesse (m/s)
    double rcs = 10;         // surface equivalente radar (m2)
    double turnRate = 0;     // vitesse de virage (rad/s), 0 = ligne droite

    double range() const;
    double azimuth() const;          // [0, 2pi)
    double radialVelocity() const;   // dR/dt, > 0 si la cible s'eloigne
    double speed() const;
};

// Fait avancer la cible de dt. Si elle sort de [minRange, maxRange], elle
// fait demi-tour progressivement.
void moveTarget(Target& t, double dt, double minRange, double maxRange);

double wrapAngle(double a);                 // -> [0, 2pi)
double angleDiff(double a, double b);       // a - b dans [-pi, pi)

} // namespace radar
