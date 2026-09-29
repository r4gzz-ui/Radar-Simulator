#include "core/target.hpp"

#include <cmath>

#include "core/params.hpp"

namespace radar {

namespace {
constexpr double kSteerRate = 0.2;   // rad/s (~11 deg/s) pour le demi-tour
} // namespace

double wrapAngle(double a) {
    a = std::fmod(a, kTwoPi);
    return a < 0 ? a + kTwoPi : a;
}

double angleDiff(double a, double b) {
    double d = std::fmod(a - b + kPi, kTwoPi);
    if (d < 0)
        d += kTwoPi;
    return d - kPi;
}

double Target::range() const { return std::hypot(x, y); }

double Target::azimuth() const { return wrapAngle(std::atan2(x, y)); }

double Target::radialVelocity() const {
    const double r = range();
    return r > 1e-9 ? (x * vx + y * vy) / r : 0.0;
}

double Target::speed() const { return std::hypot(vx, vy); }

void moveTarget(Target& t, double dt, double minRange, double maxRange) {
    // Virage : propre a la cible, plus un demi-tour progressif si elle sort
    // de la zone [minRange, maxRange] (un rebond instantane serait une
    // manoeuvre impossible, que le filtre de Kalman ne pourrait pas suivre).
    double turn = t.turnRate;
    const double r = t.range();
    const double vr = t.radialVelocity();
    if ((r > maxRange && vr > 0) || (r < minRange && vr < 0)) {
        // On tourne du cote qui ramene la vitesse vers le centre (ou vers
        // l'exterieur si la cible est trop proche du radar).
        const double cross = t.x * t.vy - t.y * t.vx;   // > 0 : la cible tourne dans le sens trigo
        const double toward = (r > maxRange) ? 1.0 : -1.0;
        turn += toward * (cross > 0 ? 1.0 : -1.0) * kSteerRate;
    }
    if (turn != 0) {
        const double a = turn * dt;
        const double c = std::cos(a), s = std::sin(a);
        const double vx = t.vx * c - t.vy * s;
        const double vy = t.vx * s + t.vy * c;
        t.vx = vx;
        t.vy = vy;
    }
    t.x += t.vx * dt;
    t.y += t.vy * dt;
}

} // namespace radar
