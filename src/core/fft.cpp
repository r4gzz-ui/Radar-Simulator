#include "core/fft.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

#include "core/params.hpp"

namespace radar {

bool isPowerOfTwo(int n) { return n > 0 && (n & (n - 1)) == 0; }

void fft(cfloat* a, int n, bool inverse) {
    if (!isPowerOfTwo(n))
        throw std::invalid_argument("fft: la taille doit etre une puissance de 2");

    // 1. Permutation "bit-reversal" : l'element d'indice i va a l'indice dont
    //    les bits sont ceux de i dans l'ordre inverse. Apres ca, on peut
    //    combiner les sous-FFT par paires adjacentes.
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            std::swap(a[i], a[j]);
    }

    // 2. Papillons : a chaque etape, on fusionne des FFT de taille len/2 en
    //    FFT de taille len. X[k] = E[k] + w^k O[k], X[k+len/2] = E[k] - w^k O[k].
    for (int len = 2; len <= n; len <<= 1) {
        const double angle = kTwoPi / len * (inverse ? 1.0 : -1.0);
        const std::complex<double> wStep(std::cos(angle), std::sin(angle));
        const int half = len / 2;
        for (int i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (int k = 0; k < half; ++k) {
                const cfloat u = a[i + k];
                const cfloat v = a[i + k + half] * cfloat(w);
                a[i + k] = u + v;
                a[i + k + half] = u - v;
                w *= wStep;
            }
        }
    }

    if (inverse) {
        const float scale = 1.0f / static_cast<float>(n);
        for (int i = 0; i < n; ++i)
            a[i] *= scale;
    }
}

void fft(std::vector<cfloat>& data, bool inverse) {
    fft(data.data(), static_cast<int>(data.size()), inverse);
}

std::vector<float> hannWindow(int n) {
    std::vector<float> w(n);
    for (int i = 0; i < n; ++i)
        w[i] = static_cast<float>(0.5 - 0.5 * std::cos(kTwoPi * i / n));
    return w;
}

} // namespace radar
