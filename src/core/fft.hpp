#pragma once

#include <complex>
#include <vector>

namespace radar {

using cfloat = std::complex<float>;

bool isPowerOfTwo(int n);

// FFT radix-2 de Cooley-Tukey, iterative et en place.
// n doit etre une puissance de 2. inverse = true calcule la FFT inverse
// (normalisee par 1/n).
void fft(cfloat* data, int n, bool inverse = false);
void fft(std::vector<cfloat>& data, bool inverse = false);

// Fenetre de Hann : attenue les bords du signal avant la FFT pour limiter
// les lobes secondaires (une cible forte "baverait" sinon sur toutes les cases).
std::vector<float> hannWindow(int n);

} // namespace radar
