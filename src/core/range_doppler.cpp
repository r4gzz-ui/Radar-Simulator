#include "core/range_doppler.hpp"

#include <algorithm>
#include <complex>

namespace radar {

RangeDopplerProcessor::RangeDopplerProcessor(const Params& params) : p_(params) {}

void RangeDopplerProcessor::process(const Frame& frame, RangeDopplerMap& out) {
    const int ns = frame.ns;
    const int nc = frame.nc;
    if (static_cast<int>(winRange.size()) != ns)
        winRange = hannWindow(ns);
    if (static_cast<int>(winDoppler.size()) != nc)
        winDoppler = hannWindow(nc);

    cube.resize(static_cast<size_t>(ns) * nc);
    column.resize(nc);
    out.nRange = ns;
    out.nDoppler = nc;
    out.power.assign(static_cast<size_t>(ns) * nc, 0.0f);
    out.rangeProfile.assign(ns, 0.0f);

    // Normalisation : le bruit (sigma = 1) donne ~1 en moyenne dans les deux sorties.
    const float normRange = 1.0f / (ns * 0.375f * 2.0f);
    const float normRd = 1.0f / (static_cast<float>(ns) * nc * 0.375f * 0.375f * 2.0f);

    // 1. FFT distance sur chaque chirp.
    for (int m = 0; m < nc; ++m) {
        const cfloat* src = frame.chirp(m);
        cfloat* dst = cube.data() + static_cast<size_t>(m) * ns;
        for (int n = 0; n < ns; ++n)
            dst[n] = src[n] * winRange[n];
        fft(dst, ns);
        for (int r = 0; r < ns; ++r)
            out.rangeProfile[r] += std::norm(dst[r]) * normRange / nc;
    }

    // 2. FFT Doppler sur chaque case distance, puis recentrage (fftshift).
    for (int r = 0; r < ns; ++r) {
        for (int m = 0; m < nc; ++m)
            column[m] = cube[static_cast<size_t>(m) * ns + r] * winDoppler[m];
        fft(column.data(), nc);
        for (int k = 0; k < nc; ++k) {
            const int d = (k + nc / 2) % nc;
            const float pw = std::norm(column[k]) * normRd;
            out.power[static_cast<size_t>(d) * ns + r] = pw;
        }
    }
    // 3. Niveau de bruit : la mediane est insensible aux quelques cases
    //    occupees par des cibles. Pour un bruit exponentiel, moyenne = mediane / ln 2.
    median = out.power;
    std::nth_element(median.begin(), median.begin() + median.size() / 2, median.end());
    out.noiseLevel = std::max(median[median.size() / 2] / 0.693147f, 1e-12f);
}

} // namespace radar
