#pragma once

#include "raylib.h"

// Palette "phosphore vert" des ecrans radar.
namespace app::theme {

inline constexpr Color background{4, 10, 7, 255};
inline constexpr Color panel{8, 18, 12, 255};
inline constexpr Color border{25, 70, 40, 255};
inline constexpr Color grid{20, 60, 32, 255};
inline constexpr Color gridText{50, 150, 80, 255};
inline constexpr Color text{150, 230, 170, 255};
inline constexpr Color dimText{80, 140, 100, 255};
inline constexpr Color phosphor{70, 255, 120, 255};
inline constexpr Color tx{70, 255, 120, 255};
inline constexpr Color rx{255, 150, 60, 255};
inline constexpr Color otherRx{220, 200, 90, 255};
inline constexpr Color track{255, 196, 60, 255};
inline constexpr Color tentative{150, 125, 60, 255};
inline constexpr Color truth{190, 200, 195, 255};
inline constexpr Color selected{255, 110, 60, 255};
inline constexpr Color detection{255, 255, 255, 255};

} // namespace app::theme
