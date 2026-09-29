cmake_minimum_required(VERSION 3.16)
project(radar_sim LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(BUILD_APP "Compiler l'interface graphique (raylib)" ON)
option(USE_SYSTEM_RAYLIB "Utiliser raylib deja installe au lieu de le telecharger" OFF)

# --- Coeur de calcul : aucune dependance, testable sans ecran ---------------
add_library(radar_core STATIC
    src/core/fft.cpp
    src/core/target.cpp
    src/core/simulator.cpp
    src/core/range_doppler.cpp
    src/core/cfar.cpp
    src/core/tracker.cpp
)
target_include_directories(radar_core PUBLIC src)
target_compile_options(radar_core PRIVATE -Wall -Wextra)

add_executable(test_core tests/test_core.cpp)
target_link_libraries(test_core PRIVATE radar_core)

# --- Interface graphique -----------------------------------------------------
if(BUILD_APP)
    if(USE_SYSTEM_RAYLIB)
        find_package(raylib 5.0 REQUIRED)
    else()
        include(FetchContent)
        set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(BUILD_GAMES OFF CACHE BOOL "" FORCE)
        FetchContent_Declare(raylib
            GIT_REPOSITORY https://github.com/raysan5/raylib.git
            GIT_TAG 5.5
            GIT_SHALLOW TRUE)
        FetchContent_MakeAvailable(raylib)
    endif()

    add_executable(radar_sim
        src/app/main.cpp
        src/app/ppi.cpp
        src/app/plots.cpp
    )
    target_link_libraries(radar_sim PRIVATE radar_core raylib)
    target_compile_options(radar_sim PRIVATE -Wall -Wextra)
endif()
