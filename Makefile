BUILD      := build
BUILD_TEST := build-test
JOBS       ?= $(shell nproc 2>/dev/null || echo 4)

# raylib est telecharge et compile automatiquement au premier "make"
# (pas besoin de sudo). Pour utiliser un raylib deja installe :
#   make CMAKE_FLAGS=-DUSE_SYSTEM_RAYLIB=ON
CMAKE_FLAGS ?=

.PHONY: all run test clean re

all:
	cmake -S . -B $(BUILD) -DCMAKE_BUILD_TYPE=Release $(CMAKE_FLAGS)
	cmake --build $(BUILD) -j $(JOBS)

run: all
	./$(BUILD)/radar_sim

# Tests du traitement du signal et du Kalman, sans raylib ni ecran.
test:
	cmake -S . -B $(BUILD_TEST) -DCMAKE_BUILD_TYPE=Release -DBUILD_APP=OFF
	cmake --build $(BUILD_TEST) -j $(JOBS) --target test_core
	./$(BUILD_TEST)/test_core

clean:
	rm -rf $(BUILD) $(BUILD_TEST)

re: clean all
