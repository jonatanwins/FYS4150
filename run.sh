#!/bin/bash
g++ code_p3/src/simulate_penning.cpp code_p3/src/utils.cpp code_p3/src/Particle.cpp code_p3/src/PenningTrap.cpp \
-I code_p3/include -o out -std=c++14 \
-I/opt/homebrew/Cellar/armadillo/14.0.2_1/include \
-L/opt/homebrew/Cellar/armadillo/14.0.2_1/lib \
-larmadillo && ./out
