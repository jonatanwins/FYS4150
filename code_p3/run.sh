#!/bin/bash

# Detect the operating system
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    INCLUDE_PATH="-I/opt/homebrew/Cellar/armadillo/14.0.2_1/include"
    LIBRARY_PATH="-L/opt/homebrew/Cellar/armadillo/14.0.2_1/lib"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Ubuntu (or other Linux)
    INCLUDE_PATH=""
    LIBRARY_PATH=""
else
    echo "Unsupported operating system"
    exit 1
fi

g++ src/simulate_penning.cpp src/utils.cpp src/Particle.cpp src/PenningTrap.cpp \
-I include -o out -std=c++17 \
$INCLUDE_PATH $LIBRARY_PATH \
-larmadillo && ./out