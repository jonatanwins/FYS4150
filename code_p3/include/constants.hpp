#pragma once

// -------------------------------------- Constants ------------------------------------
// Coulomb constant 
constexpr double k_e = 1.38935333e5; // (μm)^3 / (μs)^2 e^2

// Magnetic field strength (Tesla) and electric potential (Volt)
constexpr double T_unit = 9.64852558e1; // μm / (μs * e)
constexpr double V_unit = 9.64852558e7; // μm^2 / (μs^2 * e)

// Penning trap configuration
constexpr double B0 = 1.00; // Tesla
constexpr double B0_converted = B0 * T_unit; // μm / (μs * e)

constexpr double V0 = 25.0e-3; // 25.0 mV = 25.0e-3 V
constexpr double V0_converted = V0 * V_unit; // μm^2 / (μs^2 * e)

constexpr double d_const = 500.0; // μm
