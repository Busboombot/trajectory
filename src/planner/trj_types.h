#pragma once
#include <vector>
#include <cmath>

// Teensy 4.1 has very fast 32 float, but 64 bit floats  less so:
// "The FPU performs 32 bit float and 64 bit double precision math in hardware.
// 32 bit float speed is approximately the same speed as integer math. 64 bit
// double precision runs at half the speed of 32 bit float."

using trj_float_t = double;

using velocity_t = trj_float_t;
using VelocityVector = std::vector<velocity_t>;

using MoveArray = std::vector<int32_t>;

enum class MoveType {
    none,
    relative,
    absolute,
    jog,
    home
};

constexpr trj_float_t BV_PRIOR = -1.;
constexpr trj_float_t BV_NEXT = -2.;
constexpr trj_float_t BV_V_MAX = -3.;
constexpr trj_float_t BV_NAN = NAN;
