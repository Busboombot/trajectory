#pragma once
#include <cstdint>
#include <deque>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <cmath> // rint
#include "trj_planner_const.h"
#include "trj_types.h"
#include "json.hpp"

using json = nlohmann::json;

class Joint {

public:

    Joint(int n, float v_max, float a_max): n(n), v_max(v_max), a_max(a_max){
        small_x = pow(v_max,2) / (2 * a_max);

        max_discontinuity = a_max / v_max;  // Max vel change in 1 step
        max_at = v_max / a_max; // Max time for acell or decel to cover full velocity range
    }

    Joint(): Joint(0,0,0){}

    explicit Joint(int n): Joint(n,0,0){}

    int n;
    trj_float_t v_max;
    trj_float_t a_max;
    trj_float_t small_x;
    trj_float_t max_discontinuity;
    trj_float_t max_at;

    friend std::__1::ostream &operator<<(std::__1::ostream &output, const Joint &j );

    json dump() const;

};

