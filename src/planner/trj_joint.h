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

class Joint {

public:

    Joint(int n, float v_max, float a_max): n(n), v_max(v_max), a_max(a_max){
        small_x = pow(v_max,2) / (2 * a_max);
    }

    Joint(): Joint(0,0,0){}

    explicit Joint(int n): Joint(n,0,0){}

    int n;
    trj_float_t v_max;
    trj_float_t a_max;
    trj_float_t small_x;

    friend std::__1::ostream &operator<<(std::__1::ostream &output, const Joint &j );

};