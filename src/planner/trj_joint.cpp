#include <cstdint>
#include <deque>
#include <vector>
#include <array>
#include <iostream>
#include "trj_joint.h"


std::ostream &operator<<(std::ostream &output, const Joint &j) {

    output << "[J " << j.n << " v=" << j.v_max << " a=" << j.a_max << " d=" << j.d_max << " ]";
    return output;

}