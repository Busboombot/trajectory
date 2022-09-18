#include <cstdint>
#include <deque>
#include <vector>
#include <array>
#include <iostream>
#include "trj_joint.h"


std::ostream &operator<<(std::ostream &output, const Joint &j) {

    output << "[J " << j.n << " v=" << j.v_max << " a=" << j.a_max  << " ]";
    return output;

}

json Joint::dump() const{

    json j;
    j["_type"] = "Joint";
    j["n"] = n;
    j["v_max"] = v_max;
    j["a_max"] = a_max;

    return j;

}