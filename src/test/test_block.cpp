#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <ctime>

#include <catch2/catch_test_macros.hpp>

#include "trj_segment.h"
#include "trj_joint.h"


TEST_CASE("Basic Segment and Block test", "[planner]")
{
    Joint j(0, 5e3, 50e3);

    Segment s(0,{j, j}, {5000, 1000});

    std::cout << s << std::endl;

    s.plan(NAN, NAN);

    std::cout << s << std::endl;

    s.setBv(0, 0);

    std::cout << s << std::endl;

    s.setBv(j.v_max, j.v_max);

    std::cout << s << std::endl;

}




