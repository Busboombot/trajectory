#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <ctime>

#include <catch2/catch_test_macros.hpp>

#include "trj_segment.h"
#include "trj_joint.h"
#include "trj_types.h"

using namespace std;

TEST_CASE("Basic Segment and Block test", "[block]")
{
    Joint j(0, 5e3, 50e3);

    Segment s(0,{j, j}, {5000, 1000});

    std::cout << s << endl<< endl;

    s.plan(0, BV_NAN, BV_NAN);

    std::cout << s << endl<< endl;

    s.setBv(0, 0);

    std::cout << s << endl<< endl;

    s.setBv(j.v_max, j.v_max);

    std::cout << s << endl<< endl;

}




