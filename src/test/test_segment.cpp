#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <time.h>

#include <catch2/catch_test_macros.hpp>
#include "trj_segment.h"
#include "trj_joint.h"

using json = nlohmann::json;
using namespace std;

void dump_seg(const Segment& seg, string test_name, string tag){
    json jout;
    jout["test"] = test_name;
    jout["output"] = seg.dump(tag, true);
    cout << "JSON" << jout << endl;
}

TEST_CASE("Basic Segment Test", "[segment]")
{
    Joint j(0, 5e3, 50e3);
    std::vector<Joint> joints = {j, j, j};

    {
        Segment s(0, joints, {1000, 400, 240});
        s.plan();

        std::cout << s << endl;

        dump_seg(s, "basic_segment_1", "A");
    }

    {
        Segment s(0, joints, {1000, 1, 499});
        s.plan();

        std::cout << s << endl;

        dump_seg(s, "basic_segment_2", "B");
    }



}