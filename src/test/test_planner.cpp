#include <iostream>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "trj_segment.h"
#include "trj_joint.h"
#include "trj_planner.h"

vector<Move> get2Moves(){
    return std::vector<Move>{
            Move( 0, {10000,100}),
            Move( 0, {10000,100})
    };
}

vector<Joint> get2Joints(){
    return std::vector<Joint>{
            Joint(0, 5e3, 50e3),
            Joint(1, 5e3, 50e3)
    };
}

TEST_CASE("Basic Planner Test", "[planner]"){

    vector<Joint> joints = get2Joints();

    Planner p(joints);

    p.move({1000,1000});
    p.move({1000,1000});
    p.move({1000,1000});


    cout <<" ============ " << endl;
    cout << p << endl;

}




