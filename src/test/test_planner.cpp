#include <iostream>
#include <vector>

#include "trj_segment.h"
#include "trj_joint.h"
#include "trj_planner.h"


#include <catch2/catch_test_macros.hpp>


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

vector<Joint> get4Joints(){
    return std::vector<Joint>{
            Joint(0, 5e3, 50e3),
            Joint(1, 5e3, 50e3),
            Joint(2, 5e3, 50e3),
            Joint(3, 5e3, 50e3)
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

TEST_CASE("Large Small Planner Test", "[planner]"){

    vector<Joint> joints = get2Joints();

    Planner p(joints);

    p.move({1000,1});
    p.move({1,1000});

    cout <<" ============ " << endl;
    cout << p << endl;

}




