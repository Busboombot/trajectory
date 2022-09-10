#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <time.h>

#include <catch2/catch_test_macros.hpp>
#include "trj_planner.h"
#include "trj_segment.h"
#include "trj_joint.h"

vector<Move> get2Moves(){
    return std::vector<Move>{
        Move( 0, {10000,100}),
        Move( 0, {10000,100})
    };
}

vector<Joint> getJoints(){
    return std::vector<Joint>{
        Joint(0, 10e3, 300e3),
        Joint(1, 10e3, 300e3)
     
    };
}

TEST_CASE("Basic planner test", "[planner]")
{
    Planner p = Planner( {Joint(0, 10e3, 300e3)} );
    
    int x = 1000;
    
    p.push(Move(0,1e5, Move::MoveType::relative, {x}));
    REQUIRE(1000 == p.getPosition()[0]);
    
    cout << p <<  endl;

    p.push(Move(0,1e5, Move::MoveType::relative, {x}));
    REQUIRE(2000 == p.getPosition()[0]);
    
    cout << p <<  endl;

    p.push(Move(0,1e5, Move::MoveType::relative, {x}));
    REQUIRE(3000 == p.getPosition()[0]);
    
    cout << p <<  endl;
}




