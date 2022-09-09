#pragma once

#include <vector>
#include <limits.h>
#include <cstdint> 
#include <iostream>
#include <iomanip>

/// A Move vector, which describes distances for all axes
/**
 * The Move Vector describes the distances to move for all
 * axes, plus the maximum velocity for the whole vector. 
 * Note: The max velocity parameter is not currently used. 
*/

using std::array;
using std::cout;
using std::endl;
using std::setw;
using std::left;
using std::right;
using std::ostream;

using MoveArray = std::vector<int32_t> ;

struct Move {

    enum class MoveType {
        none, 
        relative, 
        absolute, 
        jog,
        home
    };

    uint32_t seq = 0; 

    MoveType move_type = MoveType::relative;

    // Total Vector Time, in microseconds.
    uint32_t t = 0; 

    // Distances
    MoveArray x;

    Move(int n_joints):seq(0), move_type(MoveType::relative), t(0), x(){
        x.resize(n_joints);
    }

    Move(int n_joints, uint32_t seq, uint32_t t, int v): seq(seq), t(t), x(){
        x.resize(n_joints);
    }

    Move(uint32_t seq, uint32_t t, MoveType move_type, MoveArray x): seq(seq), move_type(move_type), t(t), x(x){}

    Move(uint32_t seq, uint32_t t, MoveType move_type, std::initializer_list<int> il): 
        seq(seq), move_type(move_type), t(t), x(MoveArray(il.begin(), il.end())){}
    
    Move(uint32_t t, std::initializer_list<int> il): Move(0, t, MoveType::relative, il){}
    

private: 
    friend std::ostream &operator<<( ostream &output, const Move &p );
};

