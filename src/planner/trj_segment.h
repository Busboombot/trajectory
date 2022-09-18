#pragma once

#include <utility>
#include <vector>
#include <iostream>

#include "trj_util.h"
#include "trj_block.h"
#include "trj_move.h"
#include "trj_joint.h"
#include "trj_types.h"
#include "json.hpp"


using namespace std;
using json = nlohmann::json;

class Planner;

/* Segment: One move for all joints, with Accel, Cruise and Decel phases. 
 * 
 *
 */

class Segment {

private:

    uint32_t n;
    trj_float_t t;
    MoveType moveType = MoveType::none;

    vector<Block> blocks;
    const vector<Joint>& joints;
    MoveArray moves;

    u_long n_joints;

public:

    Segment(uint32_t n, const std::vector<Joint>&  joints, MoveArray moveS );
    Segment(uint32_t n, const std::vector<Joint>&  joints, const Move& move );

    void plan(trj_float_t t_, int v_0_, int v_1_, Segment *prior = nullptr, Segment *next = nullptr);

    void setBv(int v_0, int v_1);

    void setBv(vector<int> v_0_, vector<int> v_1_);

public:

    MoveType getMoveType() const;

    const MoveArray& getMoves() const {return moves;}

    VelocityVector getV0();
    VelocityVector getV1();

    int getN() const { return n; }

    trj_float_t min_time();
    trj_float_t time();

    friend Planner;

    friend ostream &operator<<( ostream &output, const Segment &s );

    json dump() const;


};