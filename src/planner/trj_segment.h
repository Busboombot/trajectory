#pragma once

#include <utility>
#include <vector>
#include <iostream>

#include "trj_util.h"
#include "trj_block.h"
#include "trj_move.h"
#include "trj_joint.h"
#include "trj_types.h"


using namespace std;

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

    void plan();

    void plan(double v_i, double v_f);

    void plan(VelocityVector v_0_, VelocityVector v_1_);

    void plan_ramp();

    void setBv(double v_0, double v_1);

    void setBv(VelocityVector v_0_, VelocityVector v_1_);

public:

    MoveType getMoveType() const;

    const MoveArray& getMoves() const {return moves;}

    VelocityVector getV0();
    VelocityVector getV1();

    int getN() const { return n; }



    friend Planner;

    friend ostream &operator<<( ostream &output, const Segment &s );
};