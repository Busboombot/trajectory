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
class SegmentStepper;

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
    const vector<Joint> &getJoints() const;

    Segment(uint32_t n, const std::vector<Joint>&  joints_);
    Segment(uint32_t n, const std::vector<Joint>&  joints_, MoveArray moves );
    Segment(uint32_t n, const std::vector<Joint>&  joints_, const Move& move );

    void plan(trj_float_t t_=NAN, int v_0_=0, int v_1_=0, Segment *prior = nullptr, Segment *next = nullptr);

    void setBv(int v_0, int v_1);

    void setBv(vector<int> v_0_, vector<int> v_1_);

public:

    MoveType getMoveType() const;

    const MoveArray& getMoves() const {return moves;}

    VelocityVector getV0();
    VelocityVector getV1();

    int getN() const { return n; }

    trj_float_t minTime();
    trj_float_t time();

    trj_float_t timeErr(); // RMS difference in times of blocks

    friend Planner;
    friend SegmentStepper;

    friend ostream &operator<<( ostream &output, const Segment &s );

    json dump(std::string tag="", bool dump_joints=false) const;

    static trj_float_t boundaryError(const Segment& prior, const Segment &next);

};