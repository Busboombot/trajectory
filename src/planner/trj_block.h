#pragma once
#include <tuple>
#include <stdexcept>
#include <iostream>
#include "trj_util.h"
#include "trj_types.h"

using namespace std;
using std::ostream;

class Segment;
class Joint;
class Planner;


class Block {

public:

    Block(trj_float_t x, Joint *joint, Segment *segment) :
            x(fabs(static_cast<double>(x))), joint(joint), segment(segment) {
        this->d = sign(x);
    }

    Block(trj_float_t x, trj_float_t v_0, trj_float_t v_1, Joint *joint, Segment *segment) :
        Block(static_cast<double>(x), joint, segment){
        this->v_0 = v_0;
        this->v_1 = v_1;
    }

private:
    trj_float_t x;
    trj_float_t d=0;
    trj_float_t t=0;

    trj_float_t t_a=0;
    trj_float_t t_c=0;
    trj_float_t t_d=0;

    trj_float_t x_a=0;
    trj_float_t x_c=0;
    trj_float_t x_d=0;

    trj_float_t v_0=0;
    trj_float_t v_c=0;
    trj_float_t v_1=0;

    Joint *joint;
    Segment *segment;

    int recalcs=0;

    int step_period=0;

public:
    void init();

    void plan(trj_float_t t_);

    void plan_ramp(float t_);

    trj_float_t area();

    struct ACDBlockParams params();

    friend ostream &operator<<( ostream &output, const Block &s );

    trj_float_t getT() const;

    void setBv(trj_float_t v_0_, trj_float_t v_1_); // Clip the boundary values based on the distance

    void limitBv();

    void limitBv(Block& next);

    trj_float_t getV0() const;

    trj_float_t getV1() const;

    friend Planner;

private:

    trj_float_t consistantize();

    void set_zero();

    std::tuple<trj_float_t, trj_float_t> accel_xt(trj_float_t v_i, trj_float_t v_1); // Compute trapezoid for acceleration from v_i to v_1

    std::tuple<trj_float_t, trj_float_t> accel_acd(trj_float_t v_0_, trj_float_t v_c_, trj_float_t v_1_); // Compute both accel and decel trapezoids


};