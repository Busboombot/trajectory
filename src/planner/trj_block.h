#pragma once
#include <tuple>
#include <stdexcept>
#include <iostream>
#include "trj_util.h"
#include "trj_types.h"
#include <map>
#include "json.hpp"

#include "json.hpp"

using json = nlohmann::json;

using namespace std;
using std::ostream;

class Segment;
class Joint;
class Planner;


class Block {

public:

    Block(trj_float_t x, const Joint& joint, Segment *segment) :
            x(fabs(static_cast<double>(x))), joint(joint), segment(segment) {
        this->d = sign(x);
    }

    Block(trj_float_t x, trj_float_t v_0, trj_float_t v_1, const Joint& joint, Segment *segment) :
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

    const Joint &joint;
    Segment *segment;



public:
    void init();

    void plan(trj_float_t t_, int v_0_, int v_1_, Block *prior = nullptr, Block *next = nullptr);

    trj_float_t area();

    struct ACDBlockParams params();

    friend ostream &operator<<( ostream &output, const Block &s );

    trj_float_t getT() const;

    trj_float_t getMinTime() const;

    void setBv(int v_0_, int v_1_, Block *prior = nullptr, Block *next = nullptr) ; // Clip the boundary values based on the distance

    void limitBv();

    trj_float_t getV0() const;

    trj_float_t getV1() const;

    json dump() const;

    friend Planner;

private:


    void set_zero();

    std::tuple<trj_float_t, trj_float_t> accel_xt(trj_float_t v_i, trj_float_t v_1) const; // Compute trapezoid for acceleration from v_i to v_1

    std::tuple<trj_float_t, trj_float_t> accel_acd(trj_float_t v_0_, trj_float_t v_c_, trj_float_t v_1_) const; // Compute both accel and decel trapezoids


};