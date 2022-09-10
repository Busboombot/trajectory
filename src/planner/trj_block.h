#pragma once
#include <tuple>
#include "trj_joint.h"
#include "trj_segment.h"
#include "trj_util.h"

class ACDBlock {


public:
    ACDBlock(double x, const Joint &joint, const Segment &segment) : x(x), joint(joint), segment(segment) {
        this->d = sign(x);
    }

    ACDBlock(double x, double v_0, double v_1, const Joint &joint, const Segment &segment) : ACDBlock(x, joint, segment){
        this->v_0 = v_0;
        this->v_1 = v_1;
    }

private:
    double x;
    double t=0;

    double t_a=0;
    double t_c=0;
    double t_d=0;

    double x_a=0;
    double x_c=0;
    double x_d=0;

    double v_0=0;
    double v_c=0;
    double v_1=0;
    int d=0;

    const Joint &joint;
    const Segment &segment;

    int recalcs=0;

    int step_period=0;

private:
    void init();

    void set_bv(); // Clip the boundary values based on the distance

    tuple<double, double> accel_xt(double v_0, double v_1); // Compute trapezoid for acceleration from v_0 to v_1

    tuple<double, double> accel_acd(double v_0_, double v_c_, double v_1_); // Compute both accel and decel trapezoids
};