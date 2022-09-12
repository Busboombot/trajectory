

#include <cmath> // abs

#include <iostream>
#include <vector>
#include "trj_planner.h"
#include "trj_segment.h"
#include "trj_joint.h"
#include "trj_types.h"
#include "trj_util.h"

using namespace std;

int here_count = 0; // for the HERE macro, in trj_util


Planner::Planner(std::vector<Joint> joints_) {

    joints.resize(joints_.size());

    current_position.resize(joints.size());
    current_position = {0};

    int i = 0;
    for (Joint &j: joints_) {
        j.n = i++;
        setJoint(j);
    }
}


void Planner::move(const Move &move) {
    this->move(move.x);
}

void Planner::move(const MoveArray &move) {
    segments.emplace_back(joints, move);

    plan();

}

void Planner::plan(){
    // Plan the newest segment, and the one prior. If there are boundary
    // velocity discontinuities, also plan earlier segments


    u_long i = segments.size() - 1;
    Segment& current = segments[i];

    cout << "Plan " << i << endl;

    if (i == 0) {
        current.plan();
        return;
    }


    Segment& prior = segments[i - 1];
    prior.plan_ramp();  // Uncap the v_1 for the prior tail


    return ;


    for(int n = 0; n < 20; n++){
        // For random moves, only about 10% of the segments will have more than 2
        // iterations, and 1% more than 10.

        if (i >= segments.size() ){
            break;
        }

        current = segments[i];
        prior = segments[i - 1];

        i += plan_at_boundary(prior, current);

        if (i == 0){
            // Hit the beginning, so need to go forward
            i += 2;
        }
    }

}
int Planner::plan_at_boundary(Segment &a, Segment &b) {
    // Plan two segments, the current and immediately prior, to
    // get a consistent set of boundary velocities between them

    VelocityVector a_v0 = a.getV0();

    for(int i = 0; i < joints.size(); i++){
        a.blocks[i].limitBv(b.blocks[i]);
    }

    a.plan();
    b.setBv(a.getV1(), V_NAN);
    b.plan() ;

    if (b.getV0() != a.getV1()) {
        // This means that the current could not handle the commanded v_0,
        // so prior will have to yield.
        a.setBv(V_NAN, b.getV0());
        a.plan();
    }

    if (a_v0 != a.getV0() || boundary_error(a, b) > 1) {
        // Planning segment a changed it's v_0, so we need to o back one earlier
        return -1;
    } else {
        return 1;
    }

}

/**
 * RMS difference between velocities of blocks in two segments, at their boundary
 * @param p
 * @param c
 * @return
 */
double Planner::boundary_error(Segment &p, Segment &c){

    double sq_err = 0;

    for(int i = 0; i < joints.size(); i++ ){
        sq_err =pow( (p.blocks[i].v_1 - c.blocks[i].v_0) ,2);
    }

    return sqrt(sq_err);
}


void Planner::setNJoints(int n_joints) {
    joints.resize(n_joints);
}


void Planner::setJoint(Joint &joint) {
    if ((int) joint.n < (int) joints.size()) {
        joints[joint.n] = joint;
    }
}


bool Planner::isEmpty() { return segments.size() == 0; }


ostream &operator<<(ostream &output, const Planner &p) {


    cout << blue_bg << "════  Joints ════" << creset << endl;
    for (const Joint j: p.joints) {
        output << j ;
    }

    cout << endl << blue_bg << "════ Segments ════" << creset << endl;

    for (const Segment& s: p.segments) {
        output << s;
    }




    return output;
}


