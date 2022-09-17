

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
    segments.emplace_back(seg_num, joints, move);


    plan();

    seg_num++;
}

void Planner::plan(){
    // Plan the newest segment, and the one prior. If there are boundary
    // velocity discontinuities, also plan earlier segments

    u_long i = segments.size() - 1;
    Segment* current = &segments[i];

    if (i == 0) {
        //current->plan();
        return;
    }

    Segment* prior = &segments[i - 1];

    for(int p_iter = 0; p_iter < 20; p_iter++){
        // For random moves, only about 10% of the segments will have more than 2
        // iterations, and 1% more than 10.

        cout << "Plan " << i << " " << p_iter << " seg=" <<  current->getN() << endl;

        if (i >= segments.size() ){
            break;
        }

        current = &segments[i];
        prior = &segments[i - 1];

        i += plan_at_boundary(*prior, *current);

        if (i == 0){
            // Hit the beginning, so need to go forward
            i += 2;
        }
    }

}
int Planner::plan_at_boundary(Segment &prior, Segment &current) {
    // Plan two segments, the current and immediately prior, to
    // get prior consistent set of boundary velocities between them

    VelocityVector a_v0 = prior.getV0();


    cout <<endl<< "Plan@B " << prior.getN() << " " << current.getN() << endl;

    cout << "A  P " << prior << endl;
    cout << "A  C " << current << endl;

    //prior.plan();

    //current.setBv(prior.getV1(), V_NAN);
    //current.plan() ;


    cout << "B  P " << prior << endl;
    cout << "B  C " << current << endl;


    if (current.getV0() != prior.getV1()) {
        // This means that the current could not handle the commanded v_0,
        // so prior will have to yield.
        //prior.setBv(V_NAN, current.getV0());
        //prior.plan();
    }


    cout << "C  P " << prior << endl;
    cout << "C  C " << current << endl;



    if (a_v0 != prior.getV0() || boundary_error(prior, current) > 1) {
        // Planning segment prior changed it's v_0, so we need to o back one earlier
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
        output << s << endl;
    }




    return output;
}


