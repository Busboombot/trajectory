

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

Planner::Planner()  {

}
Planner::Planner(std::vector<Joint> joints_) {

    setJoints(joints_);
}

void Planner::setJoints(std::vector<Joint> joints_){

    joints.erase (joints.begin(),joints.end());

    plannerPosition.resize(joints.size());
    plannerPosition = {0};

    int i = 0;
    for (Joint &j: joints_) {
        j.n = i++;
        joints.push_back(j);
    }


}

void Planner::move(const Move &move) {
    this->move(move.x);
}

void Planner::move(const MoveArray &move) {

    // Add the move into the planner position
    auto mi = move.begin();
    auto ppi = plannerPosition.begin();
    for (; mi != move.end() && ppi != plannerPosition.end(); mi++, ppi++){
        *ppi += *mi;
    }

    segments.emplace_back(seg_num, joints, move);
    seg_num++;

    auto last_idx = segments.size() -1;
    Segment *pre_prior = segments.size() >= 3 ? &segments[last_idx-2] : nullptr;
    Segment *prior =     segments.size() >= 2 ? &segments[last_idx-1] : nullptr;
    Segment *current   = &segments[last_idx];


    if(pre_prior != nullptr){
        prior->plan(NAN, BV_NAN, BV_V_MAX, pre_prior);
        current->plan(NAN, BV_PRIOR, BV_NAN, prior);
        plan();
    } else if (prior != nullptr){
        prior->plan(NAN, BV_NAN, BV_V_MAX);
        current->plan(NAN, BV_PRIOR, BV_NAN, prior);
        plan();
    } else {
        current->plan(NAN, 0, 0);
    }

    queue_size = segments.size();
}

trj_float_t vLimit(int p_iter, trj_float_t v_max){
    if (p_iter < 2){
        return v_max;
    } else if (p_iter < 4){
        return v_max /2 ;
    } else {
        return 0;
    }
}

void Planner::plan(){
    // Plan the newest segment, and the one prior. If there are boundary
    // velocity discontinuities, also plan earlier segments
    Segment *current,  *prior, *pre_prior;
    trj_float_t  diff, mean_bv;
    int bends = 0;

    u_long seg_idx = segments.size() - 1;

    for(int p_iter = 0; p_iter < 15; p_iter++){
        current = &segments[seg_idx];
        prior   = &segments[seg_idx - 1];
        pre_prior = seg_idx >= 2 ? &segments[seg_idx - 2] : nullptr;

        prior->plan(NAN, BV_NAN, BV_NEXT, pre_prior, current);
        current->plan(NAN, BV_PRIOR, BV_NAN, prior);

        bends = 0;
        for (int i = 0; i < joints.size(); i++){
            Block &cb = current->blocks[i];
            Block &pb = prior->blocks[i];

            if (Block::bent(pb, cb)){
                mean_bv = Block::meanBv(pb, cb);
                diff = fabs(pb.v_1 - mean_bv);
                if ( diff < vLimit(p_iter, pb.joint.v_max)){
                    pb.v_1 = cb.v_0 = mean_bv;
                    bends++;
                }
            }
        }

        if ( bends > 0 || (pre_prior != nullptr && Segment::boundaryError(*pre_prior, *prior))){
            seg_idx += -1; // Re run on one earlier
        } else if (Segment::boundaryError(*prior, *current)) {
            seg_idx += 0; // re-run at this boundary
        } else {
            seg_idx += 1; // Advance to the next segment
        }


        seg_idx = max(1UL, seg_idx);

        if (seg_idx >= segments.size()){
            break;
        }

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




bool Planner::isEmpty() { return segments.size() == 0; }


ostream &operator<<(ostream &output, const Planner &p) {


    output << blue_bg << "════  Joints ════" << creset << endl;

    cout << "N Joints:  " << p.joints.size() << endl;

    for (const Joint j: p.joints) {
        output << j ;
    }


    output << endl << blue_bg << "════ Segments ════" << creset << endl;

    for (const Segment& s: p.segments) {
        output << s << endl;
    }

    return output;
}

json Planner::dump(const std::string& tag) const{

    json j;

    if(tag.size() > 0){
        j["_tag"] = tag;
    }

    j["_type"] = "Planner";
    for (const Joint joint: joints) {
        j["joints"].push_back(joint.dump());
    }

    for (const Segment& s: segments) {
        j["segments"].push_back(s.dump());
    }


    return j;
}


