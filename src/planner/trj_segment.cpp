#include <cmath> // abs

#include "trj_segment.h"
#include "trj_block.h"

using json = nlohmann::json;

ostream &operator<<(ostream &output, const Segment &s) {

    output << setw(6) << setprecision(4) << s.t << " ";

    for (const Block &b: s.blocks) {
        output << b;
    }

    return output;
}

// This one is  just for testing.
Segment::Segment(uint32_t n, const std::vector<Joint>&  joints_) : n(n), joints(joints_), t(0){

}

Segment::Segment(uint32_t n, const std::vector<Joint>&  joints_, MoveArray move ) :
        n(n), joints(joints_), t(0), moves(std::move(move)) {

    int axis_n = 0;
    for (const Joint &joint: joints) {
        blocks.emplace_back(this->moves[axis_n], joint, this);

        axis_n++;
    }
    n_joints = joints.size();

}

Segment::Segment(uint32_t n, const std::vector<Joint>&  joints_, const Move& move ) :
    Segment(n, joints_, move.x ){

}


void Segment::setBv(vector<int> v_0_, vector<int> v_1_) {

    for (int i = 0; i < blocks.size(); i++) {
        blocks[i].setBv(v_0_[i], v_1_[i], nullptr, nullptr);
    }
}

void Segment::setBv(int v_0_, int v_1_) {
    for (Block &b: blocks) {
        b.setBv(v_0_, v_1_);
    }
}

trj_float_t Segment::boundaryError(const Segment& prior, const Segment &next) {

    trj_float_t sumsq = 0;

    for( int i = 0; i < prior.joints.size(); i++ ){
        sumsq += pow(prior.blocks[i].getV1() - next.blocks[i].getV0(),2);
    }

    return sqrt(sumsq);

}

void Segment::plan(trj_float_t t_, int v_0_, int v_1_, Segment *prior, Segment *next) {

    Block *prior_block = nullptr;
    Block *next_block = nullptr;

    trj_float_t largest_at = 0, mt=0, lower_bound_time = 0;

    // Longest time requirest for any block to accelerate or decelerate.
    for (const Joint &j: this->joints) {
        largest_at = fmax(largest_at, j.max_at);
    }

    lower_bound_time = largest_at*2;

    for(int p_iter=0; p_iter<10; p_iter++){
        if ( !isnan(t_)){
            mt = t_;
        } else if (p_iter <2) {
            mt = minTime();
        } else if (p_iter <4) {
            mt = fmax(lower_bound_time, minTime());
        } else {
            mt = fmax(lower_bound_time, time());
        }

        for(int i=0; i<blocks.size(); i++) {
            Block &b = blocks[i];
            if (prior != nullptr) prior_block = &prior->blocks[i];
            if (next != nullptr) next_block = &next->blocks[i];

            b.plan(mt, v_0_, v_1_, prior_block, next_block);
        }

        if (timeErr() < .001) {
            break;
        } else{
            for(Block &b : blocks){
                if (b.getT() < mt){
                    b.limitBv();
                }
            }
        }

    }
}

trj_float_t  Segment::minTime(){

    trj_float_t mt=0;

    for(Block &b : blocks){
        mt = fmax(mt,b.getMinTime());
    }

    return mt;

}

trj_float_t  Segment::time(){
    trj_float_t mt=0;

    for(Block &b : blocks){
        mt = fmax(mt,b.getT());
    }

    return mt;
}

trj_float_t Segment::timeErr(){
    trj_float_t mean_time = 0, rms_time=0;

    for(Block &b : blocks){
        mean_time += b.getT();
    }

    mean_time /= blocks.size();

    for(Block &b : blocks){
        rms_time += pow(b.getT()-mean_time, 2);
    }

    return sqrt(rms_time);


}

MoveType Segment::getMoveType() const {
    return moveType;
}

VelocityVector Segment::getV0() {

    VelocityVector vv;

    for (Block &b: this->blocks) {
        vv.push_back(b.getV0());
    }

    return vv;
}

VelocityVector Segment::getV1() {

    VelocityVector vv;

    for (Block &b: this->blocks) {
        vv.push_back(b.getV1());
    }

    return vv;
}

json Segment::dump(std::string tag, bool dump_joints) const{

    vector<json> o;
    json j;

    if(tag.size() > 0){
        j["_tag"] = tag;
    }

    j["_type"]="Segment";

    j["move"] = moves;

    if(dump_joints){
        for(const Joint &joint: joints){
            j["joints"].push_back(joint.dump());
        }
    }

    for(const Block &b: blocks){
       j["blocks"].push_back(b.dump());
    }


    return j;
}

const vector<Joint> &Segment::getJoints() const {
    return joints;
}

