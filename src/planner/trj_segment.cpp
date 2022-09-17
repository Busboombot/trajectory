#include <cmath> // abs

#include "trj_segment.h"
#include "trj_block.h"

ostream &operator<<(ostream &output, const Segment &s) {

    output << setw(6) << setprecision(4) << s.t << " ";

    for (const Block &b: s.blocks) {
        output << b;
    }

    return output;
}

Segment::Segment(uint32_t n, const std::vector<Joint>&  joints, MoveArray move ) :
        n(n), t(0), joints(joints), moves(std::move(move)) {

    int axis_n = 0;
    for (const Joint &joint: this->joints) {
        blocks.emplace_back(this->moves[axis_n], joint, this);
        blocks.back().init();
        axis_n++;
    }

    n_joints = joints.size();

}

Segment::Segment(uint32_t n, const std::vector<Joint>&  joints, const Move& move ) : Segment(n, joints, move.x ){

}




void Segment::setBv(vector<double> v_0_, vector<double> v_1_) {

    for (int i = 0; i < blocks.size(); i++) {
        blocks[i].setBv(v_0_[i], v_1_[i], nullptr, nullptr);
    }
}

void Segment::setBv(double v_0_, double v_1_) {
    for (Block &b: blocks) {
        b.setBv(v_0_, v_1_);
    }
}


void Segment::plan(trj_float_t t_, int v_0_, int v_1_, Segment *prior, Segment *next) {

    Block *prior_block = nullptr;
    Block *next_block = nullptr;

    trj_float_t largest_at = 0, mt=0;

    for (const Joint &j: this->joints) {
        largest_at = fmax(largest_at, j.max_at);
    }

    for(int p_iter=0; p_iter<15; p_iter++){
        if (t != 0){
            mt = t_;
        } else if (p_iter <2) {
            mt = fmax(largest_at*2, min_time());
        } else {
            mt = fmax(largest_at*2, time());
        }

        for(int i=0; i<blocks.size(); i++) {
            Block &b = blocks[i];
            if (prior != nullptr) prior_block = &prior->blocks[i];
            if (next != nullptr) next_block = &next->blocks[i];

            b.plan(mt, v_0_, v_1_, prior_block, next_block);
        }
    }
}

trj_float_t  Segment::min_time(){

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
