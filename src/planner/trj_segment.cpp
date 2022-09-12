#include <cmath> // abs

#include "trj_segment.h"
#include "trj_block.h"

ostream &operator<<(ostream &output, const Segment &s) {

    output << setw(6) << setprecision(4) << s.t << " ";

    for (const Block &b: s.blocks) {
        output << b;
    }

    output << endl;
    return output;
}

Segment::Segment(const std::vector<Joint>&  joints, MoveArray move ) :
        n(0), t(0), joints(joints), moves(std::move(move)) {

    int axis_n = 0;
    for (Joint &joint: this->joints) {
        blocks.emplace_back(this->moves[axis_n], &joint, this);
        blocks.back().init();
        axis_n++;
    }

    n_joints = joints.size();

}

Segment::Segment(const std::vector<Joint>&  joints, const Move& move ) : Segment(joints, move.x ){

}


void Segment::setBv(vector<double> v_0_, vector<double> v_1_) {

    for (int i = 0; i < blocks.size(); i++) {
        blocks[i].setBv(v_0_[i], v_1_[i]);
    }
}

void Segment::setBv(double v_0_, double v_1_) {
    for (Block &b: blocks) {
        b.setBv(v_0_, v_1_);
    }
}

void Segment::plan() {

    for (Block &b: this->blocks) {
        t = fmax(t, b.getT());
    }

    for (Block &b: this->blocks) {
        b.plan(t);
    }

}

void Segment::plan(double v_0_, double v_1_) {
    setBv(v_0_, v_1_);
    plan();
}

void Segment::plan(vector<double> v_0_, vector<double> v_1_) {
    setBv(std::move(v_0_), std::move(v_1_));
    plan();
}

void Segment::plan_ramp() {

    for (Block &b: this->blocks) {
        t = fmax(t, b.getT());
    }

    for (Block &b: this->blocks) {
        b.plan_ramp(t);
    }

}

MoveType Segment::getMoveType() const {
    return moveType;
}

VelocityVector Segment::getV0() {

    VelocityVector vv(n_joints);

    for (Block &b: this->blocks) {
        vv.push_back(b.getV0());
    }

    return vv;
}

VelocityVector Segment::getV1() {

    VelocityVector vv(n_joints);

    for (Block &b: this->blocks) {
        vv.push_back(b.getV1());
    }

    return vv;
}
