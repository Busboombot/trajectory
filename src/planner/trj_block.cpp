

#include <functional>   // std::function
#include <cmath>       // rint,  abs, roundf
#include <algorithm>    // std::min
#include <strstream>
#include <sstream>
#include <assert.h>

#include "trj_block.h"
#include "trj_joint.h"
#include "trj_segment.h"
#include "trj_stepper.h"

using json = nlohmann::json;

trj_float_t plan_err_f(trj_float_t x, trj_float_t t, trj_float_t v_0,
                       trj_float_t v_c, trj_float_t v_1, trj_float_t a) {

    trj_float_t t_a, t_d, x_a, x_d, t_c, x_c, t_ad, x_ad;

    t_ad = (fabs(v_c - v_0) + fabs(v_c - v_1)) / a;
    x_ad = fabs((pow(v_0, 2) - pow(v_c, 2)) / (2. * a)) +
           abs((pow(v_1, 2) - pow(v_c, 2)) / (2. * a));

    t_c = fmax(t - (t_ad), 0);
    x_c = fmax(v_c, 0) * t_c;

    return x - (x_ad + x_c);
}

trj_float_t plan_ramp_err_f(trj_float_t x, trj_float_t v_0, trj_float_t v_c, trj_float_t v_1, trj_float_t a) {

    trj_float_t t_a, t_d, x_a, x_d, t_c, x_c, t, x_err;

    t_a = fabs(v_0 - v_1) / a;  // Time to change from v0 to v1 at max acceleration
    x_a = (v_0 + v_1) / 2 * t_a;

    t_c = t - t_a;

    x_c = v_c * t_c;
    x_err = x - (x_a + x_c);

    return x_err;
}

/**
 * @brief Use a binary search to find the value of v_c that produces an value for
 * f(v_c) of close to 0. 
 * 
 * @param f  Error functon
 * @param v_min minimum v_c value
 * @param v_guess minital v_c value
 * @param v_max max v_c value. 
 * @return float 
 */
trj_float_t binary_search(const std::function<trj_float_t(trj_float_t)> &f, trj_float_t v_min, trj_float_t v_guess,
                          trj_float_t v_max) {

    trj_float_t old_guess;

    for (int i = 0; i < 20; i++) {

        trj_float_t x = f(v_guess);

        if (roundf(x) > 0) {
            old_guess = v_guess;
            v_guess = (v_max + v_guess) / 2.;
            v_min = old_guess;

        } else if (roundf(x) < 0) {
            old_guess = v_guess;
            v_guess = (v_min + v_guess) / 2.;
            v_max = old_guess;

        } else {
            return v_guess;
        }

        if (fabs(v_min - v_max) < 1)
            return v_guess;
    }

    return NAN;
}

void Block::set_zero() {

    x_a = x_d = x_c = 0;
    t_a = t_d = t_c = 0;
    v_0 = v_c = v_1 = 0;

}


tuple<trj_float_t, trj_float_t> Block::accel_xt(trj_float_t v_i, trj_float_t v_f) const {
    // Distance and time required to accelerate from v0 to v1 at
    // acceleration a

    if (v_i == v_f) {
        return {0, 0};
    }


    trj_float_t t_ = fabs((v_f - v_i) / joint.a_max);  // Time to change from v0 to v1 at max acceleration
    trj_float_t x_ = fabs((v_i + v_f) / 2 * t_);


    return {x_, t_};

}

/**
 * COmpute the x and t for both the accel and decel phases
 * @param v_0_ Initial velocity
 * @param v_c_ Cruise velocity
 * @param v_1_ Final velocity
 * @return Tuple of x and t values
 */
tuple<trj_float_t, trj_float_t> Block::accel_acd(trj_float_t v_0_, trj_float_t v_c_, trj_float_t v_1_) const {
    // Distance and time required to accelerate from v0 to v1 at
    // acceleration a


    trj_float_t t1, t2, x1, x2;

    tie(x1, t1) = accel_xt(v_0_, v_c_);
    tie(x2, t2) = accel_xt(v_c_, v_1_);

    return {x1 + x2, t1 + t2};

}

/**
 * Calculate the area (distance) of the block.
 * @return
 */
trj_float_t Block::area() {

    trj_float_t x_ad, t_ad, t_c_;

    tie(x_ad, t_ad) = accel_acd(v_0, v_c, v_1);

    if (t_ad > t && t - t_ad < 1e-7) { /// Avoid very small negatives
        t_c_ = 0;
    }

    x_c = v_c * t_c;

    if (round(x_c) < 0 or t_c_ < 0) {
        throw std::runtime_error("Negative x_c");
    }

    return x_ad + x_c;
}

trj_float_t Block::getMinTime() const {

    trj_float_t x_ad_, t_ad_, t_c_, v_c_;

    if (x == 0) {
        v_c_ = 0;

    } else if (x < 2.* joint.small_x) {
        v_c_ = (sqrt(4. * joint.a_max * x +
                     2. * pow(v_0, 2) +
                     2. * pow(v_1, 2.)
                     ) / 2);

    } else {
        // If there is more area than the triangular profile for these boundary
        // velocities, the v_c must be v_max. In this case, it must also be true that:
        //    x_ad, t_ad = accel_acd(self.v_0, v_max, self.v_1, a_max)
        //    assert self.x > x_ad
        v_c_ = joint.v_max;
    }


    tie(x_ad_, t_ad_) = accel_acd(v_0, v_c_, v_1);

    if (v_c_ != 0) {
        t_c_ = (x - x_ad_) / v_c_;
    } else {
        t_c_ = 0;
    }

    t_c_ = max(t_c_, t_ad_ / 2);  // Enforce 1/3 rule, each a,c,d is 1/3 of total time

    return t_c_ + t_ad_;
}

void Block::plan(trj_float_t t_, int v_0_, int v_1_, Block *prior, Block *next) {

    trj_float_t v_c_;

    if (isnan(t_)) {
        t = getMinTime();
    } else {
        t = t_;
    }

    setBv((int) v_0_, (int) v_1_, prior, next);


    if (x == 0 || t == 0) {
        set_zero();
        t_c = t_;
        return;
    }

    // Find v_c with a binary search, then patch it up if the selection changes the segment time.

    auto err = [this](float v_c_) {

        trj_float_t err = plan_err_f(x, t, v_0, v_c_, v_1, joint.a_max);
        return err;
    };

    v_c_ = fmin(binary_search(err, 0, x / t, joint.v_max), joint.v_max);


    v_c = fmin(v_c_, joint.v_max);

    tie(x_a, t_a) = accel_xt(v_0, v_c);
    tie(x_d, t_d) = accel_xt(v_c, v_1);

    // consistantize will make the values consistent with each other,
    // and if anything is wrong, it will show up in a negative x_c

    x_c =  x - (x_a + x_d);

    t_a = abs((v_c - v_0) / joint.a_max);
    t_d = abs((v_c - v_1) / joint.a_max);

    if (round(x_c) == 0 and x_c < 0) { // get rid of small negatives
        x_c = 0;
    }

    if (v_c != 0) {
        t_c = abs(x_c / v_c);
    } else {
        t_c = 0;
    }

    t = t_a + t_c + t_d;

}


void Block::setBv(int v_0_, int v_1_, Block *prior, Block *next) {

    trj_float_t x_a_, t_a_;

    assert(v_0_ != BV_NEXT);
    assert(v_1_ != BV_PRIOR);

    if (v_0_ == (int) BV_PRIOR) {
        assert(prior != nullptr);
        v_0 = prior->v_1;
    } else if (v_0_ != (int) BV_NAN) {
        v_0 = v_0_;
    }

    if (v_1_ == (int) BV_NEXT) {
        assert(next != nullptr);
        v_1 = next->v_0;
    } else if (v_1_ == (int) BV_V_MAX) {
        v_1 = this->joint.v_max;
    } else if (v_1_ != (int) BV_NAN) {
        v_1 = v_1_;
    }

    if (prior != nullptr) {
        // If the current block has a different sign -- changes direction --
        // then the boundary velocity must be zero.
        if (!same_sign(prior->d, d) || prior->x == 0 or x == 0) {
            v_0 = 0;
        }
    }

    tie(x_a_, t_a_) = accel_xt(v_0, 0);
    trj_float_t x_d_ = x - x_a_;

    if (x_d_ < 0) {
        v_0 = int(min(v_0, sqrt(2.0 * joint.a_max * static_cast< double >(x))));
        v_1 = 0;

    } else if (x == 0) {
        v_0 = 0;
        v_1 = 0;

    } else {
        v_1 = int(min(v_1, sqrt(2.0 * joint.a_max * x_d_)));

    }

    v_0 = min(v_0, joint.v_max);
    v_1 = min(v_1, joint.v_max);

}

void Block::limitBv() {
    if (v_1 > joint.v_max / 2) {
        v_1 /= 2;
        return;
    }

    if (v_0 > joint.v_max / 2) {
        v_0 /= 2;
        return;
    }

    if (v_1 > 1) {
        v_1 /= 2;
        return;
    }

    if (v_0 > 1) {
        v_0 /= 2;
        return;
    }
}

array<StepperPhase,3> Block::getStepperPhases() const{

    auto dd = double(d);
    return array<StepperPhase,3>{
            StepperPhase{int(d)*int(round(x_a)),dd*v_0,dd*v_c},
            StepperPhase{int(d)*int(round(x_c)),dd*v_c,dd*v_c},
            StepperPhase{int(d)*int(round(x_d)),dd*v_c,dd*v_1} };

}

bool Block::bent(Block& prior, Block &current) {
    int cd = (int)current.d;
    int pd = (int)prior.d;

    int s1 = sign(pd * prior.v_c - pd * prior.v_1);
    int s2 = sign(cd * current.v_0 - cd * current.v_c);

    return (s1 * s2) < 0;
}


trj_float_t Block::meanBv(Block& prior, Block &next) {

    trj_float_t a;
    trj_float_t  mv;

    if (prior.t_d + next.t_a != 0){
        a = (next.v_c - prior.v_c) / (prior.t_d + next.t_a);
        mv = prior.v_c + a*prior.t_d;
    } else  {
        mv =  ( next.v_c +  prior.v_c)/2.;
    }

    return mv;
}




string fs_a(trj_float_t x, trj_float_t v, trj_float_t t) {

    std::stringstream ss;
    ss <<
       blue << setw(5) << (int) (round(v)) << " " <<
       green << setw(5) << (int) (round(x)) << " " <<
       yellow << setw(5) << (int) (t * 1000) << creset;
    return ss.str();
}

string fs_d(trj_float_t x, trj_float_t v, trj_float_t t) {

    std::stringstream ss;
    ss <<
       green << setw(5) << (int) (round(x)) << " " <<
       yellow << setw(5) << (int) (t * 1000) << " " <<
       blue << setw(5) << (int) (round(v)) << " " << creset;
    return ss.str();
}

ostream &operator<<(ostream &output, const Block &b) {

    output << "["
           << setw(48) << fs_a(b.x_a, b.v_0, b.t_a) << "|"
           << setw(48) << fs_a(b.x_c, b.v_c, b.t_c) << "|"
           << setw(48) << fs_d(b.x_d, b.v_1, b.t_d) <<
           "] ";

    return output;
}

trj_float_t Block::getT() const {
    return t;
}


trj_float_t Block::getV0() const {
    return v_0;
}

trj_float_t Block::getV1() const {
    return v_1;
}

json Block::dump(std::string tag) const {

    json m;

    if(tag.size() > 0){
        m["_tag"] = tag;
    }

    m["_type"] = "Block";
    m["x"] = x;
    m["d"] = d;
    m["t"] = t;
    m["t_a"] = t_a;
    m["t_c"] = t_c;
    m["t_d"] = t_d;
    m["x_a"] = x_a;
    m["x_c"] = x_c;
    m["x_d"] = x_d;
    m["v_0"] = v_0;
    m["v_c"] = v_c;
    m["v_1"] = v_1;

    return m;

}



