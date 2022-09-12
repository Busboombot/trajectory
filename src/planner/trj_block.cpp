

#include <functional>   // std::function
#include <cmath>       // rint,  abs, roundf
#include <algorithm>    // std::min
#include <strstream>
#include <sstream>

#include "trj_block.h"
#include "trj_joint.h"
#include "trj_segment.h"

trj_float_t plan_err_f(trj_float_t x, trj_float_t t, trj_float_t v_0, trj_float_t v_c, trj_float_t v_1, trj_float_t a) {

    trj_float_t t_a, t_d, x_a, x_d, t_c, x_c;

    t_a = fabs(v_0 - v_1) / a;  // Time to change from v0 to v1 at max acceleration
    x_a = (v_0 + v_1) / 2 * t_a;

    t_d = fabs(v_1 - v_c) / a;  // Time to change from v0 to v1 at max acceleration
    x_d = (v_c + v_1) / 2 * t_d;

    t_c = fmax(t - (t_a + t_d), 0);
    x_c = fmax(v_c, 0) * t_c;

    return x - (x_a + x_c + x_d);
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
trj_float_t binary_search(const std::function<trj_float_t(trj_float_t)> &f, trj_float_t v_min, trj_float_t v_guess, trj_float_t v_max) {

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

trj_float_t Block::consistantize() {
    // Recalculate t to make x and v values integers, and everything more
    // consistent This operation will maintain the original value for x,
    // but may change t

    x_c = (int)x - (int)(x_a + x_d);

    t_a = abs((v_c - v_0) / joint->a_max);
    t_d = abs((v_c - v_1) / joint->a_max);

    if (round(x_c) == 0) { // get rid of small negatives
        x_c = 0;
    }

    if (v_c != 0) {
        t_c = abs(x_c / v_c);
    } else {
        t_c = 0;
    }


    t = t_a + t_c + t_d;

    return x_c;

}

void Block::init() {




    setBv(0, 0);

    v_0 = fmin(v_0, joint->v_max);
    v_1 = fmin(v_1, joint->v_max);


    if (x == 0) {
        v_c = v_0 = v_1 = 0;
    } else if (x < 2. * joint->small_x) {

        v_c = (sqrt(4. * joint->a_max * x + 2. * pow(v_0, 2) + 2. * pow(v_1, 2.)) / 2.);
    } else {
        v_c = joint->v_max;
    }

    tie(x_a, t_a) = accel_xt(v_0, v_c);
    tie(x_d, t_d) = accel_xt(v_c, v_1);

    x_c = x - (x_a + x_d);

    if (v_c != 0) {
        t_c = x_c / v_c;
    } else {
        t_c = 0;
    }



    t = t_a + t_c + t_d;



}

void Block::setBv(trj_float_t v_0_, trj_float_t v_1_) {
    // Reduce v_0 and v_1 for small values of x when there is an
    // error in x after planning


    trj_float_t x_a_, t_a_, t_d_, x_d_;

    if (!isnan(v_0_)) {
        v_0 = v_0_;
    }
    if (!isnan(v_1_)) {
        v_1 = v_1_;
    }

    tie(x_a_, t_a_) = accel_xt(v_0, 0);
    x_d_ = x - x_a_;


    if (x_d_ < 0) {
        v_0 = int(min(v_0,  sqrt(2.0 *  static_cast< double >(joint->a_max) *  static_cast< double >(x))));

        v_1 = 0;
    } else if (x == 0) {
        v_0 = 0;
        v_1 = 0;
    } else {
        v_1 = int(min(sqrt(2.0 * joint->a_max * t_d_), v_1));
    }


}

tuple<trj_float_t, trj_float_t> Block::accel_xt(trj_float_t v_i, trj_float_t v_f) {
    // Distance and time required to accelerate from v0 to v1 at
    // acceleration a

    if (v_i == v_f) {
        return {0, 0};
    }


    trj_float_t t_ = fabs((v_f - v_i) / joint->a_max);  // Time to change from v0 to v1 at max acceleration
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
tuple<trj_float_t, trj_float_t> Block::accel_acd(trj_float_t v_0_, trj_float_t v_c_, trj_float_t v_1_) {
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


void Block::plan(trj_float_t t_) {

    t = t_;

    if (x == 0 || t == 0) {
        set_zero();
        t_c = t_;
        return;
    }

    // Find v_c with a binary search, then patch it up if the selection changes the segment time.

    auto err = [this](float v_c) {

        trj_float_t err = plan_err_f(x, t, v_0, v_c, v_1, joint->a_max);

        return err;
    };

    v_c = fmin(binary_search(err, 0, x / t, joint->v_max), joint->v_max);

    tie(x_a, t_a) = accel_xt(v_0, v_c);
    tie(x_d, t_d) = accel_xt(v_c, v_1);

    // consistantize will make the values consistent with each other,
    // and if anything is wrong, it will show up in a negative x_c

    x_c = consistantize();

    trj_float_t x_err = abs(round(area()) - x);

    if (round(x_c) < 0 or (x > 25 and x_err > 1)) {
        // We've probably got a  small move, and we've already tried
        // reducing boundary velocities, so get more aggressive. If that
        // doesn't work, allow the time to expand.

        trj_float_t new_t = max(t_a + t_d, t);

        if (v_1 > 0) {
            v_1 = 0;
            return plan(t);

        } else if (v_0 > 0) {
            v_0 = 0;
            return plan(t);

        } else if (abs(new_t - t) > 0.0005) {
            return plan(t);

        } else {
            throw std::runtime_error("Unsolvable profile, incorrect area: ");
        }

    }

    if (round(t * 1000) != round(t * 1000)) {
        // This block is still too long.
        if (v_1 != 0) {
            v_1 = 0;
            return plan(t);
        } else if (round(v_0) > 0) {
            v_0 = v_0 / 2;
            return plan(t);
        }
    }
}



void Block::plan_ramp(float t_) {
    // Calculate a ramp profile for a fixed time

    trj_float_t guess;

    t = t_;

    if (x == 0 || t == 0) {
        set_zero();
        t_c = t_;
    }

    trj_float_t sqrt_ = joint->a_max * (joint->a_max * pow(t, 2) + 2 * t * v_0 - 2 * x);

    if (sqrt_ < 0) {
        guess = joint->a_max * t + v_0 - sqrt(sqrt_);
    } else {
        guess = fmin(x / t, joint->v_max);
    }

    auto err = [this](float v_c) {
        return plan_ramp_err_f(x, v_0, v_c, v_1, joint->a_max);
    };


    v_1 = v_c = fmin(binary_search(err, 0, guess, joint->v_max), joint->v_max);

    tie(x_a, t_a) = accel_xt(v_0, v_c);

    t_d = 0;
    x_d = 0;

    consistantize();
}

string fs(trj_float_t x, trj_float_t v, trj_float_t t){

    std::stringstream ss;
    ss << green << (int)(round(x)) << " " <<  blue <<  (int)(round(v)) << " " <<
    yellow << std::setprecision(4) << t << creset;
    return ss.str();
}

ostream &operator<<(ostream &output, const Block &b) {

    output <<  "["
            <<  setw(30)  << fs(b.x_a, b.v_0,b.t_a) << "|"
            <<  setw(30)  << fs(b.x_c, b.v_c,b.t_c) << "|"
            <<  setw(30)  << fs(b.x_d, b.v_1,b.t_d) <<
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

void Block::limitBv() {

    trj_float_t  x_a_, t_a_, x_d_;

    tie(x_a_, t_a_) = accel_xt(v_0, v_c);

    x_d = x - x_a;

    if (x_d < 0){
        v_0 = int(fmin(v_0, sqrt(2. * joint->a_max * x)));
    } else if (x == 0) {
        v_0 = 0;
        v_1 = 0;
    } else {
        v_1 =  int(fmin(sqrt(2. * joint->a_max * x_d), v_1));
    }

}

void Block::limitBv(Block& next) {

    limitBv();

    if (next.x == 0){
        v_1 = 0;
    }
}





