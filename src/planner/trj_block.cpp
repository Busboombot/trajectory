

#include <functional>   // std::function
#include <cmath>       // rint,  abs, roundf
#include <algorithm>    // std::min
#include "trj_block.h"

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
double binary_search(const std::function<double(double)>& f, double v_min, double v_guess, double v_max){

    double old_guess;

    for(int i=0; i < 20; i++){

        double x = f(v_guess);
       
        if (roundf(x) > 0){
            old_guess = v_guess;
            v_guess = (v_max + v_guess) / 2.;
            v_min = old_guess;

        } else if (roundf(x) < 0){
            old_guess = v_guess;
            v_guess = (v_min + v_guess) / 2.;
            v_max = old_guess;

        } else {
            return v_guess;
        }

        if (fabs(v_min-v_max) < 1)
            return v_guess;
    }

    return NAN;
}


void ACDBlock::init(){
    float a_max = joint.a_max;
    float v_max = joint.v_max;

    set_bv();

    v_0 = fmin(v_0, v_max);
    v_1 = fmin(v_1, v_max);

    if (x == 0) {
        v_c = v_0 = v_1 = 0;
    } else if (x < 2 * joint.small_x) {
        v_c = (sqrt(4. * a_max * x + 2. * pow(v_0,2) + 2. * pow(v_1,2.)) / 2.);

    } else {
        v_c = v_max;

        tie(x_a, t_a) = accel_xt(v_0, v_c);
        tie(x_d, t_d) = accel_xt(v_c, v_1);
    }

    x_c = x - (x_a + x_d);
    if (v_c != 0){
        t_c = x_c / v_c;
    } else {
        t_c = 0;
    }

    t = t_a + t_c + t_d;

}

void ACDBlock::set_bv(){
    // Reduce v_0 and v_1 for small values of x when there is an
    // error in x after planning

    float a_max = this->joint.a_max;

    tie(x_a, t_a) = accel_xt(v_0, 0);
    x_d = x - x_a;

    if (x_d < 0){
        v_0 = int(sqrt(2 * a_max * x));
        v_1 = 0;
    } else if (x == 0) {
        v_0 = 0;
        v_1 = 0;
    } else {
        v_1 = int(min(sqrt(2 * a_max * x_d), v_1));
    }
}

tuple<double, double> ACDBlock::accel_xt(double v_0, double v_1){
    // Distance and time required to accelerate from v0 to v1 at
    // acceleration a
    double a = this->joint.a_max;

    if (v_0 == v_1){
        return {0, 0};
    }

    if (v_1 < v_0){
        a = -a;
    }

    t = (v_0 - v_1) / a;  // Time to change from v0 to v1 at max acceleration
    x = (v_0 + v_1) / 2 * t;

    return {x,t};

}
/**
 * COmpute the x and t for both the accel and decel phases
 * @param v_0_ Initial velocity
 * @param v_c_ Cruise velocity
 * @param v_1_ Final velocity
 * @return Tuple of x and t values
 */
tuple<double, double> ACDBlock::accel_acd(double v_0_, double v_c_, double v_1_){
    // Distance and time required to accelerate from v0 to v1 at
    // acceleration a
    double a = this->joint.a_max;

    double t1,t2,x1,x2;

    tie(x1, t1) = accel_xt(v_0_, v_c_);
    tie(x2, t2) = accel_xt(v_c_, v_1_);

    return {x1+x2,t1+t2};

}
