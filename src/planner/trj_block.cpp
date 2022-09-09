

#include <functional>   // std::function
#include <math.h>       // rint,  abs, roundf
#include <algorithm>    // std::min

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
float binary_search(std::function<float(float)> f, float v_min, float v_guess, float v_max){

    float old_guess;

    for(int i=0; i < 20; i++){

        float x = f(v_guess);
       
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