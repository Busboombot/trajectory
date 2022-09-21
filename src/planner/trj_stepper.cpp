#include <tuple>
#include <utility>
#include "trj_stepper.h"

StepperState::StepperState(int period, int timebase) : period(period), timebase(timebase) {

    delay_inc = (float)period / (float)timebase;
}

StepperState::StepperState() : StepperState(4,1000000) { }


void StepperState::loadPhases(vector<StepperPhase> phases_) {

    phases = std::move(phases_);
    phases_left = phases.size();
    phase_n = 0;
    done = false;
}

void StepperState::loadPhases(array<StepperPhase,3> phases_) {

    phases = vector<StepperPhase>(phases_.begin(), phases_.end());
    phases_left = phases.size();
    phase_n = 0;
    done = false;
}


void StepperState::next_phase() {

    phase  = &phases[phase_n];

    direction = sign(phase->x);
    steps_left = abs(phase->x);

    t_f = (phase->vi + phase->vf) != 0 ? fabsf((2.f * (float)steps_left) / (phase->vi + phase->vf)) : 0;
    a =  t_f != 0 ? (phase->vf - phase->vi) / t_f : 0;

    phase_t = 0;

    float v = a * delay_inc + phase->vi;
    delay = v!=0 ? fabs(1/v) : 0;
    delay_counter += delay_inc;

    periods_left = int(round(t_f / delay_inc));
    done = false;

    phase_n += 1;
    phases_left -= 1;
}


int StepperState::next() {

    if (steps_left <= 0 || periods_left <= 0){
        if (done or phases_left==0){
            done = true;
            return 0;
        } else {
            next_phase();
        }
    }

    int r = 0;
    if (delay_counter > delay) {
        delay_counter -= delay;
        steps_left -= 1;
        steps_stepped += 1;
        r = direction;
    }

    periods_left -= 1;

    float v = phase->vi + a * phase_t;

    delay = v!=0 ? abs(1 / v) : 1;
    delay_counter += delay_inc;

    t += delay_inc;
    phase_t += delay_inc;

    float calc_x = fabsf((a * (float)pow(phase_t,2)) / 2.f + phase->vi * phase_t);
    float x_err = (float)steps_stepped - calc_x;

    if (abs(x_err) > .5 && phase_n == 1) {
        float s = (x_err / abs(x_err)); // Just want the sign
        delay_counter += -s * delay_inc * .1f; //  Only Make the adjustment slowly, 10% at a step
    }

    return r;
}
