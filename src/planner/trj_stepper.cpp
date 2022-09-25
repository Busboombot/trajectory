#include <tuple>
#include <utility>
#include "trj_stepper.h"

StepperState::StepperState(int period, int timebase) : period(period), timebase(timebase) {

    delay_inc = (double) period / (double) timebase;
}

StepperState::StepperState() : StepperState(4, 1000000) {}


void StepperState::loadPhases(vector<StepperPhase> phases_) {

    phases = std::move(phases_);
    phases_left = phases.size();
    phase_n = 0;
    done = false;
}

void StepperState::loadPhases(array<StepperPhase, 3> phases_) {

    loadPhases(vector<StepperPhase>(phases_.begin(), phases_.end()));

}


void StepperState::next_phase() {
    phase = &phases[phase_n];

    direction = sign(phase->x);
    steps_left = abs(phase->x);

    t_f = (phase->vi + phase->vf) != 0 ? fabs((2.f * (double) steps_left) / (phase->vi + phase->vf)) : 0;
    a = t_f != 0 ? (phase->vf - phase->vi) / t_f : 0;

    phase_t = 0;

    double v = a * delay_inc + phase->vi;
    delay = v != 0 ? fabs(1 / v) : 0;
    delay_counter += delay_inc;

    periods_left = int(round(t_f / delay_inc));
    done = false;

    phase_n += 1;
    phases_left -= 1;

    if (stepper != nullptr) {
        stepper->setDirection(direction);
    }
}

int StepperState::next(double dtime) {

    if (steps_left <= 0) { //} || periods_left <= 0) {
        if (done or phases_left == 0) {
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

        if (stepper != nullptr) {
            stepper->writeStep();
        }
    } else  if (stepper != nullptr){
        stepper->clearStep();
    }

    periods_left -= 1;

    double v = phase->vi + a * phase_t;

    delay = v != 0 ? abs(1 / v) : 1;
    delay_counter += dtime;

    t += dtime;
    phase_t += dtime;

    return 1;
}


ostream &operator<<(ostream &output, const Stepper &s) {
    output << "[Stp " << (int) s.axis << " ]";
    return output;
}

