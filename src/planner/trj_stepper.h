#pragma once
#include <vector>
#include <array>
#include <trj_util.h>

using namespace std;

struct StepperPhase{
    int x;
    double vi;
    double vf;
};

class StepperState {
private:
    // """A step segment for simulating the step interval algorithm. """

    int steps_left = 0;
    int steps_stepped = 0;
    int direction = 0;


    double t= 0;
    double t_f = 0;
    double phase_t = 0;
    double delay = 0;
    double delay_inc = 0;
    double delay_counter= 0;

    double a;

    int period = 0;
    int timebase = 0;
    int periods_left= 0;

    bool done = false;

    int phase_n;
    int phases_left = 0;
    vector<StepperPhase> phases;
    const StepperPhase *phase; // Current pahse.

public:
    StepperState(int period, int timebase);
    StepperState();

    void loadPhases(vector<StepperPhase> phases);
    void loadPhases(array<StepperPhase,3> phases);
    void next_phase();

    int next();

    inline const int isDone() { return done ? 1 : 0; };

};


class SegmentStepper {


};

