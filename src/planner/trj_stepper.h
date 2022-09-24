#pragma once
#include <utility>
#include <vector>
#include <array>
#include "trj_util.h"
#include "trj_planner.h"

using namespace std;

class StepInterface {
public:
    StepInterface()= default;
    virtual int step(double t, int step) = 0;
    virtual void setDirection(int direction) = 0;
};

struct StepperPhase{
    int x;
    double vi;
    double vf;
};

class StepperState {
private:

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

    StepInterface *stepper=nullptr;

public:
    StepperState(int period, int timebase);
    StepperState();

    void loadPhases(vector<StepperPhase> phases);
    void loadPhases(array<StepperPhase,3> phases);
    void setStepper(StepInterface* stepper_){ stepper = stepper_;}
    void next_phase();
    int next();

    inline int isDone() const { return done ? 1 : 0; };

};


class SegmentStepper {

public:

    explicit SegmentStepper(Planner &planner);
    int next();


private:

    Planner & planner;
    vector<StepperState> stepperStates;
    vector<StepInterface*> steppers;
public:
    void setSteppers(const vector<StepInterface *> &steppers);

private:

    unsigned long activeAxes = 0;

};

