#pragma once
#include "trj_planner.h"
#include "trj_types.h"

class Stepper;

using StepperPtr = shared_ptr<Stepper>;

class SegmentStepper {

public:

    explicit SegmentStepper(Planner &planner);

    int next(double dtime);

    void setSteppers( vector<StepperPtr> steppers);

    void reloadJoints();

    unsigned long getTotalPeriods() const { return totalPeriods; }

    unsigned long getActiveAxes() const { return activeAxes;}

    double getTime() const { return time; }

private:

    Planner & planner;
    vector<StepperState> stepperStates;
    vector<StepperPtr> steppers;

    unsigned long totalPeriods = 0;

    unsigned long activeAxes = 0;
    double time;

};

