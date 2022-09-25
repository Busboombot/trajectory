#include "trj_segstepper.h"
#include "trj_planner.h"

SegmentStepper::SegmentStepper(Planner &planner) : planner(planner) {
    reloadJoints();
}

void SegmentStepper::reloadJoints() {
    stepperStates.erase(stepperStates.begin(), stepperStates.end());

    for (const Joint &j: planner.getJoints()) {
        stepperStates.emplace_back();
    }
}


int SegmentStepper::next(double dtime) {

    time += dtime;
    totalPeriods += 1;

    if (activeAxes == 0 && !planner.segments.empty()) {

        Segment &seg = planner.segments.front();

        auto ssi = stepperStates.begin();
        for (const Block &b: seg.blocks) {
            (*ssi++).loadPhases(b.getStepperPhases());
        }
    }

    activeAxes = 0;
    for (StepperState &s: stepperStates) {
        activeAxes += s.next(dtime);
    }

    if (activeAxes == 0 && !planner.segments.empty()) {
        planner.segments.pop_front();
    }

    return (int) activeAxes;

}

void SegmentStepper::setSteppers(vector<StepperPtr> steppers_ ){

    this->steppers = steppers_;

    auto ssi = stepperStates.begin();

    for (StepperPtr &si: steppers) {
        si->setDirection(0);
        (*ssi++).setStepper(si);
    }
}
