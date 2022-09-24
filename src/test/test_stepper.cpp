#include <iostream>
#include <memory>
#include <vector>
#include <numeric>

#include "trj_segment.h"
#include "trj_joint.h"
#include "trj_planner.h"
#include "trj_stepper.h"
#include "catch2/catch_all.hpp"

extern vector<Move> get2Moves();

extern vector<Joint> get2Joints();

class CoutStepper : public Stepper {

public:

    CoutStepper(int axis) : Stepper(axis) { }

    ~CoutStepper() override {}

    void writeStep() override {
        Stepper::writeStep();
        count += direction;
        lastStep = 1;
    }

    void clearStep() override {
        Stepper::clearStep();
        lastStep = 0;
    }

    void setDirection(Direction direction_) override {
        Stepper::setDirection(direction_);
    }


public:
    int lastStep = 0;
    int count = 0;
};

TEST_CASE("Basic Stepper Test", "[stepper]") {

    vector<Joint> joints = get2Joints();

    Planner p(joints);

    p.move({-1000, 5000});
    p.move({-500, 10000});
    p.move({1000, -15000});

    cout << " ============ " << endl;
    cout << p << endl;


    SegmentStepper ss(p);

    vector<StepperPtr> steppers;
    steppers.push_back(std::make_shared<CoutStepper>(0 ));
    steppers.push_back(std::make_shared<CoutStepper>(1));

    ss.setSteppers(steppers);

    do {
        ss.next();
    } while (!p.empty());


    // Check that it doesn't crash after segments are exhausted.
    ss.next();
    ss.next();
    ss.next();

    auto stp = std::dynamic_pointer_cast<CoutStepper>(steppers[0]);
    REQUIRE(stp->count == -500);
    stp = std::dynamic_pointer_cast<CoutStepper>(steppers[1]);
    REQUIRE(stp->count == 0);

    cout << endl;

}
