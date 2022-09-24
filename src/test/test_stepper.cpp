#include <iostream>
#include <vector>
#include <numeric>

#include "trj_segment.h"
#include "trj_joint.h"
#include "trj_planner.h"
#include "trj_stepper.h"
#include "catch2/catch_all.hpp"

extern vector<Move> get2Moves();

extern vector<Joint> get2Joints();

class CoutStepper : public StepInterface {

public:

    CoutStepper(int axis, vector<int> &acc, double &time) : axis(axis), acc(acc), last_time(time) {}

    int step(double time, int step) override {
        last_time = time;
        acc[axis] = step;
        return step;
    }

    void setDirection(int direction_) override {
        direction = direction_;
    }

public:
    int direction;
    vector<int> &acc;
    int axis;
    double &last_time;

};

TEST_CASE("Basic Stepper Test", "[stepper]") {

    vector<Joint> joints = get2Joints();

    Planner p(joints);

    p.move({-1000, 5000});
    p.move({-500, 10000});
    p.move({1000, -15000});

    cout << " ============ " << endl;
    cout << p << endl;

    vector<int> acc(joints.size());
    vector<int> counts(joints.size());


    double last_time;
    int stepSum;
    SegmentStepper ss(p);

    vector<StepInterface *> steppers{
            new CoutStepper(0, acc, last_time),
            new CoutStepper(1, acc, last_time)
    };
    ss.setSteppers(steppers);

    //for(int i =0 ;i < 100'000; i++)
    do {

        ss.next();

        stepSum = accumulate(acc.begin(), acc.end(), 0);

        auto ai = acc.begin();
        for (int &c: counts) {
            c += *ai++;
        }

        if (stepSum > 0) {
            // cout << last_time << " ";
            for (auto s: acc) {
                //cout << s << " ";
            }
            //cout << endl;
        }
    } while (p.getNSegments() > 0);

    // Check that it doesn't crash after segments are exhausted.
    ss.next();
    ss.next();
    ss.next();

    REQUIRE(counts[0] == -500);
    REQUIRE(counts[1] == 0);
    REQUIRE( round(last_time*10)/10 == 6.2 );

    cout << last_time << " ";
    for (int &c: counts) cout << c << " ";
    cout << endl;


}
