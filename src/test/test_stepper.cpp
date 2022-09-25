#include <iostream>
#include <memory>
#include <vector>
#include <numeric>
#include <fstream>


#include "trj_segment.h"
#include "trj_joint.h"
#include "trj_planner.h"
#include "trj_stepper.h"
#include "catch2/catch_all.hpp"


#include "boost/filesystem.hpp"   // includes all needed Boost.Filesystem declarations
#include <iostream>               // for std::cout

using namespace boost::filesystem;

extern vector<Move> get2Moves();
extern vector<Joint> get2Joints();
extern vector<Joint> get4Joints();
extern std::vector<int> extractIntegerWords(const string& str);
using Ints = vector<int>;

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

#define csdc(p) ( std::dynamic_pointer_cast<CoutStepper>(p))


TEST_CASE("Basic Stepper Test", "[stepper]") {

    double dtime = 5./1e6; // 5 us

    vector<Joint> joints = get2Joints();
    Planner p(joints);

    p.move({-1000, 5000});
    p.move({-500, 10000});
    p.move({1000, -15000});

    cout << " ============ " << endl;
    cout << p << endl;

    SegmentStepper ss(p);
    array<int,2> acc;

    vector<StepperPtr> steppers;
    steppers.push_back(std::make_shared<CoutStepper>(0 ));
    steppers.push_back(std::make_shared<CoutStepper>(1));

    ss.setSteppers(steppers);

    do {
        ss.next(dtime);
    } while (!p.empty());


    // Check that it doesn't crash after segments are exhausted.
    ss.next(dtime);
    ss.next(dtime);
    ss.next(dtime);

    cout << "Final: 1:"<< csdc(steppers[0])->count<<" 2: "<< csdc(steppers[1])->count <<endl;

    REQUIRE(csdc(steppers[0])->count == -501);
    REQUIRE(csdc(steppers[1])->count == 0);

    cout << endl;

}

/* Read moves from a file */
TEST_CASE("Stepper File Test", "[stepper]") {

    double dtime = 5. / 1e6; // 5 us

    vector<Joint> joints{
            Joint(0, 5e3, 50e3),
            Joint(1, 5e3, 50e3),
            Joint(2, 5e3, 50e3)
    };
    Planner p(joints);
    fstream inputFile;
    array<int, 3> counts = {0};
    path inputFilePath = current_path().parent_path().parent_path() / "test" / "stepper_file_test.txt";
    inputFile.open( inputFilePath.string(), ios::in);

    if (inputFile.is_open()) {
        cout << "Loading  "<< inputFilePath << endl;
        for (std::string line; std::getline(inputFile, line);) {
            if (line[0] == ' ' || line[0] == '#') continue;
            Ints ints = extractIntegerWords(line); // Get all integers on a line
            p.move({ints[0], ints[1], ints[2]});
            for(int i=0; i < 3; i++) counts[i]+=ints[i];
        }
        inputFile.close(); //close the file object

        cout << "Loaded " << p.getQueueSize() << " moves,  Counts: "; for(int i=0; i < 3; i++) cout <<counts[i]<<" "; cout << endl;
        //cout << p << endl;

    } else {
        cout << "Err: not opened:  "<< inputFilePath << endl;
    }

    //
    // Run the steppers
    //

    SegmentStepper ss(p);

    vector<StepperPtr> steppers;
    steppers.push_back(std::make_shared<CoutStepper>(0 ));
    steppers.push_back(std::make_shared<CoutStepper>(1));
    steppers.push_back(std::make_shared<CoutStepper>(2));

    ss.setSteppers(steppers);

    auto start = chrono::steady_clock::now();

    int n_iter = 0;
    do {
        ss.next(dtime);
        if (++n_iter % 5'000'000 == 0){
            auto end  = chrono::steady_clock::now();
            auto diff = chrono::duration_cast<chrono::microseconds>(end - start);
            cout << "Periods: "<< ss.getTotalPeriods()<<" ("<<double(diff.count()/double(ss.getTotalPeriods()))<<")" <<
            "us/p Time: "<<ss.getTime() << " sec "<< endl;
        }
    } while (!p.empty());


    // Check that it doesn't crash after segments are exhausted.
    ss.next(dtime);
    ss.next(dtime);
    ss.next(dtime);

    cout << "Final: 1:"<< csdc(steppers[0])->count<<
            " 2: "<< csdc(steppers[1])->count <<
            " 3: "<< csdc(steppers[2])->count << endl;

    cout << "Total Periods: "<< ss.getTotalPeriods()<<" Time: "<<ss.getTime() << " sec "<< endl;

    cout << endl;


}

