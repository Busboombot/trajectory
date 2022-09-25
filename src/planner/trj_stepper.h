#pragma once
#include <utility>
#include <vector>
#include <array>
#include "trj_util.h"

using namespace std;


class Stepper {

protected:

    int8_t axis;
    bool enabled = false;
    Direction direction;
    int stepState = 0;

public:

    Stepper() : axis(0), enabled(false){};
    Stepper(int8_t axis) : axis(axis), enabled(false){ };
    virtual ~Stepper(){}

    virtual void writeStep(){ stepState = 1; }
    virtual void clearStep(){ stepState = 0; };
    virtual void enable(){enabled = true;};
    virtual void enable(Direction dir){ setDirection(dir);enable();}
    virtual void disable() { setDirection(STOP); enabled = false;}
    virtual void setDirection(Direction dir){  direction = dir; }
    virtual void setDirection(int dir){  this->setDirection(static_cast<Direction>(dir)); };

private:
    friend ostream &operator<<( ostream &output, const Stepper &s );

};

using StepperPtr = shared_ptr<Stepper>;


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
    const StepperPhase *phase; // Current phase.

    StepperPtr stepper;

public:
    StepperState(int period, int timebase);
    StepperState();

    void loadPhases(vector<StepperPhase> phases);
    void loadPhases(array<StepperPhase,3> phases);
    void setStepper(StepperPtr stepper_){
        stepper = std::move(stepper_);

    }
    void next_phase();

    int next(double dtime);

    inline int isDone() const { return done ? 1 : 0; };

};

