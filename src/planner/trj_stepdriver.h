#pragma once
#include <Arduino.h>
#include <limits.h>

#include "trj_jointss.h"
#include "trj_planner_const.h" // For N_AXES
#include "trj_config.h"
#include "trj_planner.h"
#include "trj_util.h"
#include "trj_debug.h"

typedef enum
{
    CCW = -1,  ///< Clockwise
    STOP = 0,  ///< Clockwise
    CW  = 1   ///< Counter-Clockwise

} Direction;


// 2,000 rpm for a 1.8deg stepper is 400,000 steps per min, 7K steps 
// per sec. For a 10 ustep driver, 70KHz step pulse. 


/**
 * @brief Track the current conditions for the queue and positions, as of the time
 * of ACKs and DONE messages. 
 * 
 */
struct CurrentState {
  int32_t queue_length = 0;
  uint32_t queue_time = 0;
  int32_t positions[N_AXES] = {0};
  int32_t planner_positions[N_AXES] = {0};
}; 


class Stepper {

protected:


    AxisConfig config;
    int8_t axis;
    bool enabled = false;
    Direction direction = Direction::STOP;


public:

    Stepper() : axis(0), enabled(false){};
    Stepper(int8_t axis) : axis(axis), enabled(false){ config.axis=axis;};
    Stepper(AxisConfig config) : config(config), axis(config.axis) {};
    virtual ~Stepper(){}

    virtual void writeStep(){ }
    virtual void clearStep(){};
    virtual void enable(){enabled = true;};
    virtual void enable(Direction dir){ setDirection(dir);enable();}
    virtual void disable() { setDirection(STOP); enabled = false;}
    void setDirection(Direction dir){direction = dir;};
    void setConfig(AxisConfig config){ this->config = config;}
    AxisConfig getConfig(){ return config; }

private: 
    friend ostream &operator<<( ostream &output, const Stepper &s );
  

};

class StepperState {
    
protected: 

    Direction direction=STOP; 
    uint32_t  stepsLeft=0;
    long position = 0;

    int period = 0; 

    float delay_counter=0;
    float delay=0;
    float delay_inc=0;
    float v=0;
    float a=0;
    float t=0;   // Running time
    float t_s=0; // segment time, in sections
    float v_i=0; // Initial velocity


public:

    StepperState() {}
    
    
    inline uint32_t getStepsLeft(){ return stepsLeft;}
    
    inline int getVelocity(){ return v;}
    
    inline Direction getDirection() { return direction; }

    inline long getPosition(){ return position; }

    inline void setPosition(long p){ position = p; }
    
    void setParams(uint32_t segment_time, uint32_t v0, uint32_t v1, int32_t x, int period);

    int step(Stepper *stepper);

 private:
        friend std::ostream & operator<<(std::ostream &os, const StepperState& sd);
   

};

class SubSegment {
public:
    uint16_t seq = 0;
    uint16_t code = 0;
    uint32_t segment_time = 0; // total segment time, in microseconds 
    uint32_t lastTime;
    StepperState axes[N_AXES];
public:
    SubSegment() {}
};


class StepDriver {

public:

    StepDriver() : period(10) {
        
    }
    ~StepDriver(){}

    void stop() { running = false; }

    void start() { running = true; }

    void enable();

    inline void enable(uint8_t axis){
        if (steppers[axis] != nullptr){
            steppers[axis]->enable(state[axis].getDirection());
        }
    }

    void disable();

    inline void disable(uint8_t axis){
        if (steppers[axis] != nullptr){
            steppers[axis]->disable();
        }
    }

    void setAxisConfig(uint8_t axis, unsigned int v_max, unsigned int a_max);

    void setStepper(uint8_t axis, Stepper* stepper){steppers[axis] = stepper;}

    inline StepperState& getState(uint8_t n){ return state[n]; }
    inline Planner& getPlanner(){ return planner;}

    inline int step(uint8_t n){ 
        return state[n].step(steppers[n]); 
    }

    inline void clear(uint8_t n){
        if (steppers[n] != 0){
            steppers[n]->clearStep(); 
        }
    }

    int update(); // Run all of the stepping and state updates


    int loadNextPhase();

    bool isEmpty(){ return planner.isEmpty() & !phaseIsActive; }

    void push(Move m);


    double sincePhaseStart();

    void clear(){planner.clear();}

    void zero(){
        for (StepperState ss : state){
            ss.setPosition(0);
        }
    }

    void setNAxes(int n){ n_axes = n; }
    void setPeriod(int p){period = p;}


    int checkIsDone(){
        if(segment_is_done >= 0){
            int seq = segment_is_done;
            segment_is_done  = -1 ;
            return seq;
        } else {
            return -1;
        }
    }

    bool checkIsEmpty(){
        if(is_empty){
            is_empty = false;
            return true;
        } else{
            return false;
        }
    }

protected:

    int n_axes = 0; // Gets increased in setAxisConfig

    bool running = false;
    bool enabled = false;
    bool phaseIsActive = false;

    uint32_t nextUpdate = 0; // microseconds since start of current phase for next update to steppers
    uint32_t nextClear = 0; // microsecond until next clear of step pins. 

    // flags
    int segment_is_done = -1; // If last segment completed, holds the sequence number
    bool is_empty = false; // Latching flag signalling emptyness. 

    uint32_t period; // Inter-interrupt time

    Planner planner;
    uint32_t now;
  
    StepperState state[N_AXES];

    Stepper *steppers[N_AXES] = {nullptr};

    int lastSeq; // Most recent sequence number, from the last Phase joints loaded. 

 private:
        friend std::ostream & operator<<(std::ostream &os, const StepDriver& sd);
      

};

