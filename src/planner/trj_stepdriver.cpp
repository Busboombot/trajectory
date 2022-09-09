#include <limits>
#include <chrono>

#include "trj_stepdriver.h"
#include "trj_util.h"
#include "trj_debug.h"
#include "trj_debug.h"

#include <sstream>

#ifdef TRJ_ENV_HOST
long micros(){
  return 0;
}
#endif

CurrentState current_state;
auto lastPhaseTime = steadyClock::now();



void StepDriver::setAxisConfig(uint8_t axis, unsigned int v_max, unsigned int a_max){

  if(axis>=n_axes){
    n_axes = axis+1;
    planner.setNJoints(n_axes);
  }

  Joint joint(axis,static_cast< float >(v_max), static_cast< float >(a_max));
  planner.setJoint(joint);
 
}

void StepDriver::enable(){
  enabled = true;

  switch (n_axes){
      case 6: enable(5);
      case 5: enable(4);
      case 4: enable(3);
      case 3: enable(2);
      case 2: enable(1);
      case 1: enable(0);
  }
}

void StepDriver::disable(){
  enabled = false;
  switch (n_axes){
      case 6: disable(5);
      case 5: disable(4);
      case 4: disable(3);
      case 3: disable(2);
      case 2: disable(1);
      case 1: disable(0);
  }
}

void StepDriver::push(Move move){

    
    for (int axis = 0; axis < n_axes; axis++){
      if(move.move_type == Move::MoveType::absolute){
        current_state.planner_positions[axis] = move.x[axis];
      } else {
        current_state.planner_positions[axis] += move.x[axis];
      }
    }
    
    current_state.queue_length = planner.getQueueSize();
    current_state.queue_time = planner.getQueueTime();

    planner.push(move);

}

int StepDriver::loadNextPhase(){
  const PhaseJoints&  pj = planner.getNextPhase();

  int active_axes = 0;

  for (int axis = 0; axis < n_axes; axis++){
    const JointSubSegment &jss = pj.moves[axis];
    if(jss.x != 0){
      active_axes++;
      
      state[axis].setParams(jss.t, jss.v_0, jss.v_1, jss.x, period);

      if (steppers[axis] != nullptr){
        steppers[axis]->enable(state[axis].getDirection());
      }
   
    } else {
      state[axis].setParams(0, 0, 0, 0, 0);
    }
  }
  
  lastSeq = pj.seq;
  return active_axes;
}



int StepDriver::update(){

  static unsigned long  stepsLeft = 0;

  uint32_t t = usince();
  
  if( !phaseIsActive  & !planner.isEmpty() ){ 
    int active_axes = loadNextPhase();
    
    if( active_axes > 0){
        phaseIsActive = true;
        start_usince();
        t = nextUpdate =  usince(); // SHould be 0, or there about
    }
  }

  if(phaseIsActive and (t >= nextUpdate) ){

    nextUpdate = t+period;
    stepsLeft = 0;
  
    switch (n_axes) {
      case 6: stepsLeft += step(5);
      case 5: stepsLeft += step(4);
      case 4: stepsLeft += step(3);
      case 3: stepsLeft += step(2);
      case 2: stepsLeft += step(1);
      case 1: stepsLeft += step(0);
      case 0: ;
    }


    if(stepsLeft == 0){ 

      segment_is_done = lastSeq;
      phaseIsActive = false;

      if (isEmpty()){
        is_empty = true;
      }
    } else {
     
      if(nextClear == 0){
        nextClear = t + period - 1;
      }
    }

    
  } 

  if ( (nextClear != 0 ) & (t >= nextClear)){
    switch (n_axes) {
      case 6: clear(5);
      case 5: clear(4);
      case 4: clear(3);
      case 3: clear(2);
      case 2: clear(1);
      case 1: clear(0);
      case 0: ;
    }
    nextClear = 0;
    
  }

 
  return stepsLeft;

}

double StepDriver::sincePhaseStart(){ 
  return 0;
}


void StepperState::setParams(uint32_t segment_time, uint32_t v0, uint32_t v1, int32_t x, int period){

        if (x > 0){
            direction = CW;
        } else if (x < 0){
            direction = CCW;
        } else {
            direction = STOP;
        }

        t = 0; // cumulative time
        delay = 0;
        delay_counter = 0;
    
        if(v1 + v0 != 0 and abs(x) != 0){
            stepsLeft = abs(x);
            // This is from the Dave Austin algorithm, via the AccelStepper arduino library. 
            t_s =  fabs(2.0 * ((float)x)) / ( ((float)v1) + ((float)v0 ) );
            v_i = (float)v0;
            a = ((float)v1-(float)v0) / t_s;
        } else {
            stepsLeft = 0;
            t_s = 0;
            v_i = 0;
            a = 0;
        }
        
        delay_inc = ((float)period) / ((float)TIMEBASE);
       
    }

  int StepperState::step(Stepper *stepper){


        if (stepsLeft == 0){
            return 0;
        }
        
        if (delay_counter >= delay){

            delay_counter -= delay;
            stepsLeft--;
            
            if (stepper != nullptr){
              
              stepper->writeStep();
            }
          
            position += direction;
        }

        // Note that t is always increasing, even if there isn't a step, so 
        // v is always changing. If the first step from a v0 == 0 is very large, delay can be
        // so large that the first step may not be schedule for several seconds, but if this function
        // is called regularly, t will increase, and v will increase, so delay will decrease. After
        // a couple of calls, delay will be reasonable again. 
        //
        // Alternatively, we could set a minimum velocity, but that is probably not necessary.

        v = a * t + v_i;

        if(v != 0){
            delay = 1.0/abs(v);

        } else {
            delay = 0;
        }

        delay_counter += delay_inc;
        t += delay_inc;

        
        return stepsLeft;
  }

ostream &operator<<( ostream &output, const Stepper &s )  { 
    output << "[Stp " << (int)s.axis  << " ]";
    return output;
}
        
ostream &operator<<( ostream &output, const StepperState &s ) { 
    output << "[SS sl="<<s.stepsLeft << " rt=" << s.t << " dl="<<s.delay 
           << " dc=" << s.delay_counter << " p=" << s.position << " ]" ;
    return output;

}

ostream &operator<<( ostream &output, const StepDriver &sd ) { 

    output << "--- Steppers --" <<endl;
    for(int i =0; i < sd.n_axes; i++){
      Stepper *s = sd.steppers[i];
      if(s == nullptr){
        output << "[]" << endl;
      } else {
        output << s << " " << *s << endl;
      }
    }
    output << endl;
    output << "--- State --" <<endl;
    for(int i =0; i < sd.n_axes; i++){
      output << sd.state[i] <<endl;
    }

  
    output << "--- Segments --" << endl;
    output << sd.planner << endl;

    return output;
}


