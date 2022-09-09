
#pragma once

#include <cstdint> 
#include <deque>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <math.h> // rint
#include <initializer_list> 

#include "trj_planner_const.h" // For N_AXES
#include "trj_jointss.h"
#include "trj_util.h"
#include "trj_move.h"
#include "trj_segment.h"

using std::array;
using std::cout;
using std::endl;
using std::setw;
using std::left;
using std::right;
using std::ostream;


class Segment;
class Joint;
class PhaseJoints;


/// Trajectory Planner. Turns a sequence of moves into moves plus velocities
/*
*/

class Planner {

protected:

    // Joint configuration
    std::vector<Joint> joints;

    // Running segment number
    int seg_number;

    std::deque<Segment*> segments;

    int current_phase = 0;

    // Current phase that is being stepped out?
    PhaseJoints phase_joints;

    int32_t queue_size=0;
    int32_t queue_time=0;

    MoveArray current_position;

public:

    Planner(std::vector<Joint> joints);
    
    Planner();

    // Reset the joints
    void setNJoints(int n_joints){
        joints.resize(n_joints); 
    }
    void setJoint(Joint &joint);

    const std::vector<Joint> &getJoints(){ return joints;}

    const Joint &getJoint(int i){ return joints[i];}

    const std::deque<Segment*> & getSegments() { return segments; }

    // Add a move, processing it into a Segment
    void push(const Move& move);

    void push(int seq, int t, MoveArray x);

    Segment& peekSegment(int i){
        return *segments[i];
    }

    int getSegmentsSize(){  return segments.size();}

    uint32_t getQueueTime(){ return queue_time;  }

    uint32_t getQueueSize(){ return  queue_size; }

    MoveArray getPosition(){ return current_position; }

    // Return ref to an internal PhaseJoints, loaded with the parameters for
    // the current phase of a give axis. 
    const PhaseJoints&  getNextPhase();

    const PhaseJoints& getCurrentPhase();

    void clear(){
        while (!isEmpty()){
            getNextPhase();
        }
    }

    bool isEmpty(){ return segments.size() == 0; }

    Move::MoveType getCurrentMoveType(){ 
        if (!isEmpty()) {
            return segments.front()->move_type;
        } else {
            return Move::MoveType::none;
        }
    }

    friend ostream &operator<<( ostream &output, const Planner &p );
    
    

};
