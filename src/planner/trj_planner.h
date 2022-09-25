
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
#include "trj_util.h"
#include "trj_move.h"
#include "trj_segment.h"
#include "trj_joint.h"
#include "trj_types.h"
#include "trj_stepper.h"
#include "trj_segstepper.h"
#include "json.hpp"


using std::array;
using std::cout;
using std::endl;
using std::setw;
using std::left;
using std::right;
using std::ostream;
using json = nlohmann::json;

class Segment;
class Joint;
class SegmentStepper;

class Planner {

protected:

    std::vector<Joint> joints; // Joint configuration

    std::deque<Segment> segments;

    int32_t queue_size=0;
    int32_t queue_time=0;
    int32_t seg_num = 0;

    MoveArray plannerPosition;
    MoveArray completedPosition;

public:

    explicit Planner();
    explicit Planner(std::vector<Joint> joints);

    void setJoints(std::vector<Joint> joints);

    // Add a move, processing it into a Segment
    void move(const Move& move);
    void move(const MoveArray& move);

    void plan();

    bool isEmpty();

    //SegmentStepper& getSegmentStepper();

public:

    // Fpr passing in to set_bv for boundaries you don't want to change.
    VelocityVector V_NAN = {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN};

    unsigned long getNSegments(){  return segments.size();}
    bool empty(){  return segments.empty();}

    uint32_t getQueueTime() const{ return queue_time;  }

    uint32_t getQueueSize() const{ return  queue_size; }

    MoveArray getPosition(){ return plannerPosition; }

    MoveType getCurrentMoveType(){
        if (!isEmpty()) {
            return segments.front().getMoveType();
        } else {
            return MoveType::none;
        }
    }

    const std::vector<Joint> &getJoints(){ return joints;}

    const Joint &getJoint(int i){ return joints[i];}

    json dump(const std::string& tag="") const;

    friend ostream &operator<<( ostream &output, const Planner &p );

    friend SegmentStepper;

private:


private:

    double boundary_error(Segment &p, Segment &c);



};
