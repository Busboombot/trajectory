#pragma once

#include <cstdint> 
#include <deque>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <cmath> // rint

#include "trj_planner_const.h"
#include "trj_joint.h"

#pragma once
using std::cout;
using std::endl;
using std::setw;
using std::left;
using std::right;
using std::ostream;

class Segment;

enum class SubSegName { NONE=4, ACCEL=0, CRUISE=1, DECEL=2 };

// Parameters for one phase of a segement for one joint. 
struct JointSubSegment {
    int32_t t=0; // Time in microseconds
    int32_t x=0; // Distance in steps. 
    int32_t v_0=0;
    int32_t v_1=0;
    SubSegName ssn = SubSegName::NONE;

    JointSubSegment(int32_t t, int32_t x, int32_t v_0, int32_t v_1, SubSegName ssn):
        t(t), x(x), v_0(v_0), v_1(v_1), ssn(ssn) {}

    JointSubSegment(): t(0), x(0), v_0(0), v_1(0), ssn(SubSegName::NONE) {}

    friend ostream &operator<<( ostream &output, const Joint &pj );
};


// All the move parameters for one phase, for all joints.

using AxisMoves = std::array<JointSubSegment, N_AXES>;

/**
 * @brief Hold the parameters for the steps that are currently being stepped out
 * 
 */
struct PhaseJoints {

    int seq = 0; // sequence number
    SubSegName ssn = SubSegName::NONE;

    int32_t t; // phase time in microseconds
    
    AxisMoves moves; // One JointSubSegment per axis

    friend ostream &operator<<( ostream &output, const PhaseJoints &pj );
};

// All three phases of one move for a join. 
struct JointSubSegment3{
    JointSubSegment a;
    JointSubSegment c;
    JointSubSegment d;

    friend ostream &operator<<( ostream &output, const JointSubSegment3 &j );
};

using SubSegments3 = std::vector<JointSubSegment3>;


class JointSegment {

public:

    JointSegment(int n, Joint &joint, Segment *segment, int x );

    int n;

    float x=0;
    int sign=0;
    float x_err=0;

    float v_max=0;

    float v_0=0;
    float v_0_max=0;

    float v_c;

    float v_1=0;
    float v_1_max=0;

    // Minimum executing times for various sections
    float t_min=0;
    float t_a_min=0;
    float t_c_min=0;
    float t_d_min=0;

    float x_a=0;
    float x_c=0;
    float x_d=0;

    Joint &joint;
    Segment *segment;
    JointSegment *next_js=nullptr;
    JointSegment *prior_js=nullptr;

    JointSegment *getNextJs();
    JointSegment *getPriorJs();

    void update();
    
    float update_sub_segments();
    void update_start_velocity_limit(bool is_first, bool sign_change);
    void update_end_velocity_limit(bool is_last);
    void update_boundary_velocity(JointSegment *prior_js_, JointSegment* next_js_);

    float calc_x(float v_c);
    float search_v_c();
    void update_t_min();
    void update_v_c();
    void update_max_boundary_velocities();
    float mean_v();
    void update_boundary_velocities();

    JointSubSegment3 getSubSegments3() const;
    
    // Copy current parameters into a JointSubSegment, for one phase. 
    void loadJointSubSeg(JointSubSegment& jss, SubSegName phase);

};

typedef std::vector<JointSegment>::iterator jsIter;

