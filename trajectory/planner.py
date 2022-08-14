#
"""

THIS IS THE NEW VERSION OF THE PLANNER. The old one is in segments.py

Create a list of movements for a set of axes, then calculate
an optimum trapezoidal velocity profile.

Each movement is a segment, a period of time in which all axes
must start and stop at the same time. The segments have three periods:

    a: Acceleration
    c: Constant velocity
    d: Deceleration

( the names for the A and D periods refer to an idealized trapezoidal profile,
where a is the first part of the profile and d is the last. in actual cases, a
may have a deceleration, and d an acceleration. The c period is always 0
acceleration )

The optimization procedure attempts to execute every segment in the minimum
time, by ensuring that the longest axis is accelerate to the max velocity
as the max acceleration, and that there is a minimal change in velocity between
adjacent segments.

The algorithms for SimSegment and the c and cn parameters are based on
the article, 'Generate stepper-motor speed profiles in real time',
from https://www.embedded.com/generate-stepper-motor-speed-profiles-in-real-time/

------

Joint segment shapes



"""
from collections import deque
from copy import deepcopy
from enum import Enum
from math import sqrt
from typing import List

import pandas as pd

## Parameters for simulation/step generation

# NUmber of ticks of the step function per second
TIMEBASE = 1_000_000  # ticks per second

# Shapes
TRAPEZOID = 1
TRIANGLE = 2

N_BIG = 2 ** 32

MAX_ERR_T = 1e-4


def sign(x):
    if x == 0:
        return 0
    elif x > 0:
        return 1
    else:
        return -1


def same_sign(a, b):
    return int(a) == 0 or int(b) == 0 or sign(a) == sign(b)


class SegmentError(Exception):
    pass


class ConvergenceError(SegmentError):
    pass


class ConstraintError(SegmentError):
    pass


class ValidationError(SegmentError):
    pass


class TimeRecalc(SegmentError):
    """We need to recalc segment time"""
    pass


def accel_tx(v0, v1, a):
    """Distance required to accelerate from v0 to v1 at acceleration a"""

    dt = (v1 - v0) / a  # Time to change from v0 to v1 at max acceleration

    return abs(((v0 + v1) / 2) * dt), abs(dt)  # Area of velocity trapezoid


def t_accel(vi, vf, a):
    """Time to accelerate from vi to vf"""
    return (vf - vi) / a


def t_accel_x(x, v, a):
    """Time required to move distance x, accelerating  from velocity vi"""

    return sqrt(2 * a * x + v ** 2) / a - (v / a)


def max_v0_for_x(x, a):
    """Maximum V0 for segment to allow decelereration and not cover more than
        distance x; it sets the v_0_max"""

    return sqrt(2 * x / a) * a


def v_c_max(self, x, v_0, v_1, a):
    """Velocity after accelerating for x/2 from the start and end,
    the highest velocity for which t_c is positive.
    This is the maximum achievable velocity before deceleration must begin. """
    ta = t_accel_x(x / 2, v_0, a)
    va = v_0 + a * ta

    td = t_accel_x(x / 2, v_1, a)

    vd = v_1 + a * td

    return min(va, vd)

def nearly_equal_velocity(a, b):
    """Close enough not to trigger an update"""
    try:
        return abs(a-b)/abs(a+b) < .02
    except ZeroDivisionError:
        return a == b # Must still check b/c options are either a == b == 0 or a == -b

def nearly_equal_time(a, b):
    return abs(a-b)/abs(a+b) < .02

class Joint(object):
    """Represents a single joint, with a maximum velocity and acceleration"""

    def __init__(self, v_max: float = None, a_max: float = None, d_max: float = None):
        self.n = None
        self.v_max = v_max
        self.a_max = a_max
        self.d_max = d_max if d_max is not None else self.a_max


class SubSegment(object):
    """A sub segment is a portion of the trajectory with a constant
    acceleration,  one of the aceleration, constant (cruise) or decleration portions. """

    def __init__(self, joint: "Joint", t: float, v_i: float, v_f: float, x: float, ss: float) -> None:
        self.joint = joint
        self.t = t  # Total ssubsegment time
        self.v_i = v_i  # Initial velocity
        self.v_f = v_f  # final velocity
        self.x = x  # Segment distance
        self.ss = ss  # Section, 'a'->accleration, 'c'-> constant, 'd'->deceleration
        self.direction = 1  # 1 for forward, -1 for reverse

        assert v_i == 0 or v_f == 0 or sign(v_i) == sign(v_f), f"Inconsistent directions {v_i} {v_f}"

    def set_direction(self, sign) -> None:
        self.direction = sign
        self.v_i *= sign
        self.v_f *= sign
        self.x *= sign

    def __repr__(self):
        dir = '' if self.direction > 0 else '-'

        return f"<{self.joint.n} {self.ss} {self.t:2.5f} {int(self.x):5d} {self.v_i:5.0f}->{self.v_f:5.0f}> "


class JSClass(Enum):
    """
    ╱▔▔ U Acel. V0 fixed, V1 is maxed
      Should be a typical case.

    ╱╲  A Triangle. Fixed V0, V1, too short for vc
          x is non zero but is too short got velocity to reach v_max

    ╱▔╲ T Trapezoid.  Fixed V0, V1, max vc
          Normal case where v0 and v1 are both below Vc, and the segment isn't short.

    ▔▔▔ C Constant V. Prior and next have high V0, V1
          V0 and V1 are both at vmax, so v_c should be v_max

    ╲▁▁ D Start Decel. Very small x relative to other joints.
          X is so short relative to other joints that we hit v=0
          during accel phase.

    ▔▔╲ L Clif. Next segment is very small x.
          Next segment is very small x or zero

    ▁▁▁ Z Zero distance
          Zero x means that there must bt v=0 through whole segment

    """
    ACEL = 'U'
    TRIANGLE = 'A'
    AT = 'AT'  # A or T, pinned to Zero at start and finish
    TRAPZEZOID = 'T'
    CONSTANT = 'C'
    DECEL = 'D'
    CLIFF = 'L'
    DL = 'DL'  # D or L
    ZERO = 'Z'
    UNK = '?'


kind_icon_map = {
    JSClass.ACEL: '╱▔▔',
    JSClass.TRIANGLE: '╱╲',
    JSClass.AT: '╱▔╲ ╱╲',
    JSClass.TRAPZEZOID: '╱▔╲',
    JSClass.CONSTANT: '▔▔▔ ',
    JSClass.DECEL: '╲▁▁',
    JSClass.CLIFF: '▔▔╲',
    JSClass.ZERO: ' ▁▁▁',
    JSClass.DL: '╲▁▁▔▔╲',
    JSClass.UNK: '?'
}


class JointSegment(object):
    """One segment for one joint"""

    joint: Joint
    next_js: 'JointSegment'
    prior_js: 'JointSegment'

    x: float
    dir: int
    min_t: int  # Minimum segment time.

    v_0_max: float
    v_0: float = 0
    x_a: float = 0
    t_a: float = 0

    v_c_max: float = 0
    x_c: float = 0
    t_c: float = 0
    v_c: float = 0

    x_d: float = 0
    t_d: float = 0
    v_1: float = 0
    v_1_max: float  # Don't allow updating v_1


    params = ['x', 't', 'dir', 'v_0_max', 'v_0', 'x_a', 't_a', 'x_c', 't_c', 'v_c_max', 'v_c', 'x_d', 't_d', 'v_1',
              'v_1_max']

    def __init__(self, joint: Joint = None, x: float = None, init_v0: float = 0):
        """

        :param joint:
        :type joint:
        :param x:
        :type x:
        :param init_v0: # Initial segment velocity, for case when a new segment is being added to existing queue
        :type init_v0:
        """

        self.joint = joint
        self.segment = None  # set after adding to a segment

        self.x = abs(x)
        self.dir = sign(x)
        self.init_v0 = init_v0

        self.next_js = None
        self.prior_js = None

        self.reset()

    def reset(self):
        """Reset to default state, which is that there is no acceleration or
        deceleration, with segment running at minimum time, full speed. """

        vc_m = self.joint.v_max

        # Extra imposed limits on velocity. For instance, if the next segment
        # has zero velocity ( stopped, or zero movement for segment ),
        # then the v_1 has to be zero.
        self.v_0_max = vc_m
        self.v_1_max = vc_m
        self.v_c_max = vc_m

        self.x_a = 0  # Acceleration distance
        self.t_a = 0

        self.x_d = 0  # Deceleration distance
        self.t_d = 0

        if self.x != 0:
            self.v_0 = self.init_v0 if self.init_v0 else vc_m  # Initial velocity
            self.v_c = vc_m
            self.x_c = self.x  # cruise distance
            self.t_c = self.x_c / self.v_c
            self.v_1 = vc_m  # Final velocity

        else:
            self.t_c = 0
            self.x_c = 0
            self.v_1 = 0
            self.v_c = 0
            self.t = 0

        # initial guess at minimum time, which is the distance at a consistant
        # max velocity. This time will get bumped up as we recalc
        self.t = self.x / vc_m  # self.t_a + self.t_c + self.t_d

        self.err_x = 0

    def reset_for_append(self):
        """Perform the resets required when another segment is appended, such
        as removing the v_1_max of 0"""

        self.v_1_max = self.joint.v_max


    def update_boundary_velocities(self):
        """v_1_max may need to be specified if the next segment has a very short
        travel, or zero, because it may not be possible to decelerate during the phase.
        For instance, if next segment has x=0, then v_1 must be 0 """

        v_0_o, v_1_o = self.v_0, self.v_1
        v_0_max_o, v_1_max_o = self.v_0_max, self.v_1_max

        if self.x == 0:
            self.v_0 = 0

            self.v_0_max = 0
            self.v_1_max = 0
            self.v_c_max = 0

        else:

            if not self.prior_js:
                # First segments must maintain their init_v0, which is 0 for entirely new
                # segments, but could also be set to the last velocity of a prior segment.
                self.v_0_max = self.init_v0
            elif self.dir != self.prior_js.dir:
                # Changes direction, so boundary must be zero
                self.v_0_max = 0
            else:
                # Later segments must have a incomming velocity to match
                # the prior segments outgoing velocity
                self.v_0_max = min(self.joint.v_max, self.v_0_max) # min() reduces instability

            if self.next_js:
                if self.dir != self.next_js.dir:
                    # Changes direction, so must be zero
                    self.v_1_max = 0
                elif self.next_js.x == 0:
                    # No movement in next segment, so must be zero, because the
                    # next segment can't have any acel ( actually decel ) time
                    self.v_1_max = 0
                else:
                    # The next segment is just the normal case so just match the
                    # velocities across the boundary.
                    self.v_1_max = min(self.joint.v_max, self.v_1_max) # min() reduces instability
            else:
                # No next segment, so we have to decel to zero.
                self.v_1_max = 0

        # Restrict v_0_max if there is a small distance to travel an
        # we will have to decelerate to not overshoot it.
        self.v_0_max = min(max_v0_for_x(self.x, self.joint.a_max), self.v_0_max)

        # subtract off the accel-phase distance and re-calc for the decel phase
        # ( Note that in many cases, there will be deceleration in the accel phase,
        # and acceleration in the decel phase, because this handles the case of very
        # small movements.

        x, t = accel_tx(self.v_0_max, 0, self.joint.a_max)  # Recompute time and dist in acell phase
        self.v_1_max = min(max_v0_for_x(self.x - x, self.joint.a_max), self.v_1_max)

        if self.v_0 > self.v_0_max:
            self.v_0 = self.v_0_max

        if self.v_1 > self.v_1_max:
            self.v_1 = self.v_1_max

        ## Change reporting

        if not nearly_equal_velocity(self.v_0_max, v_0_max_o):
            self.segment.mark_for_update(f'{self.id} ubv_v0max {v_0_max_o}->{self.v_0_max}')
            self.segment.mark_prior_for_update(f'{self.id} ubv_v0max {v_0_max_o}->{self.v_0_max}')

        if not nearly_equal_velocity(self.v_1_max, v_1_max_o):
            self.segment.mark_for_update(f'{self.id} ubv_v1max {v_1_max_o}->{self.v_1_max}')
            self.segment.mark_next_for_update(f'{self.id} ubv_v1max {v_1_max_o}->{self.v_1_max}')

        if not nearly_equal_velocity(self.v_0, v_0_o):
            self.segment.mark_for_update(f'{self.id} ubv_v0 {v_0_o}->{self.v_0}')
            self.segment.mark_prior_for_update(f'{self.id} ubv_v0 {v_0_o}->{self.v_0}')

        if not nearly_equal_velocity(self.v_1, v_1_o):
            self.segment.mark_for_update(f'{self.id} ubv_v0 {v_1_o}->{self.v_1}')
            self.segment.mark_prior_for_update(f'{self.id} ubv_v0 {v_1_o}->{self.v_1}')


    def recalc(self):
        """Re-calculate v_c and v_1, based on  total distance and velocity limits"""

        v_0_o, v_1_o = self.v_0, self.v_1
        v_0_max_o, v_1_max_o = self.v_0_max, self.v_1_max

        if self.prior_js:
            self.v_0 = self.prior_js.v_1
            self.v_0_max = min(self.prior_js.v_1_max, self.v_0_max)

        if self.next_js:
            self.v_1_max = min(self.next_js.v_0_max, self.v_1_max)

        self.x_a, self.t_a = accel_tx(self.v_0, self.v_c, self.joint.a_max)
        self.x_d, self.t_d = accel_tx(self.v_c, self.v_1, self.joint.a_max)

        # Remaining time is for cruise
        self.x_c = self.x - self.x_a - self.x_d
        self.t_c = self.segment.t - self.t_a - self.t_d

        if self.t_c < 0:
            self.t_c = 0;  # This will cause x error that will get fixed on the next recalc
            raise TimeRecalc(f'{self.x_c}|{self.segment.t} {self.t_a},{self.x_a} {self.t_d},{self.x_d}')
            self.t = self.t_a + self.t_d
            self.segment.mark_for_update(f'{self.id} rc_tc_lt0')

        if self.x_c < 0:
            self.segment.mark_for_update(f'{self.id} rc_xc_lt0')
            #self.x_c = 0;  # This will cause x error that will get fixed on the next recalc
            #raise TimeRecalc(f'{self.x_c}|{self.segment.t} {self.t_a},{self.x_a} {self.t_d},{self.x_d}')
            pass

        # Re calc cruise velocity, and try to set the exit velocity to be
        # the same as the cruise velocity

        v_c = 0 if self.t_c == 0 else self.x_c / self.t_c

        self.v_c = min(max(v_c, 0), self.joint.v_max)

        v_1 = min(self.v_1_max, self.v_c)

        if not nearly_equal_velocity(self.v_1,v_1):
            self.segment.mark_for_update(f'{self.id} rc_v1_change {self.v_1}->{v_1}')
            self.segment.mark_next_for_update(f'{self.id} rc_v1_change {self.v_1}->{v_1}')
            self.v_1 = v_1

        # Recalculate segment time, which may differ from self.segment.t.
        # Later, we will collect these, re-calc the segment time ( it may get longer )
        # the recalc all the joint segments, and do that repeatedly until
        # error is minimized.

        # self.t_c ought to equal zero if self.v_c  == 0, but well get
        # that on the next recalc, because if it isn't there will be err_t
        t_c = self.x_c / self.v_c if self.v_c != 0 else 0

        t = t_c + self.t_a + self.t_d

        if self.t != t:
            self.segment.mark_for_update(f'{self.id} rc_t_change {t}')
            self.t = t

        self.err_x = self.x - self.calc_x

        if not nearly_equal_velocity(v_0_max_o,self.v_0_max):
            self.segment.mark_prior_for_update(f'{self.id} rc_v0max  {v_0_max_o}->{self.v_0_max}')

        if not nearly_equal_velocity(v_1_max_o, self.v_1_max):
            self.segment.mark_prior_for_update(f'{self.id} rc_v1max {v_1_max_o}->{self.v_1_max}')

        return True

    def update(self):

        self.update_boundary_velocities()

        self.recalc()
        if self.t != self.segment.t:
            self.segment.update_min_time()
            self.recalc()

    @property
    def id(self):
        return f"{self.segment.n}/{self.joint.n}"

    @property
    def min_t(self):
        return self.t
        # return max(self.min_transit_t, self.min_accel_time)

    @property
    def min_transit_t(self):
        """Minimum transit time based on distance and max velocity, """

        # Time to accelerate to max v
        def v_max_t():
            x_a, t_a = accel_tx(self.v_0, self.joint.v_max, self.joint.a_max)

            x_c = self.x - x_a  # distance left to run after hitting max v

            return t_a + x_c / self.joint.v_max

        return max(v_max_t(), self.t if self.t is not None else 0)

    @property
    def min_accel_time(self):
        """ Minimum time, considering accelerations
        :return:
        :rtype:
        """
        a = self.joint.a_max

        if self.v_c_max == 0:
            return 0

        # Carve out space for acceleration and deceleration.
        x_a, t_a = accel_tx(self.v_0_max, self.v_c_max, a)
        x_d, t_d = accel_tx(self.v_c_max, self.v_1_max, a)

        x = self.x - x_a - x_d

        if x < 0:
            # time to accellerate x/2, then decelerate
            ta = abs(t_accel_x(x / 2, self.v_c_max, a))
            td = abs(t_accel_x(x / 2, self.v_c_max, a))

            return ta + td
        else:
            return t_a + t_d + x / self.v_c_max

    @property
    def calc_x(self):
        """Calculate the distance traveled from the area of the velocity trapezoid.
        This is an alternate calculation of x that can be used to verify parameters"""

        x_a = (self.v_0 + self.v_c) / 2 * self.t_a
        x_d = (self.v_c + self.v_1) / 2 * self.t_d
        x_c = self.v_c * self.t_c

        return x_a + x_c + x_d

    @property
    def calc_t(self):
        """Calculate the time required for each subseg from the distance an velocity"""

        try:
            t_a = 2 * self.x_a / (self.v_0 + self.v_c)
        except ZeroDivisionError:
            t_a = 0

        try:
            t_d = 2 * self.x_d / (self.v_c + self.v_1)
        except ZeroDivisionError:
            t_d = 0

        try:
            t_c = self.x_c/self.v_c
        except ZeroDivisionError:
            t_c = 0

        return t_a + t_c + t_d

    @property
    def sum_x(self):
        """Calculate the distance from the individual components distances"""
        return self.x_a + self.x_c + self.x_d

    @property
    def subsegments(self):
        ss = []

        ss.append(
            SubSegment(self.joint, round(self.t_a, 7), round(self.v_0, 2), round(self.v_c, 2), round(self.x_a, 0),
                       'a'))

        ss.append(
            SubSegment(self.joint, round(self.t_c, 7), round(self.v_c, 2), round(self.v_c, 2), round(self.x_c, 0),
                       'c'))

        ss.append(
            SubSegment(self.joint, round(self.t_d, 7), round(self.v_c, 2), round(self.v_1, 2), round(self.x_d, 0),
                       'd'))

        return ss

    @property
    def classify(self):

        if self.x == 0:
            return JSClass.ZERO
        elif self.v_0_max == 0 and self.v_1_max == 0:
            return JSClass.AT
        elif self.v_0_max == 0 and self.v_1_max != 0:
            return JSClass.ACEL
        elif self.v_0_max != 0 and self.v_1_max == 0:
            return JSClass.DL
        elif self.v_0_max != 0 and self.v_1_max != 0:
            return JSClass.CONSTANT
        else:
            return JSClass.UNK

    def __str__(self):
        from colors import color, bold

        v0 = color(f"{int(self.v_0):<5d}", fg='green')
        xa = bold(f"{int(self.dir * self.x_a):>6d}")
        c = bold(f"{int(round(self.dir * self.x_c)):>6d}")
        vc = color(f"{int(self.v_c):<6d}", fg='blue')
        xd = bold(f"{int(self.dir * self.x_d):<6d}")
        v1 = color(f"{int(self.v_1):>5d}", fg='red')

        return f"[{v0} {xa}↗{c + '@' + vc}↘{xd} {v1}]"

    def __repr__(self):
        return self.__str__()

    @property
    def debug(self):

        from colors import color, bold

        v0_max = color(f"{int(self.v_0_max):<5d}", fg='green')
        v0 = color(f"{int(self.v_0):<5d}", fg='green')
        xa = bold(f"{int(self.dir * self.x_a):>6d}")
        ta = bold(f"{int(self.dir * self.t_a):<1.5f}")
        c = bold(f"{int(round(self.dir * self.x_c))}")
        vc = color(f"{int(self.v_c)}", fg='blue')
        tc = bold(f"{int(self.dir * self.t_c):<1.5f}")
        xd = bold(f"{int(self.dir * self.x_d):>6d}")
        td = bold(f"{int(self.dir * self.t_d):<1.5f}")
        v1 = color(f"{int(self.v_1):>5d}", fg='red')
        v1_max = color(f"{int(self.v_1_max):>5d}", fg='red')

        return f"[{v0_max}|{v0} {xa}%{ta} ↗ {c}@{vc}%{tc} ↘ {xd}%{td} {v1}|{v1_max}]"


class Segment(object):
    """One segment, for all joints"""

    n: int = None;
    next_seg: "Segment" = None
    prior_seg: "Segment" = None
    t: float = 0

    needs_update: bool = False

    n_updates: int = 0

    updates: List[str]

    def __init__(self, n, joint_segments: List[JointSegment]):

        self.n = n
        self.joint_segments = joint_segments

        for js in self.joint_segments:
            js.segment = self

        self.n_updates = 0

        self.updates = []

    def mark_for_update(self, reason=None):
        self.needs_update = True

        if reason:
            self.updates.append(reason)


    def mark_next_for_update(self, reason=None):
        if self.next_seg:
            self.next_seg.mark_for_update(reason)

    def mark_prior_for_update(self, reason=None):
        if self.prior_seg:
            self.prior_seg.mark_for_update(reason)


    def clear_needs_update(self):
        self.needs_update = False

    def reset(self):
        for js in self.joint_segments:
            js.reset()

    def reset_for_append(self):
        for js in self.joint_segments:
            js.reset_for_append()

    def recalc(self):
        need_tc = False;

        for s in self.joint_segments:
            need_tc = not s.recalc() or need_tc

    def link(self, prior_segment: "Segment"):
        """Link segments together"""

        prior_segment.next_seg = self
        self.prior_seg = prior_segment

        for prior_js, next_js in zip(self.prior_seg.joint_segments, self.joint_segments):
            prior_js.next_js = next_js
            next_js.prior_js = prior_js

    def update_boundary_velocities(self):
        for js in self.joint_segments:
            if js.prior_js:
                js.prior_js.update_boundary_velocities()
            js.update_boundary_velocities()

    def update_segment_time(self):
        t = [js.min_t for js in self.joint_segments]
        self.t = max(t)

    @property
    def err_t(self):
        """RMS error between joint segment calculated times and segment time"""
        err = sqrt(sum([(js.t - self.t) ** 2 for js in self.joint_segments]))

        return err/len(self.joint_segments)/self.t

    @property
    def err_x(self):
        """RMS error between joint segment calculated times and segment time"""
        dist = sum(js.x for js in self.joint_segments)
        err =  sqrt(sum([(js.err_x) ** 2 for js in self.joint_segments]))

        return err/dist

    def update_min_time(self):
        self.t = self.min_t

    @property
    def min_t(self):
        return max([js.min_t for js in self.joint_segments])

    @property
    def subsegments(self):
        for js in self.joint_segments:
            yield from js.subsegments

    def __str__(self):

        return f"{self.t:>3.4f}|" + ' '.join(str(js) for js in self.joint_segments)


class SegmentList(object):
    positions: List[float]  # Positions after last movement addition
    segments: deque[Segment]
    all_segments: List[Segment]  # All of the segments, unprocessed

    def __init__(self, joints: List[Joint]):

        self.joints = deepcopy(joints)

        self.directions = [0] * len(self.joints)

        for i, j in enumerate(self.joints):
            j.n = i

        self.segments = deque()
        self.positions = [0] * len(self.joints)

    def rmove(self, joint_distances: List[int], update=True):
        """Add a new segment, with joints expressing joint distance
        :type joint_distances: object
        """

        assert len(joint_distances) == len(self.joints)

        next_seg = Segment(len(self.segments),
                           [JointSegment(j, x=x) for j, x in zip(self.joints, joint_distances)])

        if len(self.segments) > 0:
            self.segments[-1].reset_for_append()
            self.segments[-1].next_segment = next_seg
            next_seg.prior_seg = self.segments[-1]
            next_seg.link(self.segments[-1])


        self.segments.append(next_seg)


        if update:
            self.update(next_seg)
            # Possibly run several rounds of updates to get the total error down.

        return next_seg

    def reset(self):
        for s in self.segments:
            s.reset()

    def _lim_segments(self, limit=4):
        if limit:
            return list(self.segments)[-limit:]
        else:
            return list(self.segments)

    def recalc(self, limit=3):
        for s in self._lim_segments(limit):
            s.recalc()

    def update_min_time(self, limit=4):
        for s in self._lim_segments(limit):
            s.update_min_time()

    def update(self, s=None):

        if s:
            s.update_boundary_velocities()
            s.update_segment_time()

        for i in range(4):
            for s in self.segments:
                if s.needs_update:
                    s.update_boundary_velocities()
                    s.update_segment_time()
                    s.recalc()
                    s.clear_needs_update()

            if self.count_needs_update() == 0:
                break

    def err_t(self, limit=4):

        l = self._lim_segments(limit)

        et = sum([s.err_t for s in l]) / len(l)

        if et is None:
            return -1

        return et

    def err_x(self, limit=4):

        l = self._lim_segments(limit)

        et = sum([s.err_x for s in l]) / len(l)

        if et is None:
            return -1

        return et

    def count_needs_update(self):
        return sum( s.needs_update for s in self.segments)

    @property
    def dataframe(self):

        rows = []

        for ss in self.subsegments:
            rows.append([None, ss.joint.n, ss.x, ss.v_i, ss.v_f, ss.ss, ss.t])

        df = pd.DataFrame(rows, columns=' t axis x v_i v_f ss del_t'.split())

        df['t'] = df.groupby('axis').del_t.cumsum()

        return df

    @property
    def dataframe_stacked(self):
        return self.dataframe.set_index(['t', 'axis']).stack().unstack(-2).unstack(-1)

    @property
    def params(self):
        """
        :return: a dataframe of paramaters
        :rtype:
        """

        from operator import attrgetter

        ag = attrgetter(*JointSegment.params)

        rows = []
        for i, s in enumerate(self.segments):
            for j, js in enumerate(s.joint_segments):
                d = (i, j, js.segment.t) + ag(js) + (js.calc_x, js.sum_x, js.calc_t)
                rows.append(d)

        df =  pd.DataFrame(rows, columns=['seg', 'js', 'seg_t'] + JointSegment.params + ['calc_x', 'sum_x', 'calc_t'])

        for c in ['v_c','v_0_max','v_1_max','v_0','v_1', 'calc_x','sum_x']:
            df[c] = df[c].astype('int')

        for c in ['seg_t','t','t_a','t_d','t_c']:
            df[c] = df[c].round(5)

        return df

    @property
    def subsegments(self):

        for s in self.segments:
            yield from s.subsegments

    def __getitem__(self, item):
        """Return a joint segment by the id"""
        try:
            s, js = item
            return self.segments[s].joint_segments[js]
        except TypeError:
            return self.segments[item]

    def __iter__(self):
        return iter(self.segments)

    def __len__(self):
        return len(self.segments)

    def __str__(self):
        return '\n'.join(str(s) for s in self.segments)

    def _repr_pretty_(self, p, cycle):
        p.text(str(self) if not cycle else '...')

    def debug_str(self):
        return '\n'.join(s.debug_str() for s in self.segments)


__all__ = ['SegmentList', 'Segment', 'Joint', 'JointSegment', 'SegmentError']
