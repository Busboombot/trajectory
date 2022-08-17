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
from collections import namedtuple
from copy import deepcopy
from enum import Enum
from math import sqrt
from typing import List
from warnings import warn

import pandas as pd

from .trapmath import *

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

class BoundaryError(SegmentError):
    pass

class ConstraintError(SegmentError):
    pass


class ValidationError(SegmentError):
    pass

class ShortSegment(SegmentError):
    """Segment is too short"""
    pass

class TimeRecalc(SegmentError):
    """We need to recalc segment time"""
    pass


def nearly_equal_velocity(a, b):
    """Close enough not to trigger an update"""
    try:
        return abs(a - b) / abs(a + b) < .02
    except ZeroDivisionError:
        return a == b  # Must still check b/c options are either a == b == 0 or a == -b


def nearly_equal_time(a, b):
    return abs(a - b) / abs(a + b) < .02


Params = namedtuple('Params', 't x t_a t_c t_d x_a x_c x_d v_0 v_c v_1 is_triangle input'.split())
InputParams = namedtuple('InputParams', 'x v_0 v_c v_1 a_max t'.split())


def inital_parameters(x, v_0, v_c, v_1, a_max, t=None):
    """Find the  lowest time for a segment, and
    report any parameters required to get it"""

    t_inp = t

    x_a, t_a = accel_tx(v_0, v_c, a_max)
    x_d, t_d = accel_tx(v_c, v_1, a_max)

    if (x < x_a + x_d):
        # Triangle segment, so there is no cruise time

        # Calculate new top velocity, the velocity we
        # get to accelerating from v_0 when we hit the center, where
        # we have to turn around and decelerate.

        v_0 = min(max_v0_for_x(x, a_max), v_0)
        v_1 = min(max_v0_for_x(x, a_max), v_1)

        v_c = min(v_c_max(x, v_0, v_1, a_max), v_c)

        x_a, t_a = accel_tx(v_0, v_c, a_max)
        x_d, t_d = accel_tx(v_c, v_1, a_max)

        t_c = 0
        x_c = 0
        is_triangle = True
    else:
        x_c = x - x_a - x_d
        t_c = x_c / v_c
        is_triangle = False

    t = t_a + t_c + t_d
    x = x_a + x_c + x_d

    return Params(t, x, t_a, t_c, t_d, x_a, x_c, x_d, v_0, v_c, v_1, is_triangle,
                  InputParams(x, v_0, v_c, v_1, a_max, t_inp))


def max_v1(x, t, v_0, v_c, a_max, **kwds):
    """ Calculate initial parameters, but with a variable v_1 and a fixed t"""

    ip = InputParams(x, v_0, v_c, 0, a_max, t)

    # v_c = min(a_max * (t/2), v_max) # Velocity after accelerating to midpoint

    t_a = t_accel(v_0, v_c, a_max)  # Time after accelerating to v_max
    x_a = (v_0 + v_c) / 2 * t_a  # Distance covered to get to max v

    t_r = t - t_a  # Remaining time

    v_mean = ((x - x_a) / 2) / t_r  # Mean velocity for remaining half

    # mean v == (v_c+v_1)/2
    v_1 = 2 * v_mean - v_c

    if v_1 > v_c:
        pass
        #warn(f"{kwds.get('id')} v_1 {v_1} exceeds v_c {v_c}")
        #raise ValidationError()

    if v_1 < 0:
        # negative v_1 means that we have some space for t_c

        v_1 = 0
        t_d = t_accel(v_1, v_c, a_max)
        assert t_d >= 0, (t_d, ip)
        x_d = (v_c + v_1) / 2 * t_d
        t_c = t - t_a - t_d
        x_c = x - x_a - x_d

        is_triangle = False

    else:
        t_c = 0
        x_c = 0
        t_d = 0
        x_d = 0
        is_triangle = True

    return Params(t, x, t_a, t_c, t_d, x_a, x_c, x_d, v_0, v_c, v_1, is_triangle, ip)


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

    def __init__(self, id: str, seg_n: int, joint: "Joint", js: "JointSegment",
                 t: float, v_i: float, v_f: float, x: float, ss: float) -> None:
        self.id = id
        self.segment_n = seg_n
        self.joint = joint
        self.js = js
        self.t = t  # Total ssubsegment time
        self.v_i = v_i  # Initial velocity
        self.v_f = v_f  # final velocity
        self.x = x  # Segment distance
        self.ss = ss  # Section, 'a'->accleration, 'c'-> constant, 'd'->deceleration
        self.direction = 1  # 1 for forward, -1 for reverse

        assert v_i == 0 or v_f == 0 or sign(v_i) == sign(
            v_f), f"Inconsistent directions {v_i} {v_f} for {self.id}{self.ss} "

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

    ╱-╲ T Trapezoid.  V0 & V1 == 0, max vc
          Normal case where v0 and v1 are both below Vc, and the segment isn't short.

    ╱▔╲ P pentagon.  Fixed V0, V1, but not zero, max vc
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
    PENTAGON = 'P'
    AT = 'AT'  # A or T, pinned to Zero at start and finish
    TRAPZEZOID = 'T'
    TROUGH = 'R'
    CONSTANT = 'C'
    DECEL = 'D'
    CLIFF = 'L'
    DL = 'DL'  # D or L
    ZERO = 'Z'

    UNK = '?'


kind_icon_map = {
    JSClass.ACEL: '╱▔▔',
    JSClass.TRIANGLE: '╱╲',
    JSClass.TRAPZEZOID: '╱-╲',
    JSClass.PENTAGON: '╱▔╲',
    JSClass.TROUGH: '╲▁╱',
    JSClass.CONSTANT: '▔▔▔ ',
    JSClass.DECEL: '╲▁▁',
    JSClass.CLIFF: '▔▔╲',
    JSClass.ZERO: '▁▁▁',
    JSClass.UNK: '?'
}

def classify(x, v_0, v_1, v_max, a_max):

    x_a, t_a = accel_tx(v_0, v_max, a_max)
    x_d, t_d = accel_tx(v_max, v_1, a_max)

    if x == 0:
        return JSClass.ZERO
    elif v_0 == 0 and v_1 == 0:
        if x_a + x_d > x:
            return JSClass.TRIANGLE
        else:
            return JSClass.TRAPZEZOID
    elif v_0 != 0 and v_1 != 0:
        x_break = (v_0 + v_1) / 2 * (t_a + t_d)
        return JSClass.PENTAGON if x > x_break else JSClass.TROUGH
    elif v_0 == 0 and v_1 != 0:
        return JSClass.ACEL
    elif v_0 != 0 and v_1 == 0:
        x_break = (v_0 + v_1) / 2 * (t_a+t_d)
        return JSClass.CLIFF if x > x_break else JSClass.DECEL
    elif v_0 == v_1:
        return JSClass.CONSTANT

    else:
        return JSClass.UNK

class JointSegment(object):
    """One segment for one joint"""

    joint: Joint
    segment: 'Segment'
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

        self.needs_update_boundaries = False
        self.needs_recalc = False

        self.reset()

    def mark_need_recalc(self):
        self.segment.needs_recalc = True

    def set_v_0(self, v):
        v_0 = min(self.v_0_max, v)
        if v_0 != self.v_0:
            self.v_0 = v_0

            # Mark the prior segment for updating v_1
            if self.prior_js and self.prior_js.v_1 != self.v_0:
                self.prior_js.segment.needs_update_boundaries = True

    def set_v_0_max(self, v):
        assert v<=self.v_0_max, f'Should not be increasing v_0_max ({self.v_0_max}->{v})'
        if self.v_0_max != v:
            self.v_0_max = v
            self.set_v_0(min(self.v_0, self.v_0_max))

    def set_v_1(self, v):
        v_1 = min(self.v_1_max, v)
        if v_1 != self.v_1:
            self.v_1 = v_1

    def set_v_1_max(self, v):
        assert v <= self.v_1_max, 'Should not be increasing v_1_max'
        if self.v_1_max != v:
            self.v_1_max = v
            self.set_v_1(min(self.v_1, self.v_1_max))

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

    def initialize(self):
        self.update_boundary_velocities()
        self.set_initial_parameters()

    def initialize_for_append(self):
        """Perform the resets required when another segment is appended, such
        as removing the v_1_max of 0"""

        self.v_1_max = self.joint.v_max
        self.update_boundary_velocities()
        #self.reset()
        #self.initialize()


    def set_initial_parameters(self):
        """Set initial values, which will get adjusted later as we link up adjacent
        segments. But these are pretty good guesses"""

        if self.x == 0:
            self.set_v_0(0)
            self.v_0_max = 0
            self.v_1_max = 0
            self.v_c_max = 0
            self.t = 0

        else:

            # If this is a short segment, then there is a restriction on the
            # starting velocity, which imposes a backwards restriction on the
            # boundary velocities. In the limit case, x=0 for this segment, the
            # v_0 must be zero, so the v_1 of the prior segment must be zer.

            v_0_max = min(max_v0_for_x(self.x, self.joint.a_max), self.joint.v_max, self.v_0_max)

            p = inital_parameters(self.x, v_0_max, self.joint.v_max , self.v_1_max, self.joint.a_max)
            v_0 = v_0_max = p.v_0

            self.x_a = p.x_a
            self.t_a = p.t_a

            self.v_c = self.v_c_max = p.v_c
            self.t_c = p.t_c
            self.x_c = p.x_c

            self.x_d = p.x_d
            self.t_d = p.t_d

            self.v_1 = self.v_1_max = p.v_1

            self.t = self.t_a + self.t_c + self.t_d

            assert self.x == 0 or self.t > 0, 'Time should not be zero'

            self.set_v_0_max(v_0_max)
            self.set_v_0(v_0)

    def update_initial_parameters(self, t):
        """"""

        if self.x == 0:
            self.t = t
            return

        # Re-calculate the initial parameters, but with a time constraint.
        #  max_v0_for_x is the maximum v_0 we can have and not cover more than
        #  distance x in time t, so it can set v_0_max

        v_0_max = min(max_v0_for_x(self.x,self.joint.a_max ), self.joint.v_max, self.v_0_max)

        v_0 = min(self.v_0, v_0_max)

        p = max_v1(self.x, t, v_0, self.joint.v_max, self.joint.a_max, id=self.id)

        self.x_a = p.x_a
        self.t_a = p.t_a

        self.v_c = self.v_c_max = p.v_c
        self.t_c = p.t_c
        self.x_c = p.x_c

        self.x_d = p.x_d
        self.t_d = p.t_d

        self.set_v_1_max(max(self.v_1, self.v_1_max))
        self.set_v_1(p.v_1)

        self.set_v_0_max(v_0_max)
        #self.set_v_0(p.v_0) # max_v1 shouldn't change this

        self.t = self.t_a + self.t_c + self.t_d

        if self.t == 0:
            self.t = self.t_c = self.segment.t

        assert self.x == 0 or self.t > 0, f'Time  should not be zero'

    def update_boundary_velocities(self):
        """v_1_max may need to be specified if the next segment has a very short
        travel, or zero, because it may not be possible to decelerate during the phase.
        For instance, if next segment has x=0, then v_1 must be 0

        Boundary velocity limits should only need to be recalculated when:
            - A segment is newly added to the queue
            - A segment gets a new next segment ( one is appended to the end)
        """

        v_0_max = self.v_0_max
        v_1_max = self.v_1_max
        v_0 = self.v_0
        v_1 = self.v_1

        if not self.prior_js:
            # First segments must maintain their init_v0, which is 0 for entirely new
            # segments, but could also be set to the last velocity of a prior segment.
            v_0_max = self.init_v0

        elif self.dir != self.prior_js.dir:
            # Changes direction, so boundary must be zero
            v_0_max = 0
        else:
            # Later segments must have a incoming velocity to match
            # the prior segments outgoing velocity
            v_0_max = min(self.joint.v_max, v_0_max)  # min() reduces instability
            #v_0 = self.prior_js.v_1

        if self.next_js:
            if self.dir != self.next_js.dir:
                # Changes direction, so must be zero
                v_1_max = 0
            elif self.next_js.x == 0:
                # No movement in next segment, so must be zero, because the
                # next segment can't have any accel ( actually decel ) time
                v_1_max = 0
            else:
                # The next segment is just the normal case so just match the
                # velocities across the boundary.
                v_1_max = max(self.joint.v_max, v_1_max)  # min() reduces instability
        else:
            # No next segment, so we have to decel to zero.
            v_1_max = 0

        # Restrict v_0_max if there is a small distance to travel an
        # we will have to decelerate to not overshoot it.

        v_0_max = min(max_v0_for_x(self.x, self.joint.a_max), v_0_max)

        # subtract off the accel-phase distance and re-calc for the decel phase
        # ( Note that in many cases, there will be deceleration in the accel phase,
        # and acceleration in the decel phase, because this handles the case of very
        # small movements.

        x, t = accel_tx(v_0_max, 0, self.joint.a_max)  # Recompute time and dist in accel phase

        if x > self.x:
            # if self.x-x <0, max_v0_for_x will fail ( sqrt(-1) )
            v_1_max = 0
        else:
           v_1_max = min(max_v0_for_x(self.x - x, self.joint.a_max), v_1_max)

        v_0 = min(v_0_max, v_0)
        v_1 = min(v_1_max, v_1)

        self.set_v_0(v_0)
        self.set_v_1(v_1)

        self.set_v_0_max(v_0_max)
        self.set_v_1_max(v_1_max)

    def link_velocities_forward(self):
        """Link the starting velocity in this segment to the end of the prior segment

        This should only be required when the v_1 of the prior segment changes.
        """

        if self.next_js:
            if self.v_1 != self.next_js.v_0:
                self.set_v_1(self.next_js.v_0)

            if self.v_1_max > self.next_js.v_0_max:
                self.set_v_1_max(self.next_js.v_0_max)
        else:
            assert False, "Should not have called this routine"

    def link_velocities_backward(self):
        """Link the starting velocity in this segment to the end of the prior segment

        This should only be required when the v_1 of the prior segment changes.
        """

        if self.prior_js and self.v_0 != self.prior_js.v_1:
                self.set_v_0(self.prior_js.v_1)
                self.v_0_max = self.prior_js.v_1

    def recalc(self):
        """Re-calculate v_c and v_1, based on  total distance and velocity limits"""

        self.set_v_0(min(self.v_0, self.v_0_max))
        self.set_v_1(min(self.v_1, self.v_1_max))

        if self.x == 0:
            self.t_c = self.t = self.segment.t
            return

        # Accel times and distances, assuming accel to self.v_c
        self.x_a, self.t_a = accel_tx(self.v_0, self.v_c, self.joint.a_max)
        self.x_d, self.t_d = accel_tx(self.v_c, self.v_1, self.joint.a_max)

        if self.x <= self.x_a + self.x_d:

            # The x break is the line, on the v vs t plot, from v_0 to v_1. If
            # x is above the break, the a phase accelerates, and if it's below, the a
            # phase decelerates. It is above the break, we have a triangle, below is a trough
            x_break = (self.v_0 + self.v_1) / 2 * self.t

            if self.x > x_break:
                # Triangle case ( will be trap after expansion )
                # A triangle segment is weird because it will have t_c = 0, but a non-zero
                # v_c. The v_c is just a target for accel and decel
                self.v_c = self.x / self.t
            else:
                # Trough case. Similar to triangle, but v_c  = 0
                self.v_c = 0

            self.x_a, self.t_a = accel_tx(self.v_0, self.v_c, self.joint.a_max)
            self.x_d, self.t_d = accel_tx(self.v_c, self.v_1, self.joint.a_max)
            self.x_c = self.x - self.x_a - self.x_d
            self.t_c = self.segment.t - self.t_a - self.t_d

        else:
            self.x_c = self.x - self.x_a - self.x_d
            self.t_c = self.segment.t - self.t_a - self.t_d

            if self.t_c <= 0:
                # Uh no ... we're about to get a divide by zero, so punt
                # this is likely because other joints in this segment changed the time
                # and this one ended up in a triangle.

                self.t = self.t_a + self.t_d
                raise ShortSegment('tc < 0')

                # Need v_c, even if t_c ==0, because it is the peak velocity
                self.v_c = self.joint.a_max * self.t_a + self.v_0

            self.v_c = min(self.x_c / self.t_c, self.joint.v_max)

        self.set_v_1(min(self.v_1_max, self.v_c))

        # Recalculate segment time, which may differ from self.segment.t.
        # Later, we will collect these, re-calc the segment time ( it may get longer )
        # the recalc all the joint segments, and do that repeatedly until
        # error is minimized.

        self.t = self.t_c + self.t_a + self.t_d

        if self.t == 0: # Need to have some time in segment, even if no movement.
            self.t_c = self.segment.t
            self.t = self.t_c + self.t_a + self.t_d

        # In this case, even if the distance is zero, we should not have a zero time;
        # it should have the segment tim.
        assert self.t > 0, 'Time should not be zero'

        if self.x > 0:
            assert self.x_a + self.x_c + self.x_d > 0


        return True

    def update(self):
        pass

    @property
    def id(self):
        return f"{self.segment.n}/{self.joint.n}"

    @property
    def classify(self):
        return classify(self.x, self.v_0, self.v_1, self.joint.v_max, self.joint.a_max)

    @property
    def is_triangle(self):
        """Return true if this will be a triangle segment, with x_c == 0"""
        x_a, t_a = accel_tx(self.v_0, self.joint.v_max, self.joint.a_max)
        x_d, t_d = accel_tx(self.joint.v_max, self.v_1, self.joint.a_max)

        return self.x > x_a + x_d

    @property
    def min_t(self):
        return self.t
        # return max(self.min_transit_t, self.min_accel_time)

    @property
    def err_x(self):
        """Return true of the subsegment error is larger than 3% of the x"""
        if self.x != 0:
            return self.subseg_error
        else:
            return False

    @property
    def calc_x(self):
        """Calculate the distance traveled from the area of the velocity trapezoid.
        This is an alternate calculation of x that can be used to verify parameters"""

        x_a = (self.v_0 + self.v_c) / 2 * self.t_a
        x_d = (self.v_c + self.v_1) / 2 * self.t_d
        x_c = self.v_c * self.t_c

        return x_a + x_c + x_d

    @property
    def subseg_error(self):
        """Calculate the x error, but by subsegments rather than the whole segment. This
        will catch errors in subsegment that would cancel out. """

        x_a = (self.v_0 + self.v_c) / 2 * self.t_a
        x_d = (self.v_c + self.v_1) / 2 * self.t_d
        x_c = self.v_c * self.t_c

        return abs(self.x_a-x_a)+abs(self.x_c-x_c)+abs(self.x_d-x_d)

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
            t_c = self.x_c / self.v_c
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
            SubSegment(self.id, self.segment.n, self.joint, self,
                       round(self.t_a, 7), round(self.v_0, 2), round(self.v_c, 2),round(self.x_a, 0),
                       'a'))

        ss.append(
            SubSegment(self.id, self.segment.n, self.joint, self,
                       round(self.t_c, 7), round(self.v_c, 2), round(self.v_c, 2), round(self.x_c, 0),
                       'c'))

        ss.append(
            SubSegment(self.id, self.segment.n, self.joint, self,
                       round(self.t_d, 7), round(self.v_c, 2), round(self.v_1, 2), round(self.x_d, 0),
                       'd'))

        return ss


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

    needs_recalc: bool = False
    needs_update_boundaries: bool = False

    n_updates: int = 0

    updates: List[str]

    def __init__(self, n, joint_segments: List[JointSegment]):

        self.n = n
        self.joint_segments = joint_segments

        for js in self.joint_segments:
            js.segment = self

        self.n_updates = 0

        self.updates = []

    @property
    def needs_update_parameters(self):
        """Return true if there is a segment with a time that is shorter that
        the segment time, which means we need to run update_initial_parameters"""
        return self.t != min([js.t for js in self])

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

    def initialize(self):
        for js in self:
            js.initialize()
        self.update_min_time()

    def initialize_for_append(self):
        for js in self:
            js.initialize_for_append()
        self.update_min_time()


    def update(self, _warn = True):

        self.recalc()
        self.update_initial_parameters()
        self.update_min_time()

        for i in range(10):
            self.recalc()

            if self.max_err_x < 3:
                return True
        else:
            if _warn:
                warn(f'Failed to converge {self.max_err_x}')
            return False

    def set_initial_parameters(self):
        for s in self.joint_segments:
            s.set_initial_parameters()

    def update_initial_parameters(self):
        self.update_min_time()
        for s in self.joint_segments:
            s.update_initial_parameters(self.t)

    def update_boundary_velocities(self):
        for js in self.joint_segments:
            js.update_boundary_velocities()

    def link_velocities_forward(self):
        for js in self.joint_segments:
            js.link_velocities_forward()

    def link_velocities_backward(self):
        for js in self.joint_segments:
            js.link_velocities_backward()

    def update_min_time(self):
        self.t = self.min_t
        assert self.t > 0

    def recalc(self):
        for js in self.joint_segments:
            js.recalc()

    def reset(self):
        for js in self.joint_segments:
            js.reset()

    def reset_for_append(self):
        for js in self.joint_segments:
            js.reset_for_append()

    def link(self, prior_segment: "Segment"):
        """Link segments together"""

        prior_segment.next_seg = self
        self.prior_seg = prior_segment

        for prior_js, next_js in zip(self.prior_seg.joint_segments, self.joint_segments):
            prior_js.next_js = next_js
            next_js.prior_js = prior_js

    def update_segment_time(self):
        t = [js.t for js in self.joint_segments]
        self.t = max(t)

    @property
    def err_t(self):
        """RMS error between joint segment calculated times and segment time"""
        err = sqrt(sum([(js.t - self.t) ** 2 for js in self.joint_segments]))

        return err / len(self.joint_segments) / self.t

    @property
    def err_x(self):
        """RMS error between joint segment calculated times and segment time"""
        dist = sum(js.x for js in self.joint_segments)
        err = sqrt(sum([(js.err_x) ** 2 for js in self.joint_segments]))

        return err / dist

    @property
    def max_err_x(self):
        return max([js.err_x for js in self])


    @property
    def times(self):
        return [self.t] + [js.t for js in self.joint_segments]

    @property
    def min_t(self):
        return max([js.t for js in self.joint_segments])

    @property
    def params(self):

        from operator import attrgetter
        ag = attrgetter(*JointSegment.params)

        columns = ['seg', 'js', 'seg_t'] + JointSegment.params + ['calc_x', 'sum_x', 'calc_t']

        rows = []
        for j, js in enumerate(self.joint_segments):
            r = (self.n, j, js.segment.t) + ag(js) + (js.calc_x, js.sum_x, js.calc_t)
            d = dict(zip(columns, r) )
            rows.append(d)

        df = pd.DataFrame(rows)

        df['sum_x'] = df.x_a + df.x_c + df.x_d
        df['err_x'] = abs(df.calc_x - df.sum_x).round(4)

        for c in ['v_c', 'v_0_max', 'v_1_max', 'v_0', 'v_1', 'calc_x', 'sum_x']:
            df[c] = df[c].round(0).astype('int')

        for c in ['seg_t', 't', 't_a', 't_d', 't_c']:
            df[c] = df[c].round(5)

        return df

    def save(self):
        from operator import itemgetter
        keys = ['x', 'dir', 't',
                'v_0_max', 'v_1_max', 'v_c_max',
                'x_a', 't_a', 'x_d', 't_d',  'v_0', 'v_c', 'x_c', 't_c', 'v_1' ]

        ag = itemgetter(*keys)
        return [ dict( zip(keys,ag(js.__dict__))) for js in self]

    def load(self, d):

        for js, e in zip(self.joint_segments,d):
            for k, v in e.items():
                setattr(js, k, v)

    @property
    def subsegments(self):
        for js in self.joint_segments:
            yield from js.subsegments

    def __iter__(self):
        return iter(self.joint_segments)

    def __getitem__(self, item):
        return self.joint_segments[item]

    def __str__(self):

        return f"{self.t:>3.4f}|" + ' '.join(str(js) for js in self.joint_segments)+f" {round(self.max_err_x,2)}"


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

    def new_segment(self, joint_distances):
        return Segment(len(self.segments),
                       [JointSegment(j, x=x) for j, x in zip(self.joints, joint_distances)])

    def rmove(self, joint_distances: List[int], update=True):
        """Add a new segment, with joints expressing joint distance
        :type joint_distances: object
        """

        assert len(joint_distances) == len(self.joints)

        s = self.new_segment((joint_distances))

        if len(self.segments) > 0:

            self.segments[-1].next_segment = s
            s.prior_seg = self.segments[-1]
            s.link(self.segments[-1])

        self.segments.append(s)

        if s.prior_seg:
            s.prior_seg.initialize_for_append()

        s.initialize()

        if update:
            self.update(s)

        return s

    def has_discontinuity(self, s1, s2):

        for js1, js2 in zip(s1, s2):
            if js1.v_1 != js2.v_0:
                return True

    def update(self, s=None):

        if s:

            if s.prior_seg:
                s.update_initial_parameters()
                s.prior_seg.update()

            s.link_velocities_backward()
            s.update()

            if s.prior_seg and self.has_discontinuity(s.prior_seg, s):
                s.prior_seg.link_velocities_forward()
                s.prior_seg.update()
                s.update()

        # One more for good measure
        for s in self:
            s.update()


    def _lim_segments(self, limit=4):
        if limit:
            return list(self.segments)[-limit:]
        else:
            return list(self.segments)

    def err_t(self, limit=4):

        l = self._lim_segments(limit)

        et = sum([s.err_t for s in l]) / len(l)

        if et is None:
            return -1

        return et

    def err_x(self, limit=4):
        l = self._lim_segments(limit)
        return sum([js.subseg_error for s in l for js in s])

    def count_needs_update(self):
        return sum(s.needs_update for s in self.segments)

    def validate(self, limit=3):
        """Check that the c error in each of the segments is less than X%"""
        for s in self:
            for js in s:
                if js.x != 0:
                    rel_error = abs(js.err_x / js.x)
                    if js.err_x > limit:
                        raise ValidationError(f"{js.id} Excessive error. err_x={js.err_x} re={rel_error}")

    @property
    def dataframe(self):
        from warnings import warn

        rows = []

        try:
            self.validate()
        except ValidationError as e:
            warn(str(e))

        for ss in self.subsegments:
            rows.append([None, ss.segment_n, ss.joint.n,
                         ss.x, ss.js.v_0_max, ss.v_i, ss.v_f, ss.js.v_1_max, ss.ss, ss.t])

        df = pd.DataFrame(rows, columns=' t seg axis x v0m v_i v_f v1m ss del_t'.split())

        df['t'] = df.groupby('axis').del_t.cumsum()
        df['calc_x'] = (df.v_i + df.v_f) / 2 * df.del_t
        df['err'] = df.x - df.calc_x

        return df

    @property
    def total_error(self):
        """Sum of the segment errors, from the dataframe"""
        return self.dataframe.err.abs().sum()

    @property
    def total_re(self):
        """Sum of the segment errors, from the dataframe, divided by total distance, of all axes"""
        t = self.dataframe
        return t.err.abs().sum()/t.x.sum()

    @property
    def dataframe_stacked(self):
        return self.dataframe.set_index(['t', 'axis']).stack().unstack(-2).unstack(-1)

    @property
    def params(self):
        """
        :return: a dataframe of paramaters
        :rtype:
        """

        frames = []
        for i, s in enumerate(self.segments):
            frames.append(s.params)

        df = pd.concat(frames)

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
