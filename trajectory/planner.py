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

"""
from copy import deepcopy
from itertools import count
from math import sqrt
from typing import List
from collections import deque

## Parameters for simulation/step generation

# NUmber of ticks of the step function per second
TIMEBASE = 1_000_000  # ticks per second

# Shapes
TRAPEZOID = 1
TRIANGLE = 2

N_BIG = 2**32

def binary_search(f, v_min, v_guess, v_max):

    for i in range(20):

        x = f(v_guess)

        if round(x) > 0:
            v_guess, v_min = (v_max + v_guess) / 2, v_guess

        elif round(x) < 0:
            v_guess, v_max = (v_min + v_guess) / 2, v_guess

        else:
            return v_guess

        if abs(v_min-v_max) <1:
            return v_guess

    else:
        return None

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


def accel_x(v0, v1, a):
    """Distance required to accelerate from v0 to v1 at acceleration a"""

    dt = (v1 - v0) / a  # Time to change from v0 to v1 at max acceleration

    return abs(((v0 + v1) / 2) * dt), abs(dt)  # Area of velocity trapezoid

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

    def __init__(self, joint: "Joint",  t: float, v_i: float, v_f: float, x: float, ss: float) -> None:
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


class JointSegment(object):
    """One segment for one joint"""

    x: float
    x_err: float  # Extra step value to add, somewhere

    x_a: float = 0
    x_c: float = 0
    x_d: float = 0

    v_0: float = 0
    v_0_max: float
    v_c: float = 0
    v_1: float = 0
    v_1_max: float  # Don't allow updating v_1

    shape: int = TRAPEZOID

    joint: Joint
    next_js: 'JointSegment'
    prior_js: 'JointSegment'

    def __init__(self, joint: Joint = None, x: float = None, v0: float = 0):

        self.joint = joint
        self.segment = None  # set after adding to a segment

        self.x = abs(x)
        self.s = sign(x)

        self.v_0 = v0 # Initial velocity

        self.x_a = 0 # Acceleration distance
        self.t_a = 0;

        self.v_c = self.joint.v_max
        self.x_c = self.x # cruise distance
        self.t_c = 0

        self.x_d = 0 # Deceleration distance
        self.t_d = 0

        self.v_1 = 0 # Final velocity

        # Extra imposed limits on velocity
        self.v_0_max = self.joint.v_max
        self.v_1_max = self.joint.v_max

        self.next_js = None
        self.prior_js = None

        self.x_err = 0

    def update_v_0(self):
        if not self.prior_js:
            self.v_0 = 0
        else:
            self.v_0 = self.prior_js.v_1

    def update(self):

        seg_t = self.segment.t

        for i in range(10):
            self.x_a, self.t_a = accel_x(self.v_0, self.v_c, self.joint.a_max)

            self.x_c = self.x - self.x_a
            self.t_c = seg_t - self.t_a

            v_c = self.x_c/self.t_c

            if abs(self.v_c - v_c) > 1:
                #print(self.id, self.v_c, v_c, self.x_a, self.t_a)
                self.v_c = v_c
            else:
                break

        self.v_1 = self.v_c


    @property
    def id(self):
        return f"{self.segment.n}/{self.joint.n}"

    @property
    def min_t(self):
        """Minimum transit time based on distance and max velocity"""

        # Time to accelerate to max v
        x_a, t_a = accel_x(self.v_0, self.joint.v_max, self.joint.a_max)

        x_c = self.x - x_a # distance left to run after hitting max v

        return t_a + x_c/self.joint.v_max

    @property
    def max_v_c(self):
        """Max VC achievable given the minimum transit time for the segment"""

        return self.x/self.segment.min_t

    @property
    def min_accel_t(self):
        """Min time to accelerate to max_v_c, given prior segment. """

        if self.prior_js:
            return 1
        else:
            return 0

    @property
    def subsegments(self):
        ss = []
        if round(self.t_a, 7) > 0:
            ss.append(SubSegment(self.joint, round(self.t_a, 7), round(self.v_0, 2), round(self.v_c, 2), round(self.x_a, 0),'a'))

        if round(self.t_c, 7) > 0:
            ss.append(SubSegment(self.joint, round(self.t_c, 7), round(self.v_c, 2), round(self.v_c, 2), round(self.x_c, 0),'c'))

        #if round(self.segment.t_d, 7) > 0:
        #    yield SubSegment(round(self.segment.t_d, 7), round(self.v_c, 2), round(self.v_1, 2), round(self.x_d, 0),
        #                     'd')

        return ss


    def __str__(self):
        from colors import color, bold

        v0 = color(f"{int(self.v_0):<5d}", fg='green')
        xa = bold(f"{int(self.s * self.x_a):>6d}")
        c =  bold(f"{int(round(self.s * self.x_c)):>6d}")
        vc = color(f"{int(self.v_c):<6d}", fg='blue')
        xd = bold(f"{int(self.s * self.x_d):<6d}")
        v1 = color(f"{int(self.v_1):>5d}", fg='red')

        return f"[{v0} {xa}↗{c+'@'+vc}↘{xd} {v1}]"

    def __repr__(self):
        return self.__str__()

    def debug_str(self):
        from colors import color

        a = f"{self.v_0_max:<6.0f}"
        ta = f"{round(self.min_t_a * 1000, 1):^5.1f}"
        c = f"{round(self.min_t_c * 1000, 3):^13.3f}"
        td = f"{round(self.min_t_d * 1000, 1):^5.1f}"
        d = f"{self.v_1_max:>6.0f}"

        return f"[{color(str(a), bg='green')} {ta} : {c} : {td} {color(str(d), bg='red')}]"


class Segment(object):
    """One segment, for all joints"""

    next_seg: "Segment" = None
    prior_seg: "Segment" = None

    def __init__(self, n, joint_segments: List[JointSegment]):

        self.n = n
        self.joint_segments = joint_segments

        for js in self.joint_segments:
            js.segment = self

    def link(self, prior_segment: "Segment"):
        """Link segments together"""

        prior_segment.next_seg  = self
        self.prior_seg = prior_segment

        for prior_js, next_js in zip(self.prior_seg.joint_segments, self.joint_segments):
            prior_js.next_js = next_js
            next_js.prior_js = prior_js

    def update(self):

        for js in self.joint_segments:
            js.update_v_0()

        self.t = self.min_t

        for js in self.joint_segments:
            js.update()


    @property
    def min_t(self):
        return max([ js.min_t for js in self.joint_segments])

    @property
    def subsegments(self):

        return [js.subsegments for js in self.joint_segments]

    def __str__(self):

        return f"{self.t:>3.4f}|" + ' '.join(str(js) for js in self.joint_segments)


class SegmentList(object):

    positions: List[float] # Positions after last movement addition
    segments: deque[Segment]
    all_segments: List[Segment] # All of the segments, unprocessed

    def __init__(self, joints: List[Joint]):

        self.joints = deepcopy(joints)

        self.directions = [0] * len(self.joints)

        for i, j in enumerate(self.joints):
            j.n = i

        self.segments = deque()
        self.positions = [0] * len(self.joints)


    def rmove(self, joint_distances: List[int]):
        """Add a new segment, with joints expressing joint distance
        :type joint_distances: object
        """

        assert len(joint_distances) == len(self.joints)

        next_seg = Segment(len(self.segments), [JointSegment(j, x=x, v0=0) for j, x in zip(self.joints, joint_distances)])

        if len(self.segments) > 0:
            self.segments[-1].next_segment = next_seg
            next_seg.prior_seg = self.segments[-1]
            next_seg.link(self.segments[-1])

        self.segments.append(next_seg)

        return next_seg

    def update(self):
        for s in self.segments:
            s.update()


    def min_joint_times(self):
        """Minimum time for each joint, given distance an max velocity"""

        acum = [0]*len(self.joints)

        for s in self.segments:
            for i, js in enumerate(s.joint_segments):
                acum[i] += js.x/js.joint.v_max

        return acum

    @property
    def dataframe(self):
        import pandas as pd

        rows = []

        tc = 0
        for e in self.subsegments:
            seg_header = []
            for i, ss in enumerate(e):

                rows.append([ss.joint.n, ss.x, ss.v_i, ss.v_f, ss.ss, ss.t])

        df = pd.DataFrame(rows, columns='axis x v_i v_f ss del_t'.split())

        df['t'] = df.groupby('axis').del_t.cumsum()

        return df

    @property
    def dataframe_stacked(self):
        return self.dataframe.set_index(['t', 'axis']).stack().unstack(-2).unstack(-1)

    @property
    def subsegments(self):
        from itertools import chain

        # I'd comment this, but even through I wrote it, I don't understand it.
        return [ list(chain(*e)) for e in zip(*[s.subsegments for s in self.segments]) ]


    def __iter__(self):
        return iter(self.segments)

    def __len__(self):
        return len(self.segments)

    def __str__(self):
        return '\n'.join(str(s) for s in self.segments)

    def debug_str(self):
        return '\n'.join(s.debug_str() for s in self.segments)

__all__ = ['SegmentList', 'Segment', 'Joint', 'JointSegment', 'SegmentError']