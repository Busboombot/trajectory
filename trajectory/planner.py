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
from dataclasses import dataclass, asdict, replace
from math import sqrt
from typing import List

import pandas as pd

from . import TrapMathError, ConvergenceError

## Parameters for simulation/step generation

# NUmber of ticks of the step function per second


TIMEBASE = 1_000_000  # ticks per second

# Shapes
TRAPEZOID = 1
TRIANGLE = 2

N_BIG = 2 ** 32

MAX_ERR_T = 1e-4


def binary_search(f, v_min, v_guess, v_max):
    for i in range(20):

        x = f(v_guess)

        if round(x) > 0:
            v_guess, v_min = (v_max + v_guess) / 2, v_guess

        elif round(x) < 0:
            v_guess, v_max = (v_min + v_guess) / 2, v_guess

        else:
            return v_guess

        if abs(v_min - v_max) < .05:
            return v_guess

    else:
        return None


def accel_xt(v_i, v_f, a):
    """Distance and time required to accelerate from v0 to v1 at acceleration a"""

    if v_f == v_i:
        return 0, 0

    if v_f < v_i:
        a = -a

    t = (v_f - v_i) / a  # Time to change from v0 to v1 at max acceleration
    x = (v_i + v_f) / 2 * t

    return x, t  # Area of velocity trapezoid


def accel_acd(v_0, v_c, v_1, a):
    """ Same result as running accel_xt, for v_0->v_c and v_c->v_1,
    and adding the x and t values.
    """
    t_ad = (abs(v_c - v_0) + abs(v_c - v_1)) / a
    x_ad = abs((v_0 ** 2 - v_c ** 2) / (2 * a)) + abs((v_1 ** 2 - v_c ** 2) / (2 * a))

    return x_ad, t_ad


def consistantize(b, return_error=False):
    """Recalculate t to make x and v values integers, and everything more consistent
     This operation will maintain the original value for x, but may change t"""

    b.x_a = (b.x_a)
    b.x_d = (b.x_d)
    b.x_c = b.x - (b.x_a + b.x_d)

    b.t_a = abs((b.v_c - b.v_0) / b.joint.a_max)
    b.t_d = abs((b.v_c - b.v_1) / b.joint.a_max)

    if round(b.x_c) == 0:  # get rid of small negatives
        b.x_c = 0

    if b.v_c != 0:
        b.t_c = abs(b.x_c / b.v_c)
        #assert b.t_c >= 0, (b.v_c, b.t_c, b.x_c)
    else:
        b.t_c = 0

    b.t = b.t_a + b.t_c + b.t_d

    # Check error against area calculation
    if return_error:
        x_ad, t_ad = accel_acd(b.v_0, b.v_c, b.v_1, b.joint.a_max)
        t_c = round(b.t - t_ad, 8)  # Avoid very small negatives
        x_c = b.v_c * t_c
        x = x_ad + x_c
        x_e = x - b.x
        return b.x_c, x_e
    else:
        return b.x_c


def sign(x):
    if x == 0:
        return 0
    elif x > 0:
        return 1
    else:
        return -1


def same_sign(a, b):
    return int(a) == 0 or int(b) == 0 or sign(a) == sign(b)


def maxmin(l, v, h):
    return max(min(v, h), l)


def make_area_error_func(b):
    def f(v_c):
        return b.x - b.replace(v_c=v_c).arear

    return f


def set_bv(x, v_0, v_1, a_max):
    """Reduce v_0 and v_1 for small values of x when there is an
    error in x after planning """

    x_a, t_a = accel_xt(v_0, 0, a_max)
    x_d = x - x_a

    if x_d < 0:
        v_0 = sqrt(2 * a_max * x)
        v_1 = 0
    else:
        v_0 = v_0
        v_1 = min(sqrt(2 * a_max * x_d), v_1)

    return (v_0, v_1)



@dataclass
class Joint:
    v_max: float
    a_max: float
    # distances below the small x limit will never hit v_max before needing to decelerate
    small_x: float = None

    def __post_init__(self):
        self.small_x = (self.v_max ** 2) / (self.a_max)

    def new_block(self, x, v_0=None, v_1=None):
        v_0 = v_0 if v_0 is not None else self.v_max
        v_1 = v_1 if v_1 is not None else self.v_max

        return Block(x=x, v_0=v_0, v_1=v_1, joint=self)


@dataclass
class Block:
    x: float = 0
    t: float = 0
    t_a: float = 0
    t_c: float = 0
    t_d: float = 0
    x_a: float = 0
    x_c: float = 0
    x_d: float = 0
    v_0: float = 0
    v_c: float = 0
    v_1: float = 0
    d: int = 0  # direction, -1 or 1

    min_t: float = 0

    v_0_max = None
    v_1_max = None

    joint: Joint = None
    segment: 'Segment' = None
    next: 'Block' = None
    prior: 'Block' = None

    flag: str = None
    recalcs: int = 0

    _param_names = ['x', 't', 'dir', 'v_0_max', 'v_0', 'x_a', 't_a', 'x_c', 't_c', 'v_c_max', 'v_c', 'x_d', 't_d',
                    'v_1',
                    'v_1_max']

    @property
    def id(self):
        return f"{self.segment.n}/{self.joint.n}"

    def __post_init__(self):
        self.d = sign(self.x)
        self.x = abs(self.x)
        self.v_0_max = self.joint.v_max
        self.v_1_max = self.joint.v_max

    def asdict(self):
        return asdict(self)

    def replace(self, **kwds):
        return replace(self, **kwds)

    @property
    def subsegments(self):

        rd = lambda v, n: round(v, n)
        rds = lambda v, n: self.d * round(v, n)

        subsegs = (
            (rd(self.t_a, 7), rds(self.v_0, 2), rds(self.v_c, 2), rds(self.x_a, 0), 'a'),
            (rd(self.t_c, 7), rds(self.v_c, 2), rds(self.v_c, 2), rds(self.x_c, 0), 'c'),
            (rd(self.t_d, 7), rds(self.v_c, 2), rds(self.v_1, 2), rds(self.x_d, 0), 'd')
        )

        return [SubSegment(self.id, self.segment.n, self.joint, self, *ss) for ss in subsegs]

    @property
    def area(self):
        """Calculate the distance x as the area of the profile. Ought to
        always match .x """

        x_ad, t_ad = accel_acd(self.v_0, self.v_c, self.v_1, self.joint.a_max)

        t_c = round(self.t - t_ad, 8)  # Avoid very small negatives

        if t_c < 0:
            raise TrapMathError(f'Negative t_c ({t_c}) ')

        x_c = self.v_c * t_c

        if round(x_c) < 0:
            raise TrapMathError(f'Negative x_c ({x_c}) t_c={t_c}')

        return x_ad + x_c

    @property
    def arear(self):
        """Robust area calculation. May return wrong results for bad parameters, but
        will not throw exceptions"""

        x_ad, t_ad = accel_acd(self.v_0, self.v_c, self.v_1, self.joint.a_max)

        t_c = max(self.t - t_ad, 0)
        x_c = max(self.v_c, 0) * t_c

        return x_ad + x_c

    def init(self):

        a_max = self.joint.a_max
        v_max = self.joint.v_max

        # Limit the boundary values for small moves
        self.v_0, self.v_1 = set_bv(self.x, self.v_0, self.v_1, self.joint.a_max)

        self.v_0 = min(self.v_0, v_max)
        self.v_1 = min(self.v_1, v_max)

        if self.x == 0:
            self.v_c = self.v_0 = self.v_1 = 0
            self.flag = 'Z'

        elif  self.x < self.joint.small_x:
            # The limit is the same one used for set_bv.
            # the equation here is the sympy solution to:
            # t_a = (v_c - v_0) / a
            # t_d = (v_1 - v_c) / -a
            # x_a = ((v_0 + v_c) / 2) * t_a
            # x_d = ((v_c + v_1) / 2) * t_d
            # x_c = x - (x_a + x_d)
            # solve(x_c, v_c)[1]

            #self.v_0 = 0
            #self.v_1 = 0
            self.v_c = (sqrt(4 * a_max * self.x + 2 * self.v_0 ** 2 + 2 * self.v_1 ** 2) / 2)
            #self.v_c = min(self.v_c, v_max)  # formula errs for short segments

            self.flag = 'S'
        else:
            # If there is more area than the triangular profile for these boundary
            # velocities, the v_c must be v_max. In this case, it must also be true that:
            #    x_ad, t_ad = accel_acd(self.v_0, v_max, self.v_1, a_max)
            #    assert self.x > x_ad
            self.v_c = v_max
            self.flag = 'M'

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, a_max)
        self.x_d, self.t_d = accel_xt(self.v_c, self.v_1, a_max)

        assert round(self.x_a + self.x_d) <= self.x

        self.x_c = self.x - (self.x_a + self.x_d)
        self.t_c = self.x_c / self.v_c if self.v_c != 0 else 0
        self.t = self.t_a + self.t_c + self.t_d

        assert round(self.x_c) >= 0, (self.x_c, self.v_c, self.flag, self)
        assert abs(self.area - self.x) < 1, (self.x, self.area, self.t, self.v_0, self.v_1)
        return self

    def plan_ramp(self, t):
        """Calculate a ramp profile for a fixed time"""

        self.t = t
        a_max = self.joint.a_max

        if self.x == 0 or self.t == 0:
            self.set_zero()
            self.t_c = self.t = t
            self.flag = 'Z'  # 5% of test cases
            return self

        # Run the binary search anyway, just in case.
        def err(v_c):
            x_a, t_a = accel_xt(self.v_0, v_c, self.joint.a_max)
            t_c = self.t - t_a

            x_c = v_c * t_c
            x_err = self.x - (x_a + x_c)

            return x_err

        try:
            guess = a_max * t + self.v_0 - sqrt(a_max * (a_max * t ** 2 + 2 * t * self.v_0 - 2 * self.x))
        except ValueError:
            guess = self.x/self.t

        self.v_1 = self.v_c = min(binary_search(err, 0, guess, self.joint.v_max), self.joint.v_max)

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, self.joint.a_max)

        self.t_d, self.x_d = 0, 0

        consistantize(self)

        self.flag = "PR"

        assert self.v_0 <= self.joint.v_max
        assert self.v_c <= self.joint.v_max, self.v_c
        assert self.v_1 <= self.joint.v_max

        return self

    def plan(self, t):

        self.t = t

        a_max = self.joint.a_max

        if self.x == 0 or self.t == 0:
            self.set_zero()
            self.t_c = self.t = t
            self.flag = 'Z'  # 5% of test cases
            return self

        # Find v_c with a binary search, then patch it up if the
        # selection changes the segment time.
        self.v_c = min(binary_search(make_area_error_func(self), 0,
                                 self.x / self.t, self.joint.v_max), self.joint.v_max)
        assert self.v_c <= self.joint.v_max
        self.flag = 'O'

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, a_max)
        self.x_d, self.t_d = accel_xt(self.v_c, self.v_1, a_max)

        # consistantize will make all of the values consistent with each other,
        # and if anything is wrong, it will show up in a negative x_c

        x_c = consistantize(self)

        def psfrt(t, flag):  # plan, set flag and return
            b = self.plan(t)
            b.flag = 'NT'
            return b

        x_err = abs(round(self.area) - self.x)

        if round(x_c) < 0 or (self.x > 25 and x_err > 1):
            # We've probably got a really small move, and we've already tried
            # reducing boundary velocities, so get more aggressive. If that
            # doesn't work, allow the time to expand.

            new_t = max(self.t_a + self.t_d, self.t)

            if self.v_1 > 0:
                self.v_1 = 0
                return psfrt(t, 'V1Z')

            elif self.v_0 > 0:
                self.v_0 = 0
                return psfrt(t, 'V0Z')

            elif abs(new_t - self.t) > 0.0005:
                print(round(new_t-self.t,4), new_t, self.t)
                return psfrt(new_t, 'NT')
            else:
                raise ConvergenceError(
                    'Unsolvable profile, incorrect area: '
                    f'x_c={x_c} x_err={x_err} x={self.x}, x_ad={self.x_a + self.x_d} t={self.t} t_ad={self.t_a + self.t_d}')

        if round(t, 3) != round(self.t, 3):
            # This block is still too long.
            if self.v_1 != 0:
                self.v_1 = 0
                return psfrt(t, 'HT1')
            elif round(self.v_0) > 0:
                self.v_0 = self.v_0 / 2
                return psfrt(t, 'HT1')
            else:
                pass
                # raise ConvergenceError(
                #    'Unsolvable profile, wrong time: '
                #    f't={self.t}, commanded t ={t} t_ad={self.t_a + self.t_d} v_c={self.v_c}')

        assert round(self.x_a + self.x_d) <= self.x
        #assert self.t_a + self.t_d <= self.t
        assert self.v_c >= 0, (self.v_c, self.flag)
        assert abs(self.area - self.x) < 2, (self.area, self)
        assert self.t > 0
        assert self.v_0 <= self.joint.v_max
        assert self.v_c <= self.joint.v_max
        assert self.v_1 <= self.joint.v_max

        return self

    def set_zero(self):

        self.flag = ''  # 5% of test cases
        self.x_a = self.x_d = self.x_c = 0
        self.t_a = self.t_d = self.t_c = 0
        self.v_0 = self.v_c = self.v_1 = 0

    def str(self):
        from colors import color, bold

        v0 = color(f"{int(self.v_0):<5d}", fg='green')
        xa = bold(f"{int(self.d * self.x_a):>6d}")
        c = bold(f"{int(round(self.d * self.x_c)):>6d}")
        vc = color(f"{int(self.v_c):<6d}", fg='blue')
        xd = bold(f"{int(self.d * self.x_d):<6d}")
        v1 = color(f"{int(self.v_1):>5d}", fg='red')

        return f"[{v0} {xa}↗{c + '@' + vc}↘{xd} {v1}]"


    def _repr_pretty_(self, p, cycle):
        p.text(self.str() if not cycle else '...')

    @property
    def debug(self):

        from colors import color, bold

        v0_max = color(f"{int(self.v_0_max):<5d}", fg='green')
        v0 = color(f"{int(self.v_0):<5d}", fg='green')
        xa = bold(f"{int(self.d * self.x_a):>6d}")
        ta = bold(f"{int(self.d * self.t_a):<1.5f}")
        c = bold(f"{int(round(self.d * self.x_c))}")
        vc = color(f"{int(self.v_c)}", fg='blue')
        tc = bold(f"{int(self.d * self.t_c):<1.5f}")
        xd = bold(f"{int(self.d * self.x_d):>6d}")
        td = bold(f"{int(self.d * self.t_d):<1.5f}")
        v1 = color(f"{int(self.v_1):>5d}", fg='red')
        v1_max = color(f"{int(self.v_1_max):>5d}", fg='red')

        return f"[{v0_max}|{v0} {xa}%{ta} ↗ {c}@{vc}%{tc} ↘ {xd}%{td} {v1}|{v1_max}]"


class Segment(object):
    """One segment, for all joints"""

    n: int = None
    next: "Segment" = None
    prior: "Segment" = None
    blocks: "List[Block]" = None
    t: float = 0

    n_updates: int = 0

    updates: List[str]

    def __init__(self, n, blocks: List[Block], prior=None):

        self.n = n
        self.blocks = blocks

        for b in self.blocks:
            b.segment = self

        if prior:
            Segment.link(prior, self)

    @classmethod
    def link(cls, prior: "Segment", current: "Segment"):

        prior.next = current
        current.prior = prior

        for p, c in zip(prior.blocks, current.blocks):
            p.next = c
            c.prior = p

    def init(self):
        """Re-initialize all of the blocks"""
        for j in self:
            j.init()

        return self

    def plan(self, style='T', update_bv=False):

        if update_bv:
            self.update_1_boundary()
            self.update_0_boundary()

        for i in range(6):
            mt = max(self.times)

            for j in self:
                if style == 'R':
                    j.plan_ramp(mt)
                else:
                    j.plan(mt)

            if self.times_e_rms < 0.001:
                break
        else:
            if self.times_e_rms < 0.003:
                from warnings import warn
                warn("Failed to converge sid={self.n} with err > 0.001, but < 0.003")
            else:
                raise ConvergenceError(
                    f'Too many planning updates for sid={self.n}. times={self.times}, rms={self.times_e_rms}')
                pass

    def update_0_boundary(self):
        if not self.prior:
            return

        for p, c in zip(self.prior.blocks, self.blocks):
            c.v_0 = min(c.v_0, p.v_1)

    def update_1_boundary(self):
        """ Give this block v_1s that match the v_0 of the next block
        :return:
        :rtype:
        """
        if not self.next:
            return

        for c, n in zip(self.blocks, self.next.blocks):
            c.v_1 = min(c.v_1, n.v_0)

    @property
    def final_velocities(self):
        return [j.v_1 for j in self.blocks]

    @property
    def boundary_velocities_0(self):
        """Get the boundary velocities for the v_0's of this segment
        and the v_1's of the prior
        :return:
        :rtype:
        """
        if self.prior:
            p = self.prior.blocks
        else:
            p = [None for _ in self.blocks]

        return [(p.v_1 if p is not None else None,
                 c.v_0 if c is not None else None)
                for p, c in zip(p, self.blocks)]

    @property
    def boundary_velocities_1(self):
        """Get the boundary velocities for the v_1's of this segment
        and the v_1's of the prior
        :return:
        :rtype:
        """
        if self.next:
            n = self.next.blocks
        else:
            n = [None for _ in self.blocks]

        return [(c.v_1 if c is not None else None,
                 n.v_0 if n is not None else None)
                for n, c in zip(n, self.blocks)]

    @property
    def has_0_discontinuity(self):
        return not all(a is None or a == b for a, b in self.boundary_velocities_0)

    def fix_boundary_bumps(self):
        """Check if the boundary velocity is larger than either of the v_c"""
        fixes = 0

        if not self.prior:
            return 0

        for p, c in zip(self.prior.blocks, self.blocks):

            v_c_m = (p.v_c + c.v_c) / 2

            if c.v_0 > v_c_m:
                c.v_0 = v_c_m
                p.v_1 = v_c_m
                fixes += 1

        return fixes

    @property
    def has_1_discontinuity(self):
        return not all(b is None or a == b for a, b in self.boundary_velocities_1)

    @property
    def params(self):
        return [j for j in self.blocks]

    @property
    def times(self):
        return [round(js.t, 6) for js in self.blocks]

    @property
    def times_e_rms(self):
        """Compute the RMS difference of the times from the mean time"""
        import numpy as np

        # Ignore really small axis for this calculation, because
        # they are very troublesome, and errors in them wont make a
        # big difference
        times = [round(b.t, 6) for b in self.blocks if b.x > 100]

        if not times:
            times = [0]

        s = sum(times)
        m = s / len(times)
        return np.sqrt(np.sum([(t - m) ** 2 for t in times]))

    @property
    def min_time(self):
        return min(self.times)

    @property
    def params_df(self):

        from operator import attrgetter
        ag = attrgetter(*Block._param_names)

        columns = ['seg', 'js', 'seg_t'] + Block._param_names + ['calc_x', 'sum_x', 'calc_t']

        rows = []
        for j, js in enumerate(self.blocks):
            r = (self.n, j, js.segment.t) + ag(js) + (js.calc_x, js.sum_x, js.calc_t)
            assert len(columns) == len(r)
            d = dict(zip(columns, r))
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
                'x_a', 't_a', 'x_d', 't_d', 'v_0', 'v_c', 'x_c', 't_c', 'v_1']

        ag = itemgetter(*keys)
        return [dict(zip(keys, ag(js.__dict__))) for js in self]

    def load(self, d):

        for js, e in zip(self.blocks, d):
            for k, v in e.items():
                setattr(js, k, v)

    @property
    def subsegments(self):
        for b in self.blocks:
            yield from b.subsegments

    def __iter__(self):
        return iter(self.blocks)

    def __getitem__(self, item):
        return self.blocks[item]

    def __str__(self):

        return f"{self.t:>3.4f}|" + ' '.join(js.str() for js in self.blocks)

    def _repr_pretty_(self, p, cycle):
        p.text(str(self) if not cycle else '...')

    def plot(self, ax=None):
        from .plot import plot_params
        plot_params(self.blocks, ax=ax)


@dataclass
class SubSegment:
    """A sub segment is a portion of the trajectory with a constant
    acceleration,  one of the aceleration, constant (cruise) or decleration portions. """

    id: str
    segment_n: int
    joint: "Joint"
    js: "Block"
    t: float
    v_i: float
    v_f: float
    x: float
    ss: float
    direction: int = 1

    def __post_init__(self):
        assert self.v_i == 0 or self.v_f == 0 or sign(self.v_i) == sign(self.v_f), \
            f"Inconsistent directions {self.v_i} {self.v_f} for {self.id}{self.ss} "

    def set_direction(self, sign) -> None:
        self.direction = sign
        self.v_i *= sign
        self.v_f *= sign
        self.x *= sign

    def __repr__(self):
        return f"<{self.joint.n} {self.ss} {self.t:2.5f} {int(self.x):5d} {self.v_i:5.0f}->{self.v_f:5.0f}> "

    @property
    def row(self):
        return [self.segment_n, self.joint.n, self.x,
                self.direction * self.v_i, self.direction * self.v_f,
                self.ss, self.t]


class SegmentList(object):
    positions: List[float]  # Positions after last movement addition
    segments: deque[Segment]
    all_segments: List[Segment]  # All of the segments, unprocessed

    def __init__(self, joints: List[Joint]):

        self.joints = [Joint(j.v_max, j.a_max, i) for i, j in enumerate(joints)]

        self.directions = [0] * len(self.joints)

        for i, j in enumerate(self.joints):
            j.n = i

        self.segments = deque()
        self.positions = [0] * len(self.joints)

    def new_segment(self, x: List[int], prior_velocities):

        blocks = [j.new_block(x=x, v_0=v_0, v_1=0) for j, x, v_0 in zip(self.joints, x, prior_velocities)]

        prior = self.segments[-1] if len(self.segments) > 0 else None

        s = Segment(len(self.segments), blocks, prior)

        return s

    def rmove(self, x: List[int], update=True, planner=None):
        """Add a new segment, with joints expressing joint distance
        :type x: object
        """

        self.uncap_tail()

        assert len(x) == len(self.joints)

        if len(self.segments) > 0:
            prior = self.segments[-1]
            prior_velocities = []

            for p, c_x in zip(prior.blocks, x):
                if not same_sign(p.d, sign(c_x)):
                    prior_velocities.append(0)  # changed direction, so must be 0
                else:
                    prior_velocities.append(p.v_1)

        else:
            prior_velocities = [0] * len(x)

        s = self.new_segment(x, prior_velocities)

        if len(self.segments) > 0:
            Segment.link(self.segments[-1], s)

        self.segments.append(s)

        s.init()

        if planner:
            planner(s)
        else:
            self.planr(s)

        return s

    def planr(self, s):
        """Replan recursively"""
        if s is None:
            return

        s.plan()

        if s.has_0_discontinuity:
            s.prior.update_1_boundary()
            self.planr(s.prior)
            s.update_0_boundary()
            s.plan()

    def uncap_tail(self):
        if len(self.segments) == 0:
            return

        for b in self.segments[-1]:
            b.v_1 = b.joint.v_max

        mt = self.segments[-1].min_time
        for b in self.segments[-1]:
            b.init().plan_ramp(mt)  # Replan it as a ramp

    @property
    def dataframe(self):

        rows = []

        for ss in self.subsegments:
            rows.append([None, *ss.row])

        df = pd.DataFrame(rows, columns=' t seg axis x v_i v_f ss del_t'.split())

        df['t'] = df.groupby('axis').del_t.cumsum()
        df['calc_x'] = (df.v_i + df.v_f) / 2 * df.del_t
        df['err'] = df.x - df.calc_x

        return df

    @property
    def discontinuities(self):
        """Yield segment pairs with velocity discontinuities"""
        for c, n in self.joint_pairs:
            if round(c.v_1, 2) != round(n.v_0, 2):
                yield c, n

    def has_discontinuity(self, s1, s2):

        for js1, js2 in zip(s1, s2):
            if js1.v_1 != js2.v_0:
                return True

    @property
    def dataframe_stacked(self):
        return self.dataframe.set_index(['t', 'axis']).stack().unstack(-2).unstack(-1)

    @property
    def params_df(self):
        """
        :return: a dataframe of paramaters
        :rtype:
        """

        frames = []
        for i, s in enumerate(self.segments):
            frames.append(s.params_df)

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
            return self.segments[s].blocks[js]
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

    def plot(self, ax=None, axis=None):
        from .plot import plot_trajectory

        df = self.dataframe

        if axis is not None:
            df = df[df.axis == axis]

        plot_trajectory(df, ax=ax)

