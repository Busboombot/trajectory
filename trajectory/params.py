from dataclasses import dataclass, asdict, replace
from enum import Enum
from math import sqrt
from operator import attrgetter

from scipy.optimize import minimize_scalar

from .exceptions import *


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




def halve_boundary_velocities(p):
    """Reduce v_1 by halves, then v_0, finally set them to 0 """

    if  p.v_1 > 1:  # Allow 6 divisions
        p.v_1 = p.v_1/2
    elif p.v_0 > 1:
        p.v_0 = p.v_0/2
    elif p.v_0 > 0 and p.v_1 > 0:
        p.v_0 = 0
        p.v_1 = 0
    else:
        raise TrapMathError("Can't reduce velocities any more")



def make_area_error_func(x, t, v_0, v_1, v_max, a_max, collar=True):
    # Compute an error for outside the boundaries, which has a slope back
    # to the valid range

    def f(v_c):
        from .params import accel_acd

        if v_c < 0:
            return v_max

        x_ad, t_ad = accel_acd(v_0, v_c, v_1, a_max)

        t_c = t - t_ad
        x_c = v_c * t_c

        if round(x_c) < 0:
            return v_max

        return abs(x - (x_c + x_ad))

    return f


def update_boundary_velocities(p: "Block", c: "Block", n: "Block"):
    """v_1_max may need to be specified if the next segment has a very short
    travel, or zero, because it may not be possible to decelerate during the phase.
    For instance, if next segment has x=0, then v_1 must be 0
    """

    # Note that routine only links velocities backwards,
    # that is c.v_0 <- p.v_1. Let the segment determine it's
    # v_1, subject to p.v_0_max, and then the next segment
    # will link backwards to this one.

    if p is None:
        # prior == None means that this is the first
        # in the list
        c.v_0_max = c.v_0 = 0
    else:
        if c.x == 0 or p.x == 0 or not same_sign(p.d, c.d):
            c.v_0 = c.v_0_max = 0
        else:
            c.v_0_max = min(p.v_1_max, c.v_0_max)

        c.v_0 = min(p.v_1, c.v_0_max)

    if n is None:
        # nxt == none means this is the last in the list
        c.v_1_max = c.v_1 = 0
    else:
        if c.x == 0 or n.x == 0 or not same_sign(n.d, c.d):
            c.v_1 = c.v_1_max = 0
        else:
            c.v_1_max = min(n.v_0_max, c.v_1_max)

    c.v_1 = maxmin(c.v_1_min, c.v_1, c.v_1_max)
    c.v_0 = maxmin(c.v_0_min, c.v_0, c.v_0_max)

    assert n is not None or (c.v_1 == 0), 'V_1 should be 0 on last segment. '

    if p and p.v_1 != c.v_0:
        return True  # signal that v_0 changed so may need to update the prior
    else:
        return False

@dataclass
class Joint:
    v_max: float
    a_max: float

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
    v_0_max: float = None
    v_0_min: float = 0
    v_1_max: float = None
    v_1_min: float = 0
    t_min: float = 0
    d: int = 0  # direction, -1 or 1
    joint: Joint = None
    flag: str = None
    recalcs: int = 0
    jsclass: "JSClass" = None

    segment = None
    next: "Block" = None
    prior: "Block" = None

    def __post_init__(self):
        self.d = sign(self.x)

    def asdict(self):
        return asdict(self)

    def replace(self, **kwds):
        return replace(self, **kwds)

    @property
    def subsegments(self):
        rd = round  # lambda v, n: v
        return (
            (rd(self.t_a, 7), rd(self.v_0, 2), rd(self.v_c, 2), rd(self.x_a, 0), 'a'),
            (rd(self.t_c, 7), rd(self.v_c, 2), rd(self.v_c, 2), rd(self.x_c, 0), 'c'),
            (rd(self.t_d, 7), rd(self.v_c, 2), rd(self.v_1, 2), rd(self.x_d, 0), 'd')
        )

    @property
    def area(self):
        """Calculate the distance x as the area of the profile. Ought to
        always match .x """

        x_ad, t_ad = accel_acd(self.v_0, self.v_c, self.v_1, self.joint.a_max)

        t_c = self.t - t_ad
        x_c = self.v_c * t_c

        if round(x_c) < 0:
            raise TrapMathError(f'Negative x_c' + str((x_c, x_ad,
                                                       (self.x, self.v_0, self.v_c, self.v_1, self.joint.a_max))))

        return x_ad + x_c

    def init(self):

        a_max = self.joint.a_max
        v_max = self.joint.v_max

        if self.x == 0:
            self.t = self.t_a = self.t_c = self.t_d = 0
            self.x_a = self.x_c = self.x_d = 0
            self.v_0 = self.v_1 = 0

        elif self.x < (480_000 / a_max):
            # boundary velocity limit. Below this limit, the reduction algorithm fails,
            # so just drive the velocities to zero. The 480K value is empirical; I don't
            # know what causes it.
            self.v_0 = 0
            self.v_1 = 0

        elif self.x <= (v_max ** 2) / (2 * a_max):
            # Small movements are really difficult to process without errors, so
            # lets try to push them down to be close to constant speed
            dt = v_max / a_max
            self.v_0 = min(self.x / dt, self.v_0)
            self.v_1 = min(self.x / dt, self.v_1)

        x_a, t_a = accel_xt(self.v_0, v_max, a_max)
        x_d, t_d = accel_xt(v_max, self.v_1, a_max)

        if self.x == 0:
            v_c = 0
            flag = 'Z'
            x_a = x_d = t_a = t_c = t_d = 0

        elif self.x >= x_a + x_d:
            # If there is more area than the triangular profile for these boundary
            # velocities, the v_c must be v_max
            # This is the complimentary case to the next one; if one is true the
            # other is not, so one could be an else, and the 'O' case never runs.
            v_c = v_max
            flag = 'M'

        elif self.x < (v_max ** 2) / (a_max):
            # The round() is important here. Without it, we get math domain errors
            # elsewhere.
            # t_a = (v_c - v_0) / a
            # t_d = (v_1 - v_c) / -a
            # x_a = ((v_0 + v_c) / 2) * t_a
            # x_d = ((v_c + v_1) / 2) * t_d
            # x_c = x - (x_a + x_d)
            # solve(x_c, v_c)[1]

            v_c = round(sqrt(4 * a_max * self.x + 2 * self.v_0 ** 2 + 2 * self.v_1 ** 2) / 2)
            flag = 'S'

        else:
            assert False, "This code should not execute"

            # What's left are hex profiles in cases where 0 < v_c < v_max.
            # These are triangular profiles, so x_c == 0 and x == (t_a+t_d)
            def f(v_c):
                x_ad, t_ad = accel_acd(self.v_0, v_c, self.v_1, a_max)
                return abs(self.x - x_ad)

            # Quick scan of the space to set the initial bracket.
            mv = min((f(v_c), v_c) for v_c in range(0, 5000, 50))[1]

            r = minimize_scalar(f, bracket=(mv - 10, mv + 10))

            v_c = round(r.x)

            flag = 'O'

        x_a, t_a = accel_xt(self.v_0, v_c, a_max)
        x_d, t_d = accel_xt(v_c, self.v_1, a_max)

        x_c = self.x - (x_a + x_d)
        t_c = x_c / v_c if v_c != 0 else 0

        assert round(x_c) >= 0, (x_c, v_c, flag, self)

        self.t = t_a + t_c + t_d
        self.v_c = v_c
        self.x_a = x_a
        self.x_c = x_c
        self.x_d = x_d
        self.t_a = t_a
        self.t_c = t_c
        self.t_d = t_d
        self.flag = flag

        return self

    def plan(self, t):

        ag = attrgetter(*'x t v_0 v_1'.split())
        assert not any([e is None for e in ag(self)]), f'Something is null: {ag(b)}'

        self.t = t

        a_max = self.joint.a_max
        v_h = max(self.v_0, self.v_1)

        ##
        ## First try to calculate the values
        ##

        if self.x == 0 or self.t == 0:
            self.v_c = 0
            self.flag = 'Z' # 5% of test cases

        elif self.v_0 == 0 and self.v_1 == 0:
            # should always be the lower root
            # Compare to equation 3.6 of Biagiotti
            self.v_c = a_max * self.t / 2 - sqrt(a_max * (a_max * self.t ** 2 - 4 * self.x)) / 2
            self.flag = 'T' # .16% of test cases

        elif self.x >= v_h * self.t and self.v_0 == self.v_1:
            # subtract off v_0, v_1 and the corresponding x and run again.
            x_ = self.v_0 * self.t

            # This is the same case as v_0 == v_1 == 0, but after subtracting off
            # the area of the rectangular base, from v =0 to v=v_0. Then we add back in
            # v_0.
            self.v_c = self.v_0 + a_max * self.t / 2 - sqrt(a_max * (a_max * self.t ** 2 - 4 * (self.x - x_))) / 2
            self.flag = 'T0' # about 12%% of cases
        else:
            ## Well, that didn't work,
            ## so let try minimization

            f = make_area_error_func(self.x, self.t, self.v_0, self.v_1, self.joint.v_max, self.joint.a_max)

            mv = self.x / self.t  # this usually performs just as well as the scan, and is cheaper.

            try:
                # The bounded method behaves badly some times, getting the wrong
                # minimum. In thse cases, the area call will fail, and we run the
                # other minimize call, which performs well in these cases.
                r = minimize_scalar(f, method='bounded',
                                    bracket=(mv - 10, mv + 10), bounds=(0, self.joint.v_max))

                a = self.replace(v_c=r.x).area  # Just to throw an exception
                if round(a) != self.x:
                    raise TrapMathError()
                self.flag = 'O1' # 78% of test cases
            except:
                r = minimize_scalar(f, bracket=(mv - 10, mv + 10))
                self.flag = 'O2' # About 1.5% of cases, or less

            self.v_c = r.x

            # Mop up errors. These cases occur because moving v_c also
            # increases the time for the block, so to we have to adjust the time;
            # this makes the blocks longer than they were commanded to be, and that
            # will require another update later.
            x_ad, t_ad = accel_acd(self.v_0, self.v_c, self.v_1, self.joint.a_max)
            if t_ad > self.t:
                self.recalcs += 1
                self.plan(t_ad)
                self.flag = 'RT' # .4% of test cases
            if round(self.area) != self.x:
                halve_boundary_velocities(self)
                self.init()
                x_ad, t_ad = accel_acd(self.v_0, self.v_c, self.v_1, self.joint.a_max)
                self.recalcs += 1
                self.plan(max(self.t, t_ad))
                self.flag = 'RV' # 1.5% of test cases

        assert self.v_c >= 0, (self.v_c, self.flag)
        assert round(self.area) == self.x, (self.area, self.x)

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, a_max)
        self.x_d, self.t_d = accel_xt(self.v_c, self.v_1, a_max)

        self.x_c = self.x-(self.x_a+self.x_d)
        self.t_c = self.t-(self.t_a + self.t_d)

        return self

class JSClass(Enum):
    """

    === Convex

    ╱▔▔ U Acel. V0 fixed, V1 is maxed
          Should be a typical case.

    ╱╲  A Triangle. Fixed V0, V1, too short for vc
          x is non zero but is too short got velocity to reach v_max

    ╱▔╲ T Trapezoid.  V0 & V1 == 0, max vc
          Normal case where v0 and v1 are both below Vc, and the segment isn't short.

    ╱▔╲ P Pentagon.  Fixed V0, V1, but not zero, max vc
          Normal case where v0 and v1 are both below Vc, and the segment isn't short.

    ▔▔╲ L Clif. Next segment is very small x.
          Next segment is very small x or zero



    === Flat

    ▔▔▔ C Constant V. Prior and next have high V0, V1
          V0 and V1 are both at vmax, so v_c should be v_max

    ▁▁▁ Z Zero distance
          Zero x means that there must bt v=0 through whole segment

    === Concave

    ╲▁╱ V Trough. Like a pentagon, but Vc is below either V0 or V1

    ╲▁▁ R Ramp down. Very small x relative to other joints.
          X is so short relative to other joints that we hit v=0
          during accel phase.


    """
    ACEL = 'U'
    TRIANGLE = 'A'
    PENTAGON = 'P'
    HEXAGON = 'H'
    AT = 'AT'  # A or T, pinned to Zero at start and finish
    TRAPZEZOID = 'T'
    TROUGH = 'V'
    CONSTANT = 'C'
    RAMP = 'R'
    DECEL = 'D'
    CLIFF = 'L'
    DL = 'DL'  # D or L
    ZERO = 'Z'
    XXXXX = 'X'
    UNK = '?'


kind_icon_map = {

    JSClass.TRIANGLE: '╱╲',
    JSClass.TRAPZEZOID: '╱-╲',
    JSClass.PENTAGON: '╱▔╲',
    JSClass.HEXAGON: '╱▔╲',
    JSClass.TROUGH: '╲▁╱',
    JSClass.CONSTANT: '▔▔▔ ',
    JSClass.ZERO: '▁▁▁',
    JSClass.ACEL: '╱▔▔',
    JSClass.DECEL: '╲',
    JSClass.RAMP: '╲▁▁',
    JSClass.CLIFF: '▔▔╲',
    JSClass.XXXXX: '▁▁╱',
    JSClass.UNK: '?'
}


def classify(p):
    """Assign a class type to a planner block"""
    from operator import attrgetter
    ag = attrgetter(*'x v_0 v_c t_c v_1 v_max a_max'.split())

    x, v_0, v_c, t_c, v_1, v_max, a_max = ag(p)

    if x == 0:
        return JSClass.ZERO
    elif v_c == v_1 and v_1 == v_0:
        return JSClass.CONSTANT

    elif v_0 == 0 and v_1 == 0:
        if t_c == 0:
            return JSClass.TRIANGLE
        else:
            return JSClass.TRAPZEZOID

    elif v_0 != 0 and v_1 == 0 and v_c == 0:
        return JSClass.DECEL
    elif v_0 == v_c and v_1 > v_0:
        return JSClass.RAMP  # Probably should have a new one, it's reversed from other ram
    elif v_0 == v_c and v_1 < v_0:
        return JSClass.CLIFF
    elif v_1 == v_c and v_0 > v_1:
        return JSClass.RAMP
    elif v_1 == v_c and v_0 < v_1:
        return JSClass.ACEL

    elif (v_0 != 0 or v_1 != 0) and v_0 < v_max and v_1 < v_max:
        if t_c == 0:
            return JSClass.PENTAGON
        elif v_c > max(v_0, v_1):
            return JSClass.HEXAGON
        elif v_c < min(v_0, v_1):
            return JSClass.TROUGH
        else:
            return JSClass.UNK
    elif v_c < v_0 and v_c < v_1:
        return JSClass.TROUGH


    else:
        return JSClass.UNK
