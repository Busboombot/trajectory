import typing
from dataclasses import dataclass, asdict, replace
from enum import Enum
from math import sqrt
from typing import List
from copy import copy
from scipy.optimize import minimize_scalar
from operator import attrgetter

from trajectory.exceptions import TrapMathError, ParameterAdjustmentError



@dataclass
class InputParams:
    x: float
    v_0: float
    v_1: float
    v_max: float
    a_max: float


@dataclass
class Params:
    x: float
    t: float = None
    t_a: float = None
    t_c: float = None
    t_d: float = None
    x_a: float = None
    x_c: float = None
    x_d: float = None
    v_0: float = 0
    v_c: float = None
    v_1: float = 0
    v_0_max: float = None
    v_0_min: float = 0
    v_1_max: float = None
    v_1_min: float = 0
    v_max: float = None
    a_max: float = None
    t_min: float = None
    d: int = None  # direction, -1 or 1
    ip: InputParams = None
    flag: str = None
    recalcs: int = 0
    saved = None

    def __post_init__(self):
        self.d=sign(self.x)

    def asdict(self):
        return asdict(self)

    def replace(self, **kwds):
        return replace(self, **kwds)

    @property
    def subsegments(self):
        rd = round #lambda v, n: v
        return (
            (rd(self.t_a, 7), rd(self.v_0, 2), rd(self.v_c, 2), rd(self.x_a, 0),'a'),
            (rd(self.t_c, 7), rd(self.v_c, 2), rd(self.v_c, 2), rd(self.x_c, 0),'c'),
            (rd(self.t_d, 7), rd(self.v_c, 2), rd(self.v_1, 2), rd(self.x_d, 0),'d')
        )

    def store(self):
        self.saved= [round(self.t,4), round(self.v_0), round(self.v_c), round(self.v_1)]

    @property
    def is_changed(self):

        if self.saved is None:
            return False

        for p in zip(self.saved, [round(self.t,4), round(self.v_0), round(self.v_c), round(self.v_1)] ):

            if p[0] != p[1]:
                return True

        return False

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

    ╲▁╱ R Trough. Like a pentagon, but Vc is below either V0 or V1

    ╲▁▁ D Start Decel. Very small x relative to other joints.
          X is so short relative to other joints that we hit v=0
          during accel phase.


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
        x_break = (v_0 + v_1) / 2 * (t_a + t_d)
        return JSClass.CLIFF if x > x_break else JSClass.DECEL
    elif v_0 == v_1:
        return JSClass.CONSTANT

    else:
        return JSClass.UNK


def sign(x):
    if x == 0:
        return 0
    elif x > 0:
        return 1
    else:
        return -1

def same_sign(a, b):
    return int(a) == 0 or int(b) == 0 or sign(a) == sign(b)


def accel_tx(v0, v1, a):
    """Distance and time required to accelerate from v0 to v1 at acceleration a"""

    t = abs((v1 - v0) / a)  # Time to change from v0 to v1 at max acceleration

    if t == 0:
        return 0, 0

    # re-calc a to get the right sign
    a = (v1 - v0) / t

    x = v0 * t + (1 / 2) * a * t ** 2

    return x, t  # Area of velocity trapezoid



def trap_v_c(x, t, a):
    """Height ( v_c) of a trapezoid of area x and base width t"""
    from math import sqrt

    # This is the equation worked out by hand
    x_f = a * t ** 2 / 4  # Area of the whole triangle. Our x is the trap at the bottom
    x_t = x_f - x  # Area of the triangle above our trapezoid

    if x_t < 0:
        # this is a triangle segment, probably, or the time is too short.
        # if x_t is clos to zero, then it's a triangle
        if abs(x_t) < 3:  # That 3 is completely arbitrary
            return x_f  # Triangle
        else:
            # Probably a decel
            raise TrapMathError(f"Probably DECEL: x_t < 0: x={x} t={t} x_f={x_f} x_t={x_t}")

    t_c = sqrt(4 * x_t / a)  # t_c is both the total time for the top triangle, and the cruise time for our trap
    t_a = t_d = (t - t_c) / 2  # acceleration and deccel times
    v_c = a * t_a  # v_c is at the end of the accel time

    # This is based on an equation from SymPy, as a double check
    v_c_ = a * t / 2 - sqrt(a * (a * t ** 2 - 4 * x)) / 2

    assert round(v_c_, 4) == round(v_c, 4), f"Result mismatch: {v_c_}!={v_c}"

    return v_c

def min_time_parameters(x, v_0, v_1, v_max, a_max):
    """Create initial Parameters with estimates for the minimum time for
    a given distance """

    ip = InputParams(x, v_0, v_1, v_max, a_max)
    for i in range(10): # range only need to be bigger than # of possible recalcs
        p = _min_time_parameters(x, v_0, v_1, v_max, a_max, ip)
        if p.flag == 'R': # need recalc:
            (x, v_0, v_1, v_max, a_max, ip) = p.x, p.v_0, p.v_1, p.v_max, p.a_max, p.ip
        else:
            p = p.replace(recalcs=i)
            break

    return p



def min_time_parameters_p(p: Params, inplace = False):
    """construct the minimum time parameters from a Params object. Return
    a new Params by default, but with inplace=True copy into the object provided"""

    if p is None:
        return p

    ag = attrgetter(*'x v_0 v_1 v_max a_max'.split())
    p_ = min_time_parameters(*ag(p))

    if inplace:
        for key, value in p_.asdict().items():
            setattr(p, key, value)
    else:
        return p_

def _min_time_parameters(x, v_0, v_1, v_max, a_max, ip):
    """Find the  lowest time for a segment joint, independent of time constraints
    imposed by other segments, and report any parameters required to get it
    """

    # x values below this value are very likely to require minimization,
    # and fail.
    # danger_zone = (v_max**2)/(2*a_max)
    # danger = x < danger_zone

    if ip is None:
        ip = InputParams(x, v_0, v_1, v_max, a_max)

    base_params = dict(
        v_0_max=v_max,
        v_1_max=v_max,
        v_max = v_max,
        a_max = a_max
    )

    if x <= 0:
        return Params(0, 0, v_0 = 0, v_c = 0, v_1 = 0,
                      t_a=0, x_a=0, t_c=0, x_c=0, t_d=0, x_d=0,
                      t_min = 0,
                      **base_params, ip=ip, flag='Z')

    # Absolute minimum time, because we can't accelerate from v0 t0 v1 any
    # faster
    x_min, t_min = accel_tx(v_0, v_1, a_max)

    if x < x_min:
        # Very small moves, usually at very high speeds, with a big
        # difference between v_0 and v_1
        # First try to change v_1. This will bring it closer to v_0
        a = a_max
        t = (-v_0 + sqrt(2 * a * x + v_0 ** 2)) / a
        v_1 = -v_0 + (2 * x / t)
        assert v_1 > 0, f'v_1 is negative'

        if v_1 > v_max:
            # Now lower v_0
            v_1 = v_max
            v_0 = -v_1 + (2 * x / t)

        x_min, t_min = accel_tx(v_0, v_1, a_max)

        assert round(x, 2) >= round(x_min, 2), f"v_1 reduction failed {x} < {x_min}"

    # compute the minimum time and distance to accelerate to v_max,
    # then back down to v_1
    x_a, t_a = accel_tx(v_0, v_max, a_max)
    x_d, t_d = accel_tx(v_max, v_1, a_max)

    if x < x_a + x_d:
        # The X value is smaller than the distance at which t_c > 0, and probably
        # v_c < v_max. This region has a really complicated function for finding
        # minimum t, with lots of local minimum and global ridges.

        # First, we will try minimizing the error in the area function to find v_c
        # that gives us a trajectory with the given x, but about half the time,
        # the minimization returns results that radically change the x value
        # ... which is very bad.

        def f(v_max):
            x_a, t_a = accel_tx(v_0, v_max, a_max)
            x_d, t_d = accel_tx(v_max, v_1, a_max)

            return abs(x - (x_a + x_d))

        r = minimize_scalar(f, bounds=(0, v_max), method='bounded')
        v_c = round(r.x)

        x_a, t_a = accel_tx(v_0, v_c, a_max)
        x_d, t_d = accel_tx(v_c, v_1, a_max)

        t_c = 0
        x_c = 0

        flag = 'M'  # Minimization

        if round(x_a + x_d) != round(ip.x):
            # This is the really bad case ... the minimization failed.
            # First try setting v_1 to v0, then to zero, then try v_0. If both are
            # zero, it must be solvable. These changes may cause re-calc in
            # adjacent segments, but that should be rare.

            p = Params(x, **base_params, ip=ip, flag='R')
            if v_0 != v_1 and v_1 != 0 and v_0 != 0:
                return p.replace(v_1=v_0, v_1_max=v_0)
            elif v_1 != 0:
                return p.replace(v_1=0, v_1_max=0)
            elif v_0 != 0:
                return p.replace(v_0=0, v_0_max=0)

            raise TrapMathError("Can't adjust parameters")

    else:
        # The acceleration hit the max, so we just need
        # to add cruise time.
        x_c = x - x_a - x_d
        t_c = x_c / v_max
        v_c = v_max

        flag = 'C'  # Calculated

    t = t_a + t_c + t_d
    x = x_a + x_c + x_d

    # Very important not to change distance
    if flag != 'X':  # We already flagged the "X" case, and there should not be any other cases of changing x
        assert round(x) == round(ip.x), f"distance changed: {ip.x}->{x} flag={flag}"

    return Params(x, t, t_a, t_c, t_d, x_a, x_c, x_d, v_0, v_c, v_1,
                  **base_params, t_min = t, ip=ip, flag=flag)

def hex_area(t, v_0, v_c, v_1, a_max):
    """Return the area for a velocity hexagon"""

    x_a, t_a = accel_tx(v_0, v_c, a_max)
    x_d, t_d = accel_tx(v_1, v_c, a_max)
    t_c = t - t_a - t_d
    x_c = v_c * t_c
    x = x_a + x_c + x_d

    return x if x > 0 else 0

def hex_area_p(p: Params):

    ag_ha = attrgetter(*'t v_0 v_c v_1 a_max'.split())

    return hex_area(*ag_ha(p))

def hex_v_c(x, t, v_0, v_1, v_max, a_max):
    """Calculate the v_c parameter of a velocity hexagon, given all other
    parameters """

    args = (x, t, v_0, v_1, v_max, a_max)
    assert not any([e is None for e in args]), f'Something is null: {args}'

    def f(v_c):
        x_ = hex_area(t, v_0, v_c, v_1, a_max)
        assert x_ is not None

        return abs(x - x_)



    # There are a lot of cases where the velocity hexagon will have two
    # solutions for v_c. We should (probably) always perfer the higher one.
    # Using Brent, because I can't get Bounded to start at the top of the interval.
    r = minimize_scalar(f, bracket=(v_max, 0), method='brent')

    # r = minimize_scalar(f, bounds=(0, v_max), method='bounded')
    v_c = round(r.x)

    return v_c

#####
# Group operations
#####

def uncap_end_segment(s: List[Params]):
    """Remove the v=0 requirement on the last segment"""

    for j in s:
        if j:
            j.v_1_max = j.v_1 = j.v_max

def segment_changed(s: List[Params]):
    return any([j.is_changed if j is not None else None for j in s])

def assert_equal_area(a: Params, b: Params):
    """ Changing the area of a segment absolutely cannot be tolerated
    :param a:
    :type a:
    :param b:
    :type b:
    :return:
    :rtype:
    """

    if a is None or b is None:
        return

    assert_consistent_area(a)
    assert_consistent_area(b)

    assert round(a.x) == round(b.x), (a.x, b.x)  # Changing the area is a serious error

    # no, really, it's bad
    ax = hex_area_p(a)
    bx = hex_area_p(b)
    assert round(ax) == round(bx), (ax, bx)

def assert_consistent_area(a: Params):
    """Assert that the internal x is the same as the re-calculated area"""

    if a is None:
        return

    ax = hex_area_p(a)
    assert round(a.x) == round(ax), (a.x, ax)

def save_for_change(l: List[Params]):
    """Save joints for memoization, to discover changes"""

    return [copy(e) if e else None for e in l]

def assert_unchanged_area(l: List[Params], memo:  List[Params]):
    for a, b in zip(l, memo):
        assert_equal_area(a,b)



def update_boundary_velocities(prior: Params, current: Params, nxt: Params):
    """v_1_max may need to be specified if the next segment has a very short
    travel, or zero, because it may not be possible to decelerate during the phase.
    For instance, if next segment has x=0, then v_1 must be 0
    """

    if prior is None:
        # prior == None means that this is the first
        # in the list
        current.v_0_max = current.v_0 = 0
    else:
        if prior.x == 0 or not same_sign(prior.d, current.d):
            current.v_0 = current.v_0_max = 0
        else:
            current.v_0_max = min(prior.v_1_max, current.v_0_max)

        current.v_0 = min(prior.v_1, current.v_0_max)

    if nxt is None:
        # nxt == none means this is the last in the list
        current.v_1 = 0
    else:
        if nxt.x == 0 or not same_sign(nxt.d, current.d):
            current.v_1 = current.v_1_max = 0
        else:
            current.v_1_max = min(nxt.v_0_max, current.v_1_max)

    current.v_1 = max(min(current.v_1, current.v_1_max), current.v_1_min)
    current.v_0 = max(min(current.v_0, current.v_0_max), current.v_0_min)

    if prior and  prior.v_1 != current.v_0:
        return True # signal that v_0 changed so may need to update the prior
    else:
        return False

def update_segment(joints: List[Params]) -> bool:

    if all(e is None for e in joints):
        return False

    for j in joints:
        min_time_parameters_p(j, inplace=True)

    min_time = max([e.t for e in joints])

    for p in joints:
        update_params(p, min_time)
        assert_consistent_area(p)

def update_params(p, t):
    """Update parameters to consider segment time, and set
    values for the ACD phases and Vc"""

    if p is None:
        return None

    if p.x == 0:
        p.t_c = p.t = t
        p.t_a = p.t_d = 0
        p.t_a= p.t_d = p.v_c = 0
        return p

    p.t = t
    p.v_c = hex_v_c(p.x, p.t, p.v_0, p.v_1, p.v_max, p.a_max)

    p.x_a, p.t_a = accel_tx(p.v_0, p.v_c, p.a_max)
    p.x_d, p.t_d = accel_tx(p.v_1, p.v_c, p.a_max)

    if p.x_a + p.x_d > p.x:
        # We're trying to increase the time for the segment, but v_0 and or
        # v_1 are too large

        # Going to need to move the velocities down. Just make them
        # symmetric for now, but it would probably be better to
        # move one down first.
        p.v_0 = p.v_0_max = sqrt(2 * p.a_max * p.x/2)
        p.v_1 = p.v_1_max = sqrt(2 * p.a_max * p.x/2)

        p.v_c = hex_v_c(p.x, p.t, p.v_0, p.v_1, p.v_max, p.a_max)
        p.x_a, p.t_a = accel_tx(p.v_0, p.v_c, p.a_max)
        p.x_d, p.t_d = accel_tx(p.v_1, p.v_c, p.a_max)

        p.flag = 'R' # Signal recalc
    else:
        p.flag = 'C'

    p.t_c = p.t - p.t_a - p.t_d

    assert p.t_c >= 0, (t, str( p))

    p.x_c = p.v_c * p.t_c

    x = p.x_a + p.x_c + p.x_d

    assert round(x) == round(p.x), (x, t, str(p))

    return p