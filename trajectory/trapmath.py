
from math import sqrt
from collections import namedtuple
from enum import Enum

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

    return sqrt(2 * x * a)


def v_c_max(x, v_0, v_1, a):
    """Velocity after accelerating for x/2 from the start and end,
    the highest velocity for which t_c is positive.
    This is the maximum achievable velocity before deceleration must begin. """
    ta = t_accel_x(x / 2, v_0, a)
    va = v_0 + a * ta

    td = t_accel_x(x / 2, v_1, a)

    vd = v_1 + a * td

    # The averaging produces some error, but it is small for reasonable x,
    # larger than 100, and with reasonable a_max

    return (va + vd) / 2


def simplify_pentagon(x, t, v_0, v_c, v_1, a):
    """ Divide a complex pentagon into three parts:
        * The top, all of the area above the velocity max(v0, v1), a simple trapezoid
        * A side, the trapezoid area on the v_0 or v_1 side, whichever is smaller
        * The rectangular area below the top.

        And adjustment to v_c required to reduce error will happen in the trapezoid area,
        which is a lot easier to calculate"""

    v_max = max(v_1, v_0)
    v_min = min(v_1, v_0)

    # side trapezoid (st_) area
    st_t = (v_max - v_min) / a
    st_x = (v_max + v_min) / 2 * st_t

    # Subtract off the st and base rectangle  (br_)to get the top trapezoid (tt_)
    br_t = t - st_t
    br_x = br_t * v_max

    tt_x = x - st_x - br_x

    # v_c = 2*br_x/(t+br_t)

    return  br_t, tt_x, v_max


def trap_v_c(x, t, a):
    """Height ( v_c) of a trapezoid of area x and base width t"""
    from math import sqrt

    v = a * t / 2 - sqrt(a * (a * t ** 2 - 4 * x)) / 2
    return v


def penta_area(t, v_0, v_c, v_1, a):
    """Return the are for a velocity pentagon"""

    # Full worked out equation
    x_a, t_a = accel_tx(v_0, v_c, a)
    x_d, t_d = accel_tx(v_1, v_c, a)
    t_c = t - t_a - t_d
    x_c = v_c * t_c
    x = x_a + x_c + x_d

    # Simplification from SymPy, but it returns negative values when v_c is zero,
    # and wrong results if v_c <
    # x = (2*v_c*(a*t + v_0 + v_1 - 2*v_c) - (v_0 - v_c)*(v_0 + v_c) - (v_1 - v_c)*(v_1 + v_c))/(2*a)

    # The x equation will return negative values when v_1 or _v0 is larger than v_c
    # x = abs(x)

    return x if x > 0 else 0


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


def penta_v_c(x, t, v_0, v_1, a):
    """ Find the v_c for a pentagon, through decompositions

    /--------------------
   /                     \
  /          T            \
 / ________________________\
|        B2               |S\
|_________________________|__\
|        B1                  |
------------------------------
    """

    x_i = x

    v_max = max(v_1, v_0)
    v_min = min(v_1, v_0)

    # Base B1, below min(v_0, v_1)
    x_b1 = v_min * t

    # Side s (st_) area
    t_st = (v_max - v_min) / a
    x_st = (v_max + v_min) / 2 * t_st

    # Base B2, between v_min and v_max
    t_b2 = t - t_st
    x_b2 = t_b2 * (v_max - v_min)

    # Area of the top trapezoid
    x_t = x - x_st - x_b1 - x_b2
    t_t = t_b2

    v_c = trap_v_c(x_t, t_t, a)

    return v_c + v_max


t, v_0, v_c, v_1, a = (1, 4000, 2000, 4000, 50_000)
x = penta_area(t, v_0, v_c, v_1, a)

x, t, v_c, penta_v_c(x, t, v_0, v_1, a)