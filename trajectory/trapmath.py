from math import sqrt
from operator import attrgetter

from scipy.optimize import minimize_scalar

from .params import Block

X_LOWER_LIMIT = 10


def maxmin(l, v, h):
    return max(min(v, h), l)


def sign(x):
    if x == 0:
        return 0
    elif x > 0:
        return 1
    else:
        return -1


def same_sign(a, b):
    return int(a) == 0 or int(b) == 0 or sign(a) == sign(b)


def accel_xt(v_i, v_f, a):
    """Distance and time required to accelerate from v0 to v1 at acceleration a"""

    if v_f == v_i:
        return 0, 0

    if v_f < v_i:
        a = -a

    t = (v_f - v_i) / a  # Time to change from v0 to v1 at max acceleration
    x = (v_i + v_f) / 2 * t

    return x, t  # Area of velocity trapezoid


def make_area_error_func(x, t, v_0, v_1, v_max, a_max, collar=True):
    # Compute an error for outside the boundaries, which has a slope back
    # to the valid range

    def f(v_c):
        from .params import accel_acd

        x_ad, t_ad = accel_acd(v_0, v_c, v_1, a_max)

        t_c = t - t_ad
        x_c = v_c * t_c

        if round(x_c) < 0:
            return v_max

        return abs(x - ( x_c + x_ad ))

    return f


def hex_v_c(b):
    """Calculate the v_c parameter of a velocity hexagon, given all other
    parameters """

    ag = attrgetter(*'x t v_0 v_1'.split())
    assert not any([e is None for e in ag(b)]), f'Something is null: {ag(b)}'

    a_max = b.joint.a_max
    v_h = max(b.v_0, b.v_1)

    ##
    ## First try to calculate the values
    ##

    if b.x == 0 or b.t == 0:
        return 0, 'Z'

    elif b.v_0 == 0 and b.v_1 == 0:
        # should always be the lower root
        v_c = a_max * b.t / 2 - sqrt(a_max * (a_max * b.t ** 2 - 4 * b.x)) / 2
        return v_c, 'T'

    elif b.x >= v_h * b.t and b.v_0 == b.v_1:
        # subtract off v_0, v_1 and the corresponding x and run again.
        x_ = b.v_0 * b.t

        # This is the same case as v_0 == v_1 == 0, but after subtracting off
        # the area of the rectangular base, from v =0 to v=v_0. Then we add back in
        # v_0.
        v_c = b.v_0 + a_max * b.t / 2 - sqrt(a_max * (a_max * b.t ** 2 - 4 * (b.x - x_))) / 2

        return v_c, 'H'

    ## Well, that didn't work,
    ## so let try minimization

    f = make_area_error_func(b.x, b.t, b.v_0, b.v_1, b.joint.v_max, b.joint.a_max)

    mv = b.x/b.t # this usually performs just as well as the scan, and is cheaper.

    try:
        # The bounded method behaves badly some times, getting the wrong
        # minimum. In thse cases, the area call will fail, and we run the
        # other minimize call, which performs well in these cases.
        r = minimize_scalar(f, method='bounded',
                            bracket=(mv - 10, mv + 10), bounds=(0, b.joint.v_max))
        b.replace(v_c=r.x).area
    except:
        r = minimize_scalar(f, bracket=(mv - 10, mv + 10) )

    return r.x, 'O'


def triangular_area(v_0, v_1, v_c, a_max):
    """Area of a triangular profile, which accelerates to v_c and
    immediately decelerates"""
    return (2 * a_max * (v_0 + v_1) + (v_0 - v_c) ** 2 + (v_1 - v_c) ** 2) / (2 * a_max)


#####
# Group operations
#####

def calc_boundary_values(x, v_0, v_1, v_max, a_max):
    v_0_max = v_max
    v_1_max = v_max

    if x == 0:
        # Zero length segments must be zero
        v_0 = v_0_max = 0
        v_1 = v_1_max = 0
    elif x < (480_000 / a_max):
        # boundary velocity limit. Below this limit, the reduction algorithm fails,
        # so just drive the velocities to zero. The 480K value is empirical; I don't
        # know what causes it, but it seems to apply for a variety of accelerations
        # This value is usually < 10
        v_0 = v_0_max = 0
        v_1 = v_1_max = 0
    elif x < (v_max ** 2) / (2 * a_max):
        # This is the distance traveled in a full speed-acceleration from
        # vmax to 0. This value is usually < 0
        dt = v_max / a_max
        v_0 = v_0_max = min(x / dt, v_0)
        v_1 = v_1_max = min(x / dt, v_1)
    else:
        # Just return the input values.
        pass

    return v_0, v_1, v_0_max, v_1_max


def update_boundary_velocities(p: Block, c: Block, n: Block):
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
