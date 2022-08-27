
from math import sqrt
from operator import attrgetter
from typing import List

from scipy.optimize import minimize_scalar

from trajectory.exceptions import TrapMathError
from .params import Params

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


def hex_area(t, v_0, v_c, v_1, a_max):

    x_a, t_a = accel_xt(v_0, v_c, a_max)
    x_d, t_d = accel_xt(v_c, v_1, a_max)

    t_c = t - (t_a + t_d )

    if t_c < 0:
        raise TrapMathError(f'Negative t_c'+str( ( t_c, (t, v_0, v_c, v_1, a_max))))

    x_c = v_c * t_c

    return x_a + x_c + x_d

def hex_area_p(p):
    from operator import attrgetter
    ag = attrgetter(*'t v_0 v_c v_1  a_max'.split())
    return hex_area(*ag(p))


def make_area_error_func(x, t, v_0, v_1, v_max, a_max, collar=True):
    # Compute an error for outside the boundaries, which has a slope back
    # to the valid range

    def f(v_c):
        #x_ = hex_area(t, v_0, v_c, v_1, a_max)

        x_a, t_a = accel_xt(v_0, v_c, a_max)
        x_d, t_d = accel_xt(v_c, v_1, a_max)

        t_c = t - (t_a + t_d)

        # The abs is not correct for finding the area, but
        # helps a lot in the minimization error function.
        # x_c = abs(v_c * t_c)

        x_c = max(v_c * t_c,0)

        x_ = x_a + x_c + x_d

        e = abs(x - x_)

        if collar and (x_ < 0 or v_c < 0 or v_c > v_max or t_c < 0):
            return e

        else:
            return e

    return f


def calc_v_c(x, t, v_0, v_1, v_max, a_max, return_qc=False, recurse=0):
    """ Calculate v_c for a profile. The results will be exact if v_0==v_1,
    and a guess if v_0 != v_1"""

    if x == 0 or t == 0:
        if return_qc:
            return [0,0], 0  # 0 == exact
        else:
            return [0,0]

    if v_0 == 0 and v_1 == 0:
        # should always be the lower root
        v_c = [a_max * t / 2 - sqrt(a_max * (a_max * t ** 2 - 4 * x)) / 2,
               a_max * t / 2 + sqrt(a_max * (a_max * t ** 2 - 4 * x)) / 2]

        v_c = [e for e in v_c if e <= v_max and e >= 0]

        assert len(v_c) > 0, ('Got no quesses', x, t, v_0, v_1)

        if return_qc:
            return v_c,0 # 0 == exact
        else:
            return v_c

    v_h = max(v_0, v_1)

    if x >= v_h * t:  # v_c must be above both v_0 and v_1, so it's a normal trap
        if v_0 == v_1:
            # substract off v_0, v_1 and the corresponding x and run again.
            x_ = v_0 * t
            v_c = calc_v_c(x - x_, t, 0, 0, v_max, a_max, recurse=recurse+1)

            v_c = [e+v_0 for e in v_c]

            if return_qc:
                return v_c, 0  # 0 == exact
            else:
                return v_c

    if x < v_h*t and v_0 == v_1 and v_0 == v_max:

        if return_qc:
            return [v_h,v_h], -1  # 0 == exact
        else:
            return [v_h,v_h]

        #assert False, ('Impossible distance for velocities', x, v_h, v_h*t)


    # Last resort, take the mean of v_0 and v_1
    v_m = (v_0 + v_1) / 2
    if recurse > 4:
        from .exceptions import ConvergenceError
        raise ConvergenceError()

    v_c = calc_v_c(x, t, v_m, v_m, v_max, a_max, recurse=recurse+1)
    v_c = [e + v_0 for e in v_c]

    if return_qc:
        return v_c, -1
    else:
        return v_c

def hex_v_c(x, t, v_0, v_1, v_max, a_max):
    """Calculate the v_c parameter of a velocity hexagon, given all other
    parameters """
    from scipy.optimize import minimize_scalar, minimize, basinhopping

    args = (x, t, v_0, v_1, v_max, a_max)

    assert not any([e is None for e in args]), f'Something is null: {args}'

    f = make_area_error_func(x, t, v_0, v_1, v_max, a_max)

    # Try to calculate the v_c, if if it is good, just return it.
    # If it is close, use it for the bracket.

    ms_args = {'bounds': (0, v_max) }
    guess = v_max / 2
    try:
        vcg, qc = calc_v_c(x, t, v_0, v_1, v_max, a_max, True)

        vcg = [ min(e, v_max) for e in vcg]

        if qc == 0:
            return vcg[0]
        else:
            if len(vcg) == 2:
                dv = (vcg[1] - vcg[0]) / 2
            else:
                dv = 200

            ms_args = {'bracket': (vcg[0] - 5, vcg[0] + 5),
                    'bounds': (vcg[0] - dv, vcg[0] + dv)}

            guess = vcg[0]

    except Exception as e:
        print(type(e), e, (x, t, v_0, v_1, v_max, a_max))

    r = minimize_scalar(f, method='bounded', bounds=(0,v_max))
    #r = minimize_scalar(f)
    v_c = round(r.x)

    #r = basinhopping(f, [guess])
    #v_c = round(r.x[0])

    return v_c


def hex_v_c_p(p):
    ag = attrgetter(*'x t v_0 v_1 v_max a_max'.split())

    return hex_v_c(*ag(p))


def reduce_v0_v1(dx, t, v_0, v_1):
    """Reduce v_1 and maybe v_0 to correct an error in x,
    usually because the boundary velocities are too high for the x"""

    dv = 2 * dx / t  # corresponding change in velocities

    if dv <= v_1:  # Just change v_1
        v_1 -= dv
    elif dv < v_0 + v_1:  # change both v_0 and v_1
        dv -= v_1
        v_1 = 0
        v_0 -= dv
    else:
        raise TrapMathError()

    return v_0, v_1

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
    # <2 to allow a little tolerance in calculations
    assert abs(a.x - ax) < 2, (a.x, ax)


def save_for_change(l: List[Params]):
    """Save joints for memoization, to discover changes"""
    from copy import copy

    return [copy(e) if e else None for e in l]


def assert_unchanged_area(l: List[Params], memo: List[Params]):
    for a, b in zip(l, memo):
        assert_equal_area(a, b)


def update_boundary_velocities(p: Params, c: Params, n: Params):
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
