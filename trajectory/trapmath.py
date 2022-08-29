from math import sqrt
from operator import attrgetter

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

    t_c = t - (t_a + t_d)

    if round(t_c, 4) < 0:
        raise TrapMathError(f'Negative t_c' + str((t_c, t_a + t_d, (t, v_0, v_c, v_1, a_max))))

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
        # x_ = hex_area(t, v_0, v_c, v_1, a_max)

        x_a, t_a = accel_xt(v_0, v_c, a_max)
        x_d, t_d = accel_xt(v_c, v_1, a_max)

        t_c = t - (t_a + t_d)

        # The abs is not correct for finding the area, but
        # helps a lot in the minimization error function.
        # x_c = abs(v_c * t_c)

        x_c = max(v_c * t_c, 0)

        x_ = x_a + x_c + x_d

        e = abs(x - x_)

        if collar and (x_ < 0 or v_c < 0 or v_c > v_max or t_c < 0):
            return e

        else:
            return e

    return f


def _calc_vc_min(x, v_0, v_1, v_max, a_max):
    """Find a combination of v_c and boundary velocities that can solve the given profile.
    the procedure is primarily looking for boundary values where the accel and decel
    distances don't exceed the total distance

    This one finds the minimum time profile
    """
    from scipy.optimize import minimize_scalar

    recalcs = 0

    if x == 0:
        return *[0] * 7, 'Z'

    # Reduce small entry and exit values

    if x <= (v_max ** 2) / (2 * a_max):
        # Drop the values immediately if they are small
        dt = v_max / a_max
        v_0 = min(x / dt, v_0)
        v_1 = min(x / dt, v_1)

    # boundary velocity limit. Below this limit, the reduction algorithm fails,
    # so just drive the velocities to zero. The 480K value is empirical; I don't
    # know what causes it.
    elif x < (480_000 / a_max):
        v_0 = 0
        v_1 = 0

    x_a, t_a = accel_xt(v_0, v_max, a_max)
    x_d, t_d = accel_xt(v_max, v_1, a_max)

    if x >= x_a + x_d:
        # If there is more area than the triangular profile for these boundary
        # velocities, the v_c must be v_max
        return v_0, v_max, v_1, x_a, t_a, x_d, t_d, 'M'

    # What's left are hex profiles in cases where 0 < v_c < v_max.
    # These are triangular profiles, so x_c == 0 and x == (t_a+t_d)
    def f(v_c):
        x_a, t_a = accel_xt(v_0, v_c, a_max)
        x_d, t_d = accel_xt(v_c, v_1, a_max)

        return abs(x-(x_a + x_d))

    # Quick scan of the spaceto set the inital bracket.
    mv = min((f(v_c), v_c) for v_c in range(0, 5000, 50))[1]

    r = minimize_scalar(f, bracket=(mv - 10, mv + 10))

    v_c = round(r.x)

    x_a, t_a = accel_xt(v_0, v_c, a_max)
    x_d, t_d = accel_xt(v_c, v_1, a_max)

    return v_0, v_c, v_1, x_a, t_a, x_d, t_d, 'O'

def calc_vc_min(x, v_0, v_1, v_max, a_max):
    from .params import InputParams

    v_0, v_c, v_1, x_a, t_a, x_d, t_d, r = _calc_vc_min(x, v_0, v_1, v_max, a_max)
    x_c = x - (x_a + x_d)
    assert x_c >= 0, (x_c, r)
    t_c = x_c / v_c if v_c != 0 else 0
    t = t_a + t_c + t_d

    return Params(x, t=t,
                  v_0=v_0, v_c=v_c, v_1=v_1,
                  t_a=t_a, t_c=t_c, t_d=t_d,
                  x_a=x_a, x_c=x_c, x_d=x_d,
                  v_max=v_max, a_max=a_max,
                  recalcs=r, ip=InputParams(x, v_0, v_1, v_max, a_max))


def calc_v_c_t(x, t, v_0, v_1, v_max, a_max, recurse=0):
    """ Calculate v_c for a profile. The results will be exact if v_0==v_1,
    and a guess if v_0 != v_1"""

    if recurse > 4:
        from .exceptions import ConvergenceError
        raise ConvergenceError(str((x, t, v_0, v_1, v_max, a_max)))

    if x == 0 or t == 0:
        return 0, 'Z'

    if v_0 == 0 and v_1 == 0:
        # should always be the lower root
        v_c = a_max * t / 2 - sqrt(a_max * (a_max * t ** 2 - 4 * x)) / 2
        return v_c, 'T'

    v_h = max(v_0, v_1)

    if x >= v_h * t:  # v_c must be above both v_0 and v_1, so it's a normal trap
        if v_0 == v_1:
            # subtract off v_0, v_1 and the corresponding x and run again.
            x_ = v_0 * t
            v_c = calc_v_c_t(x - x_, t, 0, 0, v_max, a_max, recurse=recurse + 1)[0]
            v_c +=  v_0

            return v_c, 'H'

    if x < v_h * t and v_0 == v_1 and v_0 == v_max:
        return v_h, 'L1'

    if x < v_h * t and v_0 == v_1 and v_0 != v_max:
        # This case might never happen
        v_c = x / t
        return v_c, 'L2'

    # Last resort, take the mean of v_0 and v_1
    v_m = (v_0 + v_1) / 2

    v_c = calc_v_c_t(x, t, v_m, v_m, v_max, a_max, recurse=recurse + 1)[0]
    v_c += v_0

    return v_c, 'C'


def hex_v_c(x, t, v_0, v_1, v_max, a_max):
    """Calculate the v_c parameter of a velocity hexagon, given all other
    parameters """
    from scipy.optimize import minimize_scalar

    args = (x, t, v_0, v_1, v_max, a_max)

    assert not any([e is None for e in args]), f'Something is null: {args}'

    f = make_area_error_func(x, t, v_0, v_1, v_max, a_max)

    # Try to calculate the v_c, if if it is good, just return it.
    # If it is close, use it for the bracket.

    vcg, qc = calc_v_c_t(x, t, v_0, v_1, v_max, a_max)
    return vcg

    r = minimize_scalar(f, method='bounded', bounds=(0, v_max))
    v_c = round(r.x)

    return v_c


def hex_v_c_p(p):
    ag = attrgetter(*'x t v_0 v_1 v_max a_max'.split())

    return hex_v_c(*ag(p))


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
