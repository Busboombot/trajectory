from math import sqrt

from .exceptions import *
from .params import *
from .trapmath import accel_xt

def limit_boundary_velocities(b):

    v_max = b.joint.v_max
    a_max = b.joint.a_max

    if b.x < (480_000 / b.a_max):
        # boundary velocity limit. Below this limit, the reduction algorithm fails,
        # so just drive the velocities to zero. The 480K value is empirical; I don't
        # know what causes it.
        b.v_0 = 0
        b.v_1 = 0

    elif b.x < (v_max ** 2) / (2 * a_max):
        # Decel limit. Can't decelerate to zero in the required distance
        dt = v_max / a_max
        b.v_0 = min(b.x / dt, b.v_0)
        b.v_1 = min(b.x / dt, b.v_1)



def plan_min_time(b: Block, inplace=True):
    """A minimum time profile for a given distance and boundary velocities."""

    from scipy.optimize import minimize_scalar

    a_max = b.joint.a_max
    v_max = b.joint.v_max

    if b.x == 0:
        b.t = b.t_a = b.t_c = b.t_d = 0
        b.x_a = b.x_c = b.x_d = 0
        b.v_0 = b.v_1 = 0

    elif b.x <= (v_max ** 2) / (2 * a_max):
        # Small movements are really difficult to process without errors, so
        # lets try to push them down to be close to constant speed
        dt = v_max / a_max
        b.v_0 = min(b.x / dt, b.v_0)
        b.v_1 = min(b.x / dt, b.v_1)

    elif b.x < (480_000 / a_max):
        # boundary velocity limit. Below this limit, the reduction algorithm fails,
        # so just drive the velocities to zero. The 480K value is empirical; I don't
        # know what causes it.
        b.v_0 = 0
        b.v_1 = 0

    x_a, t_a = accel_xt(b.v_0, v_max, a_max)
    x_d, t_d = accel_xt(v_max, b.v_1, a_max)

    if b.x == 0:
        b.v_c = 0
        flag = 'Z'
    if b.x >= x_a + x_d:
        # If there is more area than the triangular profile for these boundary
        # velocities, the v_c must be v_max
        v_c = v_max
        flag = 'M'

    else:
        # What's left are hex profiles in cases where 0 < v_c < v_max.
        # These are triangular profiles, so x_c == 0 and x == (t_a+t_d)
        def f(v_c):
            x_ad, t_ad = accel_acd(b.v_0, v_c, b.v_1, a_max)
            return abs(b.x - x_ad)

        # Quick scan of the space to set the initial bracket.
        mv = min((f(v_c), v_c) for v_c in range(0, 5000, 50))[1]

        r = minimize_scalar(f, bracket=(mv - 10, mv + 10))

        v_c = round(r.x)

        x_a, t_a = accel_xt(b.v_0, v_c, a_max)
        x_d, t_d = accel_xt(v_c, b.v_1, a_max)
        flag = 'O'

    x_c = b.x - (x_a + x_d)

    assert round(x_c) >= 0, (flag, b)

    t_c = x_c / v_c if v_c != 0 else 0

    if inplace:
        b.t = t_a + t_c + t_d
        b.v_c = v_c
        b.x_a = x_a; b.x_c = x_c; b.x_d = x_d
        b.t_a = t_a; b.t_c = t_c; b.t_d = t_d
        b.flag = flag
    else:
        return b.replace(
            t=t_a + t_c + t_d,
            v_c=v_c,
            x_a=x_a,  x_c=x_c, x_d=x_d,
            t_a=t_a,  t_c=t_c, t_d=t_d,
            flag=flag
        )


def halve_boundary_velocities(p):
    """Reduce v_1 by halves, then v_0, finally set them to 0 """

    if p.v_1 > p.v_max // 2 ** 6:
        return p.replace(v_1=p.v_1 // 2)
    elif p.v_0 > p.v_max // 2 ** 6:
        return p.replace(v_0=p.v_0 // 2)
    elif p.v_0 > 0 and p.v_1 > 0:
        return p.replace(v_0=0, v_1=0)
    else:
        raise TrapMathError()

def update_params_accel(p, t, v_0=None):
    """Update a block to a specific time, and into an accel
    form, with no decel phase. """

    v_0 = p.v_max if v_0 is None else v_0

    a_max = p.a_max

    # no decel: v_1 == v_c
    v_c = a_max * t + v_0 - sqrt(p.a_max * (p.a_max * t ** 2 + 2 * t * v_0 - 2 * p.x))
    v_c = min(v_c, p.v_max)
    x_a, t_a = accel_xt(v_0, v_c, p.a_max)

    if x_a > p.x:
        # It's just a simple trapezoid, no cruise
        v_c = sqrt(2 * p.a_max * p.x + v_0 ** 2)
        x_a, t_a = accel_xt(v_0, v_c, p.a_max)
        assert round(x_a) >= round(p.x), (x_a, p.x)

    x_c = p.x - x_a
    t_c = x_c / v_c

    p = p.replace(t=t_a + t_c, x=x_a + x_c,
                  x_a=x_a, x_c=x_c, x_d=0,
                  t_a=t_a, t_c=t_c, t_d=0,
                  v_0=v_0, v_c=v_c, v_1=v_c)
    return p

