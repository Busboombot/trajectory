from enum import Enum
from .exceptions import *
from math import sqrt
from .params import *
from .trapmath import accel_xt


def reduce_v0_v1(p, recurse=None):
    """Reduce v_1 and maybe v_0 to correct an error in x,
    usually because the boundary velocities are too high for the x"""

    if p.x_c >= 0:
        return p  # Nothing to fix

    dx = abs(p.x_c)

    dv = 2 * dx / p.t  # corresponding change in velocities

    # boundary velocity limit. Below this limit, the reduction algorithm fails,
    # so just drive the velocities to zero. The 480K value is empirical; I don't
    # know what causes it.

    if p.x < (480_000/p.a_max) or (recurse is not None and recurse == 0):
        v_0 = 0
        v_1 = 0
    elif (p.v_1 < 800 or dx < 20) and dv <= p.v_1:
        # For small v_1, it isn' worth the multiple rounds of recursion
        # since we'll probably drive it to zero anyway
        v_1 =0
        v_0 = p.v_0
    elif dv <= p.v_1:  # Just change v_1
        v_1 = round(p.v_1 - dv)
        v_0 = p.v_0
    elif dv <= p.v_0 + p.v_1:  # change both v_0 and v_1
        dv -= p.v_1
        v_1 = 0
        v_0 = round(p.v_0 - dv)
    else:
        assert False, 'Should not get here'
        # Probably never get here?
        v_0 = 0
        v_1 = 0

    #print("ARGS", dx, (p.x, v_0, v_1, p.v_max, p.a_max) )

    return basic_profile(p.x, v_0, v_1, p.v_max, p.a_max, recurse=recurse)

def basic_profile(x, v_0, v_1, v_max, a_max, throw=False, recurse=None):
    """Return a minimum time triangle, trapezoid, pentagon or hex profile """

    if recurse is not None and recurse == 0:
        raise ParameterError(None, 'Recursion depth exceeded')

    if x < 0:
        raise ParameterError(None, 'Distance must be positive')

    # Time to get to v_c, which is presumed to be less than v_max

    ####
    #### THIS EQUATION IS PROBABLY WRONG! I think it is actually the total time for the segment.
    #### not the interstion time, which are in hex_alt_area:
    ## time of lower intersection of acceleration lines
    #t_l = (a_max * t + v_0 - v_1) / (2 * a_max)
    ## time of upper intersection
    #t_u = (a_max * t - v_0 + v_1) / (2 * a_max)
    t = (-v_0 - v_1 + sqrt(4 * a_max * x + 2 * v_0 ** 2 + 2 * v_1 ** 2)) / a_max
    v_c = a_max * t / 2 + v_0 / 2 + v_1 / 2  # ... which is also the time we hit v_c

    if v_c > v_max:
        v_c = v_max

    x_a, t_a = accel_xt(v_0, v_c, a_max)
    x_d, t_d = accel_xt(v_c, v_1, -a_max)

    x_c = x - x_a - x_d
    x_c_r = round(x_c)

    flag = 'BP' # Normal return
    if x_c_r > 0:
        assert v_c == v_max, (x_c_r, v_c, v_max)
        t_c = x_c/v_c
    elif x_c_r == 0:
        t_c = 0
    else: # less than 0, which is an error
        t_c = 0
        flag = 'E' # Signal error.

    p = Params(x, t=t_a + t_c + t_d, flag=flag,
                  v_0=v_0, v_c=v_c, v_1=v_1,
                  t_a=t_a, t_c=t_c, t_d=t_d,
                  x_a=x_a, x_c=x_c_r, x_d=x_d,
                  v_max=v_max, a_max=a_max)

    if p.x_c < 0:
        if throw:
            raise ParameterError(p, f'x_c ({p.x_c}) is negative')
        else:
            # Try reducing v_0 and v_1
            return reduce_v0_v1(p, recurse-1 if recurse is not None else 3)
    else:
        return p


MIN_SEGMENT_TIME = 0.1

def min_profile(x, v_0, v_1, v_max, a_max):
    """A minimum time profile for a given distance and boundary velocities."""

    # Small movements are really difficult to process without errors, so
    # lets try to push them down to be close to constant speed.
    if x < (v_max**2)/(2*a_max):
        dt = v_max / a_max
        v_0 = min(x/dt, v_0)
        v_1 = min(x/dt, v_1)



    def find_v_c(x, v_0, v_1, v_max, a_max):
        """Find a combination o v_c and boundary velocities that can solve the given profile"""
        v_c = v_max
        recalcs = 0
        for i in range(8):
            while v_c > v_max / 16:

                x_a, t_a = accel_xt(v_0, v_c, a_max)
                x_d, t_d = accel_xt(v_c, v_1, a_max)

                if x >= x_a + x_d:
                    return v_0, v_c, v_1, x_a, t_a, x_d, t_d, recalcs

                v_c = v_c // 2
                recalcs += 1

            v_c = v_max
            v_0 = v_0 // 2
            v_1 = v_1 // 2

        raise TrapMathError('Failed to find v_c '+str((x, v_0, v_1, v_max, a_max)))

    v_0, v_c, v_1, x_a, t_a, x_d, t_d, recalcs = find_v_c(x, v_0, v_1, v_max, a_max)

    x_c = x - x_a - x_d
    t_c = x_c / v_c

    return Params(x, t=t_a + t_c + t_d,
                  v_0=v_0, v_c=v_c, v_1=v_1,
                  t_a=t_a, t_c=t_c, t_d=t_d,
                  x_a=x_a, x_c=x_c, x_d=x_d,
                  v_max=v_max, a_max=a_max,
                  recalcs = recalcs, ip=InputParams(x, v_0, v_1, v_max, a_max))


def update_params(p, t):
    """Update parameter block with a new time"""
    from .trapmath import hex_area_p, hex_v_c

    assert round(t - p.t, 6) >= 0, (t, p.t) # Time must get larger

    v_m = p.x / t  # mean time

    v_c = hex_v_c(p.x, t, p.v_0, p.v_1, p.v_max, p.a_max)

    x_a, t_a = accel_xt(p.v_0, v_c, p.a_max)
    x_d, t_d = accel_xt(v_c, p.v_1, p.a_max)

    x_c = p.x - (x_a + x_d)

    if v_c == 0 or x_c == 0:
        t_c = t - (t_a + t_d)
        x_c = 0
    else:
        t_c = x_c/v_c

    x_ = x_a + x_c + x_d
    t_ = t_a + t_c + t_d

    x_ha = hex_area_p(p)

    assert round(x_) == round(p.x), (p.x, x_, x_ha)
    assert round(t - t_, 2) >= 0, (t, t_)
    assert round(v_c) >= 0

    p = p.replace(t=t_,
                  x_a=x_a, x_c=x_c,x_d=x_d,
                  t_a=t_a, t_c=t_c, t_d=t_d,
                  v_c=v_c)

    assert round(hex_area_p(p)) == round(p.x), (hex_area_p(p), p.x)

    return p
