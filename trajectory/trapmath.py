
from math import sqrt

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