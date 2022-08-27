from math import sqrt



def hex_v_c_convex(x, t, v_0, v_1, v_max, a_max):
    """ Find the v_c for a pentagon, through decompositions

    For a velocity hexagon, and some other reduced shapes, this routine
    will break the pentagon into four regions, trivially calculating the
    areas of the first three ( B1, S, B2 ) then with the remaining area in T,
    calculate Vc.

            /--------------------          Vc
           /                     \
          /          T            \
         / ________________________\       Vmax ( V0 )
        |        B2               |S\
        |_________________________|__\     Vmin ( V1 )
        |        B1                  |
        ------------------------------
    """

    if x == 0 or t == 0:
        return 0

    vcm = min(v_max, v_c_max_t(t, v_0, v_1, a_max))  # highest velocity possible for time
    vcm_x = hex_area(t, v_0, vcm, v_1, a_max)  # distance traveled at max velocoty.

    if x - vcm_x > 1:
        x_a, t_a = accel_tx(v_0, vcm, a_max)
        x_d, t_d = accel_tx(vcm, v_1, a_max)

        #assert round(vcm_x, 2) == round(x_a + x_d, 2), (vcm_x, x_a + x_d, (t, v_0, v_1, a_max))

        raise ShortTimeError(t_a + t_d,
                             f'Distance ({x}) is larger than largest possible for time (max_x={round(vcm_x)})')

    v_high = max(v_1, v_0)
    v_low = min(v_1, v_0)

    # distance covered with acceleration directly from v_0 to v_1
    # break_x = (v_1+v_0)/2 *t
    # distance covered with a D or U shape
    du_x = hex_area(t, v_high, v_high, v_low, a_max)

    if round(du_x) == x:
        return v_high

    # =======================

    # Base B1, below min(v_0, v_1)
    x_b1 = v_low * t

    # Side s (st_) area
    t_st = (v_high - v_low) / a_max
    x_st = (v_high + v_low) / 2 * t_st

    # if v_high == v_low, there is not st nor b2
    assert not ((v_high == v_low ) ^ (x_st == 0))

    # Base B2, between v_low and v_high
    t_b2 = t - t_st
    x_b2 = t_b2 * (v_high - v_low)

    assert not ((v_high == v_low) ^ (x_b2 == 0))

    # Area of the top trapezoid
    x_t = x - x_st - x_b1 - x_b2
    t_t = t_b2

    if round(x_t) == 0:
        return v_high

    if t_t <= 0 and abs(t_t) > 1 and round(x) > 0:
        raise ShortTimeError(f'No time left for v_c t_t = {round(t_t, 6)} ')

    try:
        v_c = trap_v_c(x_t, t_t, a_max) + v_high
    except TrapMathError as e:
        # Guess this is a DECEL and there isn't enough time.
        x_, t_ = accel_tx(v_high, v_low, a_max)
        if t_ > t:
            raise ShortTimeError(t_, f'Need at least {t_}s to decel. t={round(t, 6)}')

        raise TrapMathError(f"{v_high} {v_low} {x} {x_} {t},{t_}")

    x_ = hex_area(t, v_0, v_c, v_1, a_max)

    return v_c


def v_c_max_x(x, v_0, v_1, a):
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


def v_c_max_t(t, v_0, v_1, a):
    """Find the max v_c, but doesn't assume decellerating 1/2 way through"""

    # Check that these arguments don't exceed acceleration limits
    if round(abs(v_0 - v_1) / t) > round(a):
        # cant get from v_0 to v_1 in time t within acceleration limits.
        raise LimitError(f'Exceeded acceleration limit: {abs(v_0 - v_1) / t} > {a}')

    t_a = (v_1 - v_0 + a * t) / (2 * a)

    v_c = v_0 + a * t_a

    # Check by solving from the other direction
    t_d = t - t_a
    v_c_ = v_1 + a * t_d

    assert abs(v_c - v_c_) < 1, (round(v_c), round(v_c_), t, v_0, v_1, a)

    return v_c


def min_v_time(x, v_0, v_1, v_max, a_max):
    """Calculate the minimum time to get from v_0 to v_1"""
    pass



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
        # warn(f"{kwds.get('id')} v_1 {v_1} exceeds v_c {v_c}")
        # raise ValidationError()

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

def penta_area(v_0, v_c, v_1, a):
    """Return the area and time for a velocity pentagon, which looks like a house,
    there is no cruise"""

    x_a, t_a = accel_tx(v_0, v_c, a)
    x_d, t_d = accel_tx(v_1, v_c, a)

    return x_a + x_d, t_a + t_d

    return x if x > 0 else 0

def max_v0_for_x(x, a):
    """Maximum V0 for segment to allow decelereration and not cover more than
        distance x; it sets the v_0_max"""

    return sqrt(2 * x * a)

def t_accel(vi, vf, a):
    """Time to accelerate from vi to vf"""
    return (vf - vi) / a


def t_accel_x(x, v, a):
    """Time required to move distance x, accelerating  from velocity vi"""

    return sqrt(2 * a * x + v ** 2) / a - (v / a)



def nearly_equal_velocity(a, b):
    """Close enough not to trigger an update"""
    try:
        return abs(a - b) / abs(a + b) < .02
    except ZeroDivisionError:
        return a == b  # Must still check b/c options are either a == b == 0 or a == -b


def nearly_equal_time(a, b):
    return abs(a - b) / abs(a + b) < .02

