from math import sqrt
from operator import attrgetter

from scipy.optimize import minimize_scalar

from . import same_sign, maxmin
from .params import *


def limit_boundary_velocities(b):
    v_max = b.joint.v_max
    a_max = b.joint.a_max

    if x == 0:
        # Zero length segments must be zero
        b.v_0 = b.v_0_max = 0
        b.v_1 = b.v_1_max = 0
    elif b.x < (480_000 / b.a_max):
        # boundary velocity limit. Below this limit, the reduction algorithm fails,
        # so just drive the velocities to zero. The 480K value is empirical; I don't
        # know what causes it.
        b.v_0 = b.v_0_max = 0
        b.v_1 = b.v_1_max = 0

    elif b.x < (v_max ** 2) / (2 * a_max):
        # Decel limit. Can't decelerate to zero in the required distance
        dt = v_max / a_max
        b.v_0 = min(b.x / dt, b.v_0)
        b.v_1 = min(b.x / dt, b.v_1)




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