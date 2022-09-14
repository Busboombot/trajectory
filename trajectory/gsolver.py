from dataclasses import dataclass, asdict, replace
from math import sqrt

import pandas as pd

from .exceptions import TrapMathError
from .stepper import DEFAULT_PERIOD, TIMEBASE, Stepper


def binary_search(f, v_min, v_guess, v_max):
    for i in range(20):

        x = f(v_guess)

        if round(x) > 0:
            v_guess, v_min = (v_max + v_guess) / 2, v_guess

        elif round(x) < 0:
            v_guess, v_max = (v_min + v_guess) / 2, v_guess

        else:
            return v_guess

        if abs(v_min - v_max) < .05:
            return v_guess

    else:
        return None


def accel_xt(v_i, v_f, a):
    """Distance and time required to accelerate from v0 to v1 at acceleration a"""

    if v_f == v_i:
        return 0, 0

    if v_f < v_i:
        a = -a

    t = (v_f - v_i) / a  # Time to change from v0 to v1 at max acceleration
    x = (v_i + v_f) / 2 * t

    return x, t  # Area of velocity trapezoid


def accel_acd(v_0, v_c, v_1, a):
    """ Same result as running accel_xt, for v_0->v_c and v_c->v_1,
    and adding the x and t values.
    """
    t_ad = (abs(v_c - v_0) + abs(v_c - v_1)) / a
    x_ad = abs((v_0 ** 2 - v_c ** 2) / (2 * a)) + abs((v_1 ** 2 - v_c ** 2) / (2 * a))

    return x_ad, t_ad


def sign(x):
    if x == 0:
        return 0
    elif x > 0:
        return 1
    else:
        return -1


def same_sign(a, b):
    return int(a) == 0 or int(b) == 0 or sign(a) == sign(b)

def mean_bv( prior, next_):

    if prior.t_d + next_.t_a != 0:
        a = (next_.v_c - prior.v_c) / (prior.t_d + next_.t_a)
    else:
        return (next_.v_c + prior.v_c) / 2

    v = prior.v_c + a * prior.t_d

    return v

def bent(prior, current):
    cd = current.d
    pd = prior.d

    s1 = sign(pd * prior.v_c - pd * prior.v_1)
    s2 = sign(cd * current.v_0 - cd * current.v_c)

    return s1 * s2 < 0


@dataclass
class Joint:
    v_max: float
    a_max: float
    # distances below the small x limit will never hit v_max before needing to decelerate
    small_x: float = None  # min distance for v_max->0->v_max
    max_discontinuity: float = None  # Max velocity difference between adjacent blocks
    max_at: float = None  # Max time to accel from 0 to v_max

    def __post_init__(self):
        self.small_x = (self.v_max ** 2) / (2 * self.a_max)
        self.max_discontinuity = self.a_max / self.v_max  # Max vel change in 1 step
        self.max_at = self.v_max / self.a_max

    def new_block(self, x, v_0=None, v_1=None):
        v_0 = v_0 if v_0 is not None else self.v_max
        v_1 = v_1 if v_1 is not None else self.v_max

        return ACDBlock(x=x, v_0=v_0, v_1=v_1, joint=self)


@dataclass
class ACDBlock:
    x: float = 0
    t: float = 0

    t_a: float = 0
    t_c: float = 0
    t_d: float = 0

    x_a: float = 0
    x_c: float = 0
    x_d: float = 0

    v_0: float = 0
    v_c: float = 0
    v_1: float = 0
    d: int = 0  # direction, -1 or 1

    joint: Joint = None
    segment: 'Segment' = None

    recalcs: int = 0

    step_period: int = DEFAULT_PERIOD

    _param_names = ['x', 't', 'dir', 'v_0_max', 'v_0', 'x_a', 't_a', 'x_c', 't_c', 'v_c_max', 'v_c', 'x_d', 't_d',
                    'v_1',
                    'v_1_max']

    @property
    def id(self):
        return f"{self.segment.n}/{self.joint.n}"

    def __post_init__(self):
        self.d = sign(self.x)
        self.x = abs(self.x)

    def init(self, v_0=None, v_1=None, prior=None):
        """Find a minimum time profile for the block"""

        # Limit the boundary values for small moves
        self.set_bv(v_0, v_1, prior)

        if self.x == 0:
            self.v_c = 0

        elif self.x < 2 * self.joint.small_x:
            # The limit is the same one used for set_bv.
            # the equation here is the sympy solution to:
            # t_a = (v_c - v_0) / a
            # t_d = (v_1 - v_c) / -a
            # x_a = ((v_0 + v_c) / 2) * t_a
            # x_d = ((v_c + v_1) / 2) * t_d
            # x_c = x - (x_a + x_d)
            # solve(x_c, v_c)[1]

            self.v_c = (sqrt(4 * self.joint.a_max * self.x + 2 * self.v_0 ** 2 + 2 * self.v_1 ** 2) / 2)

        else:
            # If there is more area than the triangular profile for these boundary
            # velocities, the v_c must be v_max. In this case, it must also be true that:
            #    x_ad, t_ad = accel_acd(self.v_0, v_max, self.v_1, a_max)
            #    assert self.x > x_ad
            self.v_c = self.joint.v_max

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, self.joint.a_max)
        self.x_d, self.t_d = accel_xt(self.v_c, self.v_1, self.joint.a_max)

        assert round(self.x_a + self.x_d) <= self.x

        self.x_c = self.x - (self.x_a + self.x_d)
        self.t_c = self.x_c / self.v_c if self.v_c != 0 else 0
        self.t = self.t_a + self.t_c + self.t_d

        assert round(self.x_c) >= 0, (self.x_c, self.v_c, self)
        assert abs(self.area - self.x) < 1, (self.x, self.area, self.t, self.v_0, self.v_1)

        return self

    def min_time(self, v_0=None, v_1=None, prior=None):
        """Return the smallest reasonable time to complete this block"""

        # self.set_bv(v_0, v_1, prior)

        if self.x == 0:
            v_c = 0

        elif self.x < 2 * self.joint.small_x:
            # The limit is the same one used for set_bv.
            # the equation here is the sympy solution to:
            # t_a = (v_c - v_0) / a
            # t_d = (v_1 - v_c) / -a
            # x_a = ((v_0 + v_c) / 2) * t_a
            # x_d = ((v_c + v_1) / 2) * t_d
            # x_c = x - (x_a + x_d)
            # solve(x_c, v_c)[1]

            v_c = (sqrt(4 * self.joint.a_max * self.x + 2 * self.v_0 ** 2 + 2 * self.v_1 ** 2) / 2)

        else:
            # If there is more area than the triangular profile for these boundary
            # velocities, the v_c must be v_max. In this case, it must also be true that:
            #    x_ad, t_ad = accel_acd(self.v_0, v_max, self.v_1, a_max)
            #    assert self.x > x_ad
            v_c = self.joint.v_max

        x_ad, t_ad = accel_acd(self.v_0, v_c, self.v_1, self.joint.a_max)

        x_c = self.x - x_ad
        t_c = x_c / v_c if v_c != 0 else 0

        self.t = t_c + t_ad

        return self.t

    def _plan(self, t, v_0=None, v_1=None, prior=None, next_=None):

        self.t = t

        self.set_bv(v_0=v_0, v_1=v_1, prior=prior, next_=next_)

        a_max = self.joint.a_max

        if self.x == 0 or self.t == 0:
            self.set_zero()
            self.t_c = self.t = t
            return self

        # Find v_c with a binary search, then patch it up if the
        # selection changes the segment time.

        def err(v_c):
            x_ad, t_ad = accel_acd(self.v_0, v_c, self.v_1, self.joint.a_max)

            t_c = max(self.t - t_ad, 0)
            x_c = max(v_c, 0) * t_c

            return self.x - (x_ad + x_c)

        v_c = binary_search(err, 0, self.x / self.t, self.joint.v_max)

        self.v_c = min(v_c, self.joint.v_max)

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, a_max)
        self.x_d, self.t_d = accel_xt(self.v_c, self.v_1, a_max)

        # consistantize will make all of the values consistent with each other,
        # and if anything is wrong, it will show up in a negative x_c
        self.consistantize()

        assert self.t > 0
        assert self.v_c <= self.joint.v_max
        assert self.v_c >= 0, self.v_c
        assert abs(self.area - self.x) < 2, (self.area, self)
        assert round(self.x_a + self.x_d) <= self.x
        assert self.v_0 <= self.joint.v_max
        assert self.v_c <= self.joint.v_max
        assert self.v_1 <= self.joint.v_max

        return self

    def plan(self, t, v_0=None, v_1=None, prior=None, next_=None):

        self._plan(t, v_0=v_0, v_1=v_1, prior=prior, next_=next_)

        def has_error(b):
            x_err = abs(round(b.area) - b.x)
            return round(b.x_c) < 0 or (b.x > 25 and x_err > 1) or round(t, 3) != round(b.t, 3)

        if not has_error(self):
            return self

        # Try various strategies to fix the errors

        # Reduce v_1 a bit
        while self.v_1 > 1500:
            self.v_1 = int(self.v_1 / 2)
            self._plan(t, prior=prior, next_=next_)
            if not has_error(self):
                return self

        # Maybe it just needs a re-calc to allow time to expand?
        self._plan(self.t, prior=prior, next_=next_)
        if not has_error(self):
            return self

        # Reduce v_0 to v_c
        if self.v_0 > self.v_c:
            self.v_0 = self.v_c
            self._plan(t, prior=prior, next_=next_)
            if not has_error(self):
                return self

        # Finish off reducing v_1
        while self.v_1 > 1:
            self.v_1 = int(self.v_1 / 2)
            self._plan(t, prior=prior, next_=next_)
            if not has_error(self):
                return self

        # Finish off reducing v_0
        while self.v_0 > 1:
            self.v_0 = int(self.v_0 / 2)
            self._plan(t, prior=prior, next_=next_)
            if not has_error(self):
                return self

    def consistantize(self, return_error=False):
        """Recalculate t to make x and v values integers, and everything more consistent
         This operation will maintain the original value for x, but may change t"""

        self.x_c = self.x - (self.x_a + self.x_d)

        self.t_a = abs((self.v_c - self.v_0) / self.joint.a_max)
        self.t_d = abs((self.v_c - self.v_1) / self.joint.a_max)

        if round(self.x_c) == 0:  # get rid of small negatives
            self.x_c = 0

        if self.v_c != 0:
            self.t_c = abs(self.x_c / self.v_c)
            # assert b.t_c >= 0, (b.v_c, b.t_c, b.x_c)
        else:
            self.t_c = 0

        self.t = self.t_a + self.t_c + self.t_d

        # Check error against area calculation
        if return_error:
            x_ad, t_ad = accel_acd(self.v_0, self.v_c, self.v_1, self.joint.a_max)
            t_c = round(self.t - t_ad, 8)  # Avoid very small negatives
            x_c = self.v_c * t_c
            x = x_ad + x_c
            x_e = x - self.x
            return self.x_c, x_e
        else:
            return self.x_c

    @property
    def area(self):
        """Calculate the distance x as the area of the profile. Ought to
        always match .x """

        x_ad, t_ad = accel_acd(self.v_0, self.v_c, self.v_1, self.joint.a_max)

        t_c = round(self.t - t_ad, 8)  # Avoid very small negatives

        if t_c < 0:
            raise TrapMathError(f'Negative t_c ({t_c}) ')

        x_c = self.v_c * t_c

        if round(x_c) < 0:
            raise TrapMathError(f'Negative x_c ({x_c}) t_c={t_c}')

        return x_ad + x_c

    @property
    def arear(self):
        """Robust area calculation. May return wrong results for bad parameters, but
        will not throw exceptions"""

        x_ad, t_ad = accel_acd(self.v_0, self.v_c, self.v_1, self.joint.a_max)

        t_c = max(self.t - t_ad, 0)
        x_c = max(self.v_c, 0) * t_c

        return x_ad + x_c

    def set_zero(self):

        self.x_a = self.x_d = self.x_c = 0
        self.t_a = self.t_d = self.t_c = 0
        self.v_0 = self.v_c = self.v_1 = 0



    @property
    def dataframe(self):

        rows = []
        d = self.d
        rows.append({'t': None, 'seg': 0, 'axis': 0,
                     'x': d * self.x_a, 'v_i': d * self.v_0, 'v_f': d * self.v_c, 'del_t': self.t_a})
        rows.append({'t': None, 'seg': 0, 'axis': 0,
                     'x': d * self.x_c, 'v_i': d * self.v_c, 'v_f': d * self.v_c, 'del_t': self.t_c})
        rows.append({'t': None, 'seg': 0, 'axis': 0,
                     'x': d * self.x_d, 'v_i': d * self.v_c, 'v_f': d * self.v_1, 'del_t': self.t_d})

        return pd.DataFrame(rows)

    def set_bv(self, v_0=None, v_1=None, prior=None, next_=None):

        assert v_1 != 'prior' and v_0 != 'next'

        if (v_0 == 'prior' or v_0 == 'mean') and prior is not None:
            self.v_0 = prior.v_1
        elif v_0 == 'v_max':
            self.v_0 = self.joint.v_max
        elif v_0 is not None:
            self.v_0 = v_0

        if (v_1 == 'next' or v_1 == 'mean') and next_ is not None:
            self.v_1 = next_.v_0
        elif v_1 == 'v_max':
            self.v_1 = self.joint.v_max
        elif v_1 is not None:
            self.v_1 = v_1

        if prior:
            # If the current block has a different sign -- changes direction --
            # then the boundary velocity must be zero.
            if not same_sign(prior.d, self.d) or prior.x == 0 or self.x == 0:
                self.v_0 = 0

        x_a, t_a = accel_xt(self.v_0, 0, self.joint.a_max)
        x_d = self.x - x_a

        if x_d < 0:
            self.v_0 = int(min(self.v_0, sqrt(2 * self.joint.a_max * self.x)))
            self.v_1 = 0
        elif self.x == 0:
            self.v_0 = 0
            self.v_1 = 0
        else:
            self.v_1 = int(min(sqrt(2 * self.joint.a_max * x_d), self.v_1))

        # Get rid of kinks at the boundary, where the v_1/v_0 values
        # are much different from the v_c values of the adjacent blocks, and
        # the v_c blocks are similar.
        #
        # So, change ___/\___ into ------
        #
        # We can reduce v_1/v_0 to match v_c without much concern, but
        # should generally not increase the boundary velocity, so
        # there is a limit on the size of increases.

        if v_0 == 'mean' and prior is not None and bent(prior, self):
            v = mean_bv(prior, self)
            if v <= self.v_0 and v <= prior.v_1:
                self.v_0 = v
            elif max(abs(self.v_0 - v), abs(prior.v_1 - v)) < 2000:
                self.v_0 = v


        if v_1 == 'mean' and next_ is not None and bent(self, next_):
            v = mean_bv(self, next_)
            if v <= self.v_1 and v <= next_.v_0:
                self.v_1 = v
            elif max(abs(self.v_1 - v), abs(next_.v_0 - v)) < 2000:
                self.v_1 = v


        self.v_0 = min(self.v_0, self.joint.v_max)
        self.v_1 = min(self.v_1, self.joint.v_max)

        return self.v_0, self.v_1

    def asdict(self):
        return asdict(self)

    def replace(self, **kwds):
        return replace(self, **kwds)

    def plot(self, ax=None):
        from .plot import plot_trajectory
        plot_trajectory(self.dataframe, ax=ax)

    def steppers(self, period=None):

        if period is None:
            period = self.step_period

        return [
            Stepper(self.d * self.x_a, self.v_0, self.v_c, period),
            Stepper(self.d * self.x_c, self.v_c, self.v_c, period),
            Stepper(self.d * self.x_d, self.v_c, self.v_1, period),
        ]

    def iter_steps(self, until_t=None, period=None):

        if period is None:
            period = self.step_period

        if until_t is None:
            until_t = int(self.segment.time * period * TIMEBASE)
        else:
            until_t *= int(period * TIMEBASE)

        t = 0
        for s in self.steppers(period):

            while not s.done:
                yield next(s)
                t += period
                if t >= until_t:
                    return

        while True:
            yield 0
            t += period
            if t >= until_t:
                return

    def str(self):
        from colors import color, bold

        v0 = color(f"{int(self.v_0):<5d}", fg='green')
        xa = bold(f"{int(self.d * self.x_a):>6d}")
        c = bold(f"{int(round(self.d * self.x_c)):>6d}")
        vc = color(f"{int(self.v_c):<6d}", fg='blue')
        xd = bold(f"{int(self.d * self.x_d):<6d}")
        v1 = color(f"{int(self.v_1):>5d}", fg='red')

        return f"[{v0} {xa}↗{c + '@' + vc}↘{xd} {v1}]"

    def _repr_pretty_(self, p, cycle):
        p.text(self.str() if not cycle else '...')
