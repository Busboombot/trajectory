from dataclasses import dataclass, asdict, replace
from math import sqrt

import pandas as pd

from .exceptions import TrapMathError, ConvergenceError

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


def set_bv(x, v_0, v_1, a_max):
    """Reduce v_0 and v_1 for small values of x when there is an
    error in x after planning """

    x_a, t_a = accel_xt(v_0, 0, a_max)
    x_d = x - x_a

    if x_d < 0:
        v_0 = int(sqrt(2 * a_max * x))
        v_1 = 0
    elif x == 0:
        v_0 = 0
        v_1 = 0
    else:
        v_0 = v_0
        v_1 = int(min(sqrt(2 * a_max * x_d), v_1))

    return (v_0, v_1)

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


def consistantize(b, return_error=False):
    """Recalculate t to make x and v values integers, and everything more consistent
     This operation will maintain the original value for x, but may change t"""

    b.x_a = b.x_a
    b.x_d = b.x_d
    b.x_c = b.x - (b.x_a + b.x_d)

    b.t_a = abs((b.v_c - b.v_0) / b.joint.a_max)
    b.t_d = abs((b.v_c - b.v_1) / b.joint.a_max)

    if round(b.x_c) == 0:  # get rid of small negatives
        b.x_c = 0

    if b.v_c != 0:
        b.t_c = abs(b.x_c / b.v_c)
        #assert b.t_c >= 0, (b.v_c, b.t_c, b.x_c)
    else:
        b.t_c = 0

    b.t = b.t_a + b.t_c + b.t_d

    # Check error against area calculation
    if return_error:
        x_ad, t_ad = accel_acd(b.v_0, b.v_c, b.v_1, b.joint.a_max)
        t_c = round(b.t - t_ad, 8)  # Avoid very small negatives
        x_c = b.v_c * t_c
        x = x_ad + x_c
        x_e = x - b.x
        return b.x_c, x_e
    else:
        return b.x_c


def sign(x):
    if x == 0:
        return 0
    elif x > 0:
        return 1
    else:
        return -1

def same_sign(a, b):
    return int(a) == 0 or int(b) == 0 or sign(a) == sign(b)

@dataclass
class Joint:
    v_max: float
    a_max: float
    # distances below the small x limit will never hit v_max before needing to decelerate
    small_x: float = None # min distance for v_max->0->v_max
    max_discontinuity: float = None # Max velocity difference between adjacent blocks

    def __post_init__(self):
        self.small_x = (self.v_max ** 2) / (2*self.a_max)
        self.max_discontinuity = self.a_max/self.v_max # Max vel change in 1 step

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
    next: 'ACDBlock' = None
    prior: 'ACDBlock' = None

    flag: "str" = None
    recalcs: int = 0

    _param_names = ['x', 't', 'dir', 'v_0_max', 'v_0', 'x_a', 't_a', 'x_c', 't_c', 'v_c_max', 'v_c', 'x_d', 't_d',
                    'v_1',
                    'v_1_max']

    @property
    def id(self):
        return f"{self.segment.n}/{self.joint.n}"

    def __post_init__(self):
        self.d = sign(self.x)
        self.x = abs(self.x)

    def init(self):

        a_max = self.joint.a_max
        v_max = self.joint.v_max

        # Limit the boundary values for small moves
        self.v_0, self.v_1 = set_bv(self.x, self.v_0, self.v_1, self.joint.a_max)

        self.v_0 = min(self.v_0, v_max)
        self.v_1 = min(self.v_1, v_max)

        if self.x == 0:
            self.v_c = self.v_0 = self.v_1 = 0
            self.flag = 'Z'

        elif self.x < 2*self.joint.small_x:
            # The limit is the same one used for set_bv.
            # the equation here is the sympy solution to:
            # t_a = (v_c - v_0) / a
            # t_d = (v_1 - v_c) / -a
            # x_a = ((v_0 + v_c) / 2) * t_a
            # x_d = ((v_c + v_1) / 2) * t_d
            # x_c = x - (x_a + x_d)
            # solve(x_c, v_c)[1]

            #self.v_0 = 0
            #self.v_1 = 0
            self.v_c = (sqrt(4 * a_max * self.x + 2 * self.v_0 ** 2 + 2 * self.v_1 ** 2) / 2)
            #self.v_c = min(self.v_c, v_max)  # formula errs for short segments

            self.flag = 'S'
        else:
            # If there is more area than the triangular profile for these boundary
            # velocities, the v_c must be v_max. In this case, it must also be true that:
            #    x_ad, t_ad = accel_acd(self.v_0, v_max, self.v_1, a_max)
            #    assert self.x > x_ad
            self.v_c = v_max
            self.flag = 'M'

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, a_max)
        self.x_d, self.t_d = accel_xt(self.v_c, self.v_1, a_max)

        assert round(self.x_a + self.x_d) <= self.x

        self.x_c = self.x - (self.x_a + self.x_d)
        self.t_c = self.x_c / self.v_c if self.v_c != 0 else 0
        self.t = self.t_a + self.t_c + self.t_d

        assert round(self.x_c) >= 0, (self.x_c, self.v_c, self.flag, self)
        assert abs(self.area - self.x) < 1, (self.x, self.area, self.t, self.v_0, self.v_1)

        return self

    def plan_ramp(self, t):
        """Calculate a ramp profile for a fixed time"""

        self.t = t
        a_max = self.joint.a_max

        if self.x == 0 or self.t == 0:
            self.set_zero()
            self.t_c = self.t = t
            self.flag = 'Z'  # 5% of test cases
            return self

        # Run the binary search anyway, just in case.
        def err(v_c):
            x_a, t_a = accel_xt(self.v_0, v_c, self.joint.a_max)
            t_c = self.t - t_a

            x_c = v_c * t_c
            x_err = self.x - (x_a + x_c)

            return x_err

        try:
            guess = a_max * t + self.v_0 - sqrt(a_max * (a_max * t ** 2 + 2 * t * self.v_0 - 2 * self.x))
            guess_flag = 'c'
        except ValueError:
            guess = min(self.x/self.t, self.joint.v_max)
            guess_flag = 'm'

        try:
            self.v_1 = self.v_c = min(binary_search(err, 0, guess, self.joint.v_max), self.joint.v_max)
        except TypeError:
            assert False, ('Binary search failed', (guess, guess_flag, self.x, self.t, self.v_0))

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, self.joint.a_max)

        self.t_d, self.x_d = 0, 0

        consistantize(self)

        self.flag = "PR"

        assert self.v_0 <= self.joint.v_max
        assert self.v_c <= self.joint.v_max, self.v_c
        assert self.v_1 <= self.joint.v_max

        return self

    def plan(self, t):

        self.t = t

        a_max = self.joint.a_max

        if self.x == 0 or self.t == 0:
            self.set_zero()
            self.t_c = self.t = t
            self.flag = 'Z'  # 5% of test cases
            return self

        # Find v_c with a binary search, then patch it up if the
        # selection changes the segment time.

        def err(v_c):
            return self.x - self.replace(v_c=v_c).arear

        self.v_c = min(binary_search(err, 0,
                                 self.x / self.t, self.joint.v_max), self.joint.v_max)
        assert self.v_c <= self.joint.v_max
        self.flag = 'O'

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, a_max)
        self.x_d, self.t_d = accel_xt(self.v_c, self.v_1, a_max)

        # consistantize will make all of the values consistent with each other,
        # and if anything is wrong, it will show up in a negative x_c

        x_c = consistantize(self)

        def psfrt(t, flag):  # plan, set flag and return
            b = self.plan(t)
            b.flag = 'NT'
            return b

        x_err = abs(round(self.area) - self.x)

        if round(x_c) < 0 or (self.x > 25 and x_err > 1):
            # We've probably got a really small move, and we've already tried
            # reducing boundary velocities, so get more aggressive. If that
            # doesn't work, allow the time to expand.

            new_t = max(self.t_a + self.t_d, self.t)

            if self.v_1 > 0:
                self.v_1 = 0
                return psfrt(t, 'V1Z')

            elif self.v_0 > 0:
                self.v_0 = 0
                return psfrt(t, 'V0Z')

            elif abs(new_t - self.t) > 0.0005:

                return psfrt(new_t, 'NT')
            else:
                raise ConvergenceError(
                    'Unsolvable profile, incorrect area: '
                    f'x_c={x_c} x_err={x_err} x={self.x}, x_ad={self.x_a + self.x_d} t={self.t} t_ad={self.t_a + self.t_d}')

        if round(t, 3) != round(self.t, 3):
            # This block is still too long.
            if self.v_1 != 0:
                self.v_1 = 0
                return psfrt(t, 'HT1')
            elif round(self.v_0) > 0:
                self.v_0 = self.v_0 / 2
                return psfrt(t, 'HT1')
            else:
                pass
                # raise ConvergenceError(
                #    'Unsolvable profile, wrong time: '
                #    f't={self.t}, commanded t ={t} t_ad={self.t_a + self.t_d} v_c={self.v_c}')

        assert round(self.x_a + self.x_d) <= self.x
        #assert self.t_a + self.t_d <= self.t
        assert self.v_c >= 0, (self.v_c, self.flag)
        assert abs(self.area - self.x) < 2, (self.area, self)
        assert self.t > 0
        assert self.v_0 <= self.joint.v_max
        assert self.v_c <= self.joint.v_max
        assert self.v_1 <= self.joint.v_max

        return self

    def set_bv(self, v_0=None, v_1=None):

        v_0 = v_0 if v_0 is not None else self.v_0
        v_1 = v_1 if v_1 is not None else self.v_1

        self.v_0, self.v_1 = set_bv(self.x, v_0, v_1, self.joint.a_max)

    def asdict(self):
        return asdict(self)

    def replace(self, **kwds):
        return replace(self, **kwds)

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

        self.flag = ''  # 5% of test cases
        self.x_a = self.x_d = self.x_c = 0
        self.t_a = self.t_d = self.t_c = 0
        self.v_0 = self.v_c = self.v_1 = 0

    @property
    def dataframe(self):

        rows = []
        d = self.d
        rows.append({'t':None, 'seg':0, 'axis':0,
             'x':d*self.x_a, 'v_i':d*self.v_0, 'v_f':d*self.v_c, 'del_t':self.t_a})
        rows.append({'t': None, 'seg': 0, 'axis': 0,
                     'x':d* self.x_c, 'v_i': d*self.v_c, 'v_f': d*self.v_c, 'del_t': self.t_c})
        rows.append({'t': None, 'seg': 0, 'axis': 0,
                     'x': d*self.x_d, 'v_i': d*self.v_c, 'v_f': d*self.v_1, 'del_t': self.t_d})

        return pd.DataFrame(rows)

    def plot(self, ax=None):
        from .plot import plot_trajectory
        plot_trajectory(self.dataframe, ax=ax)

    @property
    def sim(self):
        from .sim import SimSegment
        return [
            SimSegment(self.x_a, self.v_0, self.v_c, self.d),
            SimSegment(self.x_c, self.v_c, self.v_c, self.d),
            SimSegment(self.x_d, self.v_c, self.v_1, self.d),
        ]

    def iter_steps(self, t0=0):
        for s in self.sim:
            for t, step in s.iter_period(t0=t0):
                t0 = t
                yield t,step



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
