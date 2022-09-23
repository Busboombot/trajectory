from dataclasses import dataclass, asdict, replace
from math import sqrt

import pandas as pd

from .exceptions import TrapMathError
from .stepper import DEFAULT_PERIOD, Stepper


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
    #t_ad = (abs(v_c - v_0) + abs(v_c - v_1)) / a
    x_ad = abs((v_0 ** 2 - v_c ** 2) / (2 * a)) + abs((v_1 ** 2 - v_c ** 2) / (2 * a))
    #return x_ad, t_ad

    x_a, t_a = accel_xt(v_0, v_c, a)
    x_d, t_d = accel_xt(v_c, v_1, a)

    assert round((x_a+x_d)-(x_ad),2) == 0, (x_a+x_d, x_ad)

    return x_a+x_d, t_a+t_d




def sign(x):
    if x == 0:
        return 0
    elif x > 0:
        return 1
    else:
        return -1


def same_sign(a, b):
    return int(a) == 0 or int(b) == 0 or sign(a) == sign(b)


def mean_bv(prior, next_):
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
    n: int = 0

    def __post_init__(self):
        self.small_x = (self.v_max ** 2) / (2 * self.a_max)
        self.max_discontinuity = self.a_max / self.v_max  # Max vel change in 1 step
        self.max_at = self.v_max / self.a_max

    def new_block(self, x, v_0=None, v_1=None):
        v_0 = v_0 if v_0 is not None else self.v_max
        v_1 = v_1 if v_1 is not None else self.v_max

        return Block(x=x, v_0=v_0, v_1=v_1, joint=self)


@dataclass
class Block:
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

    v_c_max = None

    joint: Joint = None
    segment: 'Segment' = None

    replans: int = 0
    reductions: "List[str]" = None
    memo: object = None
    step_period: int = DEFAULT_PERIOD

    @property
    def id(self):
        return f"{self.segment.n}/{self.joint.n}"

    def __post_init__(self):

        if self.d == 0: # ONly set if it hasn't be set in ctor
            self.d = sign(self.x)
        self.x = abs(self.x)
        self.v_c_max = self.joint.v_max
        self.reductions = []
        self.memo = []
        self.errors = []
        self.reduction_step = 0

    def min_time(self):
        """Return the smallest reasonable time to complete this block"""

        # self.set_bv(v_0, v_1, prior)

        if self.x == 0:
            v_c = 0

        elif self.x < 2. * self.joint.small_x:
            # The limit is the same one used for set_bv.
            # the equation here is the sympy solution to:
            # t_a = (v_c - v_0) / a
            # t_d = (v_1 - v_c) / -a
            # x_a = ((v_0 + v_c) / 2) * t_a
            # x_d = ((v_c + v_1) / 2) * t_d
            # x_c = x - (x_a + x_d)
            # solve(x_c, v_c)[1]

            v_c = (sqrt(4. * self.joint.a_max * self.x +
                        2. * self.v_0 ** 2 +
                        2. * self.v_1 ** 2
                        ) / 2.)
        else:
            # If there is more area than the triangular profile for these boundary
            # velocities, the v_c must be v_max. In this case, it must also be true that:
            #    x_ad, t_ad = accel_acd(self.v_0, v_max, self.v_1, a_max)
            #    assert self.x > x_ad
            v_c = self.joint.v_max


        x_ad, t_ad = accel_acd(self.v_0, v_c, self.v_1, self.joint.a_max)

        t_c = (self.x - x_ad) / v_c if v_c != 0 else 0

        t_c = max(t_c, t_ad / 2)  # Enforce 1/3 rule, each a,c,d is 1/3 of total time

        return t_c + t_ad

    def plan(self, t=None, v_0=None, v_1=None, prior=None, next_=None, iter=None):

        if t == None:
            t = self.min_time()

        if self.segment is not None:
            assert not ((prior is None) ^ (self.segment.prior is None))

        self.set_bv(v_0=v_0, v_1=v_1, prior=prior, next_=next_)

        self.t = t
        self.replans += 1

        if self.x == 0 or self.t == 0:
            self.x_a = self.x_d = self.x_c = 0
            self.t_a = self.t_d = self.t_c = 0
            self.v_0 = self.v_c = self.v_1 = 0
            self.t_c = self.t = t
            return self

        # Find v_c with a binary search, then patch it up if the
        # selection changes the segment time.

        def err(v_c):
            x_ad, t_ad = accel_acd(self.v_0, v_c, self.v_1, self.joint.a_max)

            t_c = max(self.t - t_ad, 0)
            x_c = max(v_c, 0) * t_c

            err = self.x - (x_ad + x_c)

            return err

        v_c = binary_search(err, 0, self.x / self.t, self.joint.v_max)

        self.v_c = min(v_c, self.v_c_max)

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, self.joint.a_max)
        self.x_d, self.t_d = accel_xt(self.v_c, self.v_1, self.joint.a_max)

        self.x_c = self.x - (self.x_a + self.x_d)

        self.t_a = abs((self.v_c - self.v_0) / self.joint.a_max)
        self.t_d = abs((self.v_c - self.v_1) / self.joint.a_max)

        if round(self.x_c) == 0 and self.x_c < 0:  # get rid of small negatives
            self.x_c = 0

        self.t_c = abs(self.x_c / self.v_c) if self.v_c != 0 else 0

        self.t = self.t_a + self.t_c + self.t_d

        assert self.t > 0
        assert self.v_c <= self.joint.v_max
        assert self.v_c >= 0, self.v_c
        assert abs(self.area - self.x) < 2, (self.area, self)
        assert round(self.x_a + self.x_d) <= self.x
        assert self.v_0 <= self.joint.v_max
        assert self.v_c <= self.joint.v_max
        assert self.v_1 <= self.joint.v_max

        return self

    def set_bv(self, v_0=None, v_1=None, prior=None, next_=None):

        if v_0 == 'prior' and prior is not None:
            self.v_0 = prior.v_1
        elif v_0 is not None:
            self.v_0 = v_0

        if v_1 == 'next' and next_ is not None:
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
            # Basic trapezoid formula, in terms of x instead of t
            self.v_0 = int(min(self.v_0, sqrt(2 * self.joint.a_max * self.x)))
            self.v_1 = 0
        elif self.x == 0:
            self.v_0 = 0
            self.v_1 = 0
        else:
            self.v_1 = int(min(self.v_1, sqrt(2 * self.joint.a_max * x_d)))

        self.v_0 = min(self.v_0, self.joint.v_max)
        self.v_1 = min(self.v_1, self.joint.v_max)

        return self.v_0, self.v_1

    def limit_bv(self):

        if self.v_1 > self.joint.v_max / 2:
            self.v_1 //= 2
            self.reductions.append('V1A')
            return

        if self.v_0 > self.joint.v_max / 2:
            self.v_0 //= 2
            self.reductions.append('V0A')
            return

        if self.v_1 > 1:
            self.v_1 //= 2
            self.reductions.append('V1B')
            return

        if self.v_0 > 1:
            self.v_0 //= 2
            self.reductions.append('V0B')
            return


    def reductions(self, t):
        def has_error(b, t):
            assert round(b.x_c) >= 0
            assert not (b.x > 25 and abs(round(b.area) - b.x) > 1)
            return round(t, 3) != round(b.t, 3)

        if not has_error(self, t):
            self.reductions.append('N')
            return self

        # Set v_0 and v_1 to the mean velocity
        v_m = self.x / t
        self.v_0 = min(self.v_0, v_m)
        self.v_1 = min(self.v_1, v_m)

        self._plan(t)
        if not has_error(self, t):
            self.reductions.append('VC')
            return self

        if True:
            for t_ in range(11, 18):
                self._plan(t * t_ / 10.)
                if not has_error(self, t * t_ / 10.):
                    self.reductions.append('T1')
                    return self

        # Finish off reducing v_1
        while self.v_1 > 1:
            self.v_1 = int(self.v_1 / 2)
            self._plan(t)
            if not has_error(self, t):
                self.reductions.append('V1')
                return self

        # Finish off reducing v_0
        while self.v_0 > 1:
            self.v_0 = int(self.v_0 / 2)
            self._plan(t)
            if not has_error(self, t):
                self.reductions.append('V0')
                return self

        return self


    def zero(self):
        self.x_a = self.x_d = self.x_c = 0
        self.t_a = self.t_d = self.t_c = 0
        self.v_0 = self.v_c = self.v_1 = 0
        self.t = 0


    @property
    def is_dominant(self):
        if not self.segment:
            return False

        return self.segment.dominant == self.joint.n


    @property
    def dominant(self):
        if not self.segment:
            return False

        return self.segment.blocks[self.segment.dominant]


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


    def asdict(self):
        return asdict(self)


    def replace(self, **kwds):
        return replace(self, **kwds)


    def plot(self, ax=None):
        from .plot import plot_trajectory
        plot_trajectory(self.dataframe, ax=ax)

    def stepper_blocks(self):

        def ri(v):
            return int(round(v))

        # t is the segment time, which can differe a bit from the times of the
        # individual blocks. For stepping, we adjust these times, which
        # adjusts accelerations, so all stepped block have the same time

        return (
            (ri(self.d * self.x_a), ri(self.d * self.v_0), ri(self.d * self.v_c)),
            (ri(self.d * self.x_c), ri(self.d * self.v_c), ri(self.d * self.v_c)),
            (ri(self.d * self.x_d), ri(self.d * self.v_c), ri(self.d * self.v_1))
        )

    def step(self, period=DEFAULT_PERIOD, details=False):
        stp = Stepper(period, details=details)

        t = 0

        stp.load_phases(self.stepper_blocks(self.t))

        while not stp.done:
            if details:
                yield stp.next_details()
            else:
                yield stp.next()

            t += period

    def str(self):
        from colors import color, bold

        ri = lambda v: int(round(v))

        ta = color(f"{ri(self.t_a*1000):<5d}", fg='yellow')
        tc = color(f"{ri(self.t_c * 1000):<5d}", fg='yellow')
        td = color(f"{ri(self.t_d * 1000):<5d}", fg='yellow')


        v0 = color(f"{ri(self.v_0):<5d}", fg='blue')
        xa = color(bold(f"{ri(self.d * self.x_a):>6d}"), fg='green')
        xc = color(bold(f"{ri(self.d * self.x_c):>6d}"), fg='green')
        vc = color(f"{ri(self.v_c):<6d}", fg='blue')
        xd = color(bold(f"{ri(self.d * self.x_d):<6d}"), fg='green')
        v1 = color(f"{ri(self.v_1):>5d}", fg='blue')

        return f"[{v0} {xa} {ta}↗{vc} {xc} {tc}↘{xd} {td} {v1}]"


    def _repr_pretty_(self, p, cycle):
        p.text(self.str() if not cycle else '...')
