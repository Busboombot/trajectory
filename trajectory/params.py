from dataclasses import dataclass, asdict, replace
from math import sqrt
from operator import attrgetter

from .exceptions import *


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

def consistantize(b):
    """Recalculate t to make x and v values integers, and everything more consistent
     This operation will maintain the original value for x, but may change t"""

    b.x_a = round(b.x_a)
    b.x_d = round(b.x_d)
    b.x_c = b.x -(b.x_a+b.x_d)

    b.v_0, b.v_c, b.v_1 = [round(e) for e in (b.v_0, b.v_c, b.v_1)]

    b.t_a = abs((b.v_c - b.v_0) / b.joint.a_max)
    b.t_d = abs((b.v_c - b.v_1) / b.joint.a_max)

    b.t_c = b.x_c/b.v_c if b.v_c != 0 else 0

    b.t = b.t_a + b.t_c + b.t_d

    return b.x_c


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


def maxmin(l, v, h):
    return max(min(v, h), l)

def make_area_error_func(b):
    def f(v_c):
        return b.x - b.replace(v_c=v_c).arear

    return f


def set_bv(x, v_0, v_1, a_max):
    """Reduce v_0 and v_1 for small values of x when there is an
    error in x after planning """

    x_a, t_a = accel_xt(v_0, 0, a_max)
    x_d = x - x_a

    if x_d < 0:
        v_0 = sqrt(2 * a_max * x)
        v_1 = 0
    else:
        v_0 = v_0
        v_1 = min(sqrt(2 * a_max * x_d), v_1)

    return (v_0, v_1)





@dataclass
class Joint:
    v_max: float
    a_max: float
    # distances below the small x limit will never hit v_max before needing to decelerate
    small_x: float = None

    def __post_init__(self):
        self.small_x = (self.v_max ** 2) / (self.a_max)

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
    joint: Joint = None
    flag: str = None
    recalcs: int = 0

    segment = None
    next: "Block" = None
    prior: "Block" = None

    def __post_init__(self):
        self.d = sign(self.x)

    def asdict(self):
        return asdict(self)

    def replace(self, **kwds):
        return replace(self, **kwds)

    @property
    def subsegments(self):
        rd = round  # lambda v, n: v
        return (
            (rd(self.t_a, 7), rd(self.v_0, 2), rd(self.v_c, 2), rd(self.x_a, 0), 'a'),
            (rd(self.t_c, 7), rd(self.v_c, 2), rd(self.v_c, 2), rd(self.x_c, 0), 'c'),
            (rd(self.t_d, 7), rd(self.v_c, 2), rd(self.v_1, 2), rd(self.x_d, 0), 'd')
        )

    @property
    def area(self):
        """Calculate the distance x as the area of the profile. Ought to
        always match .x """

        x_ad, t_ad = accel_acd(self.v_0, self.v_c, self.v_1, self.joint.a_max)

        t_c = round(self.t - t_ad,8) # Avoid very small negatives

        if t_c < 0:
            raise TrapMathError(f'Negative t_c ({t_c}) ' )

        x_c = self.v_c * t_c

        if round(x_c) < 0:
            raise TrapMathError(f'Negative x_c ({x_c}) t_c={t_c}' )

        return x_ad + x_c

    @property
    def arear(self):
        """Robust area calculation. May return wrong results for bad parameters, but
        will not throw exceptions"""

        x_ad, t_ad = accel_acd(self.v_0, self.v_c, self.v_1, self.joint.a_max)

        t_c = max(self.t - t_ad, 0)
        x_c = max(self.v_c, 0) * t_c

        return x_ad + x_c

    def init(self):

        a_max = self.joint.a_max
        v_max = self.joint.v_max

        # Limit the boundary values for small moves
        self.v_0, self.v_1 = set_bv(self.x, self.v_0, self.v_1,self.joint.a_max)

        self.v_0 = min(self.v_0, v_max)
        self.v_1 = min(self.v_1, v_max)

        if self.x == 0:
            self.v_c = self.v_0 = self.v_1 = 0
            self.flag = 'Z'

        elif self.x < self.joint.small_x:
            # The limit is the same one used for set_bv.
            # the equation here is the sympy solution to:
            # t_a = (v_c - v_0) / a
            # t_d = (v_1 - v_c) / -a
            # x_a = ((v_0 + v_c) / 2) * t_a
            # x_d = ((v_c + v_1) / 2) * t_d
            # x_c = x - (x_a + x_d)
            # solve(x_c, v_c)[1]

            self.v_c = (sqrt(4 * a_max * self.x + 2 * self.v_0 ** 2 + 2 * self.v_1 ** 2) / 2)
            self.v_c = min(self.v_c, v_max) # formula errs for short segments
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
        assert abs(self.area-self.x) <1, (self.x, self.area, self.t, self.v_0, self.v_1)
        return self

    def plan_ramp(self, t):
        """Calculate a ramp profile for a fixed time"""

        self.t = t

        if self.x == 0 or self.t == 0:
            self.set_zero()
            self.t_c = self.t = t
            self.flag = 'Z'  # 5% of test cases
            return self

        def err(v_c):
            """Robust area calculation. May return wrong results for bad parameters, but
            will not throw exceptions"""

            x_a, t_a = accel_xt(self.v_0, self.v_c, self.joint.a_max)

            t_c = max(self.t - t_a, 0)
            x_c = max(v_c, 0) * t_c

            return self.x - (x_a + x_c)

        self.v_1 = self.v_c = binary_search(err, 0, self.x / self.t, self.joint.v_max)

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, self.joint.a_max)
        self.x_c = self.x - self.x_a
        self.t_c = self.t - self.t_a
        self.t_d, self.x_d = 0, 0

        self.t = self.t_a + self.t_c

        self.flag = "PR"

        return self

    def set_zero(self):

        self.flag = ''  # 5% of test cases
        self.x_a = self.x_d = self.x_c = 0
        self.t_a = self.t_d = self.t_c = 0
        self.v_0 = self.v_c = self.v_1 = 0

    def plan(self, t):

        self.t = t

        a_max = self.joint.a_max
        v_h = max(self.v_0, self.v_1)

        if self.x == 0 or self.t == 0:
            self.set_zero()
            self.t_c = self.t = t
            self.flag = 'Z'  # 5% of test cases
            return self

        elif self.v_0 == 0 and self.v_1 == 0:
            # should always be the lower root
            # Compare to equation 3.6 of Biagiotti
            # Sympy:
            # t_ad = v_c / a  # acel and decel are symmetric
            # x_ad = a * t_ad ** 2 / 2
            # x_c = x - 2 * x_ad
            # t_c = t - 2 * t_ad
            # solve(Eq(v_c, simplify(x_c / t_c)), v_c)
            try:
                self.v_c = a_max * self.t / 2 - sqrt(a_max * (a_max * self.t ** 2 - 4 * self.x)) / 2
            except ValueError:
                self.v_c = sqrt(self.x*a_max)

            self.flag = 'T'  # .16% of test cases

        elif self.x >= v_h * self.t and self.v_0 == self.v_1:
            # subtract off v_0, v_1 and the corresponding x and run again.
            x_ = self.v_0 * self.t

            # This is the same case as v_0 == v_1 == 0, but after subtracting off
            # the area of the rectangular base, from v =0 to v=v_0. Then we add back in
            # v_0.
            self.v_c = self.v_0 + a_max * self.t / 2 - sqrt(a_max * (a_max * self.t ** 2 - 4 * (self.x - x_))) / 2
            self.flag = 'T0'
        else:

            # Find v_c with a binary search, then patch it up if the
            # selection changes the segment time.
            self.v_c = binary_search(make_area_error_func(self), 0,
                                     self.x / self.t, self.joint.v_max)
            self.flag = 'O'

        self.x_a, self.t_a = accel_xt(self.v_0, self.v_c, a_max)
        self.x_d, self.t_d = accel_xt(self.v_c, self.v_1, a_max)

        # consistantize will make all of the values consistent with each other,
        # and if anthing is wrong, it will show up in a negative x_c

        if consistantize(self) < 0 or ( self.x > 25 and abs(round(self.area) - self.x)) > 1 :
            # We've probably got a really small move, and we've already tried
            # reducing boundary velocities, so get more aggressive.

            new_t = max(self.t_a + self.t_d, self.t)

            if self.v_1 > 0:
                self.v_1 = 0
                return self.plan(t)
            elif self.v_0 > 0:
                self.v_0 = 0
                return self.plan(t)
            elif new_t != self.t:
                return self.plan(new_t)
            else:
                raise ConvergenceError(
                    f'Unsolvable profile: x={self.x}, x_ad={self.x_a+self.x_d} t={self.t} t_ad={self.t_a + self.t_d}')


        assert round(self.x_a+self.x_d) <= self.x
        assert self.t_a + self.t_d <= self.t
        assert self.v_c >= 0, (self.v_c, self.flag)
        assert abs(self.area - self.x) <2, (self.area, self)
        assert self.t > 0

        return self
