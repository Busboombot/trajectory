
from dataclasses import dataclass, asdict, replace
from enum import Enum

from .exceptions import *

def accel_acd(v_0, v_c, v_1, a):
    """ Same result as running accel_xt, for v_0->v_c and v_c->v_1,
    and adding the x and t values.
    """
    t_ad = (abs(v_c - v_0) + abs(v_c - v_1)) / a
    x_ad = abs((v_0 ** 2 - v_c ** 2) / (2 * a)) + abs((v_1 ** 2 - v_c ** 2) / (2 * a))

    return x_ad, t_ad

@dataclass
class Joint:
    v_max: float
    a_max: float

    def new_move(self, x, v_0=None, v_1=None):

        v_0 = v_0 if v_0 is not None else self.v_max
        v_1 = v_1 if v_1 is not None else self.v_max

        return Move(x, v_0, v_1, self)

    def new_block(self, x, v_0=None, v_1=None):
        return self.new_move(x,v_0, v_1).new_block()

@dataclass
class Move:
    x: float
    v_0: float
    v_1: float
    joint: Joint

    def new_block(self):
        return Block(x=self.x, v_0=self.v_0, v_1=self.v_1,
                     move=self, joint=self.joint)

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
    v_0_max: float = None
    v_0_min: float = 0
    v_1_max: float = None
    v_1_min: float = 0
    t_min: float = 0
    d: int = 0  # direction, -1 or 1
    joint: Joint = None
    move: Move = None
    flag: str = None
    recalcs: int = 0
    jsclass: "JSClass" = None
    saved = [None] * 4  # Initial save forces status as having changed.

    def __post_init__(self):
        from .trapmath import sign
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
        t_c = self.t - t_ad
        x_c = self.v_c * t_c

        if round(t_c, 4) < 0:
            raise TrapMathError(f'Negative t_c' + str((t_c, t_ad,
                    (self.t, self.v_0, self.v_c, self.v_1, self.joint.a_max))))

        return x_ad + x_c


    def plan(self, t=None):
        from .profiles import plan_min_time

        if t is None:
            plan_min_time(self)
        else:
            assert False



class JSClass(Enum):
    """

    === Convex

    ╱▔▔ U Acel. V0 fixed, V1 is maxed
          Should be a typical case.

    ╱╲  A Triangle. Fixed V0, V1, too short for vc
          x is non zero but is too short got velocity to reach v_max

    ╱▔╲ T Trapezoid.  V0 & V1 == 0, max vc
          Normal case where v0 and v1 are both below Vc, and the segment isn't short.

    ╱▔╲ P Pentagon.  Fixed V0, V1, but not zero, max vc
          Normal case where v0 and v1 are both below Vc, and the segment isn't short.

    ▔▔╲ L Clif. Next segment is very small x.
          Next segment is very small x or zero



    === Flat

    ▔▔▔ C Constant V. Prior and next have high V0, V1
          V0 and V1 are both at vmax, so v_c should be v_max

    ▁▁▁ Z Zero distance
          Zero x means that there must bt v=0 through whole segment

    === Concave

    ╲▁╱ V Trough. Like a pentagon, but Vc is below either V0 or V1

    ╲▁▁ R Ramp down. Very small x relative to other joints.
          X is so short relative to other joints that we hit v=0
          during accel phase.


    """
    ACEL = 'U'
    TRIANGLE = 'A'
    PENTAGON = 'P'
    HEXAGON = 'H'
    AT = 'AT'  # A or T, pinned to Zero at start and finish
    TRAPZEZOID = 'T'
    TROUGH = 'V'
    CONSTANT = 'C'
    RAMP = 'R'
    DECEL = 'D'
    CLIFF = 'L'
    DL = 'DL'  # D or L
    ZERO = 'Z'
    XXXXX = 'X'
    UNK = '?'


kind_icon_map = {

    JSClass.TRIANGLE: '╱╲',
    JSClass.TRAPZEZOID: '╱-╲',
    JSClass.PENTAGON: '╱▔╲',
    JSClass.HEXAGON: '╱▔╲',
    JSClass.TROUGH: '╲▁╱',
    JSClass.CONSTANT: '▔▔▔ ',
    JSClass.ZERO: '▁▁▁',
    JSClass.ACEL: '╱▔▔',
    JSClass.DECEL: '╲',
    JSClass.RAMP: '╲▁▁',
    JSClass.CLIFF: '▔▔╲',
    JSClass.XXXXX: '▁▁╱',
    JSClass.UNK: '?'
}


def classify(p):
    """Assign a class type to a planner block"""
    from operator import attrgetter
    ag = attrgetter(*'x v_0 v_c t_c v_1 v_max a_max'.split())

    x, v_0, v_c, t_c, v_1, v_max, a_max = ag(p)

    if x == 0:
        return JSClass.ZERO
    elif v_c == v_1 and v_1 == v_0:
        return JSClass.CONSTANT

    elif v_0 == 0 and v_1 == 0:
        if t_c == 0:
            return JSClass.TRIANGLE
        else:
            return JSClass.TRAPZEZOID

    elif v_0 != 0 and v_1 == 0 and v_c == 0:
        return JSClass.DECEL
    elif v_0 == v_c and v_1 > v_0:
        return JSClass.RAMP  # Probably should have a new one, it's reversed from other ram
    elif v_0 == v_c and v_1 < v_0:
        return JSClass.CLIFF
    elif v_1 == v_c and v_0 > v_1:
        return JSClass.RAMP
    elif v_1 == v_c and v_0 < v_1:
        return JSClass.ACEL

    elif (v_0 != 0 or v_1 != 0) and v_0 < v_max and v_1 < v_max:
        if t_c == 0:
            return JSClass.PENTAGON
        elif v_c > max(v_0, v_1):
            return JSClass.HEXAGON
        elif v_c < min(v_0, v_1):
            return JSClass.TROUGH
        else:
            return JSClass.UNK
    elif v_c < v_0 and v_c < v_1:
        return JSClass.TROUGH


    else:
        return JSClass.UNK

