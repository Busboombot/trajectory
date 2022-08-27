
from dataclasses import dataclass, asdict, replace
from enum import Enum


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



@dataclass
class InputParams:
    x: float
    v_0: float
    v_1: float
    v_max: float
    a_max: float


@dataclass
class Params:
    x: float
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
    v_max: float = 0
    a_max: float = 0
    t_min: float = 0
    d: int = 0  # direction, -1 or 1
    ip: InputParams = None
    flag: str = None
    recalcs: int = 0
    jsclass: JSClass = JSClass.UNK
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

    def store(self):
        self.saved = [round(self.t, 4), round(self.v_0), round(self.v_c), round(self.v_1)]

    @property
    def is_changed(self):

        if self.saved is None:
            return False

        for p in zip(self.saved, [round(self.t, 4), round(self.v_0), round(self.v_c), round(self.v_1)]):

            if p[0] != p[1]:
                return True

        return False

    @property
    def v_0_changed(self):

        if self.saved is None:
            return False

        return self.saved[1] != self.v_0

    @property
    def changes(self):

        if self.saved is None:
            return [None] * len(self.saved)

        return [p[0] != p[1]
                for p in zip(self.saved, [round(self.t, 4), round(self.v_0), round(self.v_c), round(self.v_1)])]

