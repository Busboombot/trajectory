#
"""

THIS IS THE NEW VERSION OF THE PLANNER. The old one is in segments.py

Create a list of movements for a set of axes, then calculate
an optimum trapezoidal velocity profile.

Each movement is a segment, a period of time in which all axes
must start and stop at the same time. The segments have three periods:

    a: Acceleration
    c: Constant velocity
    d: Deceleration

( the names for the A and D periods refer to an idealized trapezoidal profile,
where a is the first part of the profile and d is the last. in actual cases, a
may have a deceleration, and d an acceleration. The c period is always 0
acceleration )

The optimization procedure attempts to execute every segment in the minimum
time, by ensuring that the longest axis is accelerate to the max velocity
as the max acceleration, and that there is a minimal change in velocity between
adjacent segments.

The algorithms for SimSegment and the c and cn parameters are based on
the article, 'Generate stepper-motor speed profiles in real time',
from https://www.embedded.com/generate-stepper-motor-speed-profiles-in-real-time/

------

Joint segment shapes



"""
from collections import deque

import pandas as pd

from .exceptions import ValidationError
from.params import *
from typing import List
from dataclasses import dataclass

## Parameters for simulation/step generation

# NUmber of ticks of the step function per second
TIMEBASE = 1_000_000  # ticks per second

# Shapes
TRAPEZOID = 1
TRIANGLE = 2

N_BIG = 2 ** 32

MAX_ERR_T = 1e-4


@dataclass
class SubSegment:
    """A sub segment is a portion of the trajectory with a constant
    acceleration,  one of the aceleration, constant (cruise) or decleration portions. """

    id: str
    segment_n: int
    joint: "Joint"
    js: "JointSegment"
    t: float
    v_i: float
    v_f: float
    x: float
    ss: float
    v_0_max: float
    v_1_max: float
    direction: int = 1

    def __post_init__(self):
        assert self.v_i == 0 or self.v_f == 0 or sign(self.v_i) == sign(self.v_f), \
            f"Inconsistent directions {self.v_i} {self.v_f} for {self.id}{self.ss} "

    def set_direction(self, sign) -> None:
        self.direction = sign
        self.v_i *= sign
        self.v_f *= sign
        self.x *= sign

    def __repr__(self):
        return f"<{self.joint.n} {self.ss} {self.t:2.5f} {int(self.x):5d} {self.v_i:5.0f}->{self.v_f:5.0f}> "

    @property
    def row(self):
        return [self.segment_n, self.joint.n, self.x, self.v_i, self.v_f, self.ss, self.t, self.v_1_max, self.v_0_max]


class JointSegment(object):
    """One segment for one joint"""

    joint: Joint
    segment: 'Segment'
    next_js: 'JointSegment'
    prior_js: 'JointSegment'

    b: Block

    _param_names = ['x', 't', 'dir', 'v_0_max', 'v_0', 'x_a', 't_a', 'x_c', 't_c', 'v_c_max', 'v_c', 'x_d', 't_d',
                    'v_1',
                    'v_1_max']

    def __init__(self, joint: Joint = None, x: float = None, v_0: float = 0, v_1: float = 0):
        """

        :param joint:
        :type joint:
        :param x:
        :type x:
        :param init_v0: # Initial segment velocity, for case when a new segment is being added to existing queue
        :type init_v0:
        """

        self.joint = joint
        self.segment = None  # set after adding to a segment

        self.next_js = None
        self.prior_js = None

        self.b = self.joint.new_block(x, v_0, v_1).init()

    @property
    def id(self):
        return f"{self.segment.n}/{self.joint.n}"

    @property
    def err_x(self):
        """Return true of the subsegment error is larger than 3% of the x"""
        if self.b.x != 0:
            return 0
        else:
            return 0

    @property
    def subsegments(self):
        return [SubSegment(self.id, self.segment.n, self.joint, self, *ss, self.b.v_0_max, self.b.v_1_max)
                for ss in self.b.subsegments]

    def __str__(self):
        from colors import color, bold

        v0 = color(f"{int(self.b.v_0):<5d}", fg='green')
        xa = bold(f"{int(self.b.d * self.b.x_a):>6d}")
        c = bold(f"{int(round(self.b.d * self.b.x_c)):>6d}")
        vc = color(f"{int(self.b.v_c):<6d}", fg='blue')
        xd = bold(f"{int(self.b.d * self.b.x_d):<6d}")
        v1 = color(f"{int(self.b.v_1):>5d}", fg='red')

        return f"[{v0} {xa}↗{c + '@' + vc}↘{xd} {v1}]"

    def __repr__(self):
        return self.__str__()

    def _repr_pretty_(self, p, cycle):
        p.text(str(self) if not cycle else '...')

    @property
    def debug(self):

        from colors import color, bold

        v0_max = color(f"{int(self.b.v_0_max):<5d}", fg='green')
        v0 = color(f"{int(self.b.v_0):<5d}", fg='green')
        xa = bold(f"{int(self.b.d * self.b.x_a):>6d}")
        ta = bold(f"{int(self.b.d * self.b.t_a):<1.5f}")
        c = bold(f"{int(round(self.b.d * self.b.x_c))}")
        vc = color(f"{int(self.b.v_c)}", fg='blue')
        tc = bold(f"{int(self.b.d * self.b.t_c):<1.5f}")
        xd = bold(f"{int(self.b.d * self.b.x_d):>6d}")
        td = bold(f"{int(self.b.d * self.b.t_d):<1.5f}")
        v1 = color(f"{int(self.b.v_1):>5d}", fg='red')
        v1_max = color(f"{int(self.b.v_1_max):>5d}", fg='red')

        return f"[{v0_max}|{v0} {xa}%{ta} ↗ {c}@{vc}%{tc} ↘ {xd}%{td} {v1}|{v1_max}]"

class Segment(object):
    """One segment, for all joints"""

    n: int = None
    next_seg: "Segment" = None
    prior_seg: "Segment" = None
    joint_segments: "List[JointSegment]" = None
    t: float = 0

    n_updates: int = 0

    updates: List[str]

    def __init__(self, n, joint_segments: List[JointSegment], prior=None):

        self.n = n
        self.joint_segments = joint_segments

        for js in self.joint_segments:
            js.segment = self

        self.t = max(js.b.t_min for js in self.joint_segments)

        if prior:
            Segment.link(prior,self)

    @classmethod
    def link(cls, prior: "Segment", current: "Segment"):

        prior.next_seg = current
        current.prior_seg = prior

        for p, c in zip(prior.joint_segments, current.joint_segments):
            p.next_js = c
            c.prior_js = p

    @property
    def final_velocities(self):
        return [j.b.v_1 for j in self.joint_segments]

    @property
    def boundary_velocities(self):
        """Get the boundary velocities for the v_0's of this segment
        and the v_1's of the prior
        :return:
        :rtype:
        """
        if self.prior_seg:
            p = self.prior_seg.joint_segments
        else:
            p = [None for _ in self.joint_segments]

        return [ (p.b.v_1, c.b.v_0) for p, c in zip(p, self.joint_segments) ]

    @property
    def params(self):
        return [j.b for j in self.joint_segments]

    @property
    def err_t(self):
        """RMS error between joint segment calculated times and segment time"""
        err = sqrt(sum([(js.t - self.t) ** 2 for js in self.joint_segments]))

        return err / len(self.joint_segments) / self.t

    @property
    def err_x(self):
        """RMS error between joint segment calculated times and segment time"""
        dist = sum(js.x for js in self.joint_segments)
        err = sqrt(sum([(js.err_x) ** 2 for js in self.joint_segments]))

        return err / dist

    @property
    def max_err_x(self):
        return max([js.err_x for js in self])

    @property
    def times(self):
        return [self.t] + [js.t for js in self.joint_segments]

    @property
    def min_t(self):
        return max([js.t for js in self.joint_segments])

    @property
    def params_df(self):

        from operator import attrgetter
        ag = attrgetter(*JointSegment._param_names)

        columns = ['seg', 'js', 'seg_t'] + JointSegment._param_names + ['calc_x', 'sum_x', 'calc_t']

        rows = []
        for j, js in enumerate(self.joint_segments):
            r = (self.n, j, js.segment.t) + ag(js) + (js.calc_x, js.sum_x, js.calc_t)
            assert len(columns) == len(r)
            d = dict(zip(columns, r))
            rows.append(d)

        df = pd.DataFrame(rows)

        df['sum_x'] = df.x_a + df.x_c + df.x_d
        df['err_x'] = abs(df.calc_x - df.sum_x).round(4)

        for c in ['v_c', 'v_0_max', 'v_1_max', 'v_0', 'v_1', 'calc_x', 'sum_x']:
            df[c] = df[c].round(0).astype('int')

        for c in ['seg_t', 't', 't_a', 't_d', 't_c']:
            df[c] = df[c].round(5)

        return df

    def save(self):
        from operator import itemgetter
        keys = ['x', 'dir', 't',
                'v_0_max', 'v_1_max', 'v_c_max',
                'x_a', 't_a', 'x_d', 't_d', 'v_0', 'v_c', 'x_c', 't_c', 'v_1']

        ag = itemgetter(*keys)
        return [dict(zip(keys, ag(js.__dict__))) for js in self]

    def load(self, d):

        for js, e in zip(self.joint_segments, d):
            for k, v in e.items():
                setattr(js, k, v)

    @property
    def subsegments(self):
        for js in self.joint_segments:
            yield from js.subsegments

    def __iter__(self):
        return iter(self.joint_segments)

    def __getitem__(self, item):
        return self.joint_segments[item]

    def __str__(self):

        return f"{self.t:>3.4f}|" + ' '.join(str(js) for js in self.joint_segments) + f" {round(self.max_err_x, 2)}"

    def _repr_pretty_(self, p, cycle):
        p.text(str(self) if not cycle else '...')


class SegmentList(object):
    positions: List[float]  # Positions after last movement addition
    segments: deque[Segment]
    all_segments: List[Segment]  # All of the segments, unprocessed

    def __init__(self, joints: List[Joint]):

        self.joints = [Joint(j.v_max, j.a_max, i) for i, j in enumerate(joints)]

        self.directions = [0] * len(self.joints)

        for i, j in enumerate(self.joints):
            j.n = i

        self.segments = deque()
        self.positions = [0] * len(self.joints)

    def new_segment(self, joint_distances, prior_velocities):

        js = [JointSegment(j, x=x, v_0=v_0) for j, x, v_0
              in zip(self.joints, joint_distances, prior_velocities)]

        prior = self.segments[-1] if len(self.segments) > 0 else None

        s = Segment(len(self.segments), js, prior)

        return s

    def rmove(self, joint_distances: List[int], update=True):
        """Add a new segment, with joints expressing joint distance
        :type joint_distances: object
        """

        self.uncap_tail()

        assert len(joint_distances) == len(self.joints)

        if len(self.segments) > 0:
            last = self.segments[-1]

            prior_velocities = last.final_velocities
        else:
            prior_velocities = [0] * len(joint_distances)

        s = self.new_segment(joint_distances, prior_velocities)

        if len(self.segments) > 0:
            Segment.link(self.segments[-1], s)

        self.segments.append(s)

        self.update(-1)

        return s

    def uncap_tail(self):
        if len(self.segments) == 0:
            return

        for j in self.segments[-1]:
            j.b.v_1 = j.b.joint.v_max
            j.b.init()

        self.plan(-1)

    def update(self, index):

        idx = len(self.segments)-1
        for _ in range(6):
            print("Updating", idx)
            u = self.update_boundary(idx)
            if u > 0:
                print("u plan")
                self.plan(idx)

            idx -= 1
            u = self.update_boundary(idx)
            if u == 0:
                break;
            else:
                self.plan(idx)

    def plan(self, index):
        s = self.segments[index]

        times = [j.b.t for j in s]

        for i in range(4):
            mt = max(times)

            for j in s:
                j.b.plan(mt)

            times = [j.b.t for j in s]

            if all ([t == times[0] for t in times]):
                break
        else:
            raise ConvergenceError('Too many planning updates')


    def update_boundary(self, index):

        s = self.segments[index]
        if not s.prior_seg:
            return 0
        else:
            p = s.prior_seg

        updates = 0
        for b1, b2 in zip(p, s):
            b1 = b1.b
            b2 = b2.b

            if  b1.v_1 != b2.v_0 or b1.v_1 != b1.v_c:
                v = min(b1.v_1, b2.v_0)
                v_m = (b1.v_c + b2.v_c) / 2

                b1.v_1 = b1.v_c
                b2.v_0 = b1.v_c

                updates += 1

            if updates > 0:
                b1.init()
                b2.init()

        return updates


    @property
    def windows(self):
        """Yield all windows of parameters"""
        for i in range(len(self) - 1):
            yield self.get_window(i)

    @property
    def segment_pairs(self):
        """Yield all valid adjacent segment"""

        for s in self.segments:
            if s.next_seg is not None:
                yield (s, s.next_seg)

    @property
    def joint_pairs(self):
        """Yield all valid adjacent segment"""

        for s in self.segments:
            for j in s.joint_segments:
                if j.next_js:
                    yield (j, j.next_js)

    def get_window(self, i: int):

        prior = []
        current = []
        nxt = []

        for js in self.segments[i]:
            prior.append(js.prior_js.b if js.prior_js else None)
            current.append(js.b)
            nxt.append(js.next_js.b if js.next_js else None)

        return prior, current, nxt

    def _lim_segments(self, limit=4):
        if limit:
            return list(self.segments)[-limit:]
        else:
            return list(self.segments)

    def err_t(self, limit=4):

        l = self._lim_segments(limit)

        et = sum([s.err_t for s in l]) / len(l)

        if et is None:
            return -1

        return et

    def err_x(self, limit=4):
        l = self._lim_segments(limit)
        return sum([js.subseg_error for s in l for js in s])

    def count_needs_update(self):
        return sum(s.needs_update for s in self.segments)

    def validate(self, limit=3):
        """Check that the c error in each of the segments is less than X%"""
        for s in self:
            for js in s:
                if js.b.x != 0:
                    rel_error = abs(js.err_x / js.b.x)
                    if js.err_x > limit:
                        raise ValidationError(f"{js.id} Excessive error. err_x={js.err_x} re={rel_error}")

    @property
    def dataframe(self):
        from warnings import warn

        rows = []

        try:
            self.validate()
        except ValidationError as e:
            warn(str(e))

        for ss in self.subsegments:
            rows.append([None, *ss.row])

        df = pd.DataFrame(rows, columns=' t seg axis x v_i v_f ss del_t v0m v1m'.split())

        df['t'] = df.groupby('axis').del_t.cumsum()
        df['calc_x'] = (df.v_i + df.v_f) / 2 * df.del_t
        df['err'] = df.x - df.calc_x

        return df

    @property
    def total_error(self):
        """Sum of the segment errors, from the dataframe"""
        return self.dataframe.err.abs().sum()

    @property
    def total_re(self):
        """Sum of the segment errors, from the dataframe, divided by total distance, of all axes"""
        t = self.dataframe
        return t.err.abs().sum() / t.x.sum()

    @property
    def discontinuities(self):
        """Yield segment pairs with velocity discontinuities"""
        for c, n in self.joint_pairs:
            if round(c.b.v_1, 2) != round(n.b.v_0, 2):
                yield c, n

    def has_discontinuity(self, s1, s2):

        for js1, js2 in zip(s1, s2):
            if js1.b.v_1 != js2.b.v_0:
                return True

    @property
    def dataframe_stacked(self):
        return self.dataframe.set_index(['t', 'axis']).stack().unstack(-2).unstack(-1)

    @property
    def params_df(self):
        """
        :return: a dataframe of paramaters
        :rtype:
        """

        frames = []
        for i, s in enumerate(self.segments):
            frames.append(s.params_df)

        df = pd.concat(frames)

        return df

    @property
    def subsegments(self):

        for s in self.segments:
            yield from s.subsegments

    def __getitem__(self, item):
        """Return a joint segment by the id"""
        try:
            s, js = item
            return self.segments[s].joint_segments[js]
        except TypeError:
            return self.segments[item]

    def __iter__(self):
        return iter(self.segments)

    def __len__(self):
        return len(self.segments)

    def __str__(self):
        return '\n'.join(str(s) for s in self.segments)

    def _repr_pretty_(self, p, cycle):
        p.text(str(self) if not cycle else '...')

    def debug_str(self):
        return '\n'.join(s.debug_str() for s in self.segments)

    def plot(self, ax=None, axis = None):
        from .plot import plot_segment_list

        df = self.dataframe

        if axis is not None:
            df = df[df.axis == axis]

        plot_segment_list(df, ax=ax)


__all__ = ['SegmentList', 'Segment', 'Joint', 'JointSegment']
