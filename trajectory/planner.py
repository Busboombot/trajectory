#
"""

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
from typing import List
from warnings import warn

import pandas as pd

from . import ConvergenceError
from .gsolver import Joint, ACDBlock


def index_clip(n, l):
    """Clip the indexer to an list to a valid range"""
    return max(-len(l), min(n, len(l) - 1))


class Segment(object):
    """One segment, for all joints"""

    n: int = None
    next: "Segment" = None
    prior: "Segment" = None
    blocks: "List[ACDBlock]" = None
    t: float = 0

    replans: int = 0

    def __init__(self, n, blocks: List[ACDBlock]):

        self.n = n
        self.blocks = blocks

        for b in self.blocks:
            b.segment = self

    def init(self):
        """Re-initialize all of the blocks"""
        for b in self:
            b.init()

        return self

    def plan(self, style='T', v_0=None, v_1=None):

        self.replans += 1

        assert self.replans < 50

        bv_none = [None] * len(self.blocks)

        if v_0 is not None and not hasattr(v_0, "__len__"):
            v_0 = [v_0 for _ in self.blocks]
        elif v_0 is None:
            v_0 = bv_none

        if v_1 is not None and not hasattr(v_1, "__len__"):
            v_1 = [v_1 for _ in self.blocks]
        elif v_1 is None:
            v_1 = bv_none

        for b, v_0, v_1 in zip(self.blocks, v_0, v_1):
            b.set_bv(v_0=v_0, v_1=v_1)

        for i in range(6):

            mt = max(self.times)

            for b in self:
                if style == 'R':
                    b.plan_ramp(mt)
                else:
                    b.plan(mt)

            if self.times_e_rms < 0.001:
                break
        else:
            if self.times_e_rms < 0.003:
                warn(f"Failed to converge sid={self.n} with err > 0.001, but < 0.003")
            else:
                warn(f'Failed to converge for sid={self.n}, rms={self.times_e_rms}, times={self.times}')

    @property
    def v_0(self):
        return [b.v_0 for b in self.blocks]

    @property
    def v_1(self):
        return [b.v_1 for b in self.blocks]

    @property
    def times(self):
        return [round(js.t, 6) for js in self.blocks]

    @property
    def time(self):
        return max(self.times)

    @property
    def times_e_rms(self):
        """Compute the RMS difference of the times from the mean time"""
        import numpy as np

        # Ignore really small axis for this calculation, because
        # they are very troublesome, and errors in them wont make a
        # big difference
        times = [round(b.t, 6) for b in self.blocks if b.x > 100]

        if not times:
            times = [0]

        s = sum(times)
        m = s / len(times)
        return np.sqrt(np.sum([(t - m) ** 2 for t in times]))

    def __iter__(self):
        return iter(self.blocks)

    def __getitem__(self, item):
        return self.blocks[item]

    def __str__(self):

        return f"{self.t:>3.4f}|" + ' '.join(js.str() for js in self.blocks)

    def _repr_pretty_(self, p, cycle):
        p.text(str(self) if not cycle else '...')

    def stepper(self):
        """Run steppers for all of the blocks in this segment"""
        from .stepper import DEFAULT_PERIOD, TIMEBASE

        steppers = [b.iter_steps(self.time, DEFAULT_PERIOD) for b in self.blocks]

        t = 0
        end_time = self.time
        time_step = DEFAULT_PERIOD / TIMEBASE

        for steps in zip(*steppers):
            yield (round(t, 6),) + steps
            t += time_step
            if t > end_time:
                break

    def plot(self, ax=None):
        from .plot import plot_params
        plot_params(self.blocks, ax=ax)


class SegmentList(object):
    segments: List[Segment]

    def __init__(self, joints: List[Joint]):

        self.joints = [Joint(j.v_max, j.a_max, i) for i, j in enumerate(joints)]

        self.directions = [0] * len(self.joints)

        for i, j in enumerate(self.joints):
            j.n = i

        self.segments = deque()

    def move(self, x: List[int], update=True, planner=None):
        """Add a new segment, with joints expressing joint distance
        :type x: object
        """

        assert len(x) == len(self.joints)

        blocks = [j.new_block(x=x, v_0=0, v_1=0) for j, x in zip(self.joints, x, )]

        prior = self.segments[-1] if len(self.segments) > 0 else None

        s = Segment(len(self.segments), blocks)

        s.init()

        self.plan(s, prior)

        if prior:
            prior.next = s
            s.prior = prior

        self.segments.append(s)

        return s

    def plan(self, s, prior, _recurse=True):

        if s.replans > 8:
            raise ConvergenceError(f"Too many replans in sid={s.n} ")

        if prior:

            prior.plan('R')
            bv = self.boundary_velocities(prior, s, prior.v_1)
            prior.plan(v_1=bv)
            s.plan(v_0=prior.v_1)

            if s.v_0 != prior.v_1:
                # This means that s could not handle the commanded v_0,
                # so prior will have to yield.
                prior.plan(v_1=s.v_0)

        else:
            s.plan()

        # Mop up remaining discontinuities
        if _recurse:
            for i in range(4):
                dis = self.discontinuities(-2)

                if len(dis) > 0:
                    for i, (p, c) in enumerate(dis):
                        self.plan(c.segment, p.segment, False)
                else:
                    break

    def boundary_velocities(self, p, c, target_v=None):
        from .gsolver import set_bv, same_sign

        if target_v is not None and not hasattr(target_v, "__len__"):
            target_v = [target_v for _ in self.joints]
        elif target_v is None:
            target_v = [j.v_max for j in self.joints]

        bv = []
        for pb, cb, v in zip(p.blocks, c.blocks, target_v):
            _, pv = set_bv(pb.x, pb.v_0, v, pb.joint.a_max)
            cv, _ = set_bv(pb.x, v, pb.v_0, pb.joint.a_max)

            if not same_sign(pb.d, cb.d) or pb.x == 0 or cb.x == 0:
                pv = cv = 0

            bv.append(min(pv, cv))

        return bv

    def discontinuities(self, index=0):
        """Yield segment pairs with velocity discontinuities"""

        index = index_clip(index, list(self.segments))

        o = []
        for s in list(self.segments)[index:]:
            if not s.next:
                continue

            for c,n in zip(s.blocks, s.next.blocks):
                if abs(c.v_1 - n.v_0) > 2:
                    o.append((c, n))
        return o

    @property
    def dataframe(self):

        frames = []
        for s_i, s in enumerate(self.segments):
            for a_i, b in enumerate(s.blocks):
                frames.append(b.dataframe.assign(axis=a_i, seg=s_i))

        df = pd.concat(frames)

        df['t'] = df.groupby('axis').del_t.cumsum()
        df['calc_x'] = (df.v_i + df.v_f) / 2 * df.del_t
        df['err'] = df.x - df.calc_x

        return df

    def __getitem__(self, item):
        """Return a joint segment by the id"""
        try:
            s, js = item
            return self.segments[s].blocks[js]
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

    @property
    def time(self):
        return sum([s.time for s in self.segments])

    def stepper(self):
        """Step through **and consume** the segments, producing steps. The SegmentList
        will be empty when this generator completes. """

        from .stepper import DEFAULT_PERIOD, TIMEBASE

        while len(self.segments) > 0:
            steppers = [b.iter_steps(self.time, DEFAULT_PERIOD) for b in self.segments[0].blocks]
            t = 0
            end_time = self.segments[0].time
            time_step = DEFAULT_PERIOD / TIMEBASE
            for steps in zip(*steppers):
                yield (round(t, 6),) + steps
                t += time_step
                if t > end_time:
                    break

            self.segments.pop()

    def plot(self, ax=None, axis=None):
        from .plot import plot_trajectory

        if len(self.segments) == 0:
            return

        df = self.dataframe

        if axis is not None:
            df = df[df.axis == axis]

        plot_trajectory(df, ax=ax)
