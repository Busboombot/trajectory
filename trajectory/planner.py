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
import math
from collections import deque
from typing import List

import pandas as pd

from .gsolver import Joint, ACDBlock, same_sign

def index_clip(n, l):
    """Clip the indexer to an list to a valid range"""
    return max(-len(l), min(n, len(l) - 1))


class Segment(object):
    """One segment, for all joints"""

    n: int = None
    t: float = 0
    blocks: "List[ACDBlock]" = None
    prior: "Segment" = None

    replans: int = 0

    def __init__(self, n, blocks: List[ACDBlock], prior: "Segment" = None):

        self.n = n
        self.blocks = blocks

        for b in self.blocks:
            b.segment = self

        if prior is not None:
            self.prior = prior

    def init(self):
        """Re-initialize all of the blocks"""
        for b in self:
            b.init()

        return self

    def set_bv(self, v_0=None, v_1=None):
        if hasattr(v_0, "__len__"):
            for b, v in zip(self.blocks, v_0):
                b.v_0 = v
        elif v_0 is not None:
            for b in self.blocks:
                b.v_0 = v_0

        if hasattr(v_1, "__len__"):
            for b, v in zip(self.blocks, v_1):
                b.v_1 = v
        elif v_1 is not None:
            for b in self.blocks:
                b.v_1 = v_1

        for b in self.blocks:
            b.limit_bv()

    def plan(self):

        mt = max(self.times)

        for b in self:
            b.plan(mt)

        return self.times_e_rms;

    def plan_ramp(self, v_0=None):

        for b in self.blocks:
            b.v_1 = b.joint.v_max

        mt = max(self.times)

        for b in self:
            b.plan_ramp(mt)

        return self.times_e_rms;

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
    replans: List[int] = None

    def __init__(self, joints: List[Joint]):

        self.joints = [Joint(j.v_max, j.a_max, i) for i, j in enumerate(joints)]

        self.directions = [0] * len(self.joints)

        for i, j in enumerate(self.joints):
            j.n = i

        self.segments = deque()

        self.replans = []

    def move(self, x: List[int], update=True, planner=None):
        """Add a new segment, with joints expressing joint distance
        :type x: object
        """

        assert len(x) == len(self.joints)

        blocks = [j.new_block(x=x, v_0=0, v_1=0) for j, x in zip(self.joints, x, )]

        prior = self.segments[-1] if len(self.segments) > 0 else None

        s = Segment(len(self.segments), blocks, prior)

        s.init()

        self.segments.append(s)

        self.plan()

        return s

    def plan_at_boundary(self, a, b):
        """Plan two segments, the current and immediately prior, to
        get a consistent set of boundary velocities between them """

        a_v0 = a.v_0


        # If the current block has a different sing -- changes direction --
        # then the boundary velocity must be zero.
        for ab, bb in zip(a.blocks, b.blocks):
            if not same_sign(ab.d, bb.d) or ab.x == 0 or bb.x == 0:
                ab.v_1 = 0

        a.plan()  # Plan a first
        b.set_bv(v_0=a.v_1)
        b.plan()  # Plan b with maybe changed velocities from a

        if b.v_0 != a.v_1:
            # This means that the current could not handle the commanded v_0,
            # so prior will have to yield.
            a.set_bv(v_1=b.v_0)
            a.plan()

        be = self.boundary_error(a, b)

        if a_v0 != a.v_0 or be > 1:
            # Planning segment a changed it's v_0, so we need to o back one earlier
            return -1
        else:
            return 1

    def plan(self):
        """Plan the newest segment, and the one prior. If there are boundary
        velocity discontinuities, also plan earlier segments"""

        i = len(self.segments) - 1
        current = self.segments[i]

        if i == 0:
            current.plan()
            return

        prior = self.segments[i - 1]
        prior.plan_ramp()  # Uncap the v_1 for the prior tail

        for _ in range(20):
            # For random moves, only about 10% of the segments will have more than 2
            # iterations, and 1% more than 10.

            if i >= len(self.segments):
                break

            current = self.segments[i]
            prior = self.segments[i - 1]

            inc_index = self.plan_at_boundary(prior, current)

            i += inc_index

            if i == 0:
                # Hit the beginning, so need to go forward
                i += 2

        return

    def boundary_error(self, p, c):
        sq_err = 0

        for pb, cb in zip(p.blocks, c.blocks):
            sq_err += (pb.v_1 - cb.v_0) ** 2

        return math.sqrt(sq_err)

    @property
    def pairs(self):

        for c, n in zip(self.segments, list(self.segments)[1:]):
            yield c, n

    def discontinuities(self, index=0):
        """Yield segment pairs with velocity discontinuities"""

        index = index_clip(index, list(self.segments))

        o = []
        for c, n in list(self.pairs)[index:]:
            for cb, nb in zip(c.blocks, n.blocks):
                if abs(cb.v_1 - nb.v_0) > 2:
                    o.append((cb, nb))
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
