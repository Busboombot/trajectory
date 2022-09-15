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
import numpy as np

from .gsolver import Joint, ACDBlock, bent, mean_bv


def index_clip(n, l):
    """Clip the indexer to an list to a valid range"""
    n = max(-len(l), min(n, len(l) - 1))

    if n < 0:
        n = len(l) + n

    return n


class Segment(object):
    """One segment, for all joints"""

    n: int = None
    t: float = 0
    blocks: "List[ACDBlock]" = None
    joints: List[Joint] = None
    prior: "Segment" = None

    replans: int = 0

    def __init__(self, n, move: List[int], joints: List[Joint], prior: "Segment" = None):

        self.n = n
        self.joints = joints

        if prior is not None:
            self.prior = prior
            self.blocks = [j.new_block(x=x, v_0=pb.v_1, v_1=0) for j, x, pb in zip(self.joints, move, prior.blocks)]
        else:
            self.blocks = [j.new_block(x=x, v_0=0, v_1=0) for j, x in zip(self.joints, move)]

        for b in self.blocks:
            b.segment = self

    def plan(self, v_0=None, v_1=None, prior=None, next_=None):

        # Planning can change the time for a block, so planning multiple
        # will ( should ) converge on a singe segment time.

        for p_iter in range(10):  # Rarely more than 1 iteration
            if p_iter < 2:
                mt = max(0.04, self.min_time)
            else:
                mt = max(0.04, self.time)

            for i, b in enumerate(self.blocks):
                pb = prior.blocks[i] if prior is not None else None
                nb = next_.blocks[i] if next_ is not None else None

                b.plan(mt, v_0, v_1, pb, nb)

            if self.times_e_rms < .001:
                break

        return self

    def zero(self):
        for b in self.blocks:
            b.zero();

    @property
    def times(self):
        return [round(js.t, 6) for js in self.blocks]

    @property
    def time(self):
        return max(self.times)

    @property
    def min_time(self):
        """Calculated maximum of the  minimum time for each block in the segment"""
        return max([b.min_time() for b in self.blocks])

    @property
    def dominant(self):
        x = list(reversed(sorted([ (b.min_time(), i) for i, b in enumerate(self.blocks) ])))
        return x[0][1]


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

        return f"{self.time:>3.4f}|" + ' '.join(js.str() for js in self.blocks)

    def _repr_pretty_(self, p, cycle):
        p.text(str(self) if not cycle else '...')

    def stepper(self, t0=0, period=None, timebase=None, details=False):
        """Run steppers for all of the blocks in this segment"""
        from .stepper import DEFAULT_PERIOD, TIMEBASE
        from itertools import zip_longest

        period = period if period else DEFAULT_PERIOD
        timebase = timebase if timebase else TIMEBASE

        steppers = [b.iter_steps(self.time, period, details=details) for b in self.blocks]

        t = t0
        end_time = self.time
        time_step = period / timebase

        if details:
            # Only return one axis  because details have ultiple columns
            for steps in steppers[0]:
                yield (round(t, 6),) + steps
                t += time_step
                if t > end_time:
                    break

        else:

            for steps in zip(*steppers):
                yield (round(t, 6),) + steps
                t += time_step
                if t > end_time:
                    break

    @property
    def dataframe(self):
        frames = []
        for i, b in enumerate(self.blocks):
            frames.append(b.dataframe.assign(seg=self.n, axis=i))

        return pd.concat(frames)

    def plot(self, ax=None):
        from .plot import plot_trajectory
        plot_trajectory(self.dataframe, ax=ax)


class SegmentList(object):
    segments: List[Segment]
    replans: List[int] = None
    move_positions: List[int] = None
    step_positions: List[int] = None
    def __init__(self, joints: List[Joint]):

        self.joints = [Joint(j.v_max, j.a_max, i) for i, j in enumerate(joints)]

        self.directions = [0] * len(self.joints)

        for i, j in enumerate(self.joints):
            j.n = i

        self.segments = deque()

        self.move_positions = [0]*len(joints)
        self.step_positions = np.array([0] * len(joints))

        self.replans = []

    def move(self, x: List[int]):
        """Add a new segment, with joints expressing joint distance
        :type x: object
        """

        for i, x_ in enumerate(x):
            self.move_positions[i] += x_

        assert len(x) == len(self.joints)

        prior = self.segments[-1] if len(self.segments) > 0 else None

        s = Segment(len(self.segments), x, self.joints, prior)

        self.segments.append(s)

        if prior is None:
            s.plan(v_0=0, v_1=0)
        else:
            prior.plan(v_1='v_max')
            s.plan(v_0='prior', v_1=0, prior=prior)
            self.plan(len(self.segments) - 1)

        # Linear ( hopefully ) re-plan of the last few segments. This will
        # finish up straightening bumps and removing discontinuities.
        if False and len(self.segments) >= 2:
            n = index_clip(-5, self.segments)  # Replan at most last 4 elements all the way through
            n = max(1, n)  # Replnning the first one unmoors v_0
            self.plan(n)

    def plan(self, i: int = None, remove_bumps=False):

        if i is None:
            i = len(self.segments) - 1

        for p_idx in range(15):

            current = self.segments[i]
            prior = self.segments[i - 1]
            pre_prior = self.segments[i - 2] if i >= 2 else None

            prior.plan(v_1='next', next_=current)  # Plan a first
            current.plan(v_0='prior', prior=prior, )  # Plan b with maybe changed velocities from a

            # Smooth out boundary bumps between segments.
            def v_limit(p_idx, v_max):
                if p_idx == 0:
                    return v_max
                elif p_idx < 2:
                    return v_max / 4
                else:
                    return 0;

            bends = 0
            for p, n in zip(prior.blocks, current.blocks):
                if True and bent(p, n):
                    diff = abs(p.v_1 - mean_bv(p, n))
                    if diff < v_limit(p_idx, p.joint.v_max):
                        p.v_1 = n.v_0 = mean_bv(p, n)

                        bends += 1

            if bends or (pre_prior is not None and self.boundary_error(pre_prior, prior)):
                i += -1  # Run it again one segment earlier
            elif self.boundary_error(prior, current):
                # This means that the current could not handle the commanded v_0,
                # so prior will have to yield.
                i += 0  # Re-run planning on this boundary
            else:
                i += 1  # Advance to the next segment

            i = max(0, i)

            if i >= len(self.segments):
                break

        self.replans.append(p_idx)

    def boundary_error(self, p, c):
        return math.sqrt(sum([(pb.v_1 - cb.v_0) ** 2 for pb, cb in zip(p.blocks, c.blocks)]))

    @property
    def pairs(self):

        for c, n in zip(self.segments, list(self.segments)[1:]):
            yield c, n

    @property
    def blocks(self):
        for s in self.segments:
            for b in s.blocks:
                yield b

    @property
    def block_pairs(self):
        for p, n in self.pairs:
            for pb, nb in zip(p.blocks, n.blocks):
                yield (pb, nb)

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

    def step(self, details=False):
        from .stepper import DEFAULT_PERIOD, TIMEBASE
        period = DEFAULT_PERIOD

        dt = period/TIMEBASE
        t = 0
        t0 = 0
        delay_counter =[0]*len(self.joints)
        for seg_n, s in enumerate(self.segments):

            steppers = [b.stepper(details=details, delay_counter=dc) for dc, b in zip(delay_counter,s.blocks)]
            if details:
                # Can only do one axis with details.
                while True:
                    steps = [stp.next_details() for stp in steppers]
                    d = steps[0]
                    d['t'] = d['t'] + t0
                    d['sg'] = seg_n
                    yield d
                    if steppers[0].done and all([stp.done for stp in steppers]):
                        t0 = d['t']
                        delay_counter = [s.delay_counter for s in steppers]

                        break
            else:
                while True:
                    steps = [next(stp) for stp in steppers]
                    self.step_positions += np.array(steps)
                    yield [t]+steps
                    t += dt
                    if steppers[0].done and all([stp.done for stp in steppers]):
                        delay_counter = [s.delay_counter for s in steppers]
                        t0 = t
                        break

    def plot(self, ax=None, axis=None):
        from .plot import plot_trajectory

        if len(self.segments) == 0:
            return

        df = self.dataframe

        if axis is not None:
            df = df[df.axis == axis]

        plot_trajectory(df, ax=ax)
