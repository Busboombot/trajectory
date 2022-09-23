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

import numpy as np
import pandas as pd

from .gsolver import Joint, Block, bent, mean_bv


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
    blocks: "List[Block]" = None
    joints: List[Joint] = None
    prior: "Segment" = None
    move: "List[int]" = None

    replans: int = 0

    def __init__(self, n, joints: List[Joint], move: List[int] = None, prior: "Segment" = None):

        self.n = n
        self.joints = joints

        if move is not None:
            if isinstance(move[0], Block):
                self.blocks = move
                self.move = [b.x for b in self.blocks]
            else:
                self.move = move
                if prior is not None:
                    self.prior = prior
                    self.blocks = [j.new_block(x=x, v_0=pb.v_1, v_1=0) for j, x, pb in
                                   zip(self.joints, move, prior.blocks)]
                else:
                    self.blocks = [j.new_block(x=x, v_0=0, v_1=0) for j, x in zip(self.joints, move)]

                for b in self.blocks:
                    b.segment = self

    def plan(self, v_0=None, v_1=None, prior=None, next_=None, t=None, iter=10):

        # Planning can change the time for a block, so planning multiple
        # will ( should ) converge on a singe segment time.

        largest_at = max([j.max_at for j in self.joints])
        lower_bound_time = largest_at * 2

        for p_iter in range(iter):  # Rarely more than 1 iteration
            if t is not None:
                mt = t
            elif p_iter < 2:
                mt = self.min_time
            elif p_iter < 4:
                mt = max(lower_bound_time, self.min_time)
            else:
                mt = max(lower_bound_time, self.time)

            for i, b in enumerate(self.blocks):
                pb = prior.blocks[i] if prior is not None else None
                nb = next_.blocks[i] if next_ is not None else None

                b.plan(mt, v_0, v_1, pb, nb, iter=p_iter)

            if self.times_e_rms < .001:
                break
            else:
                for b in self.blocks:
                    if b.t < mt:
                        b.limit_bv()

        self.t = self.time
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
    def times_e_rms(self):
        """Compute the RMS difference of the times from the mean time"""
        import numpy as np

        times = [b.t for b in self.blocks] if len(self.blocks) else [0]

        m = sum(times) / len(times)

        return np.sqrt(np.sum([(t - m) ** 2 for t in times]))

    @property
    def v_1(self):
        return [b.v_1 for b in self.blocks]

    def __iter__(self):
        return iter(self.blocks)

    def __getitem__(self, item):
        return self.blocks[item]

    def __str__(self):

        return f"{self.time:>3.4f}|" + ' '.join(js.str() for js in self.blocks)

    def _repr_pretty_(self, p, cycle):
        p.text(str(self) if not cycle else '...')

    @property
    def stepper_blocks(self):

        return [b.stepper_blocks() for b in self.blocks]

    def stepper(self, t0=0, period=None, timebase=None, details=False):
        """Run steppers for all of the blocks in this segment"""
        from .stepper import DEFAULT_PERIOD, TIMEBASE

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
    seg_num: int = 0
    planner_position: List[int] = None
    step_position: List[int] = None

    def __init__(self, joints: List[Joint]):

        self.joints = [Joint(j.v_max, j.a_max, i) for i, j in enumerate(joints)]

        self.directions = [0] * len(self.joints)

        for i, j in enumerate(self.joints):
            j.n = i

        self.segments = deque()

        self.planner_position = [0] * len(joints)
        self.distance = [0] * len(joints)
        self.step_position = np.array([0] * len(joints))
        self.seg_num = 0;
        self.replans = []

        self.queue_length = 0
        self.queue_time = 0

    def set_position(self, pos):
        assert len(pos) == len(self.joints)

        self.planner_position = pos

    def move(self, x: List[int], v_max: List[int] = None):
        """Add a new segment, with joints expressing joint distance
        :type x: object
        """

        for i, x_ in enumerate(x):
            self.planner_position[i] += x_
            self.distance[i] += abs(x_)

        assert len(x) == len(self.joints)

        prior = self.segments[-1] if len(self.segments) > 0 else None

        s = Segment(self.seg_num, self.joints, x, prior)
        self.seg_num += 1;

        if v_max is not None:
            for b, vm in zip(s.blocks, v_max):
                b.v_c_max = vm

        self.segments.append(s)

        if prior is None:
            s.plan(v_0=0, v_1=0)
        else:
            prior.plan(v_1='v_max', prior=prior.prior)
            s.plan(v_0='prior', v_1=0, prior=prior)
            self.plan(len(self.segments) - 1)

        self.queue_length +=1
        self.queue_time = sum(s.t for s in self.segments) # Because times will change with replanning


    def amove(self, x: List[int]):
        """Move to an absolute planner position"""

        rmove = [x_ - p_ for x_, p_ in zip(x, self.planner_position)]

        return self.move(rmove)

    def vmove(self, t, v: List[int]):
        """Move by running at a target velocity for a given time"""

        v_max = max(v)
        x_max = v_max * t

        move = [x_max * (v_ / v_max) for v_ in v]

        return self.move(move, v_max=v)

    def jmove(self, t, v: List[int]):
        """Like a vmove, but also removes all but the last two segments"""

        self.segments = deque(list(self.segments)[:2])
        self.queue_length = len(self.segments)

        return self.vmove(t, v)

    def plan(self, seg_idx: int = None):

        if seg_idx is None:
            seg_idx = len(self.segments) - 1

        for p_idx in range(15):

            current = self.segments[seg_idx]
            prior = self.segments[seg_idx - 1]
            pre_prior = self.segments[seg_idx - 2] if seg_idx >= 2 else None

            assert prior == current.prior, (seg_idx, prior.n, current.prior.n)

            prior.plan(v_1='next', prior=pre_prior, next_=current)  # Plan a first
            current.plan(v_0='prior', prior=prior)  # Plan b with maybe changed velocities from a

            # Smooth out boundary bumps between segments.
            def v_limit(p_idx, v_max):
                if p_idx < 2:
                    return v_max
                elif p_idx < 4:
                    return v_max / 2
                else:
                    return 0;

            bends = 0
            for pb, cb in zip(prior.blocks, current.blocks):
                if bent(pb, cb):
                    diff = abs(pb.v_1 - mean_bv(pb, cb))
                    if diff < v_limit(p_idx, pb.joint.v_max):
                        pb.v_1 = cb.v_0 = mean_bv(pb, cb)
                        bends += 1

            if bends or (pre_prior is not None and self.boundary_error(pre_prior, prior)):
                seg_idx += -1  # Run it again one segment earlier
            elif self.boundary_error(prior, current):
                # This means that the current could not handle the commanded v_0,
                # so prior will have to yield.
                seg_idx += 0  # Re-run planning on this boundary
            else:
                seg_idx += 1  # Advance to the next segment

            seg_idx = max(1, seg_idx)

            if seg_idx >= len(self.segments):
                break

        self.replans.append(p_idx)

    @property
    def edist(self):
        """Euclidean distances"""
        return [np.sqrt(np.sum(np.square(s.move))) for s in self.segments]

    def boundary_error(self, p, c):
        return math.sqrt(sum([(pb.v_1 - cb.v_0) ** 2 for pb, cb in zip(p.blocks, c.blocks)]))

    @property
    def pairs(self):

        for c, n in zip(self.segments, list(self.segments)[1:]):
            yield c, n

    @property
    def triplets(self):
        for p, c, n in zip(self.segments, list(self.segments)[1:], list(self.segments)[2:]):
            yield p, c, n

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

    @property
    def front(self):
        """Reference to the front of the segments"""
        return self.segments[0]

    def pop(self):
        """Remove the front of the segments"""
        s = self.segments.popleft()
        self.queue_time -= s.t
        self.queue_length -= 1

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
        t = 0
        o = ''
        for i, s in enumerate(self.segments):
            o += f'{i:3d} {t:>2.4f} ' + str(s) + "\n"
            t += s.time

        return o

    def _repr_pretty_(self, p, cycle):
        p.text(str(self) if not cycle else '...')

    def debug_str(self):
        return '\n'.join(s.debug_str() for s in self.segments)

    @property
    def time(self):
        return sum([s.time for s in self.segments])

    @property
    def stepper_blocks(self):

        for seg_n, s in enumerate(self.segments):
            yield seg_n, [b.stepper_blocks() for b in s.blocks]

    def step(self, details=False):
        from .stepper import DEFAULT_PERIOD, TIMEBASE, Stepper
        period = DEFAULT_PERIOD

        dt = period / TIMEBASE
        t = 0

        steppers = [Stepper(details=details) for _ in self.joints]

        for seg_n, stepper_blocks in self.stepper_blocks:

            for stp, sb in zip(steppers, stepper_blocks):
                stp.load_phases(sb)

            if details:
                while not steppers[0].done and not all([stp.done for stp in steppers]):
                    d = steppers[0].next_details()
                    d['sg'] = seg_n
                    yield d
            else:
                while not steppers[0].done and not all([stp.done for stp in steppers]):
                    steps = [stp.next() for stp in steppers]
                    self.step_position += np.array(steps)
                    yield [t] + steps
                    t += dt

    def plot(self, ax=None, axis=None):
        from .plot import plot_trajectory

        if len(self.segments) == 0:
            return

        df = self.dataframe

        if axis is not None:
            df = df[df.axis == axis]

        plot_trajectory(df, ax=ax)


