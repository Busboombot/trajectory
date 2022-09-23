from collections import namedtuple
from math import sqrt

StepParams = namedtuple('StepParams', 'n tn  cn ca1 xn vn dir tf x')

TIMEBASE = 1_000_000  # ticks per second
DEFAULT_PERIOD = 4


def sign(x):
    if x == 0:
        return 0
    elif x > 0:
        return 1
    else:
        return -1


class Stepper(object):
    """A step segment for simulating the step interval algorithm. """

    steps_left: int = 0
    direction: int = 0
    x: int = 0

    t: float = 0
    t_f: int = 0
    delay: float = 0
    delay_counter: float = 0
    periods_left: int = 0

    done: bool = False

    def __init__(self, period=DEFAULT_PERIOD, details=False):
        """Return segment parameters given the initial velocity, final velocity, and distance traveled. """

        self.period = period
        self.delay_inc = period / TIMEBASE
        self.details = details
        self.delay_inc = self.period / TIMEBASE

        self.t = 0
        self.t0 = self.t
        self.phase_t = 0
        self.phase = 0
        self.phases_left = 0
        self.steps_left = 0
        self.steps_stepped = 0
        self.nst = 0

        self.done = False

        self.iter = 0

    def load_phases(self, phases):
        """ Phases is a list of tuples, (x, vi, vf)
        :param phases:
        :type phases:
        :return:
        :rtype:
        """
        self.phases = phases
        self.phase = 0
        self.phases_left = len(phases)
        self.done = False

        return self

    def init_next_phase(self):

        self.t0 += self.phase_t

        x, self.vi, self.vf = self.phases[self.phase]

        self.direction = sign(x)
        self.x = abs(x)
        self.steps_left = int(round(self.x))

        self.steps_stepped = 0
        self.phase_t = 0

        self.phases_left -= 1
        self.phase += 1
        self.done = False

        self.nst = self.next_step_time(1) if self.steps_left > 0 else 0

    def next_step_time(self, step):

        # If they are the same, a==0, we get divide by zero, so
        # make them slightly different and it works.
        if self.vi == self.vf:
            vi = self.vi + .1
            vf = self.vf - .1
        else:
            vi = self.vi
            vf = self.vf

        t = 2 * self.x / (vi + vf)
        a = (vf - vi) / t if t != 0 else 0

        try:
            return abs((-self.vi + sqrt(2 * a * step + self.vi ** 2)) / a)
        except ValueError:
            return abs( self.vi  / a)

    def next(self):

        if self.done:
            return 0

        if self.steps_left <= 0:

            if self.phases_left != 0:
                self.init_next_phase()
            else:
                self.done = True
                return 0

        if self.phase_t >= self.nst:

            self.steps_left -= 1
            self.steps_stepped += 1

            r = self.direction

            if self.steps_left > 0:
                self.nst = self.next_step_time(self.steps_stepped + 1)

        else:
            r = 0

        self.t += self.delay_inc
        self.phase_t += self.delay_inc

        return r

    def next_details(self):
        step = self.next()
        d = dict(s=step, dr=self.direction, t=self.t, pt=self.phase_t, tf=self.t_f, v=self.v,
                 a=self.a,
                 # vi=self.vi, vf=self.vf,
                 sl=self.steps_left, pl=self.periods_left,
                 sg=None, ph=self.phase,
                 dl=self.delay, dc=self.delay_counter,
                 xc=self.calc_x,
                 # xe=self.x_err,
                 # isdone=self.done
                 )

        return d

    def __iter__(self):
        while not self.done:
            if self.details:
                yield self.next_details()
            else:
                yield self.next()


class SegmentStepper:

    def __init__(self, sl: "SegmentList", step_if: list[object] = None, details: bool = False):
        from .stepper import DEFAULT_PERIOD, TIMEBASE, Stepper

        self.sl = sl
        self.step_if = step_if  # Step interface
        self.details = details
        self.t = 0

        self.steppers = [Stepper(details=self.details) for _ in self.sl.joints]

        self.step_position = [0] * len(self.sl.joints)

        period = DEFAULT_PERIOD
        self.dt = period / TIMEBASE

        self.seg = None

    def step(self):
        import numpy as np

        if len(self.sl.segments) == 0:
            return None

        # Load the next segment
        if self.seg is None:

            self.seg = self.sl.front

            for b, stp in zip(self.seg.blocks, self.steppers):
                stp.load_phases(b.stepper_blocks())

        # Get the next set of steps
        if self.details:
            r = self.steppers[0].next_details()
            r['sg'] = self.seg.n
            r['qt'] = self.sl.queue_time
            r['ql'] = self.sl.queue_length

            done = self.steppers[0].done

        else:
            steps = [stp.next() for stp in self.steppers]
            self.step_position += np.array(steps)
            r = [self.t] + steps
            self.t += self.dt

            done = self.steppers[0].done and all([stp.done for stp in self.steppers])

        # When done, invalidate the segment and remove it.
        if done:
            self.sl.pop()
            self.seg = None

        return r

    @property
    def done(self):
        return [stp.done for stp in self.steppers]

    def __next__(self):
        s = self.step()
        if s is None:
            raise StopIteration
        else:
            return s

    def __iter__(self):
        return self
