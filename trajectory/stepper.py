from collections import namedtuple

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
        self.delay_counter = 0

        self.t = 0
        self.phase_t = 0
        self.phase = 0
        self.steps_left = 0
        self.steps_stepped = 0
        self.calc_x = 0
        self.details = details
        self.v = 0
        self.step_inc = 1
        self.done = False

    def load_phases(self, phases):
        """ Phases is a list of tuples, (x, vi, vf)
        :param phases:
        :type phases:
        :return:
        :rtype:
        """
        self.phases = phases
        self.phase = 0
        self.done = False

    def init_next_phase(self):

        x, vi, vf = self.phases[self.phase]

        self.direction = sign(x)

        self.x = abs(x)
        self.vi = vi
        self.vf = vf
        self.delay_inc = self.period / TIMEBASE
        self.step_inc = 1

        self.t_f = abs((2. * self.x) / (vi + vf)) if (vi + vf) != 0 else 0
        self.a = (self.vf - self.vi) / self.t_f if self.t_f != 0 else 0

        self.steps_left = int(round(self.x))
        self.steps_stepped = 0
        self.phase_t = 0

        self.calc_x = 0

        # Slow start. Without this, we get a step early in the phase,
        # which results in very high instantaneous velocity.

        self.v = self.a * self.delay_inc + self.vi
        self.delay = abs(1 / self.v) if self.v else 0
        self.delay_counter += self.delay_inc

        self.periods_left = int(round(self.t_f / self.delay_inc))

        self.done = False

    def next(self):

        if self.steps_left <= 0 or self.periods_left <= 0:
            if self.done or self.phase == 3:
                self.done = True
                return 0
            else:
                self.init_next_phase()
                self.phase += 1

        if self.delay_counter > self.delay:
            self.delay_counter -= self.delay
            self.steps_left -= self.step_inc
            self.steps_stepped += self.step_inc
            r = self.direction
        else:
            r = 0

        self.periods_left -= 1

        self.v = self.vi + self.a * self.phase_t

        self.delay = abs(1 / self.v) if self.v else 1
        self.delay_counter += self.delay_inc

        self.t += self.delay_inc
        self.phase_t += self.delay_inc

        self.calc_x = abs((self.a * self.phase_t ** 2) / 2 + self.vi * self.phase_t)
        self.x_err = self.steps_stepped - self.calc_x

        if abs(self.x_err) > .5 and self.phase == 1:
            s = self.x_err / abs(self.x_err)
            self.delay_counter += -s * self.delay_inc * .1

        return r

    def next_details(self):
        step = self.next()
        d = dict(s=step, dr=self.direction, t=self.t, pt=self.phase_t, tf=self.t_f, v=self.v,
                 a=self.a,
                 # vi=self.vi, vf=self.vf,
                 sl=self.steps_left, pl=self.periods_left,
                 sg=None, ph=self.phase,
                 dl=self.delay, dc=self.delay_counter,
                 xc=self.calc_x, xe=self.x_err,
                 # isdone=self.done
                 )

        return d

    def __iter__(self):
        if self.details:
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
        s =  self.step()
        if s is None:
            raise StopIteration
        else:
            return s

    def __iter__(self):
        return self
