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
    steps_left: int = 0
    periods_left: int = 0

    done: bool = False
    blocks: "List[tuple]" = []

    def __init__(self, blocks, period=DEFAULT_PERIOD, delay_counter=0, details=False):
        """Return segment parameters given the initial velocity, final velocity, and distance traveled. """

        self.blocks = blocks
        self.period = period
        self.delay_inc = period / TIMEBASE

        self.delay_counter=delay_counter

        self.t = 0
        self.phase_t = 0
        self.phase = 0
        self.steps_left = 0
        self.steps_stepped = 0
        self.calc_x = 0
        self.details = details
        self.v = 0
        self.step_inc = 1

        self.total_steps = abs(sum([e[0] for e in blocks]))

    def init_next_block(self):

        x, vi, vf = self.blocks[self.phase]

        self.direction = sign(x)

        self.x = abs(x)
        self.vi = vi
        self.vf = vf
        self.delay_inc = self.period / TIMEBASE
        self.step_inc = 1

        self.t_f = abs((2. * self.x) / (vi + vf))
        #self.a = (vf ** 2 - vi ** 2) / (2 * x) if x != 0 else 0
        self.a = (self.vf-self.vi)/ self.t_f if self.t_f != 0 else 0

        self.steps_left = int(round(self.x))
        self.steps_stepped = 0
        self.phase_t = 0

        self.calc_x = 0
        self.delay = 0

        # Slow start. Without this, we get a step early in the phase,
        # which results in very high instantaneous velocity.

        self.v = self.a * self.delay_inc + self.vi
        self.delay = abs(1 / self.v) if self.v else 0
        self.delay_counter += self.delay_inc

        self.periods_left = int(round(self.t_f / self.delay_inc))

    def __iter__(self):
        if self.details:
            while True:
                s = self.next_details()
                if s['isdone']:
                    return
                yield s

        else:
            while True:
                s = next(self)
                yield s
                if s == -2:
                    return

    def __next__(self):

        if self.steps_left == 0 or self.periods_left == 0:
            if self.phase == 3:
                self.done = True
                self.delay_counter += self.delay_inc
                return 0
            else:
                self.init_next_block()
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

        self.delay = abs( 1 / self.v) if self.v else 0
        self.delay_counter += self.delay_inc

        self.t += self.delay_inc
        self.phase_t += self.delay_inc

        self.calc_x = abs((self.a * self.phase_t ** 2) / 2 + self.vi * self.phase_t)
        self.x_err = self.steps_stepped - self.calc_x

        if abs(self.x_err) > .5 and self.phase == 1:
            s = self.x_err / abs(self.x_err)
            self.delay_counter += -s * self.delay_inc

        return r

    def next_details(self):
        step = next(self)
        d = dict(s=step, dr=self.direction, t=self.t, pt=self.phase_t, tf=self.t_f, v=self.v,
                 a=self.a,
                 #vi=self.vi, vf=self.vf,
                 sl=self.steps_left, pl=self.periods_left,
                 sg=None, ph=self.phase,
                 dl=self.delay, dc=self.delay_counter,
                  xc=self.calc_x, xe=self.x_err,
                 #isdone=self.done
                 )

        return d
