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
    vi: float = 0
    done: bool = False

    def __init__(self, x, vi, vf, period=DEFAULT_PERIOD):
        """Return segment parameters given the initial velocity, final velocity, and distance traveled. """

        vi = abs(vi)
        vf = abs(vf)

        self.direction = sign(x)
        self.x = abs(x)

        self.vi = vi
        self.steps_left = self.x
        self.t = 0
        self.delay = 0
        self.delay_counter = 0
        self.delay_inc = period / TIMEBASE

        self.t_f = (2.*self.x)/(vi+vf) if (vi+vf) > 0 else 0# final time
        self.periods_left = round(self.t_f/self.delay_inc)

        self.a = (vf ** 2 - vi ** 2) / (2 * x) if x != 0 else 0

    def __iter__(self):
        t = 0
        while True:
            t += self.delay_inc

            if self.steps_left <=0:
                return

            yield next(self)

    def __next__(self):

        if self.steps_left <= 0:
            self.done = True

        if self.delay_counter > self.delay:
            # Using -=delay is has some variation over =0, but is more accurate.
            self.delay_counter -= self.delay
            self.steps_left -= 1
            r = self.direction
        else:
            r = 0

        v = self.a * self.t + self.vi
        self.delay = 1 / v if v else 0

        self.delay_counter += self.delay_inc
        self.t += self.delay_inc

        return r
