from itertools import count
from math import sqrt
from collections import namedtuple

StepParams = namedtuple('StepParams', 'n tn  cn ca1 xn vn dir tf x')


TIMEBASE = 1_000_000  # ticks per second

def sign(x):
    if x == 0:
        return 0
    elif x > 0:
        return 1
    else:
        return -1

class SimSegment(object):
    """A step segment for simulating the step interval algorithm. """

    fp_bits = 11  # Number of bits for the fixed point fraction

    def __init__(self, v0, v1, x=None, t=None, t0=0, timebase=None):
        """Return segment parameters given the initial velocity, final velocity, and distance traveled. """

        # Each segment must have a consistent, single-valued accceleration
        if abs(sign(v0) + sign(v1)) == 0 and (abs(v0) + abs(v1)) > 0:  # both signs are nonzero and opposite
            raise Exception("The velocity trajectory cannot cross zero. Got {}->{}".format(v0, v1))

        assert t != 0.0

        self.dir = sign(v0 + v1)
        self.x = None
        self.v0 = abs(v0)
        self.v1 = abs(v1)  # Velocity at last step
        self.t = None

        if x is None:
            self.x = int(round(.5 * (v1 + v0) * t, 0))
            self.t = t

        elif t is None:
            self.x = abs(x)
            if self.v1 + self.v0 == 0:
                self.t = 0
                self.a = 0
            else:
                self.t = abs(2. * float(self.x) / float(self.v1 + self.v0))
                self.a = (self.v1 - self.v0) / self.t

        self.x = int(round(self.x if self.x is not None else x, 0)) ;

        if x == 0 or t == 0:
            self.a = 0
            self.tf = 0  # time at last step
            self.tn = 0  # Total transit time for step n
            self.vn = 0  # running velocity
            self.xn = 0  # running position. Done when xn = x
            self.v1 = 0
            self.v0 = 0
        else:
            self.tf = t  # time at last step
            self.tn = t0  # Total transit time for step n
            self.vn = abs(v0) * self.dir  # running velocity
            self.xn = 0  # running position. Done when xn = x

        self.timebase = timebase if timebase is not None else TIMEBASE

        self.constV = False

        self.sim_time = 0;  # Time to completion in simulation iteration/

        assert self.t is not None

    def initial_params(self):
        """Set initial parameters, which are a bit different for the first step,
        to account for low torque. """

        # If n is positive, there is a non zero velocity and we are accelerating
        # If it is negative, there is a non zero velocity and we are decelerating

        v0 = float(self.v0)
        v1 = float(self.v1)

        if v0 == 0 and v1 == 0:
            # Going nowhere.
            n = 0
            cn = 0
            self.constV = True

        elif v0 == 0:
            # Starting from a stop, so need special
            # first delay.
            a = abs(v1) / self.t

            n = 0  # n will always be positive, so accelerating
            cn = 0.676 * sqrt(2.0 / a) * self.timebase  # c0 in Equation 15

            assert cn > 0

        elif v0 == v1:
            # Constant velocity.
            n = 0
            cn = self.timebase / abs(v0)
            self.constV = True

        else:
            # Normal case, between two non-zero velocities
            a = abs(v1 - v0) / self.t

            n = (v0 * v0) / (2.0 * a)  # Equation 16
            cn = self.timebase / abs(v0)

            assert n >= 0

            # Need to put the sign back on n; n must be negative for deceleration
            if abs(v1) < abs(v0):
                n = -n

        # Fixed point version of cn
        ca = int(cn * (1 << SimSegment.fp_bits))
        n = round(n)
        self.ca1 = 0

        return n, cn, ca

    @staticmethod
    def next_params(n, cn):

        n += 1

        ca1 = (2.0 * cn) / ((4.0 * n) + 1)

        cn = cn - ca1  # Equation 13

        return n, cn, ca1

    @staticmethod
    def next_params_xp(n, ca):
        """Fixed point version of the next_params function"""

        n += 1

        if n != 0:
            ca1 = (ca << 1) / ((n << 2) + 1);

            ca = ca - ca1,  # Equation 13

        return n, ca, ca1

    def next_delay(self):
        """Call this after each step to compute the delay to the next step"""

        self.n, self.cn, self.ca1 = self.next_params(self.n, self.cn)

        self.xn += 1
        self.tn += abs(self.cn)
        try:
            self.vn = abs(TIMEBASE / self.cn) * self.dir
        except ZeroDivisionError:
            self.vn = 0

        return self.n, self.cn, self.ca1

    row_header = 'n tn cn xn vn dir tf x'.split()

    def __iter__(self):

        self.n, self.cn, self.ca = self.initial_params()

        while True:
            yield StepParams(self.n, round(self.tn / TIMEBASE, 6), round(self.cn, 2), round(self.ca1, 2), self.xn,
                             int(round(self.vn, 0)), self.dir, self.tf, self.x)
            self.next_delay()
            if self.x - self.xn <= 0:
                break

    def iter_discrete(self, period):
        """Iterate discretely, in units of time of the period"""

        n, cn, ca = self.initial_params()

        t = 0
        stepsLeft = self.x
        lastTime = 0

        self.sim_time = 0

        for i in count():

            lastTime = lastTime + (period << SimSegment.fp_bits);

            self.sim_time += period / self.timebase

            if lastTime > ca:
                lastTime -= ca
                if not self.constV:
                    n += 1;
                    ca1 = int(((2 * ca) / ((4 * n) + 1)))
                    ca = abs(ca - ca1)

                stepsLeft -= 1;
                yield 1
            else:
                yield 0

            if stepsLeft == 0:
                return i

    def iter_time(self, period):

        t = 0
        steps_left = self.x

        delay_counter = 0
        v = self.v0
        delay = 0

        delay_inc = period / self.timebase

        for i in count():

            if delay_counter >= delay:
                delay_counter -= delay

                steps_left -= 1
                yield (i, round(t, 6) * self.timebase, steps_left, v, delay_counter, delay)

            v = self.a * t + self.v0
            delay = 1 / v if v else 0

            if steps_left == 0:
                break;

            delay_counter += delay_inc
            t += delay_inc

        self.sim_time = t

    @property
    def dataframe(self):
        import pandas as pd

        def _y():
            for e in self:
                yield self.n, self.tn / TIMEBASE, self.cn, self.xn, self.vn, self.dir, self.tf, self.x

        return pd.DataFrame(list(_y()), columns='n tn cn xn vn dir tf x'.split())

    # def __repr__(self):
    #    return '{:<10.6} t={:<8.4f} x={:<8.0f} n={:<6d} cn={:<8.4f} xn={:<6.0f} vn={:<6.0f}' \
    #        .format(self.tn / TIMEBASE, self.tf, self.x, self.n, self.cn, self.xn, self.vn)