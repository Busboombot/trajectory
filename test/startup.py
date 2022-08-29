
from trajectory.trapmath import *
from trajectory.profiles import *
from trajectory.params import *
from trajectory.plot import *
from trajectory import ParameterError
from itertools import product 
from more_itertools import chunked
from operator import attrgetter
from copy import copy
import pandas as pd
import numpy as np
import pickle
import seaborn as sns

from IPython.display import display
import pandas as pd

from math import sqrt
from random import randint, random, choices

pd.set_option('display.max_columns', None)

a_max = 50_000
v_max = 5_000

from functools import cache

@cache
def limits():
    """Return a list of (x,v_0, v_1) tuples for use in testing. """
    distances =  [0] + list(range(50,2000, 25))
    velocities = list(range(0,5000, 50))+[249,251, 4999, v_max]
    
    return list(enumerate(product(distances, velocities, velocities)))
