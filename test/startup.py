
from trajectory.planner import *
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
from tqdm.auto import tqdm

from IPython.display import display, HTML, Markdown
import pandas as pd

from math import sqrt
from random import randint, random, choices

pd.set_option('display.max_columns', None)

a_max = 50_000
v_max = 5_000

from functools import cache

distances =  [0,10,20,40,45, 50, 100, 200, 249,250,251, 300,400, 450, 499,500,501, 1000, 2000, 5000] 
velocities = [0,1,50]+ list(range(250,v_max, 250))+[v_max-1, v_max]

@cache
def limits():
    """Return a list of (x,v_0, v_1) tuples for use in testing. """

    
    return list(enumerate(product(distances, velocities, velocities)))