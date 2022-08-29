
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


a_max = 50_000
v_max = 5_000