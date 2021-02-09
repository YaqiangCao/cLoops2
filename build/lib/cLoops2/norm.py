#!/usr/bin/env python
#--coding:utf-8--
"""
norm.py
cLoops2 data normalization module.
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import warnings
warnings.filterwarnings("ignore")
import os
import json
import random
from glob import glob

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm

#cLoops2
from cLoops2.io import parseIxy
from cLoops2.cmat import getObsMat, getExpMat
from cLoops2.utils import isTool, callSys
from cLoops2.settings import *

