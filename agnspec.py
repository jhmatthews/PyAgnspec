#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
import dispar3_sub as dispar 


mass = 1e9
mdot = 5
angm = 0
alpha = 0.01
rout = 0.
tmin = 9000.
deltar=0.1
nre=0 

# ----------------------
# Basic input parameters
# ----------------------
#
#     MSTAR      - M(star), either in M(Sun), or in grams;
#     MDOT       - M(dot), either in M(Sun)/year; or in g/s
#     AA         - angular momentum of black hole
#     ALPHA      - Shakura-Sunyaev viscosity parameter
#     NRE        - number of radial annuli

# run dispar3 
dispar.dispar3 (mass, mdot, angm, alpha, rout, tmin, deltar, nre)


