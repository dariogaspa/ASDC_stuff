#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 13:24:42 2021

@author: dariogasparrini
"""
import string
import os
#import pyfits
from astropy.io import fits #replace pyfits
from astropy.io import ascii
from numpy import * 
import sys
from math import *
#import numarray
import math
from numpy import *
import numpy as np
from matplotlib import pylab
from matplotlib import pyplot
from pylab import load
from scipy import optimize

f=ascii.read('/Users/dariogasparrini/Documents/Likelihood_ratio/4LAC-DR3/RASS/reliability.qdp')
    
# load input variables from a file
x_values = f.columns[0]
y_values = f.columns[1]


# objective function
def objective(x, p0, p1):
	return (1-p0*(1/exp(x * p1)))

pyplot.scatter(x_values, y_values)
startpoint = 20
xfit=x_values[startpoint:]
yfit=y_values[startpoint:]
print "start at ", x_values[startpoint]
pyplot.scatter(xfit, yfit)

# fit curve
popt, pcov = optimize.curve_fit(objective, xfit, yfit)
a, b = popt
print a, b
print pcov
# define a sequence of inputs between the smallest and largest known inputs
x_line = arange(min(xfit)-0.5, max(xfit), 1)
# calculate the output for the range
y_line = objective(x_line, a, b)
print objective(3, a, b)
# create a line plot for the mapping function
pyplot.plot(x_line, y_line, '--', color='red', label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
pyplot.legend()
pyplot.show()