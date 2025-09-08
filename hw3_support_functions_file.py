"""
This module contains support functions for Homework 3.
Includes square() and squareplot() for mathematical operations.
"""
import numpy as np

#Remember that all libraries needed from your functions should be
#imported at the start:

import os # standard library imports first

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt  # NOTE: pyplot, *not* pylab!
import astropy.io.fits as fits


#writing function
def square(x):
    '''
    Returns the square of the input
    x can be int/float or a numpy array
    '''
    return np.square(x)

def squareplot(low, high, npts, saveplot=False):
    """
    Plots the square of numbers between the listed bounds
    low and high are float values
    npts is the number of points
    saveplot saves the pdf if a string is provided
    """
    x = np.linspace(low,high,npts)
    y = square(x)

    plt.figure()
    plt.plot(x,y)
    plt.xlabel("Input")
    plt.ylabel("Output")
    plt.title("Square")

    if saveplot:
        plt.savefig(saveplot, format="pdf")
    else:
        plt.show()
