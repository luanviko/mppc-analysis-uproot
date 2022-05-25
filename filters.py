import ROOT
from array import array
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from numpy.fft import fft, ifft, fftfreq
import os, sys

def runningAverage(x, y, nbins):

    # Store filtered y
    yfiltered = []

    # Run through samples in waveform
    for i in range (nbins//2, len(y)-nbins//2):
        
        # Sum and average samples
        sumsamples = 0.
        for j in range(i-nbins//2, i+nbins//2):
            sumsamples += y[j]
        yfiltered.append(sumsamples/nbins)
    
    # Return array
    return np.asarray(yfiltered)

def convolution():
    pass