import matplotlib.pyplot as plt
import numpy as np
import uproot4 as uproot
from rootlib import *
from array import array
from scipy.signal import find_peaks
import sys
from libs import *
from filters import *                                    

# p.e. height
## For run 39 only
pe_area = 164.87 # mV.ns
pe_height = 6.981 # mV

# Test number of arguments passed in from terminal
if (len(sys.argv) != 3):
    print("Usage: python mppc-analysis.py _output.root_ _01_")
    sys.exit(0)

# ROOT file and branch
output_file = sys.argv[1]
dtbranch = "dt5743_wave"+sys.argv[2]

# Digitizer specs
mppc_DIGITIZER_FULL_SCALE_RANGE = 2.5 # Vpp
mppc_DIGITIZER_RESOLUTION       = 12 # bits
mppc_DIGITIZER_SAMPLE_RATE      = 3200000000 # S/s
digiCounts = 2.**mppc_DIGITIZER_RESOLUTION
nSamples = int(1024)

# Filter parameter
nbins = 14

# Scale factors
verticalScaleFactor = 1.e3*mppc_DIGITIZER_FULL_SCALE_RANGE/digiCounts; # V/bank
horizontalScaleFactor = 1.e9/mppc_DIGITIZER_SAMPLE_RATE #ns/sample

# Find run_number from input file
run_number = findIndex(output_file)

waveforms = []

# Use uproot to retrieve waveforms
with uproot.open(output_file) as fin:

    key_counter = 0
    for key in fin.keys():
        tree = fin[key]
        table   = tree.arrays([dtbranch], library="np")

        if key_counter == 0:
            waveforms = table[dtbranch]
            key_counter += 1
        else:
            waveforms = np.concatenate((waveforms, table[dtbranch]))
        
# Find baseline over 1% of total waveforms.
# Use 20 initial samples of waveform.
# Convert baseline to mV
print("Finding baseline for this set of waveforms...")
baseline, baseline_error = findBaseline_uproot(waveforms, 0.01, 20)
baseline = baseline*verticalScaleFactor
baselineError = baseline_error*verticalScaleFactor
print("Baseline: {0:.2f} mV".format(baseline))
print("Baseline error: {0} mV".format(baselineError))

# Convert units
waveforms = waveforms*verticalScaleFactor-baseline

# Number of entries
nEntries = len(waveforms)

# Sample space
x = np.asarray(range(int(0), nSamples))

# Convert to ns and mV
x_ns = x*horizontalScaleFactor

# Store information
nPulses, pulseIndices, pulseHeights, pulseWidths, pulseAreas = [], [], [], [], []
globalHeightsError, pulseWidthsError, pulseAreasError = [], [], []
# globalIndices, globalHeights = [], []

global_imax = np.array([], "int")
global_ymax = np.array([], "float")
pulseAreas = np.array([], "float")
pulseAreasError = np.array([], "float")

print("No. of events in tree: ", nEntries)

for i in range(0, nEntries):

    print("Finding maxima and areas. Progress: {0:4.2f}%.".format((i+1.)/nEntries*100.), end="\r")

    # Find global maxima
    imax = np.argmax(waveforms[i])
    ymax = waveforms[i][imax]
    global_imax = np.append(global_imax, [imax], axis=None)
    global_ymax = np.append(global_ymax, [ymax], axis=None)

    # Find area and area uncertainty with numpy 
    A, B = find_width(waveforms[i], 4., horizontalScaleFactor) # 4 mV
    if ((A == None) or (B == None)):
        A = 40
        B = 300
    area = np.sum(horizontalScaleFactor*waveforms[i][A:B])
    deltaArea = np.asarray([horizontalScaleFactor*y for y in waveforms[i][A:B]])
    sumAreaErrorSquared = np.sum((1./deltaArea**2)*( (np.sqrt((0.05*waveforms[i][A:B])**2+baselineError**2)/waveforms[i][A:B])**2 + 1. ) )
    pulseAreas = np.append(pulseAreas, [area], axis=None)
    pulseAreasError = np.append(pulseAreasError, [np.sqrt(sumAreaErrorSquared)], axis=None)

print("Finding maxima and areas. Progress: 100.00%.")

# Output file name
output_name = "./uproot_test.root"

# Save new data to ROOT file
saveGlobalMaximaToROOT_uprootinfo(run_number, dtbranch, baseline, baselineError, horizontalScaleFactor, verticalScaleFactor, global_imax, global_ymax, pulseAreas, pulseAreasError)