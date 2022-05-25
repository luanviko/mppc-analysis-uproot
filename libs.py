import ROOT
from array import array
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from numpy.fft import fft, ifft, fftfreq
import os, sys

def find_width(waveform, height_threshold, horizontalScaleFactor):


    x1 = []
    x2 = []
    for k in range(0, len(waveform)-1 ):
        if ( (waveform[k] <= height_threshold) and (waveform[k+1] > height_threshold) ):
            x1.append(k)
        if ( (waveform[k] > height_threshold) and (waveform[k+1] <= height_threshold) and ((k+1)*horizontalScaleFactor < 300) ):
            x2.append(k)
    try:
        A = min(x1)
        B = max(x2)
    except:
        A = None 
        B = None

    return A, B


def findIndex(output_file):

    # Series of split and replace to get run_number from original file.
    name1 = output_file.split("/")
    name2 = name1[-1].replace(".root","")
    name3 = name2.replace("output","")
    run_number = name3.replace("0","")

    return run_number

def applyFilter(y):

    # Fourir filter
    Y = fft(y)
    sampling_rate = 3200000000
    
    f = fftfreq(1024, sampling_rate)

    x = range(0,1024)
    print(np.abs(Y)**2)
    plt.plot(x,y)
    # plt.plot(f,np.abs(Y)**2)
    plt.show()

def findBaseline(output_file, nWaveforms, sampleCap):

    ## Sum over the *sampleCap* first samples
    #  of the waveform to find baseline.
    #  Average this value over nWaveforms percentage
    #  of the total waveforms.

    # Open file and get tree
    file = ROOT.TFile(output_file, "READ")
    tree = file.Get("waveform_tree")
    
    # Array to store waveform information
    waveform = array('d', [0.0]*1024 )

    # Get branch and set address
    dataBranch = tree.GetBranch("dt5743_wave00")
    tree.SetBranchAddress("dt5743_wave00", waveform) 

    # Initialize baseline variable
    baseline = 0.
    sumBaselineErrorSquared = 0.

    # Average over a sampleWaveforms number.
    sampleWaveforms = int(nWaveforms*tree.GetEntries() )
    for i in range(0, sampleWaveforms ):

        # Get entries
        tree.GetEntry(i)

        # Initialize sum of samples as zero
        sumSamples = 0.
        sumSamplesErrorSquared = 0.

        # Run throught the sampleCap samples
        for j in range(0,sampleCap):
            sumSamples += waveform[j]
            sumSamplesErrorSquared += (0.05*waveform[j])**2
        waveformBaselineError = np.sqrt(sumSamplesErrorSquared)/sampleCap

        # Average sample out 
        baseline += sumSamples/sampleCap
        sumBaselineErrorSquared += waveformBaselineError**2

    # Close file
    file.Close()

    # And average over the number of sample waveforms
    return baseline/sampleWaveforms, np.sqrt(sumBaselineErrorSquared)/sampleWaveforms

def findBaseline_uproot(table, nWaveforms, sampleCap):

    ## Sum over the *sampleCap* first samples
    #  of the waveform to find baseline.
    #  Average this value over nWaveforms percentage
    #  of the total waveforms.

    # Initialize baseline variable
    baseline = 0.
    sumBaselineErrorSquared = 0.

    # Average over a sampleWaveforms number.
    sampleWaveforms = int(nWaveforms*len(table) )
    for i in range(0, sampleWaveforms ):

        # Get waveform from table
        waveform = table[i]

        # Initialize sum of samples as zero
        sumSamples = 0.
        sumSamplesErrorSquared = 0.

        # Run throught the sampleCap samples
        for j in range(0,sampleCap):
            sumSamples += waveform[j]
            sumSamplesErrorSquared += (0.05*waveform[j])**2
        waveformBaselineError = np.sqrt(sumSamplesErrorSquared)/sampleCap

        # Average sample out 
        baseline += sumSamples/sampleCap
        sumBaselineErrorSquared += waveformBaselineError**2

    # And average over the number of sample waveforms
    return baseline/sampleWaveforms, np.sqrt(sumBaselineErrorSquared)/sampleWaveforms


def checkFolder(yes, output_folder):

    # Do we want to save plots?
    if (yes == False):
        return False
    
    # If yes:
    elif (yes == True):

        # Check if folder already exists
        if os.path.exists(output_folder):
            return True 
        else:

            # Try to create the folder
            try:
                os.mkdir(output_folder)
                return True
            
            # Skip figure creation.
            except:
                print("Folder with figures could not be created. Saving figures will be skipped.")
                return False

def plotSimple(x, y, xlabel, ylabel, output_folder, N):

    ## Plot waveforms
    
    # Clear plot
    plt.cla()

    # Plot curve and pulses
    plt.plot(x, y, color='gray')

    # Set label
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    # Change range
    plt.xlim([x[0], x[-1]] )

    # Save file
    plt.savefig("{0}/waveform-nopeaks-{1}.png".format(output_folder, N))

def plotDouble(x1, y1, x2, y2, xlabel, ylabel, output_folder, N):

    plt.cla()

    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, sharey=True, constrained_layout=True, sharex=False, figsize=(8,4))

    ax0.set_ylabel(ylabel)
    ax0.set_xlabel(xlabel)
    ax1.set_xlabel(xlabel)

    ax0.set_xlim(x1[0], x1[-1])
    ax1.set_xlim(x1[0], x1[-1])

    ax0.plot(x1, y1, color='gray')
    ax1.plot(x2, y2, color='gray')

    ax0.set_title("Original pulse")
    ax1.set_title("Moving average")

    plt.savefig("{0}/waveform-nopeaks-{1}.png".format(output_folder, N)) 

def plotOver(x1, y1, x2, y2, nbins, pulses, xlabel, ylabel, output_folder, N):

    plt.cla()

    fig, (ax0) = plt.subplots(nrows=1, ncols=1, sharey=True, constrained_layout=True, sharex=False, figsize=(8,4))

    ax0.set_ylabel(ylabel)
    ax0.set_xlabel(xlabel)

    ax0.set_xlim(x1[0], x1[-1])
    ax0.plot(x1, y1, color='gray')
    ax0.plot(x2, y2, color='black')
    ax0.plot(x2[pulses], y2[pulses], 'o', mfc='red', mec='black')

    plt.title("Moving average ({0} samples)".format(nbins))

    plt.savefig("{0}/waveform-{1}.png".format(output_folder, N)) 

def plotWithPulses(x, y, peaks, N, xlabel, ylabel, output_folder):

    ## Plot waveforms with respective peaks
    
    # Clear plot
    plt.cla()

    # Plot curve and pulses
    plt.plot(x, y, color='gray')
    plt.plot(x[peaks], y[peaks], 'x', color='red')
    
    # Set label
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    # Change range
    plt.xlim([x[0], x[-1]] )

    # Save file
    plt.savefig("{0}/waveform-{1}.png".format(output_folder, N))

    

def plotWithWidthBar(x, y, peaks, height_threshold, N, xlabel, ylabel, output_folder):

    ## Plot waveforms with respective peaks
    
    # Clear plot
    plt.cla()

    # Plot curve and pulses
    plt.plot(x, y, color='gray')
    # plt.plot(x[peaks], y[peaks], 'x', color='red')
    plt.hlines(y=height_threshold, xmin=x[peaks[0]], xmax=x[peaks[-1]], linestyles='solid', label='{0}'.format(x[peaks[-1]]-x[peaks[0]]) )

    # Set label
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    # Change range
    plt.xlim([x[0], x[-1]] )

    # Save file
    plt.savefig("{0}/waveform-{1}.png".format(output_folder, N))

def plotWithAnnotations(x, y, peaks1, peaks, prominences, N, xlabel, ylabel, output_folder):

    ## Plot waveforms with respective peaks
    
    # Clear plot
    plt.cla()

    # Plot curve and pulses
    plt.plot(x, y, color='gray')
    plt.plot(x[peaks1], y[peaks1], 'x', color='blue')
    plt.plot(x[peaks], y[peaks], 'x', color='red')
    for k in range(0, len(prominences)):
        txt = "{0}".format(prominences[k])
        plt.annotate(txt, (x[peaks[k]], y[peaks[k]]) )
    
    # Set label
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    # Change range
    plt.xlim([x[0], x[-1]] )

    # Save file
    plt.savefig("{0}/waveform-{1}.png".format(output_folder, N))

def saveMultipleToROOT(run_number, bLine, hFactor, vFactor, nPulses, pulseIndices, pulseHeights):

    ## Save pulse information into a root file.
    #  Two branches: digitizer information and
    #  pulse information.

    # Write final 
    nameFinal = "./local-{0}.root".format(run_number)
    
    # Let user know
    print("Saving pulse information to: ", nameFinal)

    # Open ROOT file
    file  = ROOT.TFile(nameFinal,"RECREATE")

    # Digitizer information
    tree1 = ROOT.TTree("digitizer_info", "Digitizer information")
    
    # Digitizer information branches: baseline
    baseline = np.empty((1), dtype="float32")
    tree1.Branch("baseline", baseline, "baseline/F")
    
    # Digitizer information branches: vertical scale factor
    vfactor = np.empty((1), dtype="float32")
    tree1.Branch("verticalScaleFactor", vfactor, "verticalScaleFactor/F")
    
    # Digitizer information branches: horizontal scale factor
    hfactor = np.empty((1), dtype="float32")
    tree1.Branch("horizontalScaleFactor", hfactor, "horizontalScaleFactor/F")
    
    # Setting values
    baseline[0] = bLine
    vfactor[0]  = vFactor
    hfactor[0]  = hFactor

    #Filling tree
    tree1.Fill()

    # For each event, store array with pulse information
    tree2 = ROOT.TTree("pulse_info", "Pulse information")
    
    # Number of pulses
    pnumber = np.empty((1),dtype="float32")
    tree2.Branch("pulse_number", pnumber, "pulse_number/F")

    # Pulse indices
    pindices = ROOT.std.vector('float')(1)
    tree2.Branch("pulse_indices", pindices, "pulse_indices/F")
    
    # Pulse heights
    pheights = ROOT.std.vector('float')(1)
    tree2.Branch("pulse_heights", pheights, "pulse_heights/F")

    # Loop over each waveform
    for i in range(0, len(nPulses)):
        
        # Prepare pulse information
        pnumber[0]  = nPulses[i]

        # Prepare pulse index information
        for j in range(0,len(pulseIndices[i])):
            pindices.push_back(pulseIndices[i][j])

        # Prepare pulse height information
        for j in range(0,len(pulseHeights[i])):
            pindices.push_back(pulseHeights[i][j])

        # Fill tree
        tree2.Fill()

    # Write and close ROOT file
    file.Write()
    file.Close()

def saveGlobalMaximaToROOT(run_number, dtbranch, bLine, bLineError, hFactor, vFactor, pulseIndices, pulseHeights, pulseWidths, pulseAreas, pulseAreasError):

    ## Save pulse information into a root file.
    #  Two branches: digitizer information and
    #  pulse information.

    # Write final 
    nameFinal = "./global-{0}-{1}.root".format(run_number, dtbranch)

    print("Saving pulse information to: {0}".format(nameFinal))

    # Open ROOT file
    file  = ROOT.TFile(nameFinal,"RECREATE")

    # Digitizer information
    tree1 = ROOT.TTree("digitizer_info", "Digitizer information")
    
    # Digitizer information branches: baseline
    baseline = np.empty((1), dtype="float32")
    tree1.Branch("baseline", baseline, "baseline/F")

    # print(bLineError)

    # Digitizer information branches: baseline
    baselineError = np.empty((1), dtype="float32")
    tree1.Branch("baselineError", baselineError, "baselineError/F")
    
    # Digitizer information branches: vertical scale factor
    vfactor = np.empty((1), dtype="float32")
    tree1.Branch("verticalScaleFactor", vfactor, "verticalScaleFactor/F")
    
    # Digitizer information branches: horizontal scale factor
    hfactor = np.empty((1), dtype="float32")
    tree1.Branch("horizontalScaleFactor", hfactor, "horizontalScaleFactor/F")
    
    # Setting values
    baseline[0] = bLine
    baselineError[0] = bLineError
    vfactor[0]  = vFactor
    hfactor[0]  = hFactor

    #Filling tree
    tree1.Fill()

    # For each event, store array with pulse information
    tree2 = ROOT.TTree("pulse_info", "Pulse information")

    # Pulse index
    pindex = np.empty((1), dtype="int8")
    tree2.Branch("pulse_index", pindex, "pulse_index/I")
    
    # Pulse height
    pheight = np.empty((1), dtype="float32")
    tree2.Branch("pulse_height", pheight, "pulse_height/F")

    # Pulse width
    pwidth = np.empty((1), dtype="float32")
    tree2.Branch("pulse_width", pwidth, "pulse_width/F")

    # Pulse area
    parea = np.empty((1), dtype="float32")
    tree2.Branch("pulse_area", parea, "pulse_area/F")

    # Pulse area error 
    pareaErr = np.empty((1), dtype="float32")
    tree2.Branch("pulse_area_error", pareaErr, "pulse_area_error/F")

    # Loop over each waveform
    for i in range(0, len(pulseHeights)):
        
        # Prepare information
        pindex[0]   = pulseIndices[i]
        pheight[0]  = pulseHeights[i]
        parea[0]    = pulseAreas[i]
        pareaErr[0] = pulseAreasError[i]
        pwidth[0]   = pulseWidths[i]

        # Fill tree
        tree2.Fill()

    # Write and close ROOT file
    file.Write()
    file.Close()    

def saveGlobalMaximaToROOT_uprootinfo(run_number, dtbranch, bLine, bLineError, hFactor, vFactor, pulseIndices, pulseHeights, pulseAreas, pulseAreasError):

    ## Save pulse information into a root file.
    #  Two branches: digitizer information and
    #  pulse information.

    # Write final 
    nameFinal = "./global_uproot-{0}-{1}.root".format(run_number, dtbranch)

    print("Saving pulse information to: {0}".format(nameFinal))

    # Open ROOT file
    file  = ROOT.TFile(nameFinal,"RECREATE")

    # Digitizer information
    tree1 = ROOT.TTree("digitizer_info", "Digitizer information")
    
    # Digitizer information branches: baseline
    baseline = np.empty((1), dtype="float32")
    tree1.Branch("baseline", baseline, "baseline/F")

    # print(bLineError)

    # Digitizer information branches: baseline
    baselineError = np.empty((1), dtype="float32")
    tree1.Branch("baselineError", baselineError, "baselineError/F")
    
    # Digitizer information branches: vertical scale factor
    vfactor = np.empty((1), dtype="float32")
    tree1.Branch("verticalScaleFactor", vfactor, "verticalScaleFactor/F")
    
    # Digitizer information branches: horizontal scale factor
    hfactor = np.empty((1), dtype="float32")
    tree1.Branch("horizontalScaleFactor", hfactor, "horizontalScaleFactor/F")
    
    # Setting values
    baseline[0] = bLine
    baselineError[0] = bLineError
    vfactor[0]  = vFactor
    hfactor[0]  = hFactor

    #Filling tree
    tree1.Fill()

    # For each event, store array with pulse information
    tree2 = ROOT.TTree("pulse_info", "Pulse information")

    # Pulse index
    pindex = np.empty((1), dtype="int8")
    tree2.Branch("pulse_index", pindex, "pulse_index/I")
    
    # Pulse height
    pheight = np.empty((1), dtype="float32")
    tree2.Branch("pulse_height", pheight, "pulse_height/F")

    # Pulse area
    parea = np.empty((1), dtype="float32")
    tree2.Branch("pulse_area", parea, "pulse_area/F")

    # Pulse area error 
    pareaErr = np.empty((1), dtype="float32")
    tree2.Branch("pulse_area_error", pareaErr, "pulse_area_error/F")

    # Loop over each waveform
    for i in range(0, len(pulseHeights)):
        
        # Prepare information
        pindex[0]   = pulseIndices[i]
        pheight[0]  = pulseHeights[i]
        parea[0]    = pulseAreas[i]
        pareaErr[0] = pulseAreasError[i]

        # Fill tree
        tree2.Fill()

    # Write and close ROOT file
    file.Write()
    file.Close()    



def saveGlobalMaximaToROOTold(run_number, dtbranch, pulseIndices, pulseHeights, pulseAreas):

    ## Save pulse information into a root file.
    #  Two branches: digitizer information and
    #  pulse information.

    # Write final 
    nameFinal = "./mppc-analysis_results-{0}-{1}.root".format(run_number, dtbranch)
    
    # Let user know
    print("Saving pulse information to: ", nameFinal)

    # Open ROOT file
    file  = ROOT.TFile(nameFinal,"RECREATE")

    # Digitizer information
    tree1 = ROOT.TTree("pulse_tree","Pulse Tree")
    
    # Waveform
    waveformIndices = range(0, len(pulseAreas))
    pwaveform = ROOT.std.vector('float')(1)
    tree1.Branch("waveform_index", pwaveform)

    # Pulse indices
    pindices = ROOT.std.vector('float')(1)
    tree1.Branch("sample_index", pindices)
    
    # Pulse heights
    pheights = ROOT.std.vector('float')(1)
    tree1.Branch("pulse_height", pheights)
    
    # Pulse areas
    pareas = ROOT.std.vector('float')(1)
    tree1.Branch("pulse_area", pareas)

    for i in range(0, len(pulseAreas)):

        # Prepare information to be stored in vectors
        pwaveform.push_back(waveformIndices[i])
        pindices.push_back(pulseIndices[i])
        pheights.push_back(pulseHeights[i])
        pareas.push_back(pulseAreas[i])

    #Filling tree
    tree1.Fill()

    # Write and close ROOT file
    file.Write()
    file.Close()    