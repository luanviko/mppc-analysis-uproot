import ROOT
import numpy as np

def import_from_root(root_file_name):

    ## Open ROOT file generated by SPTF analysis,
    ## returning the information stored in the 
    ## four branches: waveform_index,
    ## sample where maxima are,
    ## maxima values (pulse heights), and
    ## waveform area (pulse area)
    
    ## Open ROOT file
    file = ROOT.TFile(root_file_name, "READ")
    tree = file.Get("pulse_tree")
    
    ## Define vectors
    vec_waveform_indices = ROOT.std.vector("Int_t")()
    vec_sample_indices   = ROOT.std.vector("Int_t")()
    vec_maxima_values    = ROOT.std.vector("Double_t")()
    vec_maxima_areas     = ROOT.std.vector("Double_t")()
    
    ## Get branches
    branch_waveform_indices = tree.GetBranch("waveform_index")
    branch_waveform_indices = tree.GetBranch("sample_index")
    branch_pulse_heights    = tree.GetBranch("pulse_height")
    branch_pulse_areas      = tree.GetBranch("pulse_area")
    
    ## Point branches to vectors
    tree.SetBranchAddress("waveform_index", vec_waveform_indices)
    tree.SetBranchAddress("sample_index", vec_sample_indices)
    tree.SetBranchAddress("pulse_height", vec_maxima_values)
    tree.SetBranchAddress("pulse_area", vec_maxima_areas)
    
    ## Get information
    tree.GetEntry(0)
    
    # Transform arrays into numpy arrays
    waveform_indices = np.asarray(vec_waveform_indices)
    sample_indices   = np.asarray(vec_sample_indices)
    pulse_heights    = np.asarray(vec_maxima_values)
    pulse_areas      = np.asarray(vec_maxima_areas)
    
    # Return numpy arrays to 
    return waveform_indices, sample_indices, pulse_heights, pulse_areas

def load_MPPC_information(config_file):

    ## Loads up information in config_file 
    ## and stores in a dictionary called info.
    ## Input file can have as many entries as necessary, 
    ## but organized as:
    ## Entry: value

    # Empty dictionary
    info = {}

    # Open config file
    with open(config_file, "r") as config:

        # Read lines from file
        lines = config.readlines()

        # Go through each line
        for line in lines:

            # Ignore lines with '#'
            if not ("#" in line):
                
                # Split entry's name from entry's value
                entries = line.split(": ")

                # Add information to dictionary
                info[ "{0}".format(entries[0])] = "{0}".format(entries[1].replace("\n",""))
                
    # Return information
    return info

def getEntries(root_file_name):

    ## Open ROOT file generated by SPTF analysis
    ## and gets the number of entries
    
    # Open ROOT file
    file = ROOT.TFile(root_file_name, "READ")
    tree = file.Get("pulse_tree")

    vec_waveform_indices = ROOT.std.vector("Int_t")()
    
    ## Get branches
    branch_waveform_indices = tree.GetBranch("waveform_index")
    
    ## Point branches to vectors
    tree.SetBranchAddress("waveform_index", vec_waveform_indices)
    
    ## Get information
    tree.GetEntry(0)

    # Get number of entries
    return len(vec_waveform_indices)