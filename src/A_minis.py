# This script contains mini extra functions
import os

# Quick functions to access the var/cond under which the given SRR is stored
def var(SRR_path):
    folders = SRR_path.split('/') # splitting the path into its folders
    var_number = folders[-3].split('_')[1] # counting backwards from known number of '/'
    return var_number
def cond(SRR_path):
    folders = SRR_path.split('/') # splitting the path into its folders
    cond_number = folders[-2].split('_')[1] # counting backwards from known number of '/'
    return cond_number
def name(SRR_path):
    folders = SRR_path.split('/') # splitting the path into its folders
    name = folders[-1].split('_')[1] # counting backwards from known number of '/'
    return name


class Bioinf_Profile():
    def __init__(self, profile_path=".bioinf-profile"):
        self.profile_path = profile_path
        self.dict = { # using a dictionary to make writing to text file easier
            "STEP": "",
            "PATH": "",
            "SRA_TOOLKIT_PATH": ""
        }
        return
    
    def read_profile(self):
        # read path information from the .bioinf-profile file
        # assign the paths to the Bioinf_Profile object properties
        file = open(self.profile_path, 'r')
        for line in file:
            split = line.split('=',1) 
                # split on first occurence of '='
            self.dict[split[0]] = split[1].rstrip()
                # the first word will be the variable name/key
                # second word will be the variable value
                # rstrip() removes trailing '\n'
        return
    
    def update_profile(self):        
        file = open(self.profile_path, 'w')
        for key in self.dict:
            file.write(key + '=' + self.dict[key] + '\n')
        return

if os.path.isfile(".bioinf-profile"):
    profile1 = Bioinf_Profile()
    profile1.read_profile()
    
