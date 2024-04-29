# This script contains mini extra functions
import os
import pandas as pd
import subprocess

# Quick functions to access the var/cond under which the given SRR is stored
def var(SRR_path):
    folders = SRR_path.split('/') # splitting the path into its folders
    var_number = folders[-3].split('_')[1] # counting backwards from known number of '/'
    return var_number
def cond(SRR_path):
    folders = SRR_path.split('/') # splitting the path into its folders
    cond_number = folders[-2].split('_')[1] # counting backwards from known number of '/'
    return cond_number
def name(SRR_path): # returns just the SRR?
    folders = SRR_path.split('/') # splitting the path into its folders
    name = folders[-1]
    return name

def listdir_visible(directory):
    # This function only lists the visible files and directories in a directory
        # in comparison, the os.listdir() function lists both visible and hidden directories
    return [item for item in os.listdir(directory) if not item.startswith('.')]

class Bioinf_Profile():
    def __init__(self, profile_path=".bioinf-profile"):
        self.profile_path = profile_path
        self.dict = { # using a dictionary to make writing to text file easier
            "STEP": "",
            "PATH": "",
            "SRA_TOOLKIT_PATH": "",
            "FASTQC_PATH": "",
            "RNA_STAR_PATH": "",
            "SEQ2FUN_PATH": ""
        }
        self.VAR_names = []
        self.COND_names = []
        self.COND_corr = []
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


def remove_apple_quarantine(tool_path):
     # On MacOS, a folder downloaded from the internet will automatically be flagged as coming
     # from an unknown developper and this script will not be able to open the folder
    # check for the com.apple.quarantine flag and remove it if the folder has it
    xattr_list=subprocess.run(["xattr", tool_path], capture_output=True, text=True).stdout.split('\n')
    if "com.apple.quarantine" in xattr_list:
        subprocess.run(["xattr", "-r", "-d","com.apple.quarantine", tool_path])