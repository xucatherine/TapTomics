# This script contains mini extra functions

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