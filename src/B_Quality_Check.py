## Should we suggest troubleshooting or links to follow for bad results?
## main script should ask user if they want to see multiQC report before generating it
## specify var/cond for each reported successful/failed FastQC analysis

          
# For this script, we need fastqc and multiqc installed
# First function checks the quality of a FASTA file, so for several it must be implemented in a for-loop
# Second function gives a visual report summarizing FastQCs for a given condition

# Imports
import subprocess
import A_minis # to be able to call var(), cond() and name()
import os
import time


# FastQC Inputs
    # needs path to FASTQ file's SRR folder
    # needs path to FastQC (FastQC must be downloaded by user)
# Ouputs
    # outputs 'FastQC' file in appropriate SRR folder

# Running FastQC Quality Check
def run_FastQC(fastqc_path, SRR_path):
    fastq_list = ["rawF.fastq", "rawR.fastq"]

    for fastq_file in fastq_list: # iterate through the forward and reverse files
        cmd = [fastqc_path, "--extract", "--delete", "--quiet", os.path.join(SRR_path, fastq_file)] 
            # Set up the command for FastQC
            # note that FastQC outputs the result files in the same directory as the input directory
            ## FastQC produces several different files; we only need the raw data file in the SRR folder.
                ## Later, I want to make it so that the other files get deleted or moved to another folder.
        attempt = 0
        wait_time = 5
        while attempt < 3:
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE)
            except subprocess.CalledProcessError:
                print(f"Error during FastQC analysis for variable {A_minis.var(SRR_path)}, condition {A_minis.cond(SRR_path)}'s {A_minis.name(SRR_path)}'s {fastq_file}.")
                print(f"Retrying in {wait_time} seconds.")
                time.sleep(wait_time)
                attempt += 1
            else: # no error
                print(f"FastQC analysis for variable {A_minis.var(SRR_path)}, condition {A_minis.cond(SRR_path)}'s {A_minis.name(SRR_path)} {fastq_file} completed.")
                break

            '''result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # running FastQC
            if result.returncode == 0: # Check if FastQC was successful
                print(f"FastQC analysis for variable {A_minis.var(SRR_path)}, condition {A_minis.cond(SRR_path)}'s {A_minis.name(SRR_path)} {fastq_file} completed.") # specifies using MAIN's functions
                break
            else:
                print(f"Error during FastQC analysis for variable {A_minis.var(SRR_path)}, condition {A_minis.cond(SRR_path)}'s {A_minis.name(SRR_path)}'s {fastq_file}.")
                print(f"Retrying in {wait_time} seconds.")
                time.sleep(wait_time)
                attempt += 1'''
        
        if attempt == 3:
            print(f"FastQC has failed {attempt} times. There may be an issue with the fastq files or FastQC.")
            return -1
    return



# MultiQC Inputs
    # needs path to relevant condition folder
# Outputs
    # outputs link to view results in a browser

# Generating MultiQC Quality Check Visuals
def run_MultiQC(cond_path):
    folders = cond_path.split('/') # splitting the path into its folders to extract cond
    var = folders[-2].split('_')[1]
    cond = folders[-1].split('_')[1]

    cmd = ["multiqc", cond_path] # making command term, MultiQC will scan directory for compatible files
    result = subprocess.run(cmd, capture_output=True, text=True) # running MultiQC
    # check if MultiQC ran successfully
    if result.returncode == 0:
        print("MultiQC ran successfully for variable "+var+"'s condition "+cond+"! To view results, open this file in your browser:\n")
        print(result.stdout)
    else:
        print("MultiQC encountered an error for variable "+var+"'s condition "+cond+":")
        print(result.stderr)
    
    return