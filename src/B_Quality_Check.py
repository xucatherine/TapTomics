
## ChatGPT highlights that FastQC is designed for FASTQ files, not FASTA files
## What does MultiQC compare? All FastQCs for a given variable, or condition?
## Should we suggest troubleshooting or links to follow for bad results?
## main script should ask user if they want to see multiQC report before generating it

## Potential for-loop - scans Samples folder, runs run_FastQC on each SRA folder
import os
Samples_path = 'this needs to have been assigned' # replace with real path
fastqc_path = 'assign this too'
for variable in os.listdir(Samples_path): # for each variable in the Samples folder
        if os.path.isdir(os.path.join(Samples_path, variable)):
            n = 1 # starting count of variables
            var_path = Samples_path+"/var_"+str(n)
            for condition in os.listdir(var_path): # for each condition in the var_n folder
                k = 1 # starting count of conditions
                cond_path = var_path+"/cond_"+str(k)
                if os.path.isdir(os.path.join(var_path, condition)):
                    for SRA in os.listdir(cond_path): # for each SRA in the cond_k folder
                        if os.path.isdir(os.path.join(cond_path, SRA)):
                            SRA_path = cond_path+"/"+os.path.basename(os.getcwd())
                            run_FastQC(fastqc_path,SRA_path)
                k += 1
            n += 1



# For this script, we need fastqc and multiqc installed
# First function checks the quality of a FASTA file, so for several it must be implemented in a for-loop
# Second function gives a visual report summarizing FastQCs


# Imports
import subprocess

# FastQC

# Inputs
    # needs path to FASTA file's SRA folder
    # needs path to FastQC (FastQC must be downloaded by user)

# Running Quality Check
def run_FastQC(fastqc_path,SRA_path):
    FASTA_path = SRA_path+"/FASTA"
    FastQC_path = SRA_path+"/FastQC" # output directory
    cmd = [fastqc_path, FastQC_path+".fastq", '-o', FastQC_path] # making command term
    subprocess.run(cmd, check=True) 
    try:
        subprocess.run(cmd, check=True) # running FastQC
        print(f"FastQC analysis completed.")
    except subprocess.CalledProcessError as e:
        print(f"Error during FastQC analysis: {e}")


# MultiQC

# Inputs
    # needs path to relevant condition folder

# Generating Quality Check Visuals
def run_MultiQC(cond_path):
    cmd = ["multiqc", cond_path] # making command term, MultiQC will scan directory for compatible files
    result = subprocess.run(cmd, capture_output=True, text=True) # running MultiQC
    # check if MultiQC ran successfully
    if result.returncode == 0:
        print("MultiQC ran successfully! To view results, open this file in your browser:\n")
        print(result.stdout)
    else:
        print("MultiQC encountered an error:")
        print(result.stderr)
