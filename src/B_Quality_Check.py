
## ChatGPT highlights that FastQC is designed for FASTQ files, not FASTA files 

## main script should ask user if they want to see multiQC report before generating it
## Should we suggest troubleshooting or links to follow for bad results?


# First function checks the quality of a FASTA file, so for several it must be implemented in a for-loop
# Second function gives a visual report summarizing FastQCs

# For this script, we need fastqc and multiqc installed


# FastQC

# Inputs
    # needs path to FASTA file's SRA folder
    # needs path to FastQC (FastQC must be downloaded by user)

# Imports
import subprocess

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

#MultiQC
## What does MultiQC compare? All FastQCs for a given variable, or condition?

# Inputs
    # needs path to relevant condition folder

# Imports
import subprocess

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