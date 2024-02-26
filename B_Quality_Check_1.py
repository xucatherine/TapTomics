
## ChatGPT highlights that FastQC is designed fro FASTQ files, not FASTA files 
## make sure FastQC is downloaded in Download script
## main script should ask user if they want to see multiQC report - then call B_Quality_Check_2.py
## Should we suggest troubleshooting or links to follow for bad results?


# This function checks the quality of a FASTA file, so for several it must be implemented in a for-loop
# For this script, we need fastqc installed

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