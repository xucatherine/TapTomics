## Should we suggest troubleshooting or links to follow for bad results?
## main script should ask user if they want to see multiQC report before generating it
## specify var/cond for each reported successful/failed FastQC analysis


          
# For this script, we need fastqc and multiqc installed
# First function checks the quality of a FASTA file, so for several it must be implemented in a for-loop
# Second function gives a visual report summarizing FastQCs for a given condition

# Imports
import subprocess

# FastQC Inputs
    # needs path to FASTQ file's SRR folder
    # needs path to FastQC (FastQC must be downloaded by user)
# Ouputs
    # outputs 'FastQC' file in appropriate SRR folder

# Running FastQC Quality Check
def run_FastQC(fastqc_path,SRR_path):
    FASTQ_path = SRR_path+"/FASTQ"
    FastQC_path = SRR_path+"/FastQC" # output directory
    cmd = [fastqc_path, FASTQ_path, '-o', FastQC_path] # making command term
    subprocess.run(cmd, check=True) 
    try:
        subprocess.run(cmd, check=True) # running FastQC
        print(f"FastQC analysis for variable "+var(SRR_path)+", condition "+cond(SRR_path)+" completed.") # specifies using MAIN's functions
    except subprocess.CalledProcessError as e:
        print(f"Error during FastQC analysis for variable "+var(SRR_path)+", condition "+cond(SRR_path)+": {e}")

# MultiQC Inputs
    # needs path to relevant condition folder
# Outputs
    # outputs link to view results in a browser

# Generating MultiQC Quality Check Visuals
def run_MultiQC(cond_path):
    folders = cond_path.split('/') # splitting the path into its folders to extract cond
    cond = folders[-1].split('_')[1]

    cmd = ["multiqc", cond_path] # making command term, MultiQC will scan directory for compatible files
    result = subprocess.run(cmd, capture_output=True, text=True) # running MultiQC
    # check if MultiQC ran successfully
    if result.returncode == 0:
        print("MultiQC ran successfully for condition "+cond+"! To view results, open this file in your browser:\n")
        print(result.stdout)
    else:
        print("MultiQC encountered an error for condition "+cond+":")
        print(result.stderr)