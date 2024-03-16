## Should we suggest troubleshooting or links to follow for bad results?
## main script should ask user if they want to see multiQC report before generating it
## specify var/cond for each reported successful/failed FastQC analysis

          
# For this script, we need fastqc and multiqc installed
# First function checks the quality of a FASTA file, so for several it must be implemented in a for-loop
# Second function gives a visual report summarizing FastQCs for a given condition

# Imports
import subprocess
import A_minis
    # to be able to call var(), cond() and name()

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
        print(f"FastQC analysis for variable "+A_minis.var(SRR_path)+", condition "+A_minis.cond(SRR_path)+"'s "+A_minis.name(SRR_path)+" completed.") # specifies using MAIN's functions
    except subprocess.CalledProcessError as e:
        print(f"Error during FastQC analysis for variable "+A_minis.var(SRR_path)+", condition "+A_minis.cond(SRR_path)+"'s "+A_minis.name(SRR_path)+": {e}")

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