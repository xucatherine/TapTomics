# Potential FastQC for-loop - iterates through all folders in Samples and runs fastqc for each SRR folder - for main script
import os
Samples_path = 'this needs to have been assigned' # replace with real path
fastqc_path = 'assign this too'
for var in os.listdir(Samples_path): # for each variable in the Samples folder
    if os.path.isdir(os.path.join(Samples_path, var)):
        var_path = os.path.join(Samples_path, var)
        for cond in os.listdir(var_path): # for each condition in the var_n folder
            if os.path.isdir(os.path.join(var_path, cond)):
                cond_path = os.path.join(var_path, cond)
                for SRR in os.listdir(cond_path): # for each SRR in the cond_k folder
                    if os.path.isdir(os.path.join(cond_path, SRR)):
                        SRR_path = cond_path+"/"+os.path.basename(os.getcwd())
                        run_FastQC(fastqc_path,SRR_path)



## Should we suggest troubleshooting or links to follow for bad results?
## main script should ask user if they want to see multiQC report before generating it
                        

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
        print(f"FastQC analysis completed.")
    except subprocess.CalledProcessError as e:
        print(f"Error during FastQC analysis: {e}")


# MultiQC Inputs
    # needs path to relevant condition folder
# Outputs
    # outputs link to view results in a browser

# Generating MultiQC Quality Check Visuals
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