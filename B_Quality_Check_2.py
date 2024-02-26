
## What does MultiQC compare? All FastQCs for a given variable, or condition?

# This function gives a visual report summarizing FastQCs
# For this script, we need multiqc installed

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