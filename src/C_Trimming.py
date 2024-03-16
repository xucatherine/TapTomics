# For this script, we need cutadapt installed
# Runs Cutadapt on a given FASTQ file to trim suspected adapter sequences


# imports
import subprocess
import A_minis
    # to be able to call var(), cond() and name()

# Inputs
    # needs path to FASTQ file's SRR folder
# Ouputs
    # outputs new 'trimmed' fastq file in SRR folder

# Running Cutadapt
def run_cutadapt(SRR_path):
    fastq_path = SRR_path+"/FastQC"
    output_path = SRR_path+"/trimmed"

    # building list of suspected adaptamer sequences to trim, from FastQC file
    adapter_sequences = []  # list to store adapter sequences
    with open(fastq_path, 'r') as file:
        lines = file.readlines()    
    start = False # flag to start recording sequences when we're in the right section
    for line in lines:
        if '>>Overrepresented sequences' in line: # identifying start of overepresented sequences in file
            start = True
            continue
        if start and '>>END_MODULE' in line: # stop loop at end of section
            break
        if start and not line.startswith('#') and line.strip(): # record non-comment, non-empty sequences
            sequence = line.split('\t')[0]
            source = line.split('\t')[-1].strip() # checking if 'adapter' is mentioned in the source column to filter sequences
            if "adapter" in source.lower():
                adapter_sequences.append(sequence) # if so, add sequence to list

    # constructing Cutadapt command
    command = ['cutadapt', '-o', output_path] # start by specifying output file
    for adapter in adapter_sequences: # add each adapter sequence to the command with its option
        command.extend(['-a', adapter])
    command.append(fastq_path) # add the input FASTQ file to the command
    
    # executing Cutadapt command
    try:
        subprocess.run(command, check=True) 
    except subprocess.CalledProcessError as e:
        print(f"Error running Cutadapt for variable {A_minis.var(SRR_path)}, {A_minis.cond(SRR_path)}'s {A_minis.name(SRR_path)}: {e}")