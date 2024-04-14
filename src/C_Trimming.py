# For this script, we need cutadapt installed
# Runs Cutadapt on a given FASTQ file to trim suspected adapter sequences


# imports
import subprocess
from A_minis import var, cond, name

    # to be able to call var(), cond() and name()

# Inputs
    # needs path to FASTQ file's SRR folder
# Ouputs
    # outputs new 'trimmed' fastq file in SRR folder

# Running Cutadapt
def run_cutadapt(SRR_path):
    #making list of fastqc data paths
    #first item is forward fastqc and second is reverse fastqc

    #for testing
    fastqc_paths = [SRR_path+"/SRR22218921_forward fastQC.txt", SRR_path+"/SRR22218921_reverse fastQC.txt"]
    fastq_paths = [SRR_path + '/SRR22218921_1.fastq', SRR_path + '/SRR22218921_2.fastq']
    output_paths = [SRR_path+"/trimmedF.fastq", SRR_path+"/trimmedR.fastq"]

    '''
    fastqc_paths = [SRR_path+"/rawF_fastqc/fastqc_data.txt", SRR_path+"/rawR_fastqc/fastqc_data.txt"]
    fastq_paths = [SRR_path + '/rawF.fastq', SRR_path + '/rawR.fastq']
    output_paths = [SRR_path+"/trimmedF.fastq", SRR_path+"/trimmedR.fastq"]
    '''
    #here we are assuming that the adpter sequences are only found in the forward QC
    #also assuming only 3' end adapters since illumina is most common sequencing

    # building list of suspected adaptamer sequences to trim, from FastQC file
    for_adapter_sequences = []
    with open(fastqc_paths[0], 'r') as file:
        lines = file.readlines()    
    start = False # flag to start recording sequences when we're in the right section
    for line in lines:
        if '>>Overrepresented sequences' in line: # identifying start of overepresented sequences in file
            start = True
            continue
        if start and '>>END_MODULE' in line: # stop loop at end of section
            break
        #if start and not line.startswith('#') and line.strip(): # record non-comment, non-empty sequences
        #why the line.strip???? better to write 'and not line.strip()==''
        #this is better but I don't think we need the check for empty lines
        if start and not line.startswith('#'): # record non-comment sequences
            sequence = line.split('\t')[0]
            source = line.split('\t')[-1].strip() # checking if 'adapter' is mentioned in the source column to filter sequences
            if "adapter" in source.lower():
                for_adapter_sequences.append(sequence) # if so, add sequence to list
    

    # constructing Cutadapt command
    command = ['cutadapt']
    for i in for_adapter_sequences: # add each adapter sequence to the command with its option
        command.extend(['-a', i])
        #get reverse adapter
        #not sure if this is 100% right
        command.extend(['-A', i])
    command.extend(['-o', output_paths[0], '-p', output_paths[1]]) # specify output file
    command.append(' '.join(str(e) for e in fastq_paths)) # add the input FASTQ file to the command
    #converting into string
    command = ' '.join(str(e) for e in command)


    # executing Cutadapt command
    try:
        subprocess.run(command, check=True) 
    except subprocess.CalledProcessError as e:
        print(f"Error running Cutadapt for variable {var(SRR_path)}, {cond(SRR_path)}'s {name(SRR_path)}: {e}")