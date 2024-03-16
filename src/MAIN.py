
## Missing: ask the user for brief description each variable and condition, and store entries in a dictionary
## Get rid of FASTA files?


# Welcome to the MAIN script for the Bioinformatic Pipeline! 
# The user runs this script, and is guided through the steps of transcriptomic analysis

# Folder Architecture (for reference, this is how downloaded data will be stored)
'''
folder: Samples
    folder: var_x    # var_1, var_2, var_3, ...
        folder: cond_y    # cond_1, cond_2, cond_3,...
            folder: SRR (subfiles can be added later)
                file: FASTA        (suggestion: raw.fasta)
                file: FASTQ        (suggestion: raw.fastq)
                file: FastQC        
                file: trimmed    (suggestion: trimmed.fastq)
                file: BAM        (suggestion: aligned.bam)
                [...]
            [more SRR folders]
        [more cond folders]
    [more var folders]

folder: References
    [type unsure]: ref_genome
    [...]

folder: Results
    [...]
'''

# Imports
from git import Repo
    # using gitpython we can access GitHub to make new folders for user's data
from Bio import Entrez
    # 'Entrez' allows us to send requests to NCBI databases
import os
    # we need this to edit folders
import subprocess
    # we need this to call fastq-dump

# Setting up path to save folders in & user email to access NCBI
path = str(input("Please enter a path link for where you would like this pipeline's running data stored: "))
repo = Repo.clone_from("https://github.com/xucatherine/bioinformatic-pipeline/src.git", path)
Entrez.email = str(input("NCBI requires an email address to track usage of their services.\nPlease input your email address: "))
    # NCBI Entrez requires an email address (according to BioBuddy)
    ## Not sure if email address needs to be linked to an account already

# Setting up Samples, References and Results folders (within src)
Samples_path = path+"/Samples"
os.makedirs(Samples_path, exist_ok=True)
References_path = path+"/References"
os.makedirs(References_path, exist_ok=True)
Results_path = path+"/Results"
os.makedirs(Results_path, exist_ok=True)

# Entering data into Samples folder - starting with number of variables to study
print("Using differential analysis to decode metabolic pathways involves identifying experimental variables under which expression of the phenotype of interest differs.")
print("Through changing the variable's intensity, we create different experimental conditions under which transcriptomes can be collected and compared against each other.")
print("Ex: In studying yeast respiration, temperature is a variable. Low, medium and high heat would correspond to 3 conditions, between which cellular respiration varies.")
print("\nHow many variables are you considering for your pathway? (minimum 1)")
print("For each variable, you will be asked to input the number of conditions studied, and for each condition list the corresponding SRR accession numbers.")
n = int(input("number of variables studied: ")) #[user inputs answer]
## initialize a dictionary here to which var and cond names can be stored?

# For each variable, make the right number of conditions folders and add data
for i in range(n):
    var_path = Samples_path+"/var_"+str(i+1) # folders will be named var_1, var_2, var_3...
    os.makedirs(var_path, exist_ok=True)
    ## ask user to name variable and store name in a dictionary?

    # Adding conditions folders
    k = int(input("number of conditions studied for variable " + str(i+1) + ': ')) #[user inputs answer]
    for j in range(k):
        cond_path = var_path+"/cond_"+str(j+1) # folders will be named cond_1, cond_2, cond_3...
        os.makedirs(cond_path, exist_ok=True)
        ## ask user to name condition and store name in a dictionary?

    # Getting SRRs and downloading their FASTA & FASTQ files, for each variable
    for cond in os.listdir(var_path): # for condition folder in parent folder var_n
        print(f"List the SRR numbers belonging to variable {i+1}'s condition {cond.split('_')[-1]}, separated by spaces.")
        print("ex: SRR12345678 SRR91011109 SRR87654321")
        SRRs = input("SRR numbers: ") #[user inputs SRRs]
        SRRs_list = SRRs.split() # making list of SRRs
        
        # Starting the downloads
        if os.path.isdir(os.path.join(var_path, cond)):
            cond_path = os.path.join(var_path, cond)
            for SRR in SRRs_list: # loop through each SRR
                SRR_path = os.path.join(cond_path, SRR) # new folder named as SRR number
                os.makedirs(SRR_path, exist_ok=True)
                
                # retrieve FASTA file (script by ChatGPT)
                handle = Entrez.efetch(db="nucleotide", id=SRR, rettype="fasta", retmode="text")
                fasta_data = handle.read()
                handle.close()
                with open(SRR_path+"/FASTA", "w") as file: # adding FASTA file to folder
                    file.write(fasta_data)
                
                # retrieve FASTQ file (script mostly by ChatGPT)
                try:
                    subprocess.run(["fastq-dump", "--split-files", "--gzip", "--outdir", SRR_path, SRR], check=True)
                    # renaming the file as FASTQ
                    for file in os.listdir(SRR_path):
                        if file.startswith(SRR) and file.endswith(".fastq.gz"):
                            os.rename(os.path.join(SRR_path, file), os.path.join(SRR_path, "FASTQ.gz"))
                except subprocess.CalledProcessError as e:
                    print(f"An error occurred while downloading {SRR}: {e}")                
print("Thank you! Your transcriptomes are ready for processing.")

# Creating list of SRRs for easy processing - iterates through all folders in Samples
SRR_paths = []
for var in os.listdir(Samples_path): # for each variable in the Samples folder
    var_path = os.path.join(Samples_path, var)
    for cond in os.listdir(var_path): # for each condition in the var_n folder
        cond_path = os.path.join(var_path, cond)
        for SRR in os.listdir(cond_path): # for each SRR in the cond_k folder
            SRR_path = os.path.join(cond_path, SRR)
            SRR_paths.append(SRR_path)

# Running Quality Checks on SRR data
import B_Quality_Check
fastqc_path = 'resolve' # path to previously downloaded FastQC, within user's computer
    ## save this path earlier when dowloading, or have user manually input it here
for SRR in SRR_paths:
    B_Quality_Check.run_FastQC(fastqc_path,SRR)

# Optional: Running MultiQC on FastQC data to create visual for user, per condition
print("Would you like to run a visualizer for the quality checks done on each of your conditions?")
print("This will output links which you can open in a browser.")
query = input("Run MultiQC? y/n : ")
if query == 'y':
    for var in os.listdir(Samples_path): # for each variable in the Samples folder
        var_path = os.path.join(Samples_path, var)
        for cond in os.listdir(var_path): # for each condition in the var_n folder
            B_Quality_Check.run_MultiQC(cond_path)

# Trimming SRR data using Quality Check results
import C_Trimming
for SRR in SRR_paths:
    C_Trimming.run_cutadapt(SRR)

# Next...
    
    
