## Ask the user for brief description of variables and conditions?
    ## could be stored in text files to keep track of data meaning but not necessary
## Ask for a reference genome?


# For this script, we need gitpython, biopython and the SRA Toolkit installed
# This script classifies user-inputted SRA data within a folder architecture
    # it downloads both the FASTA and FASTQ files for given SRAs
    # I'm not sure if the FASTA files have any use

# Folder Architecture
'''
folder: Samples
    folder: var_x
        folder: cond_1
            folder: SRA (subfiles can be added later)
                file: FASTA
            [more SRA folders]
        folder: cond_2
            folder: SRA
                file: FASTA
            [more SRA folders]
        [more cond folders]
    folder: var_y
        [...]
    [more variable folders]
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


# Setting up path to save folders in & email to access NCBI
path = str(input("Please copy-paste a path link for where you would like this program's functioning data stored: "))
repo = Repo.clone_from("https://github.com/xucatherine/bioinformatic-pipeline/src.git", path)
Entrez.email = str(input("NCBI requires an email address to track usage of their services.\nPlease input your email address: "))
    # NCBI Entrez requires an email address (according to BioBuddy)
    ## Not sure if email address needs to be linked to an account already

# Setting up Samples folder (within src)
Samples_path = path+"/Samples"
os.makedirs(Samples_path, exist_ok=True)

# Entering data - starting with number of variables to study
print("Using differential analysis to decode metabolic pathways involves identifying experimental variables under which expression of the phenotype of interest differs.")
print("Through changing the variable's intensity, we create different experimental conditions under which transcriptomes can be collected and compared against each other.")
print("Ex: In studying yeast respiration, temperature is a variable. Low, medium and high heat would correspond to 3 conditions, between which cellular respiration varies.")
print("\nHow many variables are you considering for your pathway? (minimum 1)")
print("For each variable, you will be asked to input the number of conditions studied, and for each codnition list the corresponding SRA accession numbers.")
n = int(input("number of variables studied: ")) #[user inputs answer]

# For each variable, make the right number of conditions folders and add data
for i in range(n):
    var_path = Samples_path+"/var_"+str(i+1) # folders will be named var_1, var_2, var_3...
    os.makedirs(var_path, exist_ok=True)

    # Adding conditions folders
    k = int(input("number of conditions studied for variable " + str(i+1) + ': ')) #[user inputs answer]
    for j in range(k):
        cond_path = var_path+"/cond_"+str(j+1) # folders will be named cond_1, cond_2, cond_3...
        os.makedirs(cond_path, exist_ok=True)

    # Getting SRAs and downloading their FASTA & FASTQ files, for each variable
    for cond in os.listdir(var_path): # for condition folder in parent folder var_n
        print(f"List the SRA numbers belonging to variable {i+1}'s condition {cond.split('_')[-1]}, separated by spaces.")
        print("ex: SRR12345678 SRR91011109 SRR87654321")
        SRAs = input("SRA numbers: ") #[user inputs SRAs]
        SRAs_list = SRAs.split() # making list of SRAs
        
        # Starting the downloads
        if os.path.isdir(os.path.join(var_path, cond)):
            cond_path = os.path.join(var_path, cond)
            for SRA in SRAs_list: # loop through each SRA
                SRA_path = os.path.join(cond_path, SRA) # new folder named as SRA number
                os.makedirs(SRA_path, exist_ok=True)
                
                # retrieve FASTA file (script by ChatGPT)
                handle = Entrez.efetch(db="nucleotide", id=SRA, rettype="fasta", retmode="text")
                fasta_data = handle.read()
                handle.close()
                with open(SRA_path+"/FASTA", "w") as file: # adding FASTA file to folder
                    file.write(fasta_data)
                
                # retrieve FASTQ file (script mostly by ChatGPT)
                try:
                    subprocess.run(["fastq-dump", "--split-files", "--gzip", "--outdir", SRA_path, SRA], check=True)
                    # renaming the file as FASTQ
                    for file in os.listdir(SRA_path):
                        if file.startswith(SRA) and file.endswith(".fastq.gz"):
                            os.rename(os.path.join(SRA_path, file), os.path.join(SRA_path, "FASTQ.gz"))
                except subprocess.CalledProcessError as e:
                    print(f"An error occurred while downloading {SRA}: {e}")                
print("Thank you! Your transcriptomes are ready for processing.")