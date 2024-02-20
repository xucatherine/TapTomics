# Should we ask the user for brief description of conditions and positive/negative sets?
# could be stored in text files to keep track of data meaning but not necessary


# Folder Architecture
'''
folder: Samples
    folder: condition_x
        folder: pSRAs
            folder: SRA (subfiles can be added later)
                file: FASTA
            [more SRA folders]
        folder: nSRAs
            folder: SRA
                file: FASTA
            [more SRA folders]
    folder: condition_y
        [...]
    [more conditions folders]
'''

# Installs
pip install gitpython
from git import Repo
    # using gitpython we can access GitHub to make new folders for user's data
pip install biopython
from Bio import Entrez, SeqIO
    # 'Entrez' allows us to send requests to NCBI databases
    # SeqIO is used for reading sequences
import os
    # we need this to edit folders

# Setting up path to save folders in & email to access NCBI
path = str(input("Please copy-paste the path link where you would like this program's functioning data stored: "))
repo = Repo.clone_from("https://github.com/xucatherine/bioinformatic-pipeline/src.git", path)
Entrez.email = str(input("NCBI requires an email address to track usage of their services.\nPlease input your email address: "))
    # NCBI Entrez requires an email address (according to BioBuddy)
    ## Not sure if email address needs to be linked to an account already

# Setting up Samples folder (within src)
Samples_path = path+"/Samples"
os.makedirs(Samples_path, exist_ok=True)

# Adding condition folders
print("Using differential analysis to decode metabolic pathways involves identifying conditions under which expression of the phenotype of interest differs.")
print("Each condition thus has a positive and negative set, where positive means the phenotype of interest is increased/expressed/observed.")
print("Ex: In studying yeast respiration, temperature is a condition, where heat is positive (increased respiration) and freezing would be negative.")
print("\nHow many differential conditions are you considering for your pathway? (minimum 1)")
print("For each condition, you will be asked to list the SRA accession numbers for the positive and negative sets.")
n = input("number of conditions studied: ") #[user inputs answer]
# make a condition folder for each condition, with 'p_SRAs' and 'n_SRAs' set
for i in range(n):
    cond_path = Samples_path+"/condition_"+str(n) # folders will be named condition_1, condition_2, condition_3...
    os.makedirs(cond_path, exist_ok=True)
    p_SRAs_path = cond_path+"/p_SRAs"
    os.makedirs(p_SRAs_path, exist_ok=True)
    n_SRAs_path = cond_path+"/n_SRAs"
    os.makedirs(n_SRAs_path, exist_ok=True)

# Extracting SRA FASTA files for each condition
for i in range(n):
    print("List the SRA numbers of your positive set (phenotype of interest increased/expressed/observed), seperated by spaces.")
    print("ex: SRR12345678 SRR91011109 SRR87654321")
    p_SRAs = str(input("SRA numbers of your positive set: ")) #[user inputs pSRAs]
    n_SRAs = str(input("SRA numbers of your negative set, seperated by spaces: ")) #[user inputs nSRAs]
    # get list of SRAs
    p_SRAs_list = p_SRAs.split()
    n_SRAs_list = n_SRAs.split()
    # retrieve FASTA files for pSRAs and store in folder
    for i in len(p_SRAs_list):
        SRA = p_SRAs_list[i-1] # should give something like SRR12345678
        # new SRA folder
        SRA_path = Samples_path+"/condition_"+str(n)+"/p_SRAs"+"/"+SRA # new folder named as SRA
        os.makedirs(SRA_path, exist_ok=True)
        # retrieve FASTA file (script by ChatGPT)
        handle = Entrez.efetch(db="nucleotide", id=SRA, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        # add FASTA file to folder
        with open(SRA_path+"/FASTA", "w") as file:
            file.write(fasta_data)
    # same thing but for nSRAs
    for i in len(n_SRAs_list):
        SRA = n_SRAs_list[i-1] 
        SRA_path = Samples_path+"/condition_"+str(n)+"/n_SRAs"+"/"+SRA
        os.makedirs(SRA_path, exist_ok=True)
        handle = Entrez.efetch(db="nucleotide", id=SRA, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        with open(SRA_path+"/FASTA", "w") as file:
            file.write(fasta_data)
    # done!     

print("Thank you! Your transcriptomes are ready for processing.")