
## Ask the user for brief description of variables and conditions?
    ## could be stored in text files to keep track of data meaning but not necessary
## Ask for a refernce genome?
## Needs to be modified to download Fastq files (not just FASTA), but this is done differently


# For this script, we need gitpython and biopython installed

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
from Bio import Entrez, SeqIO
    # 'Entrez' allows us to send requests to NCBI databases
    # SeqIO is used for reading sequences
import os
    # we need this to edit folders

# Setting up path to save folders in & email to access NCBI
path = str(input("Please copy-paste a path link for where you would like this program's functioning data stored: "))
repo = Repo.clone_from("https://github.com/xucatherine/bioinformatic-pipeline/src.git", path)
Entrez.email = str(input("NCBI requires an email address to track usage of their services.\nPlease input your email address: "))
    # NCBI Entrez requires an email address (according to BioBuddy)
    ## Not sure if email address needs to be linked to an account already

# Setting up Samples folder (within src)
Samples_path = path+"/Samples"
os.makedirs(Samples_path, exist_ok=True)

# Adding condition folders
print("Using differential analysis to decode metabolic pathways involves identifying experimental variables under which expression of the phenotype of interest differs.")
print("Through changing the variable's intensity, we create different experimental conditions under which transcriptomes can be collected and compared against each other.")
print("Ex: In studying yeast respiration, temperature is a variable. Low, medium and high heat would correspond to 3 conditions, between which cellular respiration varies.")
print("\nHow many variables are you considering for your pathway? (minimum 1)")
print("For each variable, you will be asked to input the number of conditions studied, and for each codnition list the corresponding SRA accession numbers.")
n = input("number of variables studied: ") #[user inputs answer]
# make a variable folder for each variable, with 'p_SRAs' and 'n_SRAs' set
for i in range(n):
    var_path = Samples_path+"/var_"+str(n) # folders will be named var_1, var_2, var_3...
    os.makedirs(var_path, exist_ok=True)
    k = input("number of conditions studied for variable " + str(n) + ': ') #[user inputs answer]
    for j in range(k):
        cond_path = var_path+"/cond_"+str(k) # folders will be named cond_1, cond_2, cond_3...
        os.makedirs(cond_path, exist_ok=True)

# Extracting SRA FASTA files for each variable
for i in range(n):
    var_path = Samples_path+"/var_"+str(n)
    # counting conditions within var_n - there may be a more efficient way of doing this
    cond_count = 0
    for cond in os.listdir(var_path): # for condition folder in parent folder var_n
        if os.path.isdir(os.path.join(var_path, cond)):
            cond_count += 1          
    # getting SRAs
    for cond in os.listdir(var_path): # for condition folder in parent folder var_n
        folder_name = os.path.basename(os.getcwd()) # this should give "cond_k", from which we can get k
        k = int(folder_name.split("_")[1]) ## extract k from current cond
        print("List the SRA numbers belonging to variable "+str(n)+"'s condition "+k+", seperated by spaces.")
        print("ex: SRR12345678 SRR91011109 SRR87654321")
        SRAs = str(input("SRA numbers: ")) #[user inputs SRAs]
        # get list of SRAs
        SRAs_list = SRAs.split()
        # retrieve FASTA files for SRAs and store in folder
        for i in len(SRAs_list):
            SRA = SRAs_list[i-1] # should give something like SRR12345678
            # new SRA folder
            SRA_path = Samples_path+"/var_"+str(n)+"/cond_"+str(k)+"/SRAs"+"/"+SRA # new folder named as SRA number
            os.makedirs(SRA_path, exist_ok=True)
            # retrieve FASTA file (script by ChatGPT)
            handle = Entrez.efetch(db="nucleotide", id=SRA, rettype="fasta", retmode="text")
            fasta_data = handle.read()
            handle.close()
            # add FASTA file to folder
            with open(SRA_path+"/FASTA", "w") as file:
                file.write(fasta_data)
    # done!     

print("Thank you! Your transcriptomes are ready for processing.")