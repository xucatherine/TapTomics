
    ## Before making it to this script we should have a preface for the user
    ## it should explain how the pipeline works and how to select RNA data

# RNA transcriptome Selection and Import - My initial thoughts
    # asks the user for the NCBI accession numbers of transcriptomes of interest
    # uses 'Download & Extract' to retrieve said data
        # makes a new data folder
        # ensures it's not dowloading duplicate data 
    # places downloaded RNA transcripts in proper folders


# put folders under src but seperate from others

# folder: Samples
    # folder: Condition_x
        # folder: Positives
             # folder: SRA (for however many you have)
                 # data: raw FASTA file
                 # data: other things added later
        # folder: Negatives
             # folder: SRA (for however many you have)
                 # data: raw FASTA file
                 # data: other things added later
    # folder: Condition_y

# Should we ask the user for brief description of conditions and positive/negative sets?
# these could be stored in text files and help someone reading through the files understand, 
# but it's not necessary


pip install biopython
from Bio import Entrez
    # from Biopython we import 'Entrez' which allows us to send requests to NCBI databases

Entrez.email = input("NCBI requires an email address to track usage of their services."
                     "\nPlease input your email address: \n")
    # NCBI Entrez requires an email address (according to BioBuddy)
    ## Not sure if email address needs to be linked to an account already

print("Using differential analysis to decode metabolic pathways involves identifying conditions under which expression of the phenotype of interest differs.")
print("Each condition thus has a positive and negative set, where positive means the phenotype of interest is increased/expressed/observed.")
print("Ex: In studying yeast respiration, temperature is a condition, where heat is positive (increased respiration) and freezing would be negative.")
print("\nHow many differential conditions are you considering for your pathway? (minimum 1)")
print("For each condition, you will be asked to list the SRA accession numbers for the positive and negative sets.")
n = input("number of conditions studied: ") #[user inputs answer]
## create 'Samples' and condition folders
## create 'p_SRAs' and 'n_SRAs' within each condition

for i in n: # for n amount of times (per condition)
    p_SRAs = input("List the SRA numbers of your positive set (phenotype of interest increased/expressed/observed), seperated by spaces: ") #[user inputs SRAs]
    n_SRAs = input("List the SRA numbers of your negative set, seperated by spaces: ") #[user inputs SRAs]
    ## these lists get split by spaces
    ## using SRA, FASTA file is recovered
    ## data is stored within appropriate condition folder

print("Thank you!")

## should a 'break' be integrated to restart the code early incase the user makes a typo?
