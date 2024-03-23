
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
                file: counts.csv
                [...]
            [more SRR folders]
        [more cond folders]
    [more var folders]

folder: References
    file: genome.gtf
    folder: compiled_counts
        folder: var_x
            comp_counts.csv
            norm_counts.csv
            vst_counts.csv
            metadata.csv
        [more var folders]
    [...]

folder: Results
    DESeq2
        DESeq2_results_var_1.csv
        DESeq2_results_var_2.csv
        [results for each var_x]
    [...]
'''

# Imports
# from git import Repo
    # using gitpython we can access GitHub to make new folders for user's data
from Bio import Entrez
    # 'Entrez' allows us to send requests to NCBI databases (from BioPython)
import os
    # we need this to edit folders
import subprocess
    # we need this to call fastq-dump

from A_minis import Bioinf_Profile

print("\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("\t~~~ Welcome to our bioinformatic pipeline! ~~~")
print("\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

Steps = {
    "A": "Sample Setup and Download",
    "B": "Quality Check",
    "C": "Trimming",
    "D": "Reference-based Assembly/Mapping",
    "E": "Abundance Estimation",
    "F": "De novo Analysis with Seq2Fun",
    "G": "Differential Expression Analysis",
    "H": "Coexpression Analysis"
}

# Initialize Bioinf_profile object 
profile = Bioinf_Profile()
# Path variables that are assigned values within the setup functions
# and that will be used in multiple places in MAIN.py:
# These are listed here more to help remember which variables are meant to be global
SRA_toolkit_path, Samples_path, References_path, Results_path, SRR_paths = None, None, None, None, None
VAR_names, COND_names = None, None
    # VAR_names = [var_1 name, var_2 name, ...]
    # COND_names = [[cond_1 name, cond_2 name, ...],  (var_1)
    #               [cond_1 name, cond_2 name, ...],  (var_2)
    #               ...]

def paths_setup():
    # Check if the pipeline has already been run & information has been saved
    global SRA_toolkit_path, Samples_path, References_path, Results_path, SRR_paths
    if os.path.isfile(".bioinf-profile"):
        profile.read_profile()

        print("It looks likes you ran the bioinformatic pipeline previously")
        print(f"and reached Step {profile.dict["STEP"]}: {Steps[profile.dict["STEP"]]}")
        print("Would you like to...")
        print(f"\t (1) continue the pipeline starting from Step {profile.dict["STEP"]}: {Steps[profile.dict["STEP"]]}?")
        print("\t (2) or rerun the pipeline from the start?")
        m = input("Please enter 1 or 2: ")
        if m == '1': 
            start_over = False
            path = profile.dict["PATH"]
            SRA_toolkit_path = profile.dict["SRA_TOOLKIT_PATH"]
            Samples_path = path+"/Samples"
            References_path = path+"/References"
            Results_path = path+"/Results"
        elif m == '2': start_over = True
    else: # if the .bioinf-profile file doesn't exist, run the pipeline from the start
        start_over = True
    
    if start_over: 
        # only runs if user has never run the pipeline before, or if user wants to restart the pipeline
        # Setting up path to save folders in & user email to access NCBI
        # this loops until the user inputs a valid path
        while True:
            path = str(input("Please enter a path for where you would like this pipeline's running data and results stored: "))
            ## We should maybe set a default to use bioinformatic-pipeline
            if not os.path.isdir(path):
                print("Error: invalid directory path. Please enter a valid path to a directory/folder.")
            else:
                break
        #repo = Repo.clone_from("https://github.com/xucatherine/bioinformatic-pipeline/src.git", path)
        #Entrez.email = str(input("\nNCBI requires an email address to track usage of their services.\nPlease input your email address: "))
            # NCBI Entrez requires an email address (according to BioBuddy)
            ## Not sure if email address needs to be linked to an account already

        # Ask the user to install the SRA-toolkit from NCBI. Ask the user for the path of the folder.
        print("\nTo obtain FastQ files from NCBI SRRs, we will be using NCBI's SRA Toolkit.")
        print("Please download the SRA Toolkit zip file for your operating system from \nthe following page: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit")
        print("Then, extract the zip file where you would like to store the toolkit.")
        SRA_toolkit_path = input("Please input the path to the SRA Toolkit folder (e.g. /Users/..../sratoolkit.3.x.x-mac-x86_64):")

        # On MacOS, the folder will automatically be flagged as coming from an unknown developper and this script will not be able to open the folder
        # check for the com.apple.quarantine flag and remove it if the folder has it
        xattr_list=subprocess.run(["xattr", SRA_toolkit_path], capture_output=True, text=True).stdout.split('\n')
        if "com.apple.quarantine" in xattr_list:
            subprocess.run(["xattr", "-r", "-d","com.apple.quarantine", SRA_toolkit_path])
        
        # Setting up Samples, References and Results folders (within the user's inputted path)
        Samples_path = path+"/Samples"
        os.makedirs(Samples_path, exist_ok=True)
        References_path = path+"/References"
        os.makedirs(References_path, exist_ok=True)
        Results_path = path+"/Results"
        os.makedirs(Results_path, exist_ok=True)

        # update the variables in .bioinf-profile
        profile.dict["PATH"] = os.path.abspath(path)
        profile.dict["SRA_TOOLKIT_PATH"] = os.path.abspath(SRA_toolkit_path)
        profile.update_profile()
    
    ## Always Runs ##
    # Creating list of SRRs for easy processing - iterates through all folders in Samples
    SRR_paths = []
    for var in os.listdir(Samples_path): # for each variable in the Samples folder
        var_path = os.path.join(Samples_path, var)
        for cond in os.listdir(var_path): # for each condition in the var_n folder
            cond_path = os.path.join(var_path, cond)
            for SRR in os.listdir(cond_path): # for each SRR in the cond_k folder
                SRR_path = os.path.join(cond_path, SRR)
                SRR_paths.append(SRR_path)


def samples_setup():
    fasterq_dump_path = os.path.join(SRA_toolkit_path, "bin", "fasterq-dump") # path to the tool that will fetch fastQ files

    # Entering data into Samples folder - starting with number of variables to study
    print("\nUsing differential analysis to decode metabolic pathways involves identifying \nexperimental variables under which expression of the phenotype of interest differs.")
    print("Through changing the variable's intensity, we create different experimental conditions \nunder which transcriptomes can be collected and compared against each other.")
    print("\tEx: In studying yeast respiration, temperature is a variable. Low, medium and high heat \n\twould correspond to 3 conditions, between which cellular respiration varies.")
    print("\nHow many variables are you considering for your pathway? (minimum 1)")
    print("For each variable, you will be asked to input the number of conditions studied, \nand for each condition, list the corresponding SRR accession numbers.\n")
    while True: # continue asking for input until they input a valid integer
        try:
            n = int(input("Number of variables studied (min. 1): ")) #[user inputs answer]
        except ValueError:
            print("Error: Please enter a valid, positive integer.")
            continue # continue looping
        else:
            if n < 1: # if the value they enter give is less than 1
                print("Error: Please enter an integer that is greater than or equal to 1.")
                continue # continue looping
            else: break # otherwise, the number is valid; break the loop

    ## initialize a dictionary here to which var and cond names can be stored?
    # Initialize global variables to store the variable names and conditions

    # For each variable, make the right number of conditions folders and add data
    for i in range(n):
        print("\n")
        var_path = Samples_path+"/var_"+str(i+1) # folders will be named var_1, var_2, var_3...
        os.makedirs(var_path, exist_ok=True)
        ## ask user to name variable and store name in a dictionary?

        # Adding conditions folders
        while True:
            try:
                k = int(input("Number of conditions studied for variable " + str(i+1) + ' (min. 2): ')) #[user inputs answer]
            except ValueError:
                print("Error: Please enter a valid, positive integer.")
                continue # continue looping
            else:
                if k < 2: # if the value they enter give is less than 2
                    print("Error: Please enter an integer that is greater than or equal to 2.")
                    continue # continue looping
                else: break # otherwise, the number is valid; break the loop

        for j in range(k):
            cond_path = var_path+"/cond_"+str(j+1) # folders will be named cond_1, cond_2, cond_3...
            os.makedirs(cond_path, exist_ok=True)
            ## ask user to name condition and store name in a dictionary?

        # Getting SRRs and downloading their FASTQ files, for each variable
        SRRs_matrix = []
        for cond in sorted(os.listdir(var_path)): # for condition folder in parent folder var_n
            print(f"\nList the SRR numbers belonging to variable {i+1}'s condition {cond.split('_')[-1]}, separated by spaces.")
            print("ex: SRR12345678 SRR91011109 SRR87654321")
            SRRs = input("SRR numbers: ") #[user inputs SRRs]
            SRRs_list = SRRs.split() # making list of SRRs
            SRRs_matrix.append(SRRs_list)
            
        print("\nWe will now download the FastQ files for each SRR. Each SRR can take ~5-10 minutes.")

        # downloading the FASTQ files
        cond_list = sorted(os.listdir(var_path))
        for i in range(len(cond_list)): # for condition folder in parent folder var_n
            # Starting the downloads
            print(f"\nDownloading SRRs for condition {i+1}")
            cond_path = os.path.join(var_path, cond_list[i])
            SRRs_list = SRRs_matrix[i] # SRRs inputted by the user for the current condition
            for SRR in SRRs_list: # loop through each SRR
                SRR_path = os.path.join(cond_path, SRR) # new folder named as SRR number
                os.makedirs(SRR_path, exist_ok=True)

                # Check if the fasta files have already been downloaded in this folder
                if os.path.isfile(os.path.join(SRR_path, "rawF.fastq")) and os.path.isfile(os.path.join(SRR_path,"rawR.fastq")):
                    print(f"It looks like the fastq files for {SRR} have already been downloaded!")
                    print("Skipping to the next SRR...")
                    continue
                
                print(f"Downloading {SRR}...")
                # retrieve FASTQ file using SRR Toolkit
                try:
                    subprocess.run([fasterq_dump_path, SRR, "--outdir", SRR_path, "--skip-technical", "--split-3"])
                        # this will create two, or three files:
                            # SRR000000_1.fastq ;  side 1 of paired reads
                            # SRR000000_2.fastq ;  side 2 of paired reads
                            # SRR000000_x.fastq ;  any unpaired/unmatched reads
                    # renaming the files to generic name
                    for file in os.listdir(SRR_path):
                        if file.startswith(SRR):
                            if file.endswith("_1.fastq"):
                                new_name = "rawF.fastq"
                                os.rename(os.path.join(SRR_path, file), os.path.join(SRR_path, new_name))
                            elif file.endswith("_2.fastq"):
                                new_name = "rawR.fastq"
                                os.rename(os.path.join(SRR_path, file), os.path.join(SRR_path, new_name))
                            else:
                                os.remove(os.path.join(SRR_path, file)) # Delete any extra files that aren't the reverse or forward reads
                except subprocess.CalledProcessError as e:
                    print(f"An error occurred while downloading {SRR}: {e}") 
            
    print("Thank you! Your transcriptomes are ready for processing.")

paths_setup()

## Use if statements to jump ahead in the code
profile.dict["STEP"] = "A"

if profile.dict["STEP"] == "A": 
    print("\n\t################################")
    print("\t# A: Sample Setup and Download #")
    print("\t################################")
    samples_setup()
    profile.dict["STEP"] = "B"

######################################
# Running Quality Checks on SRR data #
######################################

if profile.dict["STEP"] == "B":
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
                cond_path = os.path.join(var_path, cond)
                B_Quality_Check.run_MultiQC(cond_path)
    
    profile.dict["STEP"] = "C"

# Trimming SRR data using Quality Check results
import C_Trimming
for SRR in SRR_paths:
    C_Trimming.run_cutadapt(SRR)

# Next...
## Ask if there is a genome somewhere, jump to Reference-based assembly or De novo

#################################### 
# Reference-based assembly/mapping #
####################################



#####################
# De novo - Seq2Fun #
#####################



####################################
# Differential Expression Analysis #
####################################

from I_differential_exp import setup_genecount_matrix, organize_DESeq2_genecounts, run_DESeq2_R

print("\nDifferential Expression Analysis will now be performed using DESeq2.")
print("For one variable, DESEq2 will compare pairs of conditions, using the gene counts from all of the conditions")
print("to calculate the 'degree' of differential expression for each gene.")
print("This degree of differential expression is represented by statistical numbers like log2-fold change and p-value.")
print("For more information on how to interpret the DESeq2 results:")
print("https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05b_wald_test_results.html#p-values")
print("For more information on how DESeq2 works:")
print("https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html")
    
# Folders to save information for DESeq2
    # References
        # compiled_counts
            # var_1
                # comp_counts.csv
                # metadata.csv
                # norm_counts.csv
                # vst_counts.csv
            # var_2
                # ...
    # Results
        # DESeq2
            # DESeq2_results_var_1.csv
            # DESeq2_results_var_2.csv
            # ...

# make DESeq2 results folder
deseq2_results_folder = os.path.join(Results_path, "DESeq2")
if not os.path.exists(deseq2_results_folder): os.mkdir(deseq2_results_folder) # make the folder if it doesn't exist already

# make compiled_counts folder
compiled_counts_path = os.path.join(References_path, "compiled_counts")
if not os.path.exists(compiled_counts_path): os.mkdir(compiled_counts_path)
    
# make the sub folders under compiled_counts, one for each variable
var_list = os.listdir(Samples_path) # use Samples folder to know how many variables/names of the variable folders
for var_dir in var_list:
    compiled_counts_var_path = os.path.join(compiled_counts_path, var_dir)
    if not os.path.exists(compiled_counts_var_path): os.mkdir(compiled_counts_var_path) 
    
# Setup data for DESeq2. 
for var_dir in var_list:
    samples_var_path = os.path.join(Samples_path, var_dir)
    compiled_counts_var_path = os.path.join(compiled_counts_path, var_dir)
    organize_DESeq2_genecounts(samples_var_path, results_folder=compiled_counts_var_path)
        # Compiles all of the counts from Samples/var_n into 1 table (comp_counts.csv)
        # For each sample/SRR in the table, records the corresponding condition in metadata.csv

# Run DESeq2 for each var
for var_dir in var_list:
    comp_counts_csv_path = os.path.join(compiled_counts_path, var_dir, "comp_counts.csv")
    metadata_path = os.path.join(compiled_counts_path, var_dir, "metadata.csv")
    result_file = "DESeq2_results_" + var_dir
    deseq2_result_file_path = os.path.join(deseq2_results_folder, result_file)
    compiled_counts_var_path = os.path.join(compiled_counts_path, var_dir)

    run_DESeq2_R(comp_counts_csv_path, metadata_path, deseq2_result_file_path, compiled_counts_var_path, "TRUE", "TRUE")


##########################
# Co-Expression Analysis #
##########################
    
