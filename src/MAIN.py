
## Missing: ask the user for brief description each variable and condition, and store entries in a dictionary

# Welcome to the MAIN script for the Bioinformatic Pipeline! 
# The user runs this script, and is guided through the steps of transcriptomic analysis

# Folder Architecture - for reference; this is how downloaded data is stored
'''
folder: Samples
    folder: var_x    # var_1, var_2, var_3, ...
        folder: cond_y    # cond_1, cond_2, cond_3,...
            folder: SRR
                file: rawF.fastq
                file: rawR.fastq
                file: trimmedF.fastq
                file: trimmedR.fastq
                file: aligned.bam
                file: counts.csv
                folder: rawF_fastqc
                    file: fastqc_data.txt
                folder: rawR_fastqc
                    file: fastqc_data.txt
            [more SRR folders]
        [more cond folders]
    [more var folders]

folder: References
    file: genome.fasta
    file: genome.gtf
    folder: compiled_counts
        folder: var_x
            file: comp_counts.csv
            file: norm_counts.csv
            file: vst_counts.csv
            file: metadata.csv
        [more var folders]
    file: database

folder: Results
    folder: DESeq2
        file: DESeq2_results_var_1.csv
        file: DESeq2_results_var_2.csv
        [results file for each var_x]
    folder: PyWGCNA
        folder: Resultsfigures
            file: module-traitRelationships.pdf
            file: summary_power
            file: eigengenes.pdf
            file: module_heatmap_eigengene_[module name].pdf
                [module_heatmap_eigengene file for every module]
            file: network.pdf
        folder: Results
            file: pathway_modules.csv
            file: [module name]_genes.csv
                [genes file for every pathway module]
    folder: Seq2Fun
        file: Seq2Fun_summary_all_samples.html
'''

# Imports
# from git import Repo
    # using gitpython we can access GitHub to make new folders for user's data
#from Bio import Entrez
    # 'Entrez' allows us to send requests to NCBI databases (from BioPython)
import os
    # we need this to edit folders
import subprocess
    # we need this to call fastq-dump

from A_minis import Bioinf_Profile, listdir_visible, name

print("\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("\t~~  Welcome to Taptomics!  ~~")
print("\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

Steps = {
    "A": "Sample Setup and Download",
    "B": "Quality Check",
    "C": "Trimming",
    "D": "Reference-based Assembly/Mapping",
    "E": "De novo Analysis with Seq2Fun",
    "F": "Differential Expression Analysis",
    "G": "Coexpression Analysis"
}

# Initialize Bioinf_profile object 
profile = Bioinf_Profile()

###########################################################################
# Important variables that are assigned values within the setup functions #
###########################################################################
# and that will be used in multiple places in MAIN.py:
# These are listed here more to help remember which variables are meant to be global
SRA_toolkit_path, Samples_path, References_path, Results_path, SRR_paths = None, None, None, None, None
VAR_names, COND_names, COND_quant, COND_corr = [], [], [], []
    # VAR_names = [var_1 name, var_2 name, ...]
    # COND_names = [[cond_1 name, cond_2 name, ...],  (var_1)
    #               [cond_1 name, cond_2 name, ...],  (var_2)
    #               ...]
    # COND_corr = [[cond_1 correlation, cond_2 correlation, ...],  (var_1)
    #               [cond_1 correlation, cond_2 correlation, ...],  (var_2)
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
        print("Then, unzip the zip file where you would like to store the toolkit.")
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
        profile.dict["STEP"] = "A"
        profile.update_profile()
    
    ## Always Runs ##
    # Creating list of SRRs for easy processing - iterates through all folders in Samples
    SRR_paths = []
    for var in sorted(listdir_visible(Samples_path)): # for each variable in the Samples folder
        var_path = os.path.join(Samples_path, var)
        for cond in sorted(listdir_visible(var_path)): # for each condition in the var_n folder
            cond_path = os.path.join(var_path, cond)
            for SRR in sorted(listdir_visible(cond_path)): # for each SRR in the cond_k folder
                SRR_path = os.path.join(cond_path, SRR)
                SRR_paths.append(SRR_path)


def sample_setup():
    fasterq_dump_path = os.path.join(SRA_toolkit_path, "bin", "fasterq-dump") # path to the tool that will fetch fastQ files

    # Entering data into Samples folder - starting with number of variables to study
    print("\nUsing differential analysis to decode metabolic pathways involves identifying \nexperimental variables under which expression of the phenotype of interest differs.")
    print("Through changing the variable's intensity, we create different experimental conditions \nunder which transcriptomes can be collected and compared against each other.")
    print("\tEx: In studying yeast respiration, temperature is a variable. Low, medium and high heat \n\twould correspond to 3 conditions, between which cellular respiration varies.")
    print("\nFor each variable, you will be asked to input the number of conditions studied,")
    print("and for each condition, to specify if it correlates positively or negatively with the phenotype of interest")
    print("\tand to list the corresponding SRR accession numbers.\n")
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

    # For each variable, make the right number of conditions folders and add data & ask for SRRs
    SRRs_matrix = [] # for saving all of the SRRs
            # [["SRR000001 SRR000002 SRR000003 ...", "SRR000004 SRR000005 SRR000006 ...", ...], (var_1)
            #   ["SRR000011 SRR000012 SRR000013 ...", "SRR000014 SRR000015 SRR000016 ...", ...], (var_2)
            #           (cond_1)                            (cond_2)                      (cond_n)
    for i in range(n):

        print("\n")
        var_path = Samples_path+"/var_"+str(i+1) # folders will be named var_1, var_2, var_3...
        os.makedirs(var_path, exist_ok=True)
        ## ask user to name variable and store name in a list
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print(f"Now collecting data for VARIABLE {i+1}.\n")
        print(f"You may input a label for variable {i+1} to be used in the final ouputs.") 
        print(f"(examples: light, oxidative stress. If the field is left empty, we will just use var_{i+1}.)")
        var_name = input(f"Label for variable {i+1}: ")
        if var_name.strip() == "": var_name = "var_" + str(i+1) # name variable var_i+1 if user left input empty
        VAR_names.append(var_name) # 'masterlist' for the variable names

        print("")

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

        # For each condition, make folder
        for j in range(k): 
            cond_path = var_path+"/cond_"+str(j+1) # folders will be named cond_1, cond_2, cond_3...
            os.makedirs(cond_path, exist_ok=True)

        ## Ask for name & SRRs for every condition of the current variable ##
        
        SRRs_list_temp = [] # temporary list to save SRRs for current variable
        cond_names_temp = [] # temporary list to store the condition names for the current variable 
        cond_corr_temp = [] # temporary list to store the correlation types for the current variable
                            # (these will later be appended to the SRRs_matrix, COND_names, and COND_corr)

        for cond in sorted(listdir_visible(var_path)): # for condition folder in parent folder var_n
            print('\n')
            cond_n = cond.split('_')[-1]
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(f"Now collecting data on VARIABLE {i+1}'s CONDITION {cond_n}.\n")
            print("You may input a label for the condition to be used in the final outputs.")
            print(f"(examples: high light or low light. If the field is left empty, we will just use cond_{cond_n}.")
            cond_name = input(f"Label for variable {i+1}'s condition {cond_n}: ")
            if cond_name.strip() == "": cond_name = "cond_" + str(cond_n) # name condition cond_n if input was left blank
            cond_names_temp.append(cond_name) 
            
            print("")
            while True:
                corr_type = input(f"Correlation for variable {i+1}'s condition {cond_n} (positive/negative): ")
                if corr_type.lower() == "positive" or corr_type.lower() == "negative": break
                else: print("Error: please input 'positive' or 'negative'")
            cond_corr_temp.append(corr_type.lower())

            print(f"\nList the SRR numbers belonging to variable {i+1}'s condition {cond_n}, separated by spaces.")
            print("ex: SRR12345678 SRR91011109 SRR87654321")
            SRRs = input("SRR numbers: ") #[user inputs SRRs]
            SRRs_list_temp.append(SRRs) # add the SRRs to the list for current variable
        
        SRRs_matrix.append(SRRs_list_temp) # append 

        COND_names.append(cond_names_temp) # add condition names to the 'masterlist'
        COND_corr.append(cond_corr_temp)
    
    profile.VAR_names=VAR_names
    profile.COND_names=COND_names
    profile.COND_cor=COND_corr

    #### Implement a step that saves these variables^^ to the .bioinf-profile ####
    
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("\nWe will now download the FastQ files for each SRR. Each SRR can take ~5-10 minutes.")

    # Downloading ALL of the FASTQ files
    for i in range(n): # iterate through the variables
        var_path = Samples_path+"/var_"+str(i+1)
        cond_list = sorted(listdir_visible(var_path))

        for j in range(len(cond_list)): # for condition folder in parent folder var_n
            # Starting the downloads
            print(f"\nDownloading SRRs for VARIABLE {i+1} CONDITION {j+1}")
            cond_path = os.path.join(var_path, cond_list[j])
            SRRs_list = SRRs_matrix[i][j].split() # SRRs inputted by the user for the current condition
            for SRR in SRRs_list: # loop through each SRR
                SRR_path = os.path.join(cond_path, SRR) # new folder named as SRR number
                os.makedirs(SRR_path, exist_ok=True)

                # Check if the fastq files have already been downloaded in this folder
                if os.path.isfile(os.path.join(SRR_path, "rawF.fastq")) and os.path.isfile(os.path.join(SRR_path,"rawR.fastq")):
                    print(f"It looks like the fastq files for {SRR} have already been downloaded! Skipping to the next SRR...")
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

if profile.dict["STEP"] == "A": 
    print("\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("\t~ A: Sample Setup and Download ~")
    print("\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    sample_setup()
    profile.dict["STEP"] = "B"
    profile.update_profile()

##################################
# Run Quality Checks on SRR data #
##################################

if profile.dict["STEP"] == "B":
    print("\n\t~~~~~~~~~~~~~~~~~~~")
    print("\t~ B: Quality Check ~")
    print("\t~~~~~~~~~~~~~~~~~~~~\n")

    # Add brief explanation about quality checks?

    import B_Quality_Check
    if profile.dict["FASTQC_PATH"] == "": # Check if no path has been entered before
        ## Add download instructions - or maybe do this in path_setup?
        print("From the FastQC website (https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)")
        print("download and unzip the zip file 'FastQC vX.XX.X (Win/Linux zip file)' (even if you are on MacOS).")
        fastqc_path = input("Please enter the path to the FastQC folder (e.g. User/..../FastQC): ") # path to previously downloaded FastQC, within user's computer
            ## save this path earlier when dowloading, or have user manually input it here
        # check for the com.apple.quarantine flag and remove it if the folder has it
        xattr_list=subprocess.run(["xattr", fastqc_path], capture_output=True, text=True).stdout.split('\n')
        if "com.apple.quarantine" in xattr_list:
            subprocess.run(["xattr", "-r", "-d","com.apple.quarantine", SRA_toolkit_path])

        profile.dict["FASTQC_PATH"] = os.path.abspath(fastqc_path)
        profile.update_profile()
    else: 
        fastqc_path = profile.dict["FASTQC_PATH"]

    fastqc_exe = os.path.join(fastqc_path, "fastqc") # the executable for fastqc is located inside the FastQC folder
    print("We will now run FastQC on each SRR's .fastq files.")
    for SRR in SRR_paths:
        print(f"Running FastQC on {name(SRR)}...")
        # Check if FastQC has already been run on this SRR by checking for existence of rawF_fastqc & rawR_fastqc folders
        if os.path.isdir(os.path.join(SRR, 'rawF_fastqc')) and os.path.isdir(os.path.join(SRR, 'rawR_fastqc')): 
            print(f"It looks like FastQC was already run on {name(SRR)}! Skipping to the next SRR...")
            continue
        B_Quality_Check.run_FastQC(fastqc_exe,SRR)
    print("Done analyzing all SRRs with FastQC!")

    # Optional: Running MultiQC on FastQC data to create visual for user, per condition
    print("\nWould you like to run a visualizer for the quality checks done on each of your conditions?")
    print("This will output links which you can open in a browser.")
    query = input("Run MultiQC? y/n: ")
    if query == 'y':
        for var in listdir_visible(Samples_path): # for each variable in the Samples folder
            var_path = os.path.join(Samples_path, var)
            for cond in listdir_visible(var_path): # for each condition in the var_n folder
                cond_path = os.path.join(var_path, cond)
                B_Quality_Check.run_MultiQC(cond_path)
                #### MULTIQC saves the results in the current working directory automatically ###
                ### maybe we move it to Results ####
                ### also it currently tells the user to open the file in your browser, but nothing shows up###
    
    profile.dict["STEP"] = "C"
    profile.update_profile()

#################################################
# Trimming SRR data using Quality Check results #   !Unfinished!
#################################################

if profile.dict["STEP"] == "C":
    print("\n\t~~~~~~~~~~~~~~")
    print("\t~ C: Trimming ~")
    print("\t~~~~~~~~~~~~~~~\n")
    import C_Trimming
    for SRR in SRR_paths:
        C_Trimming.run_cutadapt(SRR)

    ## Ask if there is a genome here, jump to Reference-based assembly or De novo
    
    profile.dict["STEP"] = "D"
    profile.update_profile()

# Next...
## Ask if there is a genome somewhere, jump to Reference-based assembly or De novo

#################################### 
# Reference-based assembly/mapping #    !Unfinished!
####################################
if profile.dict["STEP"] == "D":
    print("\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`")
    print("\t~ D: Reference-Based Assembly ~")
    print("\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print('For this section you must have RNA STAR installed and made on your local computer. You will also need the genome and annotations for your organism of interest downloaded in fasta and gtf files respectively. These can be found on NCBI.')
    #initializing the variables we will need
    STARpath = str(input('Please enter the path to RNA STAR: '))
    genomefolder = str(input('Please enter the path to the folder in which the genome files are stored: '))
    genomefasta = str(input('Please enter the path to the genome fasta file: '))
    genomegtf = str(input('Please enter the path to the genome gtf file: '))
    CPUcores = str(input('How many cores would you like STAR to utilize? You can check how many cores your computer has or choose how many are available in the network where you are running it: '))
    import D_ReferenceBased_Assembly
    for SRR in SRR_paths:
        resultpath = SRR
        fastqfiles = [SRR+'/rawF.fastq', SRR+'/rawR.fastq']
        #initializing object of the class
        g = D_ReferenceBased_Assembly.referencemap()
        g.star(STARpath, genomefolder,genomefasta, genomegtf, fastqfiles,resultpath,CPUcores)


    profile.dict["STEP"] = "F"
    profile.update_profile()


#####################
# De novo - Seq2Fun #   !Unfinished!
#####################
if profile.dict["STEP"] == "E":
    from E_De_Novo_Analysis import extract_info_from_paths,run_seq2fun, move_and_rename_files, print_strings_as_table
    #First introduce the tool and ask user to download Seq2Fun
    print("\n\t~~~~~~~~~~~~~~~~~~~~~~")
    print("\t~ E: De novo Analysis ~")
    print("\t~~~~~~~~~~~~~~~~~~~~~~~\n")
    print(" Welcome to the de novo branch of the pipeline! In order to proceed to the de novo analyis, a few steps are required on your end.")
    print("First please download to your working folder and make the Seq2Fun toolkit by following the directives on their github: https://github.com/xia-lab/Seq2Fun/tree/master ")
    print ("I'll let you time to do that...")
    user_input_1 = input("Tell me when you are ready for the next steps by entering - ready - ")
    while True:
            if user_input_1.lower() == 'ready':
                print("Great! Let's proceed.")
                break
            else:
                print("Please type 'ready' when you are ready.")
                user_input_1 = input ("")
    print ("Good, now please input the Path directory to Seq2Fun")
    pathSeq2Fun = input("Input the Seq2Fun path here (ex: /Users/johnnycash/bioinformatic/Seq2Fun):")
    #Maybe write a liine to confirm if the path works or not?
    print("Thanks! We can now move on to the next step! Hang tight we're almost there.")
    #The second step is to ask the user to decide which database he wants to select
    strings = ['algae', 'alveolates', 'amoebozoa', 'amphibians', 'animals', 'apicomplexans', 'arthropods', 'ascomycetes', 'basidiomycetes', 'birds', 'cnidarians', 'crustaceans', 'dothideomycetes', 'eudicots', 'euglenozoa', 'eurotiomycetes', 'fishes', 'flatworms', 'fungi', 'insects', 'leotiomycetes', 'mammals', 'mollusks', 'monocots', 'nematodes', 'plants', 'protists', 'reptiles', 'saccharomycetes', 'stramenopiles', 'vertebrates']
    print ("Please pick the database from the following table that represents the better your sampled organism")
    print_strings_as_table(strings) #Function that prints the table
    DB = input("From the list above,please pick and write in lowercase the appropriate database for Seq2FUN tool to annalyse your reference free transcriptome: ")
    print ("Now please got to https://www.expressanalyst.ca/ExpressAnalyst/docs/Databases.xhtml under the - Without a reference Transcriptome - and download the database")
    print (" See Step 2 on https://github.com/xia-lab/Seq2Fun/tree/master for more info on how to do so")
    print (" Oh yeah, and don't forget to issue the following command: tar -xzvf birds.tar.gz - as expained in Step 2")
    print (f" IMPORTANT: you MUST download the database to the Database folder within Seq2Fun. The path to it should look like: "+ pathSeq2Fun + "/Database" )
    print ("I'll let you time to do that...")
    user_input_2 = input ("Tell me when you are ready for the next steps by entering - ready - ")
    while True:
            if user_input_2.lower() == 'ready':
                print("Great! Let's proceed.")
                break
            else:
                print("Please type 'ready' when you are ready.")
                user_input_2 = input ("")
    print("We're ready to go now, thanks for your time!")
    print("Preparing files for the analysis")
    output_directory = pathSeq2Fun + "/database"
    move_and_rename_files(SRR_paths, output_directory) #This fucntion will take all the inputs and place them in the same folder. It will aslo rename them.
    print("Building Seq2Fun input table")

    #IS THAT HOW YOU DO IT?? 469, 373-378 Seq2Fun 121-124 
    extract_info_from_paths(SRR_paths, output_directory)

    print("Initiating De Novo analysis with Seq2Fun")

    sample_table_path = output_directory+ "/sample.txt" #Path to the sample table
    tfmi_path = output_directory + "/" + DB + "/" + DB + "_v2.0.fmi" #Should look like "/Users/xaviersanterre/Test/Seq2Fun/database/birds/birds_v2.0.fmi" #Path to the bird database fmi file within the Birds folder
    gene_map_path = output_directory + "/" + DB + "/" + DB + "_annotation_v2.0.txt" #Should look like "/Users/xaviersanterre/Test/Seq2Fun/database/birds/birds_annotation_v2.0.txt" #Path to the birds annotation txt file within the Birds folder
    working_directory = output_directory #This needs to be within the databse or wathever folder that has the data base as well as the sample (Maybe this si wrong and could simply be wathever we want it to be)
    threads = "8" #Should pretty much always saty 8
    seq2fun_path_real = pathSeq2Fun + "/bin/seq2fun" #Path to the seq2fun tool"
    output_dir = Results_path #Path where you want the ouput file to be stored
    # References_path --> path where the abundance table will be for further use

    run_seq2fun(output_dir,References_path, seq2fun_path_real, sample_table_path, tfmi_path, gene_map_path, working_directory, threads)

    print("Finished annalysing your data")
    print("Feel free to look at the output files present in" + output_dir)
    profile.dict["F"]
    profile.update_profile()


####################################
# Differential Expression Analysis #
####################################
if profile.dict["STEP"] == "F":
    print("\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("\t~ F: Differential Expression Analysis ~")
    print("\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    from F_Differential_Expression import organize_DESeq2_genecounts, run_DESeq2_R

    print("\nDifferential Expression Analysis will now be performed using DESeq2.")
    print("For each variable, DESeq2 will compare pairs of conditions, using the gene counts from all")
    print("of the conditions to calculate the 'degree' of differential expression for each gene.")
    print("This degree of differential expression is represented by statistical numbers like log2-fold change and p-value.")
    print("For more information on how to interpret the DESeq2 results:")
    print("https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05b_wald_test_results.html#p-values")
    print("For more information on how DESeq2 works:")
    print("https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html")

    # make DESeq2 results folder
    deseq2_results_folder = os.path.join(Results_path, "DESeq2")
    if not os.path.exists(deseq2_results_folder): os.mkdir(deseq2_results_folder) # make the folder if it doesn't exist already

    # make compiled_counts folder
    compiled_counts_path = os.path.join(References_path, "compiled_counts")
    if not os.path.exists(compiled_counts_path): os.mkdir(compiled_counts_path)
        
    # make the sub folders under compiled_counts, one for each variable
    var_list = listdir_visible(Samples_path) # use Samples folder to know how many variables/names of the variable folders
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
    
    profile.dict["STEP"] = "G"
    profile.update_profile()


##########################
# Co-Expression Analysis #  !Unfinished!
##########################
        
if profile.dict["STEP"] == "G":

    ## 

    profile.dict["STEP"] = "END"
    profile.update_profile()

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~ You have reached the end of the pipeline! ~")
print("~     Thank you for using Taptomics :)      ~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
