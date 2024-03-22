
print("i am running :)")
from Bio import Entrez
    # 'Entrez' allows us to send requests to NCBI databases (from BioPython)
import os
    # we need this to edit folders
import subprocess

print("~~~ Welcome to our bioinformatic pipeline! ~~~\n")

def samples_setup():
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
    fasterq_dump_path = os.path.join(SRA_toolkit_path, "bin", "fasterq-dump")

    # On MacOS, the folder will automatically be flagged as coming from an unknown developper and this script will not be able to open the folder
    # check for the com.apple.quarantine flag and remove it if the folder has it
    xattr_list=subprocess.run(["xattr", SRA_toolkit_path], capture_output=True, text=True).stdout.split('\n')
    if "com.apple.quarantine" in xattr_list:
        subprocess.run(["xattr", "-r", "-d","com.apple.quarantine", SRA_toolkit_path])

    # Setting up Samples, References and Results folders (within the user's inputted path)
    global Samples_path # 'global' makes the variable accessible outside the function
    Samples_path = path+"/Samples"
    os.makedirs(Samples_path, exist_ok=True)
    global References_path # 'global' makes the variable accessible outside the function"
    References_path = path+"/References"
    os.makedirs(References_path, exist_ok=True)
    global Results_path # 'global' makes the variable accessible outside the function
    Results_path = path+"/Results"
    os.makedirs(Results_path, exist_ok=True)

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
        for cond in sorted(os.listdir(var_path)): # for condition folder in parent folder var_n
            print(f"\nList the SRR numbers belonging to variable {i+1}'s condition {cond.split('_')[-1]}, separated by spaces.")
            print("ex: SRR12345678 SRR91011109 SRR87654321")
            SRRs = input("SRR numbers: ") #[user inputs SRRs]
            SRRs_list = SRRs.split() # making list of SRRs
            
            # Starting the downloads
            print("We will now download the FastQ files for each SRR. Each SRR can take ~5 minutes.")
            if os.path.isdir(os.path.join(var_path, cond)):
                cond_path = os.path.join(var_path, cond)
                for SRR in SRRs_list: # loop through each SRR
                    SRR_path = os.path.join(cond_path, SRR) # new folder named as SRR number
                    os.makedirs(SRR_path, exist_ok=True)
                    print(f"Downloading {SRR}...")
                    # retrieve FASTQ file using SRR Toolkit
                    try:
                        subprocess.run([fasterq_dump_path, SRR, "--outdir", SRR_path, "--skip-technical", "--split-3"])
                            # this will create two, or three files:
                                # SRR000000_1.fastq ;  side 1 of paired reads
                                # SRR000000_2.fastq ;  side 2 of paired reads
                                # SRR000000_3.fastq ;  any unpaired/unmatched reads
                        # renaming the files to generic name
                        for file in os.listdir(SRR_path):
                            if file.startswith(SRR) and file.endswith(".fastq"):
                                num = file[-7] # extract the number
                                new_name = "raw_" + num + ".fastq"
                                os.rename(os.path.join(SRR_path, file), os.path.join(SRR_path, new_name))
                    except subprocess.CalledProcessError as e:
                        print(f"An error occurred while downloading {SRR}: {e}")  

                        
                    

    print("Thank you! Your transcriptomes are ready for processing.")

'''
    # Creating list of SRRs for easy processing - iterates through all folders in Samples
    global SRR_paths # make it so that SRR_paths can be accessed outside the function
    SRR_paths = []
    for var in os.listdir(Samples_path): # for each variable in the Samples folder
        var_path = os.path.join(Samples_path, var)
        for cond in os.listdir(var_path): # for each condition in the var_n folder
            cond_path = os.path.join(var_path, cond)
            for SRR in os.listdir(cond_path): # for each SRR in the cond_k folder
                SRR_path = os.path.join(cond_path, SRR)
                SRR_paths.append(SRR_path)

'''
samples_setup()



