
'''Seq2FUN is an ultrafast, all-in-one functional profiling tool for RNA-seq data analysis for organisms without reference genomes.'''
#The firsts thing to do is to generate a tabular with 4 colums, 1-name of the read, 2- forward read fatsq file, 3- reverse read fastq file (optional), 4- the condition
#Depending on Karina's code, asjust this section of the code
   #We need to add a line of code for if the data is paired end or not
     # See on https://www.seq2fun.ca/case_study.xhtml for what it should look like

import pandas as pd

def extract_info_from_paths(file_paths, directory_path):
    # Initialize an empty DataFrame with specified columns
    df = pd.DataFrame(columns=['Column1', 'Column2', 'Column3'])  # Adjust based on actual requirements

    # Iterate and add rows to the DataFrame
    for file_path in file_paths:
        # Considering compatibility with Windows, use os.path.split
        parts = file_path.split(os.path.sep)
        
        # Extract the specific elements based on their position in the path
        srr = parts[-2]  # SRR is the second last element
        file_name = parts[-1]  # File name is the last element
        p_cond = parts[-3]  # p_cond is the third last element
        
        # Create a dictionary for the row
        data_row = {
            'Column1': srr,
            'Column2': file_name,
            'Column3': p_cond,
        }
        
        # Append the new row to the DataFrame
        df = df.append(data_row, ignore_index=True)
     
    # Generate a unique file name to avoid overwriting existing files
    base_filename = "sample"
    filename = f"{base_filename}.txt"
    output_path = os.path.join(directory_path, filename)
    counter = 1
    while os.path.exists(output_path):
        filename = f"{base_filename}_{counter}.txt"
        output_path = os.path.join(directory_path, filename)
        counter += 1
    new_path = os.path.join(directory_path, 'sample.txt')
    with open(new_path, 'r') as file:
        lines = file.readlines()[1:]  # Read all lines except the first one

    with open(new_path, 'w') as file:
        file.writelines(lines)
    
    print(f"File saved as {output_path}")
     

# Example usage
file_paths = ['path/to/file1', 'path/to/file2']  # Update with actual paths
output_path = '/Users/xaviersanterre/Test/Seq2Fun/database'  # Specify your output path
df = extract_info_from_paths(file_paths, output_path)

#The second step is to ask the user to decide which database he wants to select
strings = ['algae', 'alveolates', 'amoebozoa', 'amphibians', 'animals', 'apicomplexans', 'arthropods', 'ascomycetes', 'basidiomycetes', 'birds', 'cnidarians', 'crustaceans', 'dothideomycetes', 'eudicots', 'euglenozoa', 'eurotiomycetes', 'fishes', 'flatworms', 'fungi', 'insects', 'leotiomycetes', 'mammals', 'mollusks', 'monocots', 'nematodes', 'plants', 'protists', 'reptiles', 'saccharomycetes', 'stramenopiles', 'vertebrates']
#write a fucntion to make it easier to read in a table of four columns instead a long list
def print_strings_as_table(strings, num_columns=4):
    for i, string in enumerate(strings):
        print(f'{string:<20}', end='')  # Adjust 20 as needed for your string lengths
        if (i + 1) % num_columns == 0:
            print()  # Newline after every 4 items
    if len(strings) % num_columns != 0:
        print()  # Ensure ending on a newline if not divisible by num_columns

# Ask the user to pick from the list of databases
print ("please pick the database from the followign table that represents the better your sampled organism")
print_strings_as_table(strings)
print ("Now please got to https://www.expressanalyst.ca/ExpressAnalyst/docs/Databases.xhtml under the - Without a reference Transcriptome - and download the database")
databaseinput = input("Input the database file:") #I need to check how Karina did to allow inputs and put them in the appropriate folders.

DB = input("From the list above,please pick and write in lowercase the appropriate database for Seq2FUN tool to annalyse your reference free transcriptome (visit https://www.seq2fun.ca/database.xhtml for more information): ")
#Now the database can be downloaded
##The code is quite complex, so I was thinking maybe we could simply ask in the main file that the user inputs the database, just like he would input the genome if it was genome base.
# Now that the table has been generated and the database selected, the function can be run
import subprocess
import sys
import os
import shutil

def run_seq2fun(output_dir,path_seq2fun, sample_table, tfmi_path, gene_map, working_dir, threads=8):
    """
    Runs seq2fun with specified parameters.

    Args:
    - sample_table (str): Path to the sample table file.
    - tfmi_path (str): Path to the .fmi database file.
    - gene_map (str): Path to the gene map file.
    - output_dir (str): Directory where the output will be saved.
    - threads (int, optional): Number of threads to use. Defaults to 8.
    """

    # Change directory to the working directory/output directory for the code to work
    os.chdir(working_dir)

    # Constructing the seq2fun command
    command = [
        path_seq2fun, # Adjust the path according to where seq2fun is installed
        "--sampletable", sample_table,
        "--tfmi", tfmi_path,
        "--genemap", gene_map,
        "-w", str(threads),
        "--profiling",
        "-V",
        "--outputMappedCleanReads",
        "--outputReadsAnnoMap",
    ]

    try:
        # Print the command to be executed (for debugging)
        print("Executing command:", ' '.join(command))
        
        # Running the seq2fun command
        subprocess.run(command, check=True)
        
        #Now write a code to move the output file to the right destination
        candidate_output = os.path.join(working_directory, 'S2fid_abundance_table_all_samples.txt')
        new_path_output = os.path.join(working_directory, 'Seq2Fun_summary_all_samples.html')
        new_path_catherine = os.path.join(working_directory, 'S2fid_abundance_table_all_samples_submit_2_expressanalyst.txt')
        print(new_path_catherine)
       
        # Move the file
        shutil.copy(new_path_output, output_dir)
        
        #Now before moving the fiecond file, remove the second column of teh table so that it is as required for the next steps
        # Read the file
        with open(new_path_catherine, 'r') as file:
            lines = file.readlines()

        # Remove the second row (index 1 since Python lists are 0-indexed)
        if len(lines) > 1:  # Check if there's a second row to remove
            del lines[1]

        # Write the modified content back to the file
        with open(new_path_catherine, 'w') as file:
            file.writelines(lines)
        #Then move it to wherever Catherine wants it to be
        shutil.copy(new_path_catherine, output_dir) #This will need to be changed to wherever Catherine wants me to send this
        print(f"Seq2Fun completed successfully. Output saved to {output_dir}")
    except subprocess.CalledProcessError as e:
        print("Seq2Fun encountered an error:", e, file=sys.stderr)


sample_table_path = "/Users/xaviersanterre/Test/Seq2Fun/database/sample.txt" #Path to the sample table
tfmi_path = "/Users/xaviersanterre/Test/Seq2Fun/database/birds/birds_v2.0.fmi" #Path to the bird database fmi file within the Birds folder
gene_map_path = "/Users/xaviersanterre/Test/Seq2Fun/database/birds/birds_annotation_v2.0.txt" #Path to the birds annotation txt file within the Birds folder
working_directory = "/Users/xaviersanterre/Test/Seq2Fun/database" #This needs to be within the databse or wathever folder that has the data base as well as the sample (Maybe this si wrong and could simply be wathever we want it to be)
threads = "8" #Should pretty much always saty 8
seq2fun_path = "/Users/xaviersanterre/Test/Seq2Fun/bin/seq2fun" #Path to the seq2fun tool"
output_dir = "/Users/xaviersanterre/Test" #Path where you want the ouput file to be stored (maybe we could have both outputs at the same plac and then catherine fetches it from there?)

    
run_seq2fun(output_dir,seq2fun_path, sample_table_path, tfmi_path, gene_map_path, working_directory, threads)
