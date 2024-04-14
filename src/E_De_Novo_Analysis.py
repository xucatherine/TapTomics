'''Seq2FUN is an ultrafast, all-in-one functional profiling tool for RNA-seq data analysis for organisms without reference genomes.'''
#The firsts thing to do is to generate a tabular with 4 colums, 1-name of the read, 2- forward read fatsq file, 3- reverse read fastq file (optional), 4- the condition
#Depending on Karina's code, asjust this section of the code
   #We need to add a line of code for if the data is paired end or not
     # See on https://www.seq2fun.ca/case_study.xhtml for what it should look like
import shutil
import pandas as pd
import os
import subprocess
import sys
def move_and_rename_files(file_paths, output_directory):
    # Iterate over each file path
    for file_path in file_paths:
        # Extract the directory and filename from the file path
        directory, filename = os.path.split(file_path)
        
        # Generate new file paths for trimmed files
        trimmed_F_path = os.path.join(directory, filename + "/trimmed_F.fastq.gz")
    
        # Get split the path
        parts = file_path.split(os.path.sep)
        # Extract the specific elements based on their position in the path
        file_name_for = parts[-1] + "trimmed_F.fastq.gz" # File name
        
        # Determine the destination path for the file
        destination_path = os.path.join(output_directory, file_name_for)
        
        # Copy the file to the new directory and rename it
        shutil.copy(trimmed_F_path, destination_path)
        
        #Now same thing for the reverse
        # Generate new file paths for trimmed files
        trimmed_R_path = os.path.join(directory, filename + "/trimmed_R.fastq.gz")
        # Extract the specific elements based on their position in the path
        file_name_rev = parts[-1] + "trimmed_R.fastq.gz" # File name
        # Determine the destination path for the file
        destination_path = os.path.join(output_directory, file_name_rev)
        # Copy the file to the new directory and rename it
        shutil.copy(trimmed_R_path, destination_path)
#write a fucntion to make it easier to read in a table of four columns instead a long list
def print_strings_as_table(strings, num_columns=4):
    for i, string in enumerate(strings):
        print(f'{string:<20}', end='')  # Adjust 20 as needed for your string lengths
        if (i + 1) % num_columns == 0:
            print()  # Newline after every 4 items
    if len(strings) % num_columns != 0:
        print()  # Ensure ending on a newline if not divisible by num_columns


def extract_info_from_paths(file_paths, output_direct_path):
    # Initialize an empty DataFrame with specified columns
    df = pd.DataFrame(columns=['Column1', 'Column2', 'Column3', 'Column4'])  # Adjust based on actual requirements

    # Iterate and add rows to the DataFrame
    for file_path in file_paths:
        # Considering compatibility with Windows, use os.path.split
        parts = file_path.split(os.path.sep)
        
        # Extract the specific elements based on their position in the path
        srr = parts[-1]  # SRR is the second last element
        file_name_for = parts[-1] + "trimmed_F.fastq.gz" # File name
        file_name_rev = parts[-1] + "trimmed_R.fastq.gz" # File name
        p_cond = parts[-3]  # p_cond is the third last element

        
        # Create a dictionary for the row
        data_row = {
            'Column1': srr,
            'Column2': file_name_for,
            'Column3': file_name_rev,
            'Column4': p_cond
        }
        
        # Append the new row to the DataFrame
        df = df.append(data_row, ignore_index=True)

# Remove the first row of the DataFrame
    df = df.iloc[0:]

    # Generate a unique file name to avoid overwriting existing files
    base_filename = "Sample"
    csv_filename = f"{base_filename}.csv"
    csv_output_path = os.path.join(output_direct_path, csv_filename)
    txt_output_path = os.path.join(output_direct_path, f"{base_filename}.txt")
    counter = 1
    while os.path.exists(csv_output_path):
        csv_filename = f"{base_filename}_{counter}.csv"
        csv_output_path = os.path.join(output_direct_path, csv_filename)
        txt_output_path = os.path.join(output_direct_path, f"{base_filename}_{counter}.txt")
        counter += 1
    
    # Write the DataFrame to a CSV file
    df.to_csv(csv_output_path, index=False)  # Write DataFrame to CSV without row indices
    
    # Read the CSV file, replace commas with spaces, then save it as a text file
    with open(csv_output_path, 'r') as csv_file:
        lines = csv_file.read().replace(',', '\t').splitlines()[1:]  # Replace commas with spaces and remove the header line

    with open(txt_output_path, 'w') as txt_file:
        txt_file.write('\n'.join(lines))

    print(f"File saved as {txt_output_path}")
    
# Example usage
#file_paths = ['path_that_user_decides/Samples/var_1/cond_1/SRR96', 'path_that_user_decides/Samples/var_2/cond_3/SRR97']  # Update with actual paths
#output_direct_path = '/Users/xaviersanterre/Test/Seq2Fun/database'  # Specify your output path
extract_info_from_paths(file_paths, output_direct_path)

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



def run_seq2fun(output_dir,References_path, path_seq2fun, sample_table, tfmi_path, gene_map, working_dir, threads=8):
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
       
        # Move the file
        shutil.copy(new_path_output, output_dir)
        
        #Now before moving the fiecond file, remove the second column of the table so that it is as required for the next steps
        # Read the file
        with open(candidate_output, 'r') as file:
            lines = file.readlines()

        # Remove the second row (index 1 since Python lists are 0-indexed)
        if len(lines) > 1:  # Check if there's a second row to remove
            del lines[1]

        # Write the modified content back to the new file that will be sent to catherine
        with open(new_path_catherine, 'w') as file:
            file.writelines(lines)
        #Then move 
        shutil.copy(new_path_catherine, References_path)
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
