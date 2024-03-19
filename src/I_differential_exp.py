
# Take gene counts from an abundance estimator, like HTSeq or feature count, perform differential expression analysis

import pandas as pd
import os
import subprocess
#import rpy2.robjects as robjects
#from rpy2.robjects import pandas2ri
#import rpy2.robjects.packages as rpackages

def organize_DESeq2_genecounts(var_folder_path, condition_labels=-1, results_folder="./"):
    '''This function compiles individual gene count files into one table 
    and creates a metadata table that records the condition associated with each sample/replicate.
    These two tables are required inputs for DESeq2, a function that compares the gene counts 
    between pairs of expression conditions (for 1 variable). When you input more than two conditions, 
    DESeq2 performs pairwise comparisons between each possible pair.

    var_folder_path is the path to a variable folder under Samples. 
    All of the conditions and corresponding samples under this variable will be used in the analysis
    Expected folder structure:
      variable
      | condition 1
        |   SRR00000
            |   fasta.fasta
                fasta.fastq
                aligned.bam
                counts.csv
            SRR000001
            |   ...
                
        condition 2
        |   ...
    
    condition_labels should give a list of the names for each condition, in the order that they appear in the folder system. 
    If no input, the default will be to just call them 'condition 1', 'condition 2', etc.
    '''
    ### If no condition labels were given, set up default condition labels ###
    if condition_labels == -1:
        condition_labels = []
        condition_num = len(os.listdir(var_folder_path))
        for i in range(1, condition_num+1):
            label = "condition " + str(i)
            condition_labels.append(label)
            
    ### Set up the data for DESeq2 ###
    print("Setting up dataframe for DESeq2!")
    # load info from the gene count files into pandas DataFrames, store all of these DataFrames into a single list
    gene_counts_DFs_list = []
    sample_labels = []
    conditions_tracker = [] # for each sample in gene_counts_DFs_matrix, this list will store the corresponding condition at the same index

    condition_folders = sorted(os.listdir(var_folder_path))
    for i in range(condition_folders):
        # iterate through each condition
        condition_path = os.path.join(var_folder_path, condition_folders[i])
        for SRR in sorted(os.listdir(condition_path)):
            # iterate through each SRR
            counts_path = os.path.join(condition_path, SRR, "counts.csv")
            dataframe = pd.read_csv(counts_path, index_col=0, sep='\t', header=None)
                # extract the gene-count file as a pandas dataframe
                # index_col=0 tells the function to use the first column (containing gene ID) as the labels/indexes
            gene_counts_DFs_list += [dataframe] # Append dataframe to list

            # use the SRR directory name as sample label
            sample_labels += [SRR]

        # condition labels for each sample that was just added
        number_of_samples = len(os.listdir(condition_path))
        conditions_tracker += [condition_labels[i]]*number_of_samples


    # Create a metadata Dataframe (will be used by DESeq2)
        # for each sample in the list "samples", you can look at the corresponding 
        # list position in the list "conditions_tracker" to see its condition
    metadata_DF = pd.DataFrame({
        'sample': sample_labels,
        'condition': conditions_tracker
    })
    print(metadata_DF)

    # concatenate Dataframes in gene_counts_DFs_list into a single Dataframe
    counts_DF = pd.concat(gene_counts_DFs_list, axis=1) # axis=1 should tell the function to line up the rows based on the row labels/indexes
    counts_DF.columns = sample_labels     # name columns after the samples
    
    print("Concatenated Pandas Dataframe")
    print(counts_DF)

    # Write the new gene-counts and metada dataframes to csv files
    gene_counts_path = results_folder + "comp_counts.csv"
    gene_counts_metadata = results_folder + "metadata.csv"
    counts_DF.to_csv(gene_counts_path)
    metadata_DF.to_csv(gene_counts_metadata)

    return


def run_DESeq2_R(gene_counts_path, metadata_path, result_path, counts_folder, normalize="FALSE", transform="FALSE", plots="FALSE"):
    
    ## I am unable to install rpy2 on my device. 
        # As a work-around, I wrote an R script for DESeq2, and this python function will write an 
        # extra line in this R script that calls my DESeq2 function with the appropriate inputs.
        # Once the R script is done running, the last line will be deleted
    print("function is running")
    my_script = open("R_scripts/run_deseq2.R", "a")
    my_script.write('run_DESeq2("' + gene_counts_path + '","' + metadata_path + '","' + result_path + '","' + counts_folder 
                     + '",normalize=' + normalize + ',transform=' + transform + ',plots=' + plots + ')')
    my_script.close() #close the file to free up memory
  
    # Run the R script
    subprocess.run(['Rscript', 'R_scripts/run_deseq2.R'])
    
    # Delete the last line of/function call from the R script
    my_script1 = open("R_scripts/run_deseq2.R", "r")
    lines = my_script1.readlines() # saves all of the lines in the script
    my_script1.close()
    my_script1 = open("R_scripts/run_deseq2.R", "w") 
    my_script1.write(''.join(lines[:-1])) # rewrite all of the lines except the last line
    my_script1.close()
    

    return

#### FUNCTIONS THAT ARE NOT FULLY IMPLEMENTED ###

def run_DESeq2(gene_counts_matrix, condition_labels, results_path="/", save_normalized_counts="no", norm_count_results="/"):
    '''This function compares the gene counts between pairs of expression conditions (for 1 variable). 
    When you input more than two conditions, DESeq2 performs pairwise comparisons between each possible pair.

    The gene count inputs will be taken as a matrix of file names/paths where each file contains the gene counts for a sample/replicate.
    The files should be organized as follows:
        [[replicate 1, replicate 2, replicate 3, ...] (condition 1)
        [replicate 1, replicate 2, replicate 3, ...]  (condition 2)
        [replicate 1, replicate 2, replicate 3, ...]] (condition 3)
    where each row is a different condition

    ##### BUILD TABLE, INPUT VARIABLE PATH

    condition_labels should give a label for each row of this matrix/condition. 
    Otherwise, the default will be to just number them 'condition 1', 'condition 2', etc.

    Normalize counts option:
        "yes": runs the DESeq2 script AND recalculates the normalized counts and saves them
        "only": only calculates the normalized counts, does not run the DESeq2 script
        "no": only runs the DESeq2 script
    '''
    
    ### R setup stuff ###
   
    ##install DESeq2 R package## (this took forever to run - ChatGPT recommended installing stuff on R separately - maybe i'll switch to writing an R script and then just  )
    print("installing DESeq2")
    r = robjects.r # initialize R instance
    # install bioconductor if not already installed
    r('''
      if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
      ''')
    # install DESeq2 using Bioconductor's BiocManager
    r('BiocManager::install("DESeq2")')

    # activate automatic conversion between pandas DataFrames and R dataframes
    print("Activating R")
    pandas2ri.activate()

    # import DeSeq2 through rpy2
    DeSeq2 = rpackages.importr('DESeq2')

    # import R's utility package
    base = rpackages.importr('base')

    # Import R's "graphics" package
    grDevices = rpackages.importr('grDevices')

    ### Set up the data for DESeq2###
    print("Setting up dataframe for DESeq2!")
    # load info from the gene count files into pandas DataFrames, store all of these DataFrames into a single list
    gene_counts_DFs_list = []
    conditions_tracker = [] # for each sample in gene_counts_DFs_matrix, this list will store the corresponding condition at the same index

    for i in range(len(gene_counts_matrix)): # iterate through each list in the gene_counts_matrix (where each list stores samples for a different condition)
        for file in gene_counts_matrix[i]:
            dataframe = pd.read_csv(file, index_col=0, sep='\t')
            # extract gene-count file as a pandas dataframe
            # index_col=0 tells the function to use the first column (containing gene ID) as the labels/indexes
            gene_counts_DFs_list += [dataframe] # Append dataframe to list

        number_of_samples = len(gene_counts_matrix[i])
        conditions_tracker += [condition_labels[i]]*number_of_samples

    print(gene_counts_DFs_list)
    
    # Create labels for the samples (sample 1, sample 2, ... )
    samples = []
    for i in range(1,len(gene_counts_DFs_list)+1):
        samples += ["sample "+str(i)]

    # Create a metadata Dataframe (will be used by DESeq2)
        # for each sample in the list "samples", you can look at the corresponding 
        # list position in the list "conditions_tracker" to see its condition
    metadata_DF = pd.DataFrame({
        'sample': samples,
        'condition': conditions_tracker
    })

    # concatenate Dataframes in gene_counts_DFs_list into a single Dataframe
    counts_DF = pd.concat(gene_counts_DFs_list, axis=1) # axis=1 should tell the function to line up the rows based on the row labels/indexes
    counts_DF.columns = samples     # name columns after the samples
    
    print("Concatenated Pandas Dataframe")
    print(counts_DF)

    ### R Stuff ###
    print("Converting the pandas dataframe into a DESeq2 dataframe")
    # convert the pandas DataFrames into R data frames
    counts_R = pandas2ri.py2rpy(counts_DF)
    metadata_R = pandas2ri.py2rpy(metadata_DF)

    #Create a DeSeq DataSet Object
    deseqDS = deseq2.DESeqDataSetFromMatrix(countData=counts_R,
                                            colData=metadata_R,
                                            design=base.Formula('~ condition'))
    
    if save_normalized_counts=="yes" or save_normalized_counts=="no":
        print("Running the DESeq2 pipeline!")
        # Run the DeSeq pipeline
        deseqDS = deseq2.DESeq(deseqDS)

        print("DESeq2 pipeline done running!")
        # Save the results in a variable
        res = deseq2.results(deseqDS)

        print("Converting the results to a pandas dataframe")
        # Convert the results back to a Pandas Dataframe
        results_DF = pandas2ri.rpy2rpy_dataframe(res)

        print("Writing the results to results file")
        # Write the results to a text file 
        with open(results_path, 'a') as f:
            DF_string = results_DF.to_string()
            f.write(DF_string)

    if save_normalized_counts=="yes" or save_normalized_counts=="only":
        # Calculate the normalized counts
        deseqDS = deseq2.estimateSizeFactors(deseqDS)
        normalized_counts = counts(deseqDS, normalized=TRUE)
        # convert to pandas dataframe
        results_counts = pandas2ri.rpy2rpy_dataframe(normalized_counts)
        # write to file
        with open(norm_count_results, 'a') as f:
            DF_string = results_counts.to_string()
            f.write(DF_string)

    # Return the pandas dataframes (?)
    return results_DF, results_counts


## edgeR ##
def run_edgeR():
    '''May be better for small sample sizes (low number of biological replicates)
    Consider running both edgeR and DeSeq2, or maybe choose one based on number of biological replicates?'''
    return
