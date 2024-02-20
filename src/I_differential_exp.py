
# Take gene counts from an abundance estimator, like HTSeq or feature count, perform differential expression analysis

import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pands2ri
from rpy2.robjects.packages import importr

def run_DESeq2(gene_counts_positive, gene_counts_negative, results_path):
    '''This function compares the gene counts between the positive and negative expression conditions.
    The inputs will be taken as lists, to take into account duplicates, triplicates, etc. of the positive and negative samples
    gene_counts_positive/negative: list of files, where each file contains the gene counts for a sample
    '''

    # activate automatic conversion between pandas DataFrames and R dataframes
    pandas2ri.activate()

    # import DeSeq2 through rpy2
    DeSeq2 = importr('DESeq2')

    # import R's utility package
    base = importr('base')

    # Import R's "graphics" package
    grDevices = importr('grDevices')

    # Create list of labels for the samples
        # for each sample in the array "samples", you can look at the corresponding 
        # array position in the array "conditions" to see if it is a positive or negative sample
    pos_length = gene_counts_positive.length()
    neg_length = gene_counts_negative.length()
    conditions = ["positive"]*pos_length + ["negative"]*neg_length
    samples = []
    for i in range(1,pos_length+1):
        samples += ["sample"+str(i)]
    for i in range(pos_length+1, neg_length+1):
        samples += ["sample"+str(i)]

    # Create a metadata Dataframe (will be used by DESeq2)
    metadata_DF = pd.Dataframe({
        'sample': samples,
        'condition': conditions
    })

    # load info from the gene count files as pandas DataFrames
        # create 2 lists containing Dataframes, where each Dataframe contains the info from 1 gene_count file
        # index_col=0 tells the function to use the first column (containing gene ID) as the labels (indexes)
    counts_pos_DFs = [pd.read_csv(file, index_col=0) for file in gene_counts_positive] # gene counts under positive phenotype expression
    counts_neg_DFs = [pd.read_csv(file, index_col=0) for file in gene_counts_negative] # gene counts under negative phenotype expression
    
    # compile DFs into a single list, and then concatenate into a single Dataframe
    counts_DFs = counts_pos_DFs + counts_neg_DFs
    counts_DF = pd.concat(counts_DFs, axis=1) # axis = 1 should line up the rows based on the row labels (indexes)
    counts_DF.columns = samples     # name columns after the samples

    ### R Stuff ###

    # convert the pandas DataFrames into R data frames
    counts_R = pandas2ri.py2rpy(counts_DF)
    metadata_R = pandas2ri.py2rpy(metadata_DF)

    #Create a DeSeqDataSet
    deseqDS = deseq2.DESeqDataSetFromMatrix(countData=counts_R,
                                            colData=metadata_R,
                                            design=base.Formula('~ condition'))
    
    # Run the DeSeq pipeline
    deseqDS = deseq2.DESeq(deseqDS)

    # Save the results in a variable
    res = deseq2.results(deseqDS)

    # Convert the results back to a Pandas Dataframe
    results_DF = panadas2ri.rpy2rpy_dataframe(res)

    # Write the results to a text file 
    with open(results_path, 'a') as f:
        DF_string = df.to_string()
        f.write(DF_string)

    # Return the pandas dataframe too
    return results_DF


def visualize_DESeq2(results_DF):
    '''This will be a function that takes the results from DESeq2 and generates plots 
    such as MA plots, volcano plots, and heatmaps'''
    
    return
