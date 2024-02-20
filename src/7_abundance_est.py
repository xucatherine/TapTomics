
# HTSeq
    # takes BAM file from STAR

import HTSeq
import subprocess

def run_HTSeq_count(bam_file_path, gene_annotation_path, result_path):
    '''This function runs the htseq_count script from the HTSeq package.
    The script takes a BAM file of aligned RNA sequences (bam_file_path), 
    compares it to a GTF file of annotated genes (gene_annotation_path) 
    and counts the number of RNA sequences that matches with each gene ID.
    The ouput containing the gene count is written into result_path.'''

    # quick check that the inputs are correct
    if not bam_file_path.lower().endswith(".bam"):
        print("Error in run_HTSeq_count: Please input a .BAM file for the first input")
    if not gene_annotation_path.lower().endswith(".gtf"):
        print("Error in run_HTSeq_count: Please input a .GTF for the second input")
    

    # Ask the person to select stranded, unstranded or reverse stranded
    check = True # This boolean will be switched off once the person enters a valid input
    while check:
        strandedness = input(["Please input whether you're RNA-seq library is \r", 
                            "(y) forward/sense stranded \r",
                            "(n) unstranded or \r",
                            "(reverse) reverse stranded \r", 
                            "Type y, n, or reverse: "])
        # Forward stranded = the read represents the RNA sequence directly (complement of the template DNA strand)
        # Reverse stranded = the read represents the complement of the RNA sequence
        # Unstranded = there is no info on whether it is forward or reverse stranded

        strandedness = strandedness.lower() # convert person's answer to lowercase

        # Check that the input is valid
        if strandedness == "y" or strandedness == "n" or strandedness == "reverse":
            check = False
    
    # set up the command that will run the htseq-count script
    htseq_command = [
        "htseq-count",
        "-f", "bam", # Format of the input file (can be BAM, SAM, or CRAM)
        "-s",  strandedness, # Standedness. 'yes', 'no', or 'reverse'
        "-t", "exon", # feature type to count over
        "-i", "gene_id", # attribute to use as each feature's ID (in GTF files, it is gene_id)
        bam_file_path, # path to the aligned read file
        gene_annotation_path # path to the GTF annotation file
    ]

    # use subprocess to run the htseq-count script
    htseq_result = subprocess.run(htseq_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    # deal with errors
    if htseq_result.returncode == 0:
        # Command was successful, process the output
        with open(result_path, "w") as output_file:
            output_file.write(htseq_result.stdout)
    else:
        # There was an error
        print("Error in execution of htseq-count. Please check that the input paths are valid and that the file types are correct.")
        print(htseq_result.stderr)
    
    return
