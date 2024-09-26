import pandas as pd


def extract_gene_name(gtfpath, id_type = "gene_id"):
    '''function to extract just the gene id (or transcript id) and gene name from gtf file'''

    gtfFile = open(gtfpath, "r")
    # read through each line
    dict = {}

    for i in range(4): gtfFile.readline() # skip first 4 lines

    for line in gtfFile:
        line_list = line.split("\t")
        attributes = line_list[-1] # long list of attributes always stored in last column
        attr_list = attributes.split("; ")

        has_gene_name = False # use to check if the gene has a gene name already
        for attr in attr_list: # iterate through attributes, look for attribute called "gene"
            if attr.split(" ")[0] == "gene":
                has_gene_name = True
                gene_name = attr.split(" ")[1].strip('"')
        
        if has_gene_name: 
            # if the gene has a gene name, find it's gene ID (or transcript ID). Save in dictionary, using gene ID as the key
            for attr in attr_list:
                    if attr.split(" ")[0] == id_type:
                        id_name = attr.split(" ")[1].strip('"')
            dict[id_name] = gene_name

    gtfFile.close()
    df = pd.DataFrame(list(dict.items()), columns=[id_type, 'gene_name'])
    df.set_index(id_type, inplace=True)
    return df

def add_gene_name_to_results(generesultfile, gtffile, outputfile, id_type):
    '''function that adds a 'gene_name' column to any gene results files (e.g. RSEM quantification file)
        by running the 'extract_gene_name' on the gtf file to get a table of just gene_id (or transcript_id) and gene_name
        then concatenating the RSEM table with this gene_id|gene_name table, using gene_id to guide the concatenation
    '''
    # upload gene results table
    rsemDF = pd.read_csv(generesultfile, sep='\t', header=0, index_col=0)

    # get just the gene ID (or transcript ID) + gene name from gtf file
    genenameDF = extract_gene_name(gtfpath=gtffile, id_type=id_type)

    ## combine rsemDF and genenameDF using the gene_ID index as the 'guide'
    newDF = pd.concat([rsemDF, genenameDF], axis=1)

    # move the gene_name column to the front (pop it out then reinsert it)
    column_to_move = newDF.pop("gene_name")
    newDF.insert(0, "gene_name", column_to_move)

    # write to tsv file
    newDF.to_csv(outputfile, index_label=id_type, sep="\t")


## RSEM results file
file_root = "SRR27321763.isoforms"
id_type = "transcript_id"
rsempath = "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/References/rsem-output/" + file_root + ".results"

## GTF file
gtffile = "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/References/GCF_000146045.2/genomic.gtf"

destinationPath = "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/References/rsem-output/" + file_root + ".new"

add_gene_name_to_results(rsempath, gtffile, destinationPath, id_type)