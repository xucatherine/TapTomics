
run_tximport <- function(var_path, metadata_csv_path) {
    if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos='https://cloud.r-project.org/')
    }

    if (!requireNamespace("tximport", quietly = TRUE)) {
    BiocManager::install("tximport")
    }

    library(tximport)

    dir <- "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/References"
        # overarching directory containing the data, for ease of use
    list.files(dir)

    # create a vector pointing to the quantification files
        # use metadata table to get vector of the SRR's
    metadata <- read.table(file.path(dir, "compiled_counts", "metadata.csv"), header = TRUE, sep=",")
    print(metadata)

        # concatenate SRR's with the rest of the path/".genes.results" file from RSEM
    files <- file.path(dir, "rsem-output", paste0(metadata$sample, ".genes.results"))
    print(files)
    names(files) <- metadata$sample
    txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE) # preps table for DESeq2
    head(txi.rsem$counts) 
}

### NEXT IS DESEQ2 STUFF 
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#DESeq2
# dds <- DESeqDataSetFromTximport(countData=txi, colData=metadata, design=~condition)