
run_DESeq2 <- function(gene_counts_path, metadata_path, result_folder, normalize=FALSE, transform=FALSE, plots=FALSE) {
    print("R script is running")
    ### Install packages and load libraries ###
    
    # install bioconductor if not already installed
    if(!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager", repos='https://cloud.r-project.org/')
    }
    # install DESeq2 using Bioconductor's BiocManager
    if(!requireNamespace("DESeq2", quietly=TRUE)){
        BiocManager::install("DESeq2")
    }

    suppressPackageStartupMessages({
        require("DESeq2")
        require(ggplot2)    
    })


    ### "Download" the gene count data ###
    genecount_data = read.csv(gene_counts_path, header = TRUE, sep=",")
    metadata = read.csv(metadata_path, header = TRUE, sep=",",)

    ### Create DESeqDataSet Object
    dds <- DESeqDataSetFromMatrix(countData=genecount_data, colData=metadata, design=~condition, tidy=TRUE)
        # Design specifies how the counts from each gene depend on the variables in the metadata.
            # -> we care mainly about the condition associated with each gene.
            #tidy=TRUE argument tells DESeq2 to use first column of countData as row names
    

    ### Run DESeq2 ###
    dds <- DESeq(dds)
    res <- results(dds)
        ## Note:by default, the result function filters for p-adj<0.1
    
    summary(res) # prints summary of DESeq2 results

    tidy_res = results(dds, tidy=TRUE) #create a 'tidy' form of results for writing to csv
    tidy_res <- tidy_res[order(tidy_res$padj),] # sort the results by adjusted p-value (lower p-value = )
    write.csv(tidy_res, paste(result_folder, "/results.csv",sep=""), row.names=FALSE) # write this table to a csv file

    ### Normalize and Transform for Correlation analysis (optional) ###
    if(normalize){
        norm_counts = counts(dds, normalized=TRUE)
            # This DESeq2 function uses the counts and size factors from the DESeqDataSet object to normalize the counts
        write.csv(norm_counts, paste(result_folder, "/norm_counts.csv",sep=""))
    }
    if(transform){
        # Variance Stabilizing Transformation
        vst_data <- varianceStabilizingTransformation(dds, blind=FALSE)
        # Extract the matrix of transformed values
        vst_matrix <- assay(vst_data)
        write.csv(vst_matrix, paste(result_folder, "/vst_counts.csv",sep=""))
    }

    ### Make plots of the data ###
    if(plots){
        ## Volcano plot ##
        #reset par
        par(mfrow=c(1,1))
        # Make a basic volcano plot
        with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

        # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
        with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
        with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

        ## PCA plot ##
        plotPCA(vsdata, intgroup="condition")
    }

}
print("R script is running 2")

#wd_path <- getwd()
#print(wd_path)
#genecounts_path <- paste(wd_path,"/DESeq2_files/compiled_counts_for_deseq.csv",sep="")
#metadata_path <- paste(wd_path,"/DESeq2_files/metadata_for_deseq.csv",sep="")
#result_dir <- paste(wd_path,"/DESeq2-results",sep="")

#run_DESeq2(genecounts_path, metadata_path, result_dir, normalize=FALSE)
