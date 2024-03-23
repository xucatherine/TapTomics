#Reference based maping module
#inputs:trimmed and cleaned RNA transcript, assembled or annotated genome
from Bio import Entrez, SeqIO
import subprocess 
import git
from bs4 import BeautifulSoup
import os
import requests

class referencemap:

    def __init__(self):
        pass

    def star(self, STARpath, genomepath, genomefasta, genomegtf, fastqpaths, resultpath, CPUcores):
        #assuming STAR is already installed and compied

        '''
        This downloading part may not be necessary
        To download STAR use bash lines
        conda create -n star_env star -c bioconda
        conda activate star_env
        '''
        '''
        The basic STAR workflow consists of two steps
        1. Generating genome indexes from the reference genome sequences (FASTA) and annotations (GTF)
        2. Mapping reads to the genome
        '''
        # Define the STAR commands

        #first this will create a genome index
        #then it will map reads onto the generated index
        star_command = [
            #creating genome index
            STARpath,  # the full path to STAR if it's not in your PATH
            '--runThreadN', CPUcores,  # Number of threads - number of available cores on the server node
            '--runMode', 'genomeGenerate', #directs STAR to run
            '--genomeDir', genomepath,  # Path to the genome directory
            '--genomeFastaFiles', genomefasta, 
            '--sjdbGTFfile', genomegtf, #path to annotations GTF file
            '--outFileNamePrefix', os.path.join(genomepath,'/'),
            #now running a mapping job
            '--runThreadN', CPUcores,
            '--genomeDir', genomepath,
            '--readFilesIn', fastqpaths[0], fastqpaths[1], #note, fastqpaths is a list including the name and path of read1 and read2 for paired end 
            '--outFileNamePrefix', os.path.join(resultpath,'/'),
            '--outSAMtype', 'BAM Unsorted', #output in BAM format
            #now quantify reads
            #this will be saved as 'ReadsPerGene.out.tab'
            '--quantmode', 'GeneCounts',  
        ]
        #support for gzipped fastq files
        if any(fastqpath.endswith('.gz') for fastqpath in fastqpaths):
                star_command += ['--readFilesCommand', 'zcat']

        # Execute the STAR command
        try:
            # Run the command and wait for it to complete
            subprocess.run(star_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print("STAR alignment and quantification completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running STAR: {e.stderr}")
        return 
    

#in MAIN the execution will be something like this

#initializing the variables we will need
#change this part to test on different devices
downloadpath = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data'
genomefasta = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data\GCA_933822405.4_Cyc-CA1.15-REFERENCE-ANNOTATED_genomic.fna'
genomegtf = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data\genomic.gtf'
genomeurl = 'https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_933822405.4/'
STARpath = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\STAR'
#This will be the SRR folder
resultpath = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data'
#depending on how many cores available or requested on compute canada
CPUcores = input('How many cores would you like STAR to utilize: ')
fastaqfiles = [r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data\raw_1.fastq', r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data\raw_2.fastq']

#the execution should be in a for loop so that an assembly is made for each sample

#initializing object of the class
g = referencemap()
g.star(STARpath, downloadpath,genomefasta, genomegtf, fastaqfiles,resultpath,CPUcores)

