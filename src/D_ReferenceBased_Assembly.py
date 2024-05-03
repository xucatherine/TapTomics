#Reference based maping module
#inputs:trimmed and cleaned RNA transcript, assembled or annotated genome

import subprocess 
import os
import shutil

class referencemap:

    def __init__(self, STARpath, genomepath, genomefasta, genomegtf, CPUcores):
        self.STARpath = STARpath
        self.genomepath = genomepath
        self.genomefasta = genomefasta
        self.genomegtf = genomegtf
        self.CPUcores = CPUcores
        return

    def star_index(self):
        #first this will create a genome index
        #then it will map reads onto the generated index
        star_command_index = [
            #creating genome index
            self.STARpath,  # the full path to STAR if it's not in your PATH
            '--runThreadN', self.CPUcores,  # Number of threads - number of available cores on the server node
            '--runMode', 'genomeGenerate', #directs STAR to run
            '--genomeDir', self.genomepath,  # Path to the genome directory
            '--genomeFastaFiles', self.genomefasta, 
            '--sjdbGTFfile', self.genomegtf, #path to annotations GTF file
            '--outFileNamePrefix', self.genomepath,
        ]

        try:
            # Run the genome index creation command
            subprocess.run(star_command_index, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running STAR index: {e.stderr}")

        return
    
    def star_map(self, fastqpaths, resultpath, resultprefix, transcriptCoords=False):

        # setup the resultspath properly (assuming the person inputted a path to a folder)
        resultpath = os.path.join(resultpath, resultprefix)
        resultpath = resultpath + "_"
            # this will be the input for --outoutFileNamePrefix
            # STAR will use this as a prefix for all of its file outputs
        
        star_command_mapping = [
            #running a mapping job
            self.STARpath,  # the full path to STAR if it's not in your PATH
            '--runThreadN', self.CPUcores,
            '--genomeDir', self.genomepath,
            '--readFilesIn', fastqpaths[0], fastqpaths[1], #note, fastqpaths is a list including the name and path of read1 and read2 for paired end 
            '--outFileNamePrefix', resultpath,
            '--outSAMtype', 'BAM', 'Unsorted', #output in BAM format
            #now quantify reads
            #this will be saved as 'ReadsPerGene.out.tab'
            '--quantMode', 'GeneCounts',  
        ]
        if transcriptCoords: #if set to True
             star_command_mapping.append('TranscriptomeSAM') 
                # to output Aligned.toTranscriptome.out.bam, which can be used directly with RSEM to quantify transcripts
        
        # Execute the STAR command
        try:
            #Run mapping job
            subprocess.run(star_command_mapping, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print("STAR alignment and quantification completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running STAR: {e.stderr}")
        
        return


    def star(self, STARpath, genomepath, genomefasta, genomegtf, fastqpaths, resultpath, resultprefix, CPUcores, transcriptCoords=False):
        #assuming STAR is already installed and compiled (installed with conda)

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
        # setup the resultspath properly (assuming the person inputted a path to a folder)
        resultpath = os.path.join(resultpath, resultprefix)
        resultpath = resultpath + "_"
            # this will be the input for --outoutFileNamePrefix
            # STAR will use this as a prefix for all of its file outputs

        # Define the STAR commands

        #first this will create a genome index
        #then it will map reads onto the generated index
        star_command_index = [
            #creating genome index
            STARpath,  # the full path to STAR if it's not in your PATH
            '--runThreadN', CPUcores,  # Number of threads - number of available cores on the server node
            '--runMode', 'genomeGenerate', #directs STAR to run
            '--genomeDir', genomepath,  # Path to the genome directory
            '--genomeFastaFiles', genomefasta, 
            '--sjdbGTFfile', genomegtf, #path to annotations GTF file
            '--outFileNamePrefix', genomepath,
        ]

        star_command_mapping = [
            #running a mapping job
            STARpath,  # the full path to STAR if it's not in your PATH
            '--runThreadN', CPUcores,
            '--genomeDir', genomepath,
            '--readFilesIn', fastqpaths[0], fastqpaths[1], #note, fastqpaths is a list including the name and path of read1 and read2 for paired end 
            '--outFileNamePrefix', resultpath,
            '--outSAMtype', 'BAM', 'Unsorted', #output in BAM format
            #now quantify reads
            #this will be saved as 'ReadsPerGene.out.tab'
            '--quantMode', 'GeneCounts',  
        ]
        if transcriptCoords: #if set to True
             star_command_mapping.append('TranscriptomeSAM') 
                # to output Aligned.toTranscriptome.out.bam, which can be used directly with RSEM to quantify transcripts

        #support for gzipped fastq files
        if any(fastqpath.endswith('.gz') for fastqpath in fastqpaths):
                star_command_mapping += ['--readFilesCommand', 'zcat']

        # Execute the STAR command
        try:
            # Run the genome index creation command
            subprocess.run(star_command_index, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            #Run mapping job
            subprocess.run(star_command_mapping, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print("STAR alignment and quantification completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running STAR: {e.stderr}")
        return

    def organize_star_files(self, star_result_folder, prefix, SRRfolder):
        '''
        STAR files that will be kept, renamed/moved, or deleted:
            SRR00001_Log.final.out - keep
            SRR00001_SJ.out.tab - keep? (records splice junctions)
            SRR00001_Aligned.toTranscriptome.out.bam (optional) - keep
            SRR00001_Log.out - delete
            SRR00001_Log.progress.out - delete
            SRR00001_Aligned.out.bam - move to Samples/var_x/cond_x/SRRXXXXX/aligned.bam
            SRR00001_ReadsPerGene.out.tab - move to Samples/var_x/cond_x/SRRXXXXX/genecounts.tab
        
        star_result_folder: .../References/STAR_log
        SRRfolder: .../Samples/var_x/cond_x/SRRXXXXXX
        '''

        delete_names = ["Log.out", "Log.progress.out"]
        for i in range(len(delete_names)):
            path = os.path.join(star_result_folder, prefix + "_" + delete_names[i])
            if os.path.exists(path):
                 os.remove(path)

        move_names = ["Aligned.out.bam", "ReadsPerGene.out.tab"]
        new_names = ["aligned.bam", "genecounts.tab"]
        for i in range(len(move_names)):
            start_path = os.path.join(star_result_folder, prefix + "_" + move_names[i])
            end_path = os.path.join(SRRfolder, new_names[i])
            shutil.move(start_path, end_path)
        
        return
    

#in MAIN the execution will be something like this

#initializing the variables we will need
#change this part to test on different devices

#downloadpath = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data'
#genomefasta = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data\GCA_933822405.4_Cyc-CA1.15-REFERENCE-ANNOTATED_genomic.fna'
#genomegtf = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data\genomic.gtf'
#STARpath = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\STAR'
#This will be the SRR folder
#resultpath = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data'
#depending on how many cores available or requested on compute canada
#CPUcores = str(input('How many cores would you like STAR to utilize: '))
#fastaqfiles = [r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data\raw_1.fastq', r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data\raw_2.fastq']

genomefasta = "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/References/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna"
genomegtf = "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/References/GCF_000146045.2/genomic.gtf"
#where the genome indices will be stored:
genomedir = "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/References/GCF_000146045.2/indices"
os.makedirs(genomedir, exist_ok=True)

STARpath = "star"
#depending on how many cores available or requested on compute canada
CPUcores = str(input('How many cores would you like STAR to utilize: '))
fastqfiles_matrix = [["/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/Samples/var_1/cond_1/SRR27321763/rawF.fastq",
              "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/Samples/var_1/cond_1/SRR27321763/rawR.fastq"],
              ["/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/Samples/var_1/cond_1/SRR27321765/rawF.fastq",
               "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/Samples/var_1/cond_1/SRR27321765/rawR.fastq"]]
resultprefixes= ["SRR27321763", "SRR27321765"]
SRR_folders = ["/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/Samples/var_1/cond_1/SRR27321763",
              "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/Samples/var_1/cond_1/SRR27321765"]

#folder where all of the outputs will go
starLog_dir = "/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/References/STAR_log"
os.makedirs(starLog_dir, exist_ok=True)

#the execution should be in a for loop so that an assembly is made for each sample
#initializing object of the class
g = referencemap(STARpath, genomedir, genomefasta, genomegtf, CPUcores)
g.star_index() #create the genome indices for STAR

for i in range(len(resultprefixes)):
    g.star_map(fastqfiles_matrix[i], starLog_dir, resultprefixes[i])
    m = input("I am done running STAR, would you like to continue? (type any key)")
    g.organize_star_files(starLog_dir, resultprefixes[i], SRR_folders[i])