#RSEM command steps:
# install RSEM with conda, from Bioconda
# have STAR installed and accessible in PATH

## Prepare reference indices for the genomes, extract transcripts from GTF file: ##
#rsem-prepare-reference --gtf .../zimo-yeast/References/GCF_000146045.2/genomic.gtf
    #--star 
    #.../zimo-yeast/References/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna 
    #.../zimo-yeast/References/rsem/yeast

#Run on each SRR:
	#If want to use own aligner (BAM file has to be alignement to transcriptome):
#rsem-calculate-expression [options] --alignments --paired-end input reference_name sample_name

#	normal pipeline, paired-end:
#rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name 

#rsem-calculate-expression 
    #--star 
    #--paired-end 
    #.../zimo-yeast/Samples/var_1/cond_1/SRR27321765/rawF.fastq, [other fastq files for the same sample]
    #.../zimo-yeast/Samples/var_1/cond_1/SRR27321765/rawR.fastq, [other fastq files for the same sample]
    #.../Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/zimo-yeast/References/rsem/yeast 
    #.../zimo-yeast/References/rsem-output/SRR27321765
import subprocess
import shutil

class RSEMobject:
    def __init__(self, genomegtf, genomefasta):
        self.genomegtf = genomegtf
        self.genomefasta = genomefasta
        self.refdir = None
        return

    def prep_RSEM_references(self, outputprefix):
        self.refdir = outputprefix
        rsem_reference_command = [
            'rsem-prepare-reference',
            '--gtf', self.genomegtf,
            '--star',
            self.genomefasta,
            outputprefix
        ]
        subprocess(rsem_reference_command)
        return
    
    def RSEM_calc_expression(self, F_fastqs, R_fastqs, outputprefix, pairedend=True):
        # combines the fastq lists into a single string, with the fastq paths separated by commas
            #this is how the command expects the files)
        F_fastqs = ",".join(F_fastqs)
        R_fastqs = ",".join(R_fastqs)
        rsem_calc_expression = [
            'rsem-calculate-expression',
            '--star',
            '--paired-end',
            F_fastqs,
            R_fastqs,
            self.refdir,
            outputprefix
        ]
        return
    
    def organize_RSEM_folders(rsem_result_folder, prefix, SRRfolder):
        ''' RSEM output files that will be kept in the References folder 
                or renamed/moved to the Samples folder:
            
            keep:
            SRR00001.log
            SRR00001.stat

            move to SRRfolder:
            SRR00001.genes.results -> .../SRR00001/transcriptcounts.tab
            SRR00001.isoforms.results -> .../SRR00001/isoformcounts.tab
            SRR00001.transcript.bam -> .../SRR00001/transcript_aligned.bam

            Inputs:    
                rsem_result_folder: .../References/RSEM_log
                SRRfolder: .../Samples/var_x/cond_x/SRRXXXXXX
        '''

        move_names = ["genes.results", "isoforms.results", "transcript.bam"]
        new_names = ["transcriptcounts.tab", "isoformcounts.tab", "transcript_aligned.bam"]
        for i in range(len(move_names)):
            start_path = os.path.join(rsem_result_folder, prefix + "." + move_names[i])
            end_path = os.path.join(SRRfolder, new_names[i])
            shutil.move(start_path, end_path)
        return

    
from A_minis import listdir_visible
import os
# WHAT WOULD GO IN THE MAIN FILE:
def run_RSEM_on_var(var_folder): 
    genomegtf = input("Please input the path to your genome's gtf file: ")
    genomefasta = input("Please input the path to your genome's fasta file: ")
    genomename = os.path.basename(genomefasta).rsplit(".", 1)[0]
        # gets the basename of the genomefasta path, then removes the extension
    rsem = RSEMobject(genomegtf, genomefasta)
    References_path = ".../References"

    rsem_refdir = os.path.join(References_path, "RSEM_refs", genomename)
    rsem.prep_RSEM_references(rsem_refdir)
    for cond_folder in listdir_visible(var_folder):
        cond_path = os.path.join(var_folder, cond_folder)
        for SRR in listdir_visible(cond_folder):
            # raw fastq file paths
            F_fastq = os.path.join(cond_path, SRR, "rawF.fastq")
            R_fastq = os.path.join(cond_path, SRR, "rawR.fastq")
            # RSEM output path
            output_path = os.path.join(References_path, "RSEM_log", SRR)

            # run the read mapping + transcript/expression quantification 
            rsem.RSEM_calc_expression([F_fastq], [R_fastq], output_path)

            # move files out of the RSEM output path and into the Samples/var_x/cond_x/SRR00000 folders
            SRR_folder = os.path.join(cond_path, SRR)
            rsem.organize_RSEM_folders(output_path, SRR, SRR_folder)
    return
    