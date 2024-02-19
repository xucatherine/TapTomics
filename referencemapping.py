#Reference based maping module
#inputs:trimmed and cleaned RNA transcript, assembled or annotated genome
from Bio import Entrez, SeqIO
import subprocess 
import git

class referencemap:
    def __init__(self):
        
    def fetchgenome(self):
        #assuming biopython has already been installed in a previous script
        #in the main file they will need the line 'from Bio import 
        Entrez.email = input('Please enter your email to access the NCBI database: ')
        #fetching genome data in GenBank format
        query = input('Please enter the accession number of interest: ')
        handle = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='text')
        genome = SeqIO.read(handle, 'gb')
        handle.close()
        #saving genome to a file
        SeqIO.write(genome, f"{genome.id}.gb", 'gb')
        return genome

    def star(self)
        #download and install star from github repository
        #using official STAR 2.7.11b maintained by alexdobin
        repo = git.Repo.clone_from('git@github.com:alexdobin/star', '/rffiles', branch='master')
        star_command = []
        
    



 
#Step 1: Retrieve assembled or annotated genome
#search

#step 2: Map RNA transcripts onto the genome

#Step 3: Display maps or data that the tool generates if the user wants to see them

#Output: Reference-based assembly
