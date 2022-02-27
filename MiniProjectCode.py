# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 14:47:10 2022

@author: Anthony
"""

#CompBio Mini Project

#Prereq 1: SPAdes installation. Typed this into terminal:
    #wget http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Linux.tar.gz
    #tar -xzf SPAdes-3.15.4-Linux.tar.gz
    #cd SPAdes-3.15.4-Linux/bin/
    
    #NOTE on SPAdes installation above: recommended to add "SPAdes installation directory to the PATH variable." 
    #Not exactly sure what this means but if problems arise this could be why

#Prereq 2: Prokka installation. Not working for some reason. 

    
import os
from Bio import SeqIO

#import os will create the working directory for the entire script. 
#Creating new folder called "MiniProject_Script" and
#ensuring that this is my current directory
os.mkdir("MiniProject_Script")
os.chdir("MiniProject_Script")

#creating the requested output file and files listed below in a "Results" folder. Do this via an os.sytem call
Log_Results = open("miniproject.log", "w")


#1. Retrieve the Illumina reads for the resequencing of K-12 project: https://www.ncbi.nlm.nih.gov/sra/SRX5005282
    #These are single-end Illumina reads. Doing this via the backend of NCBI

os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR818/SRR8185310/SRR8185310.sra")

Illum_Reads_File = "SRR8185310.sra" #Reads stored here. Single-end in this case
#Need to uncompress the .sra file above via fastq-dump. .sra files = huge.
#fastq-dump also converts .sra file into a fastq file

os.system("fastq-dump " + Illum_Reads_File)

#create new variable to store the newly formed .fastq file after fastq-dump does its job
fastq_File = "SRR8185310.fastq"

#2. Using SPAdes, assemble the genome. Write the SPAdes command to the log file.
#Code that will call SPAdes to do the assembly

SPAdes_Assembly = "spades -k 55,77,99,127 -t 2 --only-assembler -s " + fastq_File + " -o SRR8185310_Assembly/"
Log_Results.write(SPAdes_Assembly + " \n")
os.system(SPAdes_Assembly)
## This is the handle for the assembled contig file
contig_handle = "SRR8185310_Assembly/contigs.fasta"
contigs = SeqIO.parse(contig_handle, "fasta")