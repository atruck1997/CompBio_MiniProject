# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 16:12:08 2022

@author: Anthony
"""

import os
from Bio import SeqIO

os.path.expanduser('~') #Puts you into home directory
user_path = os.getcwd() #displays the current working directory of machine. should be the home direc. or wherever python is installed
print(user_path)
directory = 'results' #This will create the desired results folder where results of pipeline will show up
path = os.path.join(user_path, directory) #will join the user_path and results

#cannot make a directory if it already exists.
#This will create directory ONLY if it doesnt already exist
if not os.path.isdir(path):
 os.mkdir(path)

#variable storing requested miniproject.log file  
#log_file = open(path + "/miniproject.log", 'w')

#1. Retrieve the Illumina reads for the resequencing of K-12 project: https://www.ncbi.nlm.nih.gov/sra/SRX5005282
 #These are single-end Illumina reads. Doing this via the backend of NCBI

os.system("prefetch SRR8185310")

#Need to uncompress the .sra file above via fastq-dump. .sra files = huge.
#fastq-dump converts .sra file of reads into a fastq file for SPAdes

os.system('fastq-dump -I --split-files SRR8185310.sra')

#2. Using SPAdes, assemble the genome. Write the SPAdes command to the log file.
#Code that will call SPAdes to do the 
#-t is the number of threads used 
#-o is the output directory. /results/SRR8185310_Assembly
#This will run SPAdes

command = 'python3 SPAdes-3.15.4-Darwin/bin/spades.py -k 55, 77, 99, 127 -t 2 --only-assembler -s SRR8185310_1.fastq -o ' + path +'/SRR8185310_Assembly/'

#Writing command to MiniProject.log
#first output that will be written to it
log_file = open(path + "/MiniProject.log","w+") #AFTER DUE DATE. But should be "/miniproject.log" here
for i in range(1):
 log_file.write("Spades command: " + command + "/n" + "/n")

#3. Code to keep contigs > 1000.
#path = '/Users/Anthony/Desktop/Comp/SRR8185310_Assembly/'

Contigs_Over_1000=[]

with open(path + "contigs.fasta", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) > 1000:
 # Add this record to our list
         Contigs_Over_1000.append(record)

long_assembly = ("There are %i contigs > 1000 in the assembly." % len(Contigs_Over_1000))
#print(long_assembly)

#RESULT of long_assembly - There are 150 contigs > 1000 in the assembly.

#Write these sequences to their own file 'long_sequences.fasta' 
#for use with Prokka. Input file to prokka - long_sequences.fasta
SeqIO.write(Contigs_Over_1000, "long_sequences.fasta", "fasta")

#IN THE SRR8185310_Assembly Folder!!.  
#Write to miniproject.log
log_file = open(path + "miniproject.log","a+")
for i in range(1):
 log_file.write(long_assembly + "\n" + "\n")

#4. Code for determining total assembly reads for contigs > 1000
length = 0
temp = 0
for record in SeqIO.parse('long_sequences.fasta', 'fasta'):
 temp = len(record)
 length += temp
final_length = str(length)
total_length = ('There are ' + str(length) + ' bp in the assembly.')
#print(total_length)

#RESULT of total assembly reads for contigs > 1000 -  
#There are 4535677 bp in the assembly

#Write total_length command to miniproject.log
#path = '/Users/Anthony/Desktop/Comp/SRR8185310_Assembly/'
log_file = open(path + "miniproject.log","a+") #AFTER DUE DATE. "/miniproject.log"
for i in range(1):
 log_file.write(total_length + "\n" + "\n")

#5. Using prokka to annotate Assembly.This command will call prokka to 
#annotate genome

# --outdir is where the results of prokka will go

os.system('prokka --outdir ' + path + '/Prokka_Output/ --prefix Ecoli --genus Escherichia long_sequences.fasta')

# Write the prokka command to miniproject.log

prokka_command = "prokka --outdir " + path + "/Prokka_Output/ --prefix Ecoli --genus Escherichia long_sequences.fasta"
log = open(path + "/miniproject.log","a+")
for i in range(1):
 log.write(prokka_command + "\n" + "\n") #AFTER DUE DATE: "log_file.write"

#6. Writing the results of prokka as txt in the log file

Prokka_Output = open(path + "/Prokka_Output/Ecoli.txt").readlines()
for line in Prokka_Output:
 log_file.append(line + "/n") #AFTER DUE DATE: wasn't sure how to add this the same way as above. Maybe should have kept it as ".write" wasnt sure if it would override previous entries

#7. Discrepancy between assmbled genome in RefSeq and prokka annotation
#First storing the number of CDS and tRNAs that Prokka found
CDS = int(path + '/Prokka_Output/'[4].strip().split(" ")[1])
tRNA = int(path + 'Prokka_Output/'[5].strip().split(" ")[1])

#since we already know how many CDS and tRNAs are in the RefSeq file..
if CDS > 4140:
 log_file.append("Prokka found an additional " + str(4140-CDS) + " fewer CDS than the RefSeq ") #".write"

if tRNA > 89:
 log_file.append(" and " + str(tRNA-89) + "additional tRNA than the RefSeq " + "/n") #".write"

else:
 log_file.append(" and " + str(89-tRNA) + " fewer tRNA than the RefSeq " + "/n") #".write"

#8.Use TopHat and Cufflinks to map the reads of the of a K-12 derivative 
#and quantify their expression

#Retrieve necessary sequence(s)
os.system("prefetch SRR1411276")

#convert into fastq format via fastq-dump
os.system("fastq-dump SRR1411276")

#Retrieve .fna file
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna")

#Building a bowtie2 index from complete genome to use with TopHat
os.system("bowtie2-build NC_000913.fna EcoliK12")

#Calling TopHat and Cufflinks
os.system("tophat2 -o " + path + "/SRR1411276_output --no-novel-juncs -p 2 EcoliK12 SRR1411276_1.fastq")

os.system("cufflinks -o " + path + "/cufflinks_output -p 2 accepted_hits.bam")

#opening the results of Cufflinks to Parse data.
data = open(path + "/cufflinks_output/transcripts.gtf").readlines()
Out = open(path + "/miniproject.fpkm", "w")
#this will iterate through the results of Cufflinks
#will write desired data () #tried using a different directory to store output of this data.

#found something similar to this deep in a rabbit hole about how to use/store cufflinks data.
for line in data:
 line_data = line.split("\t")
 Out.write(line_data[0] + ", " + line_data[3] + ", " + line_data[4] + ", " + line_data[6] + ", ")
 attributes = line_data[8].split("; ")
 for attribute in attributes:
 if "FPKM" in attribute:
 Out.write(attribute + " \n")
log_file.close()
Out.close()

#Never really got it to work all at once. I'm thinking that this was likely due to having all of the tools installed about 10 different ways. Mac versions. conda versions. linux versions on server. etc.
