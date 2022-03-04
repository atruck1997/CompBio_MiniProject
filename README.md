# CompBio MiniProject

---

## Overview

This is a python script that will automate the assembly and annotation of reads from a resequencing project of the *Escherichia coli K*-12 strain. After annotation, the python script will use Tophat2 to map reads of the *E. coli* transcriptome project of a *K*-12 derivative. Subsequently, Cufflinks will be deployed to quantify these reads' expression.

---  

## Prerequisite Software Installation

It is recommended to install all software to the users $Home directory.

```
cd ~
```

- Python 3 - (must have biopython)
  
- SRA-toolkit [GitHub - ncbi/sra-tools: SRA Tools](https://github.com/ncbi/sra-tools)
  
- SPAdes [SPAdes &#8211; Center for Algorithmic Biotechnology](http://cab.spbu.ru/software/spades/)
  
- Prokka [GitHub - tseemann/prokka: Rapid prokaryotic genome annotation](https://github.com/tseemann/prokka)
  
- SAM tools [GitHub - samtools/samtools: Tools (written in C using htslib) for manipulating next-generation sequencing data](https://github.com/samtools/samtools)
  
- Tuxedo Suite
  
  - BowTie2 [Bowtie: An ultrafast, memory-efficient short read aligner](http://bowtie-bio.sourceforge.net/index.shtml)
    
  - TopHat2 http://ccb.jhu.edu/software/tophat/index.shtml
    
  - Cufflinks [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)
    

---

## Data

The script will retrieve and download all of the required files for the resequencing project.

---

## Executing the Pipeline

Two options

1. Download the python file and run it in your Linux/Ubuntu or Mac Terminal
  
  ```
  $ python3 miniprojectcode.py
  ```
  
2. Clone this repository and execute it
  
  ```
  $ git clone https:
  ```
  

---

## Results
