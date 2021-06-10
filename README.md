# NanoWgsTE
## Description
NanoWgsTE was designed to identify insertions of TEs and T-DNA in the genome using raw Nanopore whole genome sequencing data. It takes fastq file of the Nanopore reads, fasta of genome and fasta file with the target sequence. It outputs bed file with genomic coordinates of the putative insertions and .readcount file with the number of reads supporting each insertion.

## Dependencies
Please install the following packages and programs to run NanoWgsTE:
  - python
  - pysam
  - samtools
  - minimap2
  - biopython
  - bamtools
  - pandas
  - blast

We created conda environment.yml file to create a virtual environment for NanoWgsTE: 

`conda env create -f environment.yml`  
`conda activate nanowgste`  


## How to run
*positional arguments:*  
  **genome_fasta**          path to target sequence in fasta format  
  **fastq**                 path to fastq file of reads  
  **target_fasta**          path to target sequence in fasta format  
  **outBed**                path to the output bed file  

*optional arguments:*  
  -h, --help            show this help message and exit
  -mbh, --min_blast_hit   minimum length of BLAST hit for reads selection  
  -q, --map_q  minimum mapping quality  
  -mlc, --min_len_clipped   minimum length of the clipped part  
  -bp, --blastn_path        path to BLASTn program  
  -mdbp, --makeblastdb_path   path to makeblastdb program  
  -samtp, --samtools_path     path to samtools program  
  -bamtp, --bamtools_path      path to bamtools program  
  -mm2, --minimap2_path        path to minimap2 program  
  

## NanoWgsTE citation
We are currently preparing the manuscript.
