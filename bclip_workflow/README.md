# bCLIP_workflow

This workflow is used for the analysis of the clip variation protocol, called bCLIP.
The following operations are performed on the specified samples:
- extraction of UMIs
- demultiplexing of fastq files
- removal of poly(A) tails
- creation of STAR index
- mapping of reads to genome with STAR
- sorting and indexing of `.bam` files
- collapsing reads according to UMIs
- sorting and indexing of deduplicated `.bam` files
- running FastqC for quality control
- estimating the duplication level of reads based on the mapped coordinates
- separation of bam files by the duplication level of reads
- preparation of the .bed file with exonic, intronic, intergenic etc genomic segments based on the .gtf annotation
- counting the number of reads mapped to those segments
