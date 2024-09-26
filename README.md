# SimReadTruncs

### Simulate read truncations observed in ONT direct RNA sequencing data (PCR-cDNA under dev).

This Rscript simulates read truncations characteristic of Oxford Nanopore direct RNA sequencing data. It generates reads with transcript coverage distributions similar to real reads by utilising kernel density estimates to model the lengths and 3' end sites of simulated reads.

#### Dependencies (R):
- Biostrings (Bioconductor)
- data.table
- stats
- optparse
- dplyr
- this.path

#### Download:
```
git clone https://github.com/josiegleeson/SimReadTruncs.git
```

#### Usage:
```
Rscript SimReadTruncs.R -f ref_transcriptome.fasta -c read_counts.csv -l human -o truncated_reads.fasta
```
Note that the 'kde_data' folder must be in the same directory as the Rscript.


#### Introduce errors with Badread:
```
badread simulate --reference truncated_reads.fasta --quantity 1x --start_adapter_seq "" --end_adapter_seq "" --error_model nanopore2023 --junk 0 --random 0 --length 10000,10000 > truncated_w_errors_reads.fastq
```

#### Introduce errors with NanoSim:
```
example > truncated_w_errors_reads.fastq
```

#### Map to the transcriptome with minimap2:
```
minimap2 -ax map-ont -N 10 ref_transcriptome.fasta truncated_w_errors_reads.fastq | samtools view -bh -F 2052 > reads.bam
```










