# SimReadTruncs

### Simulate read truncations observed in ONT direct RNA sequencing data.

This Rscript simulates read truncations characteristic of Oxford Nanopore direct RNA sequencing data. It generates reads with transcript coverage distributions similar to real reads by utilising a linear regression and kernel density estimates to model the lengths and 3' end sites of simulated reads.

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
Rscript SimReadTruncs.R --help
Usage: SimReadTruncs.R [options]
Options:
        -f CHARACTER, --fasta=CHARACTER
                Transcriptome reference FASTA

        -c CHARACTER, --counts_file=CHARACTER
                Read count file (txt, tsv or csv). Header must be: transcript_id, read_count

        -l CHARACTER, --read_lengths=CHARACTER
                Read lengths option (sirv, human). Leave blank for default (human).

        -o CHARACTER, --output_file=CHARACTER
                Output FASTA file (output.fasta)

        -h, --help
                Show this help message and exit

# run on provided SIRV data
Rscript SimReadTruncs.R -f sirv_data/sirv_transcriptome_c.fa -c sirv_data/sirv_counts_27k_reads.csv -l sirv -o truncated_sirv_reads.fasta

# run on custom counts for human data
Rscript SimReadTruncs.R -f ref_transcriptome.fa -c read_counts.csv -l human -o truncated_reads.fasta
```
Note that the 'models' folder must be in the same directory as the Rscript.


#### Introduce errors with Badread:
```
# Using the custom model provided (trained on RNA004 SIRVs):
badread simulate --reference truncated_reads.fasta --quantity 1x --start_adapter_seq "" --end_adapter_seq "" --error_model models/rna004_error_model --junk 0 --random 0 --length 10000,10000 > truncated_w_errors_reads.fastq
# Using Badread's ONT model:
badread simulate --reference truncated_reads.fasta --quantity 1x --start_adapter_seq "" --end_adapter_seq "" --error_model nanopore2023 --junk 0 --random 0 --length 10000,10000 > truncated_w_errors_reads.fastq
```

#### Introduce errors with NanoSim (not tested):
```

simulator.py transcriptome -rt truncated_reads.fasta --model_prefix pre-trained_models/human_NA12878_dRNA_Bham1_guppy/training -o truncated_w_errors_reads
```

#### Map to the transcriptome with minimap2:
```
minimap2 -ax map-ont -N 10 ref_transcriptome.fasta truncated_w_errors_reads.fastq | samtools view -bh -F 2052 > reads.bam
```



