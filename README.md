# SimReadTruncs

Simulate read truncations observed in ONT direct RNA and cDNA data.

### Dependencies (R):
- Biostrings
- data.table
- stats
- optparse
- dplyr

### Download:
```
git clone https://github.com/josiegleeson/SimReadTruncs.git
```

### Usage:
```
Rscript SimReadTruncs.R -f transcriptome.fasta -c read_counts.csv -l human -o truncated_reads.fasta
```

### Error model introduction with Badread:
```
badread simulate --reference truncated_reads.fasta --quantity 1x --start_adapter_seq "" --end_adapter_seq "" --error_model nanopore2023 --junk 0 --random 0 --length 3000,2000 > truncated_w_errors_reads.fasta
```

### Error model introduction with NanoSim:
```
code
```
