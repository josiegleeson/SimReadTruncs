suppressWarnings({ 
  suppressPackageStartupMessages({
    library(Biostrings)
    library(data.table)
    library(stats)
    library(optparse)
    library(dplyr)
  })
})
  
option_list = list(
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Transcriptome reference FASTA", metavar="character"),
  make_option(c("-c", "--counts_file"), type="character", default=NULL,
              help="Read count file (txt, tsv or csv). Header must be: transcript_id, read_count", metavar="character"),
  make_option(c("-l", "--read_lengths"), type="character", default=NULL,
              help="Read lengths option (sirv, human) or file (txt, tsv or csv). Header must be: lengths. Leave blank for human.", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="Directory to be created and where output FASTA will be written", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# read in from command line
txs <- opt$fasta
count_input <- opt$counts_file
read_lengths <- "~/Documents/drna_sim/custom_read_lengths.csv"
output <- opt$output_file

# txs <- "~/Documents/drna_sim/sirvs/sirv_transcriptome_c.fa"
# count_input <- "~/Documents/drna_sim/sirvs/read_counts.tsv"
# output <- "~/Documents/compare_drs_data/simulated/sirv_sim.fa"

# txs <- "~/Documents/drna_sim/sirvs/sirv_transcriptome_c.fa"
# count_input <- "~/Documents/drna_sim/sirv_counts.csv"
# output <- "~/Documents/drna_sim/sirv_read_length/sirvs.fa"

message("Importing data...")

# read in counts
counts <- fread(count_input)

# import reference fasta
tx_seqs <- readDNAStringSet(txs)

# remove 0 counts
counts <- counts %>% dplyr::filter(read_count > 0)

# remove counts for transcripts not found in reference
counts <- counts %>% dplyr::filter(transcript_id %in% names(tx_seqs))

# arrange reference by counts file transcript_id
tx_seqs <- tx_seqs[counts$transcript_id]

# replicate transcript entry as per read_count
counts_txs <- rep(tx_seqs, as.integer(counts$read_count))

# give unique ids by appending an index to transcript_id
original_names <- names(counts_txs)

# apply
unique_ids <- unlist(lapply(1:length(original_names), function(i) {
  paste0("read", i, "_", original_names[i])
}))

# set the unique ids as names
names(counts_txs) <- unique_ids

# read length kde generation
# kde <- density(bamslam$read_length, from=200, to=20000, na.rm = TRUE)
# saveRDS(kde, file="~/Documents/drna_sim/kde_sirv.rds")

# set script dir
dirpath <- dirname(rstudioapi::getSourceEditorContext()$path)

# import kde
if (is.null(read_lengths) | read_lengths == "human") {
  kde <- readRDS(file = paste0(dirpath, "/kde_human.rds"))
} else if (read_lengths == "sirv") {
  kde <- readRDS(file = paste0(dirpath, "/kde_sirv.rds"))
} else {
  length_data <- fread(read_lengths)
  length_data$nt_col <- length_data[,1] 
  length_data$nt_col <- as.numeric(length_data$nt_col)
  kde <- density(length_data$nt_col, from=200, to=30000, na.rm = TRUE)
}

# sample from kde
sample_from_kde <- function(kde) {
  # length
  sample_read_length <- sample(kde$x, size = 1, prob = kde$y)
  
  return(sample_read_length)
}

# remove nt from the 5' end of sequence
remove_kde_length <- function(seq) {
  
  # get tx length
  seq_length <- nchar(seq)
  
  # get read length from kde
  simulated_read_length <- floor(sample_from_kde(kde))
  
  random_3p_end <- sample(1:10, 1)
  
  if (simulated_read_length >= seq_length & seq_length > 50 & seq_length < 1000) {
    # if nt to remove is more than sequence length (and length is <1k), return 20 nt truncation
    return(subseq(seq, start = 20, end = -abs(random_3p_end)))
  } else if (simulated_read_length >= seq_length & seq_length > 1000) {
    # if nt to remove is more than sequence length (and length is >1k), return 200 nt truncation
    return(subseq(seq, start = 200, end = -abs(random_3p_end)))
  } else if (seq_length < 50) { 
    # no truncation
    return(seq)
  } else { 
    # remove kde based number of nt
    return(subseq(seq, start = -abs(simulated_read_length), end = -abs(random_3p_end)))
  }
}

message("Truncating reads...")

# apply to DNAStringset
trunc_txs <- DNAStringSet(sapply(counts_txs, remove_kde_length))

# write out FASTA
writeXStringSet(trunc_txs, paste0(output), append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

# messages
message(paste0("Created ", output))

# Next step:
# badread simulate --reference sim.fa --quantity 1x --start_adapter_seq "" --end_adapter_seq "" --error_model nanopore2023 --junk 0 --random 0 --length 3000,2000 > sim_badread.fastq





