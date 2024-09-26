suppressWarnings({ 
  suppressPackageStartupMessages({
    library(Biostrings)
    library(data.table)
    library(stats)
    library(optparse)
    library(dplyr)
    library(this.path)
  })
})
  
option_list = list(
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Transcriptome reference FASTA", metavar="character"),
  make_option(c("-c", "--counts_file"), type="character", default=NULL,
              help="Read count file (txt, tsv or csv). Header must be: transcript_id, read_count", metavar="character"),
  make_option(c("-l", "--read_lengths"), type="character", default="human",
              help="Read lengths option (sirv, human) or file (txt, tsv or csv). Header must be: lengths. Leave blank for human.", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default="output.fasta",
              help="Output FASTA file (output.fasta)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# read in from command line
txs <- opt$fasta
count_input <- opt$counts_file
read_input <- opt$read_lengths
output <- opt$output_file

message("Importing data...")

# read in counts
counts <- fread(count_input)

# import reference fasta
tx_seqs <- readDNAStringSet(txs)

# remove 0 counts
counts <- counts %>% dplyr::filter(read_count > 0)

# remove counts for transcripts not found in reference
counts <- counts %>% dplyr::filter(transcript_id %in% names(tx_seqs))

# arrange reference by counts file transcript_ids
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

# kde <- density(bamslam$read_length, from=200, to=20000, na.rm = TRUE)
# saveRDS(kde, file="~/Documents/drna_sim/kde_sirv.rds")

# set script dir
dirpath <- dirname(this.path())
kde_3p <- readRDS(file = paste0(dirpath, "/kde_data/kde_3_prime_end.rds"))

# import kde
if (is.null(read_input) | read_input == "human") {
  kde_read_length <- readRDS(file = paste0(dirpath, "/kde_data/kde_read_length_human_brain.rds"))
  kde_long_read_length <- readRDS(file = paste0(dirpath, "/kde_data/kde_long_read_length_human_brain.rds"))
  message("Using provided human read length model")
} else if (read_input == "sirv") {
  kde_read_length <- readRDS(file = paste0(dirpath, "/kde_data/kde_read_length_sirv.rds"))
  kde_long_read_length <- readRDS(file = paste0(dirpath, "/kde_data/kde_long_read_length_sirv.rds"))
  message("Using provided sirv read length model")
} else {
  custom_data <- fread(read_input)
  length_data$nt_col <- length_data[,1]
  length_data$nt_col <- as.numeric(length_data$nt_col)
  
  kde_read_length <- density(length_data$nt_col, from=200, to=20000, na.rm = TRUE)
}

# function to sample from kde
sample_from_kde <- function(kde) {
  
  # get length
  sample_length <- sample(kde$x, size = 1, prob = kde$y)
  return(sample_length)
  
}

# function to remove nt from the 5' end of sequence
remove_kde_length <- function(seq, max_attempts = 10) {
  
  # get tx length
  seq_length <- nchar(seq)
  
  # to add: 5% of the time, no truncation occurs
  random_number <- sample(1:20, 1)

  if (random_number == 10 & seq_length < 20000) {
    
    # no truncation
    return(seq)
    
  } else {
    
    # determine which kde to use
    if (seq_length > 3000) {
      kde_to_sample_from <- kde_long_read_length
    } else {
      kde_to_sample_from <- kde_read_length
    }
    
    # sample from kde until read length is less than tx length
    attempts <- 0
    repeat {
      # sample from read length kde defined above
      simulated_read_length <- floor(sample_from_kde(kde_to_sample_from))
      attempts <- attempts + 1
      if (simulated_read_length < seq_length | attempts >= max_attempts) break
    }
    
    # if no valid length after max_attempts, use seq_length - 10
    if (attempts >= max_attempts) {
      simulated_read_length <- seq_length - 10
    }
    
    # ensure simulated_read_length is at least 1
    simulated_read_length <- max(1, simulated_read_length)
    
    # generate 3' end truncation
    simulated_3p_truncation <- floor(sample_from_kde(kde_3p))
    
    # ensure simulated_3p_truncation is not larger than the sequence length
    simulated_3p_truncation <- min(simulated_3p_truncation, seq_length - 1)
    
    # get read start and end positions
    start <- max(1, seq_length - simulated_read_length + 1)
    end <- max(start, seq_length - simulated_3p_truncation)  # ensure end is not less than start
    
    # final check to ensure start <= end <= seq_length
    if (start > end | end > seq_length) {
      # if coords are still invalid return the full sequence
      return(seq)
    }
    
    # return truncated read
    return(subseq(seq, start = start, end = end))
  }
  
}

message("Truncating reads...")

# apply to DNAStringSet
trunc_txs <- DNAStringSet(sapply(counts_txs, remove_kde_length))

# write out FASTA
writeXStringSet(trunc_txs, paste0(output), append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

message(paste0("Created ", output))

