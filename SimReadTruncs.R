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

# fit lm
# import lm 
#lm_read_lengths <- lm(readlen ~ txlen, data = uhr_real_data)

# import lm and kde
# set script dir
dirpath <- dirname(this.path())

if (is.null(read_input) | read_input == "human") {
  lm_read_lengths <- readRDS(file = paste0(dirpath, "/kde_data/human_lm.rds"))
  kde_3p <- readRDS(file = paste0(dirpath, "/kde_data/kde_3_prime_end_uhr.rds"))
  message("Using provided human read length model")
} else if (read_input == "sirv") {
  lm_read_lengths <- readRDS(file = paste0(dirpath, "/kde_data/sirv_lm.rds"))
  kde_3p <- readRDS(file = paste0(dirpath, "/kde_data/kde_3_prime_end_sirv.rds"))
  message("Using provided sirv read length model")
} else {
  message("Custom models not yet supported., please choose either 'human' or 'sirv'.")
  break()
}

original_tx_length_data <- data.frame(txlen = width(counts_txs))

# predict expected read length for each transcript
expected_length <- predict(lm_read_lengths, original_tx_length_data)

# sample residuals from the KDE
# import lm resid kde
kde_residuals <- density(residuals(lm_read_lengths))
sampled_residuals <- sample(kde_residuals$x, size = length(counts_txs), prob = kde_residuals$y, replace = TRUE)

sampled_3p_ends <- sample(kde_3p$x, size = length(counts_txs), prob = kde_3p$y, replace = TRUE)

# simulated read lengths df
simulated_lengths_df <- data.frame(txlen = original_tx_length_data$txlen,
                                   simulated_read_lengths = as.integer(expected_length + sampled_residuals),
                                   simulated_3p_ends = as.integer(sampled_3p_ends))

# update simulated length to handle neg values and reads longer than original tx length
simulated_lengths_df <- simulated_lengths_df %>%
  mutate(updated_simulated_read_lengths = case_when(
    simulated_read_lengths > txlen ~ txlen - 20,
    simulated_read_lengths < 0 ~ txlen - 20,
    TRUE ~ simulated_read_lengths
  ))

simulated_read_lengths_int <- as.vector(simulated_lengths_df$updated_simulated_read_lengths)
simulated_3p_ends_int <- as.vector(simulated_lengths_df$simulated_3p_ends)

# rename for testing
reads_input <- counts_txs

# empty list to populate
truncated_reads <- vector("list", length(reads_input))

# create progress bar
pb <- txtProgressBar(min = 0, max = length(reads_input), style = 3)

for (i in seq_along(reads_input)) {
  
  seq <- reads_input[i]
  # get original transcript length
  seq_length <- width(seq)
  
  # ensure a read is at least 50nt long
  sim_length <- max(simulated_read_lengths_int[i], 50)
  sim_3p_end <- min(simulated_3p_ends_int[i], seq_length - 1)
  
  # get start and end position for read truncation
  start <- max(1, seq_length - sim_length + 1)
  end <- max(start, seq_length - sim_3p_end)
  
  # final check to ensure valid coordinates
  if (start > end || end > seq_length) {
    truncated_reads[[i]] <- seq  # return full sequence if invalid
  } else {
    # return the truncated sequence
    truncated_reads[[i]] <- subseq(seq, start = start, end = end)
  }
  
  # update progress bar
  setTxtProgressBar(pb, i)
  
}

close(pb)

# ensure there are no names in list
names(truncated_reads) <- NULL

# convert to one DNAStringSet object
truncated_reads_stringset <- do.call(c, truncated_reads)

# write out FASTA
writeXStringSet(truncated_reads_stringset, paste0(output), append=FALSE, compress=FALSE, compression_level=NA, format="fasta")

message(paste0("Created ", output))

