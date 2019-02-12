# This script is meant to concatenate sequence files within a directory together
# into a single file. 
# 
# Usage Rscript concat_seq_files.R {input.directory} {output.fasta}
# 
# Files within the input directory should be in fasta or fasta.gz formats.
# 

options(stringsAsFactors = FALSE, scipen = 99, width = 999)

args <- commandArgs(trailingOnly = TRUE)

input_files <- list.files(path = args[1], pattern = "fasta")

input_seqs <- Biostrings::readDNAStringSet(
  filepath = file.path(args[1], input_files), format = "fasta"
)

null <- Biostrings::writeXStringSet(
  x = input_seqs, 
  filepath = args[2], 
  compress = FALSE, 
  format = "fasta", 
  width = max(Biostrings::width(input_seqs))
)

if( file.exists(args[2]) ){
  cat("Constructed single reference. Completed.")
  q(save = "no", status = 0)
}else{
  stop(
    "\n  Could not write output file.", 
    "\n  Check input sequences in ", args[1], 
    "\n  and make sure they are all in 'fasta' or 'fasta.gz' format.\n"
  )
}
