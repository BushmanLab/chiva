# This Rscript is meant solely to expand two Read alignments, in the form of a 
# PSL table with no header given two key files for read names.
# 
# Usage: Rscript expand_key_to_ids.R {R1.psl} {K1.csv} {R2.psl} {K2.csv} {output.txt}
# 

options(stringsAsFactors = FALSE, scipen = 99, width = 999)

args <- commandArgs(trailingOnly = TRUE)

R1psl <- data.table::fread(
  args[1], header = FALSE, sep = "\t", data.table = FALSE
)

K1tbl <- data.table::fread(
  args[2], header = TRUE, sep = ",", data.table = FALSE
)

R2psl <- data.table::fread(
  args[3], header = FALSE, sep = "\t", data.table = FALSE
)

K2tbl <- data.table::fread(
  args[4], header = TRUE, sep = ",", data.table = FALSE
)

R1_positive_keys <- unique(R1psl[,10])
R2_positive_keys <- unique(R2psl[,10])

K1_filtered <- K1tbl$readNames[K1tbl$seqID %in% R1_positive_keys]
K2_filtered <- K2tbl$readNames[K2tbl$seqID %in% R2_positive_keys]

key_readNames <- unique(c(K1_filtered, K2_filtered))

write.table(
  data.frame(ID = key_readNames), 
  file = args[5], 
  quote = FALSE, 
  row.names = FALSE
)

q(save = "no", status = 0)
