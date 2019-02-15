# This Rscript is meant solely to expand ids from consolidated sequence files 
# and write all identified ids to a output .txt file.
# 
# Usage: Rscript expand_key_to_ids.R {R1.fasta.gz} {K1.csv} {R2.fasta.gz} {K2.csv} {output.txt}
# 

options(stringsAsFactors = FALSE, scipen = 99, width = 999)

args <- commandArgs(trailingOnly = TRUE)

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

source(file.path(code_dir, "supporting_scripts/readKeyFile.R"))

R1_names <- as.character(ShortRead::id(ShortRead::readFasta(args[1])))

K1_tbl <- readKeyFile(args[2], stringr::str_extract(args[2], "[\\w]+$"))

R2_names <- as.character(ShortRead::id(ShortRead::readFasta(args[3])))

K2_tbl <- readKeyFile(args[4], stringr::str_extract(args[4], "[\\w]+$"))

K1_filtered <- K1_tbl$readNames[K1_tbl$seqID %in% R1_names]
K2_filtered <- K2_tbl$readNames[K2_tbl$seqID %in% R2_names]

key_readNames <- unique(c(K1_filtered, K2_filtered))

write.table(
  data.frame(ID = key_readNames), 
  file = args[5], 
  quote = FALSE, 
  row.names = FALSE
)

q(save = "no", status = 0)
