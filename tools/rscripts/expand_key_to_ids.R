# This Rscript is meant solely to expand two Read alignments, in the form of a 
# PSL table with no header given two key files for read names.
# 
# Usage: Rscript expand_key_to_ids.R {R1.psl} {K1.csv} {R2.psl} {K2.csv} {output.txt}
# 

options(stringsAsFactors = FALSE, scipen = 99, width = 999)

args <- commandArgs(trailingOnly = TRUE)

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

source(file.path(code_dir, "supporting_scripts/readPSL.R"))
source(file.path(code_dir, "supporting_scripts/readKeyFile.R"))

R1psl <- readPSL(args[1])

K1tbl <- readKeyFile(args[2], stringr::str_extract(args[2], "[\\w]+$"))

R2psl <- readPSL(args[3])

K2tbl <- readKeyFile(args[4], stringr::str_extract(args[4], "[\\w]+$"))

R1_positive_keys <- unique(R1psl$qName)
R2_positive_keys <- unique(R2psl$qName)

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
