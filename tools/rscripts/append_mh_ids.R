# This Rscript is a simple utility script that extracts the read ids from 
# multihit output objects and appends them to the read ids alread identified by
# the unique alignments.
# 
# Usage: Rscript append_mh_ids.R {input.uniq} {input.multi} {output}
# 
# {input.uniq} will be a text file with the header "ID" and a read name in each
#   following row.
#   
# {input.multi} will be an rds object that contains the multihit data from the
#   coupling script.
#   
# {output} will be the file path to a text file output. This output will have 
#   the header of ID, and each following row will have a unique read name.
# 

options(stringsAsFactors = FALSE, scipen = 99, width = 999)

suppressMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly = TRUE)

uniq_ids <- read.delim(args[1], header = TRUE)

mh_input <- readRDS(args[2])

mh_ids <- data.frame("ID" = names(mh_input$unclustered_multihits))

all_ids <- rbind(uniq_ids, mh_ids)

write.table(all_ids, file = args[3], quote = FALSE, row.names = FALSE)

q(save = "no", status = 0)
