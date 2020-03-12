# This script is meant to concatenate sequence files within a directory together
# into a single file. 
# 
# Usage Rscript concat_seq_files.R {input.directory} {output.fasta}
# 
# Files within the input directory should be in fasta or fasta.gz formats.
# 

options(stringsAsFactors = FALSE, scipen = 99, width = 999)


code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

desc <- yaml::yaml.load_file(
  file.path(code_dir, "descriptions/trim.desc.yml")
)

#' Set up and gather command line arguments
parser <- argparse::ArgumentParser(
  description = desc$program_short_description,
  usage = "Rscript trim.R <seqFile> [-h/--help, -v/--version] [optional args]"
)


parser$add_argument(
  "-i", "--input", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", help = desc$output
)


library(tidyverse)
library(ggplot2)
library(dnaplotr)
library(dnar)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

plotOrder <- c("numSeqs","trim1P","primer2P","LTR2P","trim2P","filt1P","filt2P","consol1P","consol2P","psl1","psl2")

df <- read.table(args$input, stringsAsFactors = FALSE, sep = ",", header = TRUE)

heatmapDf <- df %>%
  select(-read1File, -read2File) %>%
  gather("Variable","Value",-1) %>%
  group_by(Variable) %>%
  mutate(Value=if_else(Value<1,Value,Value/max(Value))) %>%
  ungroup() %>%
  mutate(Variable=factor(Variable,levels=plotOrder))

pdf(args$output)

ggplot(heatmapDf, aes(x=Variable, y=SampleID, fill=Value)) +
  geom_tile()

mapply(function(r1f,r2f) {
  print(r1f)
  r1Reads <- as.character(read.fastq(r1f)$seq)
  if(length(r1Reads)>1000) { r1Reads <- sample(r1Reads, 1000) }
  plotDNA(sort(r1Reads), ylab = basename(r1f))
  print(r2f)
  r2Reads <- as.character(read.fastq(r2f)$seq)
  if(length(r2Reads)>1000) { r2Reads <- sample(r2Reads, 1000) }
  plotDNA(sort(r2Reads), ylab = basename(r2f))
}, df$read1File, df$read2File)

dev.off()

q()

