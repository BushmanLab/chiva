# This script is meant to concatenate sequence files within a directory together
# into a single file. 
# 
# Usage Rscript concat_seq_files.R {input.directory} {output.fasta}
# 
# Files within the input directory should be in fasta or fasta.gz formats.
# 

library(dplyr)

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
  "--input", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", help = desc$output
)

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))


ltrAll <- read.table(args$input, sep = ",", header = TRUE, stringsAsFactors = FALSE)

ltrCombined <- ltrAll %>%
  filter(!is.na(LTRseq) & breakpoint >= 50) %>%
  mutate(SampleID=gsub("-[1-4]","",sample)) %>%
  group_by(SampleID, primer, breakpoint, LTRseq, R1over) %>%
  summarize(count=sum(count)) %>%
  arrange(primer, SampleID, -count) %>%
  group_by(SampleID) %>%
  mutate(ord=seq(1,n())) %>%
  filter(ord==1) %>%
  select(-ord) %>%
  ungroup() %>%
  select(SampleID, primer, LTRseq, R1over)

ltrDf <- data.frame(sampleName=unique(ltrAll$sample)) %>%
  mutate(SampleID=gsub("-[1-4]","",sampleName)) %>%
  left_join(ltrCombined) %>%
  rename(LTR=LTRseq,
         R1Over=R1over) %>%
  select(-SampleID)

write.csv(ltrDf, args$output, row.names = FALSE)
