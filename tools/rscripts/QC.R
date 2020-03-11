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
  "--read1", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--read2", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--trim1", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--primer2", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--LTR2", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--trim2", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--consol1", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--consol2", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--key1", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--key2", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--psl1", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--psl2", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--ltr", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "-s", "--sampleID", nargs = 1, type = "character", help = desc$output
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

result <- data.frame(SampleID=args$sampleID)
result$read1 <- as.numeric(system(paste0("zcat ", args$read1, " | wc -l"), intern = TRUE)) / 4
result$read2 <- as.numeric(system(paste0("zcat ", args$read2, " | wc -l"), intern = TRUE)) / 4
result$trim1 <- as.numeric(system(paste0("zcat ", args$trim1, " | wc -l"), intern = TRUE)) / 4
result$primer2 <- as.numeric(system(paste0("zcat ", args$primer2, " | wc -l"), intern = TRUE)) / 4
result$LTR2 <- as.numeric(system(paste0("zcat ", args$LTR2, " | wc -l"), intern = TRUE)) / 4
result$trim2 <- as.numeric(system(paste0("zcat ", args$trim2, " | wc -l"), intern = TRUE)) / 4
result$consol1 <- as.numeric(system(paste0("cat ", args$consol1, " | wc -l"), intern = TRUE)) / 2
result$consol2 <- as.numeric(system(paste0("cat ", args$consol2, " | wc -l"), intern = TRUE)) / 2



write.csv(result, file = args$output,
          quote = FALSE, row.names = FALSE, col.names = TRUE)

q()

