# 

library(dplyr)
library(RMySQL)


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
  usage = "Rscript create_supp.R <sampleInfo.csv>"
)


parser$add_argument(
  "--sampleInfo", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--supp", nargs = 1, type = "character", help = desc$output
)

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

sampleInfo <- read.csv(args$sampleInfo)
samples <- as.character(unique(sampleInfo$sampleName))

dbConn <- dbConnect(MySQL(), group = 'specimen_management')
d <- dbGetQuery(dbConn, 'select * from gtsp')

samplesToUse <- sapply(d$SpecimenAccNum, function(specimen) {
  if_else(length(grep(specimen, samples))>0, specimen, NULL)
})

samplesToUse <- samplesToUse[!is.na(samplesToUse) & samplesToUse != ""]

if(length(samplesToUse)==0) {
  print("No samples found in GTSP database!")
  q()
}

supp <- d %>%
  filter(SpecimenAccNum %in% samplesToUse) %>%
  select(SpecimenAccNum, Patient, CellType)

colnames(supp) <- c("specimen", "patient", "celltype")

write.csv(supp, args$supp, row.names = FALSE, quote = FALSE)

q()
