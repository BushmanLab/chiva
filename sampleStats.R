library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
sampleInfoFile <- args[1]
analysisDir <- args[2]

sampleInfo <- read.table(sampleInfoFile, sep = ",", header = TRUE, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(sample=stringr::str_split(stringr::str_split(sampleName,"-")[[1]][1],"_")[[1]][1]) %>%
  ungroup() %>% as.data.frame()

summary <- do.call(rbind, lapply(unique(sampleInfo$sample), function(sample) {
  U3reads = as.integer(system(paste0("zcat ", analysisDir, sample, "*[uU]3*.R1.fastq.gz | wc -l"), intern = TRUE)) / 4
  U5reads = as.integer(system(paste0("zcat ", analysisDir, sample, "*[uU]5*.R1.fastq.gz | wc -l"), intern = TRUE)) / 4
  U3R2primerReads = as.integer(system(paste0("zcat ", analysisDir, sample, "*[uU]3*.R2.trim.primer.fastq.gz | wc -l"), intern = TRUE)) / 4
  U5R2primerReads = as.integer(system(paste0("zcat ", analysisDir, sample, "*[uU]5*.R2.trim.primer.fastq.gz | wc -l"), intern = TRUE)) / 4
  U3R2LTRreads = as.integer(system(paste0("zcat ", analysisDir, sample, "*[uU]3*.R2.trim.ltr.fastq.gz | wc -l"), intern = TRUE)) / 4
  U5R2LTRreads = as.integer(system(paste0("zcat ", analysisDir, sample, "*[uU]5*.R2.trim.ltr.fastq.gz | wc -l"), intern = TRUE)) / 4
  
  return(data.frame(sample=sample,
                    U3reads=U3reads,
                    U3R2primerReads=U3R2primerReads,
                    U3R2LTRreads=U3R2LTRreads,
                    U5reads=U5reads,
                    U5R2primerReads=U5R2primerReads,
                    U5R2LTRreads=U5R2LTRreads))
}))

write.table(summary,
            paste0(args[1],".summary.csv"), sep = ",", row.names = FALSE, quote = FALSE)

