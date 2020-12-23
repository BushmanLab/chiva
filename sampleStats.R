library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
sampleInfoFile <- args[1]
analysisDir <- args[2]

sampleInfo <- read.table(sampleInfoFile, sep = ",", header = TRUE, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(sample=str_split(str_split(sampleName,"-")[[1]][1],"_")[[1]][1]) %>%
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

write.table(,
            paste0(args[1],".summary.csv"), sep = ",", row.names = FALSE, quote = FALSE)

subsample_dir <- function(dir) {
  dataDir <- paste0(dir, "/analysis_data/")
  subsampleDir <- paste0(dir, "/subsample/")
  system(paste0("mkdir ", subsampleDir))
  fastqDir <- paste0(dir, "/subsample/fastq.gz/")
  system(paste0("mkdir ", fastqDir))
  csvDir <- paste0(dir, "/subsample/csv/")
  system(paste0("mkdir ", csvDir))
  pdfDir <- paste0(dir, "/subsample/pdf/")
  system(paste0("mkdir ", pdfDir))
  
  sampleInfoFile <- paste0(dir, ".sampleInfo.csv")
  if(!file.exists(sampleInfoFile)) sampleInfoFile <- gsub("_",".",sampleInfoFile)
  
  samples <- read.table(sampleInfoFile, sep = ",", header = TRUE, stringsAsFactors = FALSE) %>%
    filter(substr(sampleName, 1, 1) == "G") %>%
    mutate(R1=paste0(dataDir,sampleName,".R1.fastq.gz"),
           R2=paste0(dataDir,sampleName,".R2.fastq.gz"),
           GTSP=substr(sampleName, 1, 8))
  
  lapply(unique(samples$GTSP), function(GTSP) {
    print(paste0("Processing ", GTSP))
    subsample(dataDir, GTSP, fastqDir)
  })
}



