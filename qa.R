library(tidyverse)
library(ggplot2)
library(ShortRead)

visualization_report <- function(dir, sampleInfoFile) {
  samples <- read.table(sampleInfoFile, sep = ",", header = TRUE, stringsAsFactors = FALSE) %>%
    filter(substr(sampleName, 1, 1) == "G") %>%
    mutate(R1=paste0(paste0(dir,"/analysis_data/"),sampleName,".R1.fastq.gz"),
           R2=paste0(paste0(dir,"/analysis_data/"),sampleName,".R2.fastq.gz"),
           GTSP=substr(sampleName, 1, 8))
  
  pdf(paste0(dir,".vis.pdf"))
  
  do.call(rbind, lapply(unique(samples$GTSP), function(thisGTSP) {
    # u3r1 <- samples %>% filter(GTSP==thisGTSP & uniqueRegion=="U3") %>% pull(R1)
    # print(dada2::plotQualityProfile(u3r1))
    # u3r2 <- samples %>% filter(GTSP==thisGTSP & uniqueRegion=="U3") %>% pull(R2)
    # print(dada2::plotQualityProfile(u3r2))
    # u5r1 <- samples %>% filter(GTSP==thisGTSP & uniqueRegion=="U5") %>% pull(R1)
    # print(dada2::plotQualityProfile(u5r1))
    # u5r2 <- samples %>% filter(GTSP==thisGTSP & uniqueRegion=="U5") %>% pull(R2)
    # print(dada2::plotQualityProfile(u5r2))
    
    #print(dada2::plotQualityProfile(c(samples %>% filter(GTSP==thisGTSP) %>% pull(R1), samples %>% filter(GTSP==thisGTSP) %>% pull(R2))))
    
    u3 <- readFastq(paste0(dir,"/subsample/fastq.gz/",thisGTSP,"_u3.R2.fastq.gz"))
    dnaplotr::plotDNA(sort(sample(as.character(sread(u3)), min(1000,length(u3)))), xlab = paste0("U3 R2 reads for ",thisGTSP))
    
    u5 <- readFastq(paste0(dir,"/subsample/fastq.gz/",thisGTSP,"_u5.R2.fastq.gz"))
    dnaplotr::plotDNA(sort(sample(as.character(sread(u5)), min(1000,length(u5)))), xlab = paste0("U5 R2 reads for ",thisGTSP))
    
  }))
  
  dev.off()
  
}

subsample <- function(dir, GTSP, outdir, sampleSize = 100000) {
  u3fastq <- readFastq(dir, paste0(GTSP, ".*u3.R2.fastq.gz"))
  u3ss <- sample(u3fastq, min(sampleSize, length(u3fastq)), replace = FALSE)
  u3outpath <- paste0(outdir, GTSP, "_u3.R2.fastq.gz")
  writeFastq(u3ss, u3outpath)
  u5fastq <- readFastq(dir, paste0(GTSP, ".*u5.R2.fastq.gz"))
  u5ss <- sample(u5fastq, min(sampleSize, length(u5fastq)), replace = FALSE)
  u5outpath <- paste0(outdir, GTSP, "_u5.R2.fastq.gz")
  writeFastq(u5ss, u5outpath)
}

subsample_dir <- function(dir, sampleInfoFile) {
  dataDir <- paste0(dir, "/analysis_data/")
  outDir <- paste0(dir, "/subsample/")
  system(paste0("mkdir ", outDir))
  outDir <- paste0(dir, "/subsample/fastq.gz/")
  system(paste0("mkdir ", outDir))
  outDir <- paste0(dir, "/subsample/csv/")
  system(paste0("mkdir ", outDir))
  outDir <- paste0(dir, "/subsample/pdf/")
  system(paste0("mkdir ", outDir))
  samples <- read.table(sampleInfoFile, sep = ",", header = TRUE, stringsAsFactors = FALSE) %>%
    filter(substr(sampleName, 1, 1) == "G") %>%
    mutate(R1=paste0(dataDir,sampleName,".R1.fastq.gz"),
           R2=paste0(dataDir,sampleName,".R2.fastq.gz"),
           GTSP=substr(sampleName, 1, 8))
  
  lapply(unique(samples$GTSP), function(GTSP) {
    subsample(dataDir, GTSP, outDir)
  })
}

subsample_dir("/home/kevin/projects/persaud/Persaud20201120")
subsample_dir("/home/kevin/projects/persaud/Persaud20201123")




visualization_report("/home/kevin/projects/persaud/Persaud20201120")
visualization_report("/home/kevin/projects/persaud/Persaud20201123")

subsample_dir("/home/kevin/projects/jones/Jones_20200928",
              "/home/kevin/projects/jones/Jones.20200928.sampleInfo.csv")
visualization_report("/home/kevin/projects/jones/Jones_20200928",
                     "/home/kevin/projects/jones/Jones.20200928.sampleInfo.csv")

