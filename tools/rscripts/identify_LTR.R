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

parser$add_argument(
  "-a", "--outputall", nargs = 1, type = "character", help = desc$output
)

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))


seqs <- as.character(ShortRead::sread(ShortRead::readFastq(args$input)))
if(length(seqs)==0) {
  write.table(character(0),
              args$output,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
  q()
}
primer <- names(sort(table(substr(seqs,1,8)), decreasing = TRUE)[1])
filtered_seqs <- seqs[substr(seqs,1,8)==primer]
LTRplus <- names(sort(table(filtered_seqs), decreasing = TRUE)[1])
potential_breaks <- gregexpr("CA",LTRplus)[[1]]

result <-
  system(paste0("bowtie2 --local -k 1 -x /data/internal/common/genomes/hg38/bowtie2/hg38 -U ",
                args$input,
                " | samtools view -S -F4 - | cut -f 6,10 | sort | uniq | cut -f 1 |",
                " grep ^[0-9][0-9]S | cut -f 1 -d S | sort | uniq -c | sed 's/^ *//'"), intern = TRUE)

if(is.null(result) | length(result)==0) {
  write.table(data.frame(primer="", breakpoint=0, LTRseq="", R1over="", count=0, sample=gsub(".R2.fastq.gz","",basename(args$input))),
              args$output,
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
  
  q()
}

resultDf <-
  do.call(rbind, lapply(result, function(X) { 
    vars <-  strsplit(X, " ")
    return(data.frame(count=as.numeric(vars[[1]][1]),
                      length=as.numeric(vars[[1]][2]),
                      stringsAsFactors = FALSE))
           })) %>%
  rowwise() %>%
  mutate(primer=primer,
         breakpoint=max(potential_breaks[potential_breaks <= length]),
         LTRseq=substr(LTRplus, 9, breakpoint-1),
         R1over=paste0("TG",dnar::revComp(substr(LTRseq, breakpoint-31,breakpoint)))) %>%
  ungroup() %>%
  group_by(primer, breakpoint, LTRseq, R1over) %>%
  summarize(count=sum(count)) %>%
  ungroup() %>%
  arrange(-count) %>%
  mutate(sample=gsub(".R2.fastq.gz","",basename(args$input)))

write.table(resultDf,
            args$output,
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
  
q()

plotSampleOverview <- function(R1) {
  R2 <- gsub("R1","R2", R1)
  R1seqs <- as.character(ShortRead::sread(ShortRead::readFastq(R1)))
  R2seqs <- as.character(ShortRead::sread(ShortRead::readFastq(R2)))
  plotDNA(c(sort(sample(R1seqs, 1000)),
            sort(sample(R2seqs, 1000))))
}

grabSampleReads <- function(R1) {
  R2 <- gsub("R1","R2", R1)
  R1seqs <- as.character(ShortRead::sread(ShortRead::readFastq(R1)))
  R2seqs <- as.character(ShortRead::sread(ShortRead::readFastq(R2)))
  return(c(sort(sample(R1seqs, 1000)),
           sort(sample(R2seqs, 1000))))
}

plotDNA(c(grabSampleReads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample9-U3-5.R1.fastq.gz"),
          grabSampleReads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample9-U3-6.R1.fastq.gz"),
          grabSampleReads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample9-U3-7.R1.fastq.gz"),
          grabSampleReads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample9-U3-8.R1.fastq.gz")))

grabSampleR2Reads <- function(R2) {
  R2seqs <- as.character(ShortRead::sread(ShortRead::readFastq(R2)))
  return(sort(sample(R2seqs, 100)))
}

plotDNA(c(grabSampleR2Reads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample9-U3-5.R2.fastq.gz"),
          grabSampleR2Reads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample9-U3-6.R2.fastq.gz"),
          grabSampleR2Reads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample9-U3-7.R2.fastq.gz"),
          grabSampleR2Reads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample9-U3-8.R2.fastq.gz"),
          grabSampleR2Reads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample3-U3-5.R2.fastq.gz"),
          grabSampleR2Reads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample3-U3-6.R2.fastq.gz"),
          grabSampleR2Reads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample3-U3-7.R2.fastq.gz"),
          grabSampleR2Reads("/home/kevin/projects/vincent/Vincent20180305/analysis_data/sample3-U3-8.R2.fastq.gz")))




q()

library(DECIPHER)
library(dendextend)
library(grid)
library(gridExtra)
library(ggdendro)
library(latticeExtra)
library(kableExtra)
library(stringr)
library(tidyverse)
library(dnaplotr)

ltr3 <- "TGGAAGGGCTAATTCACTCCCAAA"
ltr5 <- "AGTCAGTGTGGAAAATCTCTAGCA"

generate_report <- function(patient, gtsp, fastaFile, start = 8000, end = 9100, LTR_end = "TGGAAGGGCTAATTCACTCCCAAA") {

  fasta <- ShortRead::readFasta(fastaFile)
  seqs <- as.character(ShortRead::sread(fasta))
  seq_names <- as.character(ShortRead::id(fasta))
  
  NFL.index <- grep("NFL", seq_names)
  NFL_seqs <- seqs[NFL.index]
  NFL_names <- seq_names[NFL.index]
  
  LTR_terminal <- unlist(aregexec(pattern = LTR_end, text = substr(NFL_seqs, start, end), max.distance = 2))
  LTR_region_start <- start + LTR_terminal - 235
  LTR_region_end <- LTR_region_start + 300
  
  NFL_LTR_seqs <- substr(NFL_seqs, LTR_region_start, LTR_region_end)
  NFL_LTR_seqs_metric <- str_count(NFL_LTR_seqs, "[ACGT]")
  
  NFL_LTR_filter_index <- NFL_LTR_seqs_metric > 250 & substr(NFL_LTR_seqs, 235, 236) == "TG"
  
  NFL_LTR_seqs <- NFL_LTR_seqs[NFL_LTR_filter_index]
  NFL_seqs <- NFL_seqs[NFL_LTR_filter_index]
  NFL_names <- NFL_names[NFL_LTR_filter_index]
  LTR_region_start <- LTR_region_start[NFL_LTR_filter_index]
  LTR_region_end <- LTR_region_end[NFL_LTR_filter_index]
  
  NFL_df <- data.frame(sseqid=NFL_names,
                       LTR.start=LTR_region_start,
                       LTR.end=LTR_region_end)

  sample_df <- data.frame(SampleID=c("3-1","3-2","3-3","3-4","5-1","5-2","5-3","5-4"),
                          LTR=c(rep("3",4),rep("5",4)),
                          replicate=c(rep(seq(1,4),2)))
  
  primerTrimmedR2 <- do.call(rbind, lapply(sample_df$SampleID, function(sample) {
    fastq <- paste0("/home/kevin/projects/bar_hiv/",gtsp,"/analysis_data/",gtsp,"-U",sample,".R2.trim.primer.fastq.gz")
    reads <- ShortRead::readFastq(fastq)
    if(length(reads)>0) {
      return(data.frame(sample=sample,
                        id=as.character(ShortRead::id(reads)),
                        seq=as.character(ShortRead::sread(reads)),
                        stringsAsFactors = FALSE))
    } else {
      return(NULL)
    }
  }))
  
  LTRTrimmedR2 <- do.call(rbind, lapply(sample_df$SampleID, function(sample) {
    fastq <- paste0("/home/kevin/projects/bar_hiv/",gtsp,"/analysis_data/",gtsp,"-U",sample,".R2.trim.ltr.fastq.gz")
    reads <- ShortRead::readFastq(fastq)
    if(length(reads)>0) {
      return(data.frame(sample=sample,
                        id=as.character(ShortRead::id(reads)),
                        seq=as.character(ShortRead::sread(reads)),
                        stringsAsFactors = FALSE))
    } else {
      return(NULL)
    }
  }))
  
  uniq_sites_R2 <- do.call(rbind, lapply(sample_df$SampleID, function(sample) {
    uniqSitesFile <- paste0("/home/kevin/projects/bar_hiv/",gtsp,"/analysis_data/uniqSites/",gtsp,"-U",sample,".uniq.csv")
    uniqSites <-
      read.table(file = uniqSitesFile,
                 header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
      mutate(SampleID=sample)
  }))
  
  intSitesFile <- paste0("/home/kevin/projects/bar_hiv/",gtsp,"/output_data/xofil_condensed_sites.csv")
  int_sites <- read.table(file=intSitesFile,
                          header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  blast_results_R2 <- do.call(rbind, lapply(sample_df$SampleID, function(sample) {
    blastTableFile <- paste0("/home/kevin/projects/bar_hiv/",gtsp,"/analysis_data/",gtsp,"-U",sample,".blast.R2.filtered.tbl")
    blastTable <-
      read.table(file = blastTableFile,
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if(nrow(blastTable) > 0) {
      return(blastTable %>%
        mutate(SampleID=sample,
               read="R2") %>%
        left_join(NFL_df) %>%
        filter(sstart < LTR.end & send > LTR.start)
      )
    }
    return(blastTable)
  }))
  
    


#  missing_seq_names <- NFL_names[!NFL_names %in% colnames(plot.matrix)]
#  NFL_seqs <- NFL_seqs[!NFL_names %in% missing_seq_names]
#  NFL_LTR_seqs <- NFL_LTR_seqs[!NFL_names %in% missing_seq_names]
#  NFL_names <- NFL_names[!NFL_names %in% missing_seq_names]
  
  NFL_LTR_seqs_dist <- adist(NFL_LTR_seqs)
  colnames(NFL_LTR_seqs_dist) <- NFL_names
  rownames(NFL_LTR_seqs_dist) <- NFL_names
  
  hc <- hclust(as.dist(NFL_LTR_seqs_dist))
  dend <- as.dendrogram(hc)
  
  NFL_LTR_seqs_unique <- unique(NFL_LTR_seqs)
  NFL_LTR_seqs_unique_dist <- adist(NFL_LTR_seqs_unique)
  colnames(NFL_LTR_seqs_unique_dist) <- NFL_LTR_seqs_unique
  rownames(NFL_LTR_seqs_unique_dist) <- NFL_LTR_seqs_unique
  
  NFL_LTR_seqs_unique_df <-
    data.frame(seq=NFL_LTR_seqs,
               name=NFL_names) %>%
    group_by(seq) %>%
    summarize(count=n(),
              numN=max(str_count(seq, "N")),
              short=if_else(sum(grepl("SHORT",name))>0,TRUE,FALSE),
              hyper=if_else(sum(grepl("hyper",name))>0,TRUE,FALSE),
              CL1=if_else(sum(grepl("CL1",name))>0,TRUE,FALSE),
              CL2=if_else(sum(grepl("CL2",name))>0,TRUE,FALSE),
              ambig=if_else(sum(grepl("ambig",name))>0,TRUE,FALSE),
              none=if_else(sum(grepl("SHORT",name))+sum(grepl("hypermutant",name))+sum(grepl("CL1",name))+
                             sum(grepl("CL2",name))+sum(grepl("ambig",name))==0,TRUE,FALSE)) %>%
#    left_join(tsne_km_df, by="seq") %>%
    droplevels()
  
  # WHY IS SSEQID NULL for patient A02?
  intSitesMap <-
    uniq_sites_R2 %>%
    mutate(dummy=TRUE) %>%
    left_join(int_sites %>%
                rename(int_seq=seqnames,
                       int_start=start,
                       int_end=end) %>%
                select(int_seq, int_start, int_end, posid) %>%
                mutate(dummy=TRUE), by="dummy") %>%
    filter(seqnames==int_seq & start <= int_start & end >= int_end) %>%
    left_join(blast_results_R2, by=c("ID"="qseqid", "SampleID")) %>%
    filter(sseqid %in% labels(dend)) %>%
    mutate(sseqid=factor(sseqid, levels=labels(dend)))
  

  pdf(paste0("/home/kevin/projects/bar_hiv/",patient,".pdf"))
  

  par(mar=c(3,1,1,25))
  par(cex=0.5)
  #labels_colors(dend) <- as.numeric(NFL_LTR_seqs_unique_df$kmc[match(NFL_LTR_seqs[match(labels(dend), NFL_names)],NFL_LTR_seqs_unique_df$seq)])
  plot(dend, horiz = TRUE)
  #dendro.plot <- as.grob(expression(plot(dend, horiz = TRUE)))
  
  dnaplotr::plotDNA(NFL_LTR_seqs[match(labels(dend), NFL_names)])
  
  tryToPlotDNA <- function(reads) {
    if(nrow(reads) == 0) {
      return(NULL)
    }
    if(nrow(reads) > 1000) {
      reads <- reads %>% sample_n(1000)
    }
    dnaplotr::plotDNA(reads %>% arrange(seq) %>% pull(seq))
  }
  tryToPlotDNA(primerTrimmedR2 %>% filter(grepl("3",sample)))
  tryToPlotDNA(primerTrimmedR2 %>% filter(grepl("5",sample)))
  tryToPlotDNA(LTRTrimmedR2 %>% filter(grepl("3",sample)))
  tryToPlotDNA(LTRTrimmedR2 %>% filter(grepl("5",sample)))
  
  tableTheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.5), padding=unit.c(unit(1, "mm"), unit(1, "mm"))),
    colhead = list(fg_params=list(cex = 0.5), padding=unit.c(unit(1, "mm"), unit(1, "mm"))),
    rowhead = list(fg_params=list(cex = 0.5), padding=unit.c(unit(1, "mm"), unit(1, "mm"))))
  
  if(nrow(intSitesMap) > 0) {
    intSitesMatrix <- intSitesMap %>%
      group_by(sseqid, posid) %>%
      summarize(maxscore=max(bitscore)) %>%
      spread(sseqid, maxscore) %>%
      as.data.frame()
    rownames(intSitesMatrix) <- intSitesMatrix$posid
    intSitesMatrix <- intSitesMatrix %>% select(-posid)
    intSitesMatrix <- as.matrix(intSitesMatrix)
    if(ncol(intSitesMatrix) < length(labels(dend))) {
      matrixToAppend <- matrix(NA, nrow = nrow(intSitesMatrix),
                               ncol = length(labels(dend)) - ncol(intSitesMatrix))
      colnames(matrixToAppend) <- labels(dend)[!labels(dend) %in% colnames(intSitesMatrix)]
      intSitesMatrix <- cbind(intSitesMatrix, matrixToAppend)
    }
  
    intSitesMatrixToPlot <-
      as.matrix(intSitesMatrix[,match(labels(dend), colnames(intSitesMatrix))])
    intSitesMatrixTable <- t(intSitesMatrixToPlot[,ncol(intSitesMatrixToPlot):1])
    if(ncol(intSitesMatrixToPlot)==1) {
      intSitesMatrixToPlot <-
        t(as.matrix(intSitesMatrix[,match(labels(dend), colnames(intSitesMatrix))]))
      colnames(intSitesMatrixToPlot) <- colnames(intSitesMatrix)[match(labels(dend), colnames(intSitesMatrix))]
      rownames(intSitesMatrixToPlot) <- unique(intSitesDf$posid)
      intSitesMatrixTable <- as.matrix(intSitesMatrixToPlot[,ncol(intSitesMatrixToPlot):1])
      colnames(intSitesMatrixTable) <- unique(intSitesDf$posid)
      rownames(intSitesMatrixTable) <- rev(colnames(intSitesMatrixToPlot))
    }
    grid.newpage()
    grid.table(intSitesMatrixTable, theme = tableTheme)
    print(
      levelplot(intSitesMatrixToPlot,
                aspect = "fill",
                scales = list(x = list(rot = 90), y = list(cex = 0.3)),
                xlab = "Integration site",
                ylab = "NFL",
                colorkey = list(space = "left"),
                legend =
                  list(right =
                         list(fun = dendrogramGrob,
                              args =
                                list(x = dend, 
                                     side = "right",
                                     size = 3))))
    )
    
    if(FALSE) {
      
      for(sample in unique(intSitesMap$sampleName)) {
        
        intSitesDf <- intSitesMap %>% filter(sampleName==sample)
        intSitesMatrix <- as.matrix(table(intSitesDf$posid, intSitesDf$sseqid))
        
        intSitesMatrixToPlot <-
          as.matrix(intSitesMatrix[,match(labels(dend), colnames(intSitesMatrix))])
        intSitesMatrixTable <- t(intSitesMatrixToPlot[,ncol(intSitesMatrixToPlot):1])
        if(ncol(intSitesMatrixToPlot)==1) {
          intSitesMatrixToPlot <-
            t(as.matrix(intSitesMatrix[,match(labels(dend), colnames(intSitesMatrix))]))
          colnames(intSitesMatrixToPlot) <- colnames(intSitesMatrix)[match(labels(dend), colnames(intSitesMatrix))]
          rownames(intSitesMatrixToPlot) <- unique(intSitesDf$posid)
          intSitesMatrixTable <- as.matrix(intSitesMatrixToPlot[,ncol(intSitesMatrixToPlot):1])
          colnames(intSitesMatrixTable) <- unique(intSitesDf$posid)
          rownames(intSitesMatrixTable) <- rev(colnames(intSitesMatrixToPlot))
        }
        
        grid.newpage()
        grid.table(intSitesMatrixTable, theme = tableTheme)
        
        print(
          levelplot(intSitesMatrixToPlot,
                    aspect = "fill",
                    scales = list(x = list(rot = 90), y = list(cex = 0.3)),
                    xlab = "Integration site",
                    ylab = "NFL",
                    colorkey = list(space = "left"),
                    legend =
                      list(right =
                             list(fun = dendrogramGrob,
                                  args =
                                    list(x = dend, 
                                         side = "right",
                                         size = 3))))
        )
      }
      
      
      
      blast_results_paired <- do.call(rbind, lapply(sample_df$SampleID, function(sample) {
        blastTableFile <- paste0("/home/kevin/projects/bar_hiv/",gtsp,"/analysis_data/",gtsp,"-U",sample,".blast.paired.filtered.tbl")
        blastTable <-
          read.table(file = blastTableFile,
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        if(nrow(blastTable)>0) {
          return(blastTable %>%
                   mutate(SampleID=sample,
                          read="paired") %>%
                   left_join(NFL_df) %>%
                   filter((sstart.R1 < LTR.end & send.R1 > LTR.start) | (sstart.R2 < LTR.end & send.R2 > LTR.start))
          )
        }
        return(blastTable)
      }))
      
      blast_results_paired_summary <-
        blast_results_paired %>%
        group_by(SampleID, sseqid, rawseq.R1, rawseq.R2) %>%
        summarize(count=n()) %>%
        group_by(SampleID, sseqid) %>%
        summarize(uniq_seqs=n(), total_seqs=sum(count)) %>%
        droplevels() %>%
        ungroup() %>%
        mutate(sseqid=as.character(sseqid)) %>%
        left_join(sample_df)
      
      
      
      plot.matrix.df <-
        blast_results_paired_summary %>%
        select(sseqid, SampleID, total_seqs) %>%
        spread(sseqid, total_seqs)
      plot.matrix <- as.matrix(plot.matrix.df[,-1])
      plot.matrix[is.na(plot.matrix)] <- 0
      rownames(plot.matrix) <- plot.matrix.df$SampleID
      
        
      matrixToPlot <-
        t(as.matrix(plot.matrix[substr(sample, 11, 13),match(labels(dend), colnames(plot.matrix))]))
      colnames(matrixToPlot) <- colnames(plot.matrix)[match(labels(dend), colnames(plot.matrix))]
      rownames(matrixToPlot) <- sample
      matrixTable <- as.matrix(matrixToPlot[,ncol(matrixToPlot):1])
      colnames(matrixTable) <- sample
      rownames(matrixTable) <- rev(colnames(matrixToPlot))
      
      grid.newpage()
      grid.table(matrixTable, theme = tableTheme)
      
      print(
        levelplot(matrixToPlot,
                  aspect = "fill",
                  scales = list(x = list(rot = 90), y = list(cex = 0.3)),
                  xlab = "LTR (3, 5) - replicate (1-4)",
                  ylab = "NFL",
                  colorkey = list(space = "left"),
                  legend =
                    list(right =
                           list(fun = dendrogramGrob,
                                args =
                                  list(x = dend, 
                                       side = "right",
                                       size = 3))))
      )
    }
  }  
  dev.off()
}

generate_report("A01","GTSP2715","/home/kevin/dev/chiva/genomes/viral_genomes/patients/A01.fasta")
generate_report("A02","GTSP2716","/home/kevin/dev/chiva/genomes/viral_genomes/patients/A02.fasta")
generate_report("A13","GTSP2722","/home/kevin/dev/chiva/genomes/viral_genomes/patients/A13.fasta")



q()

#  blast_results_paired_summary2 <-
plotx <-
  blast_results_paired %>%
  group_by(SampleID) %>%
  mutate(prop=1/n()) %>%
  group_by(SampleID, qseqid) %>%
  summarize(propQseqs=sum(prop)) %>%
  left_join(sample_df) %>%
  ggplot(aes(x=propQseqs)) +
  geom_histogram() +
  facet_grid(LTR ~ replicate, scales="free_y")


#tsne <- Rtsne(NFL_LTR_seqs_unique_dist, perplexity = 3)

#  km <- kmeans(NFL_LTR_seqs_unique_dist, 5)
#  tsne_km_df <-
#    data.frame(seq=NFL_LTR_seqs_unique,
#               kmc=factor(km$cluster),
#               Y1=tsne$Y[,1],
#               Y2=tsne$Y[,2])

#  blast_results_paired %>%
#    group_by(SampleID) %>%
#    summarize(uniq_seqs = sum(uniq_seqs))

#  print(NFL_LTR_seqs_unique_df %>%
#          ggplot(aes(x=Y1,y=Y2, color=kmc)) +
#          geom_point())

#NFL_seqs_tree <- DECIPHER::IdClusters(NFL_LTR_seqs_dist, showPlot = TRUE, type = "dendrogram")

#  heatmap.plot <-
#    blast_results %>%
#    arrange(-uniq_seqs) %>%
#    ggplot(aes(x=SampleID, y=sseqid, fill=uniq_seqs)) +
#    geom_tile() +
#    facet_grid(. ~ read, scales = "free") +
#    theme (
#      axis.text.y = element_blank(),
#      axis.ticks.y = element_blank()
#    ) %>%
#    print()

#grid.newpage()
#grid.arrange(heatmap.plot, dendro.plot, nrow = 1)
#print(heatmap.plot, vp=viewport(x = 0.4, y=0.5, width = 0.8, height = 1.0))
#print(dendro.plot, vp=viewport(x = 0.9, y=0.445, width = 0.2, height = 1.0))

plot.matrix.df2 <-
  blast_results_paired %>%
  group_by(sseqid, best_count, SampleID) %>%
  summarize(count=n())

plot.matrix.df2 %>%
  ggplot(aes(x=best_count, y=SampleID, fill=count)) +
  geom_tile() +
  facet_grid(sseqid ~ .) +
  theme (
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle=0)
  )










generate_random_seqs <- function(nseqs, nchar) {
  sapply(seq(1,nseqs), function(i) {
    vals <- runif(nchar)
    paste0(sapply(vals, function(val) {
      if_else(val < 0.25, "A", if_else(val < 0.5, "C", if_else(val < 0.75, "G", "T")))
      }), collapse = "")
  })
}

HXB2 <- as.character(ShortRead::sread(ShortRead::readFasta("/home/kevin/dev/chiva/genomes/viral_genomes/HIV_HXB2.fasta")))
sapply(seq(15,30), function(j) {
  print(j)
  fragLength <- j
  HXB2frags <- sapply(seq(1,nchar(HXB2)-fragLength), function(i) {
    substr(HXB2,i,i+fragLength)
  }) %>% unique()
  fragmentsFile <- tempfile()
  ShortRead::writeFasta(Biostrings::DNAStringSet(HXB2frags), fragmentsFile)
  as.numeric(system(paste0("bowtie2 -f -k 1 -x /data/internal/common/genomes/hg38/bowtie2/hg38 -U ", fragmentsFile,
                           " | samtools view -S -F4 - | wc -l"), intern = TRUE))
})

args$input <- "/home/kevin/projects/bar_hiv/GTSP2715/analysis_data/GTSP2715-U3-1.R2.fastq.gz"

args$input <- "/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3042-2_u3.R2.fastq.gz"
args$input <- "/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3042-2_u3.R2.filt.fastq.gz"

args$input <- "/home/kevin/projects/bar_hiv/GTSP2715/analysis_data/GTSP2715-U3-1.R2.fastq.gz"
seqs.3.1 <- as.character(ShortRead::sread(ShortRead::readFastq(args$input)))
primer.3.1 <- names(sort(table(substr(seqs.3.1,1,8)), decreasing = TRUE)[1])
filtered_seqs.3.1 <- seqs.3.1[substr(seqs.3.1,1,8)==primer.3.1]

args$input <- "/home/kevin/projects/bar_hiv/GTSP2715/analysis_data/GTSP2715-U3-2.R2.fastq.gz"
seqs.3.2 <- as.character(ShortRead::sread(ShortRead::readFastq(args$input)))
primer.3.2 <- names(sort(table(substr(seqs.3.2,1,8)), decreasing = TRUE)[1])
filtered_seqs.3.2 <- seqs.3.2[substr(seqs.3.2,1,8)==primer.3.2]

plotDNA(c(sort(sample(filtered_seqs.3.1,1000)),
          sort(sample(filtered_seqs.3.2,1000))))


args$input <- "/home/kevin/projects/bar_hiv/GTSP2715/analysis_data/GTSP2715-U5-1.R2.fastq.gz"
seqs.5.1 <- as.character(ShortRead::sread(ShortRead::readFastq(args$input)))
primer.5.1 <- names(sort(table(substr(seqs.5.1,1,8)), decreasing = TRUE)[1])
filtered_seqs.5.1 <- seqs.5.1[substr(seqs.5.1,1,8)==primer.5.1]

args$input <- "/home/kevin/projects/bar_hiv/GTSP2715/analysis_data/GTSP2715-U5-2.R2.fastq.gz"
seqs.5.2 <- as.character(ShortRead::sread(ShortRead::readFastq(args$input)))
primer.5.2 <- names(sort(table(substr(seqs.5.2,1,8)), decreasing = TRUE)[1])
filtered_seqs.5.2 <- seqs.5.2[substr(seqs.5.2,1,8)==primer.5.2]

plotDNA(c(sort(sample(filtered_seqs.5.1,1000)),
          sort(sample(filtered_seqs.5.2,1000))))





args$input <- "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U3-1.R2.fastq.gz"
seqs.3.1 <- as.character(ShortRead::sread(ShortRead::readFastq(args$input)))
primer.3.1 <- names(sort(table(substr(seqs.3.1,1,8)), decreasing = TRUE)[1])
filtered_seqs.3.1 <- seqs.3.1[substr(seqs.3.1,1,8)==primer.3.1]

args$input <- "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U3-2.R2.fastq.gz"
seqs.3.2 <- as.character(ShortRead::sread(ShortRead::readFastq(args$input)))
primer.3.2 <- names(sort(table(substr(seqs.3.2,1,8)), decreasing = TRUE)[1])
filtered_seqs.5.2 <- seqs.3.2[substr(seqs.3.2,1,8)==primer.3.2]

args$input <- "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U3-3.R2.fastq.gz"
seqs.3.3 <- as.character(ShortRead::sread(ShortRead::readFastq(args$input)))
primer.3.3 <- names(sort(table(substr(seqs.3.3,1,8)), decreasing = TRUE)[1])
filtered_seqs.3.3 <- seqs.3.3[substr(seqs.3.3,1,8)==primer.3.3]

args$input <- "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U3-4.R2.fastq.gz"
seqs.3.4 <- as.character(ShortRead::sread(ShortRead::readFastq(args$input)))
primer.3.4 <- names(sort(table(substr(seqs.3.4,1,8)), decreasing = TRUE)[1])
filtered_seqs.3.4 <- seqs.3.4[substr(seqs.3.4,1,8)==primer.3.4]

plotDNA(c(sort(sample(filtered_seqs.3.1,1000)),
          sort(sample(filtered_seqs.3.2,1000)),
          sort(sample(filtered_seqs.3.3,1000)),
          sort(sample(filtered_seqs.3.4,1000))))



plotDNA(unique(sort(filtered_seqs)))

plotDNA(c(unique(sort(silicano.allseqs.filtered2)),
          unique(sort(silicano.allseqs.filtered))))

silicano.sites.reads <- read.csv("/home/kevin/projects/silicano/Silicano20190905/analysis_data/uniqSites/GTSP3038-3_u3.uniq.csv",
                                 header = TRUE)$ID
silicano.allseqs.sites <- as.character(ShortRead::sread(silicano.allseqs2[ShortRead::id(silicano.allseqs2) %in% silicano.sites.reads]))
silicano.allseqs.sitesr1 <- as.character(ShortRead::sread(silicano.allseqs2r1[ShortRead::id(silicano.allseqs2r1) %in% silicano.sites.reads]))
plotDNA(unique(sort(silicano.allseqs.filtered2)))
plotDNA(unique(sort(silicano.allseqs.sites)))
plotDNA(sort(substr(silicano.allseqs.sites, 105, 160)))
plotDNA(unique(sort(silicano.allseqs.filtered2r1)))
plotDNA(unique(sort(silicano.allseqs.sitesr1)))
plotDNA(unique(sort(substr(silicano.allseqs.sitesr1, 49,151))))

# FIX CODE
# Needs to handle variable length inputs

get_consensus <- function(seqs) {
  paste0(sapply(seq(1,nchar(seqs[1])), function(i) {
    charTable <- names(sort(table(substr(seqs, i, i)), decreasing = TRUE))
    return(charTable[1])
    #return(if_else(length(charTable)==1, charTable[1], if_else(charTable[1]=="-",charTable[2],charTable[1])))
  }), collapse = "")
}

get_consensus_percentage <- function(seqs, consensus) {
  sapply(seq(1,nchar(consensus)), function(i) {
    sum(substr(seqs, i, i) == substr(consensus, i, i)) / length(seqs)
  })
}

A01.allSeqs <- unique(as.character(ShortRead::readFasta("/home/kevin/projects/bar_hiv/data/A01.NFL.SEQ.10.24.19.FM.fasta")@sread))

A01.NFL <- A01.allSeqs[stringr::str_count(A01.allSeqs, "-") < 1000]

A01.NFL.consensus <- get_consensus(A01.NFL)
names(A01.NFL.consensus) <- "A01.consensus"

ShortRead::writeFasta(DNAStringSet(A01.NFL.consensus), "/home/kevin/projects/bar_hiv/data/A01.consensus.fasta")

A01.NFL.consensus_percentage <- get_consensus_percentage(A01.NFL, A01.NFL.consensus)


args$read2 <- "/home/kevin/projects/bar_hiv/GTSP2715/analysis_data/GTSP2715-U5-3.R2.trim.primer.fastq.gz"
R2 <- sample(as.character(ShortRead::readFastq(args$read2)@sread), 200)

plotDNA(c(rep(substr(A01.NFL.consensus, 8300, 8700),200),
          sort(sample(R2,200)),
          rep(substr(A01.NFL.consensus, 400, 800),200)))

#





get_common_read_lengths <- function(reads) {
  rl <- as.numeric(names(sort(table(nchar(reads)), decreasing = TRUE)))
  return(head(rl))
}


R2.files <- c(
  "/home/kevin/projects/bar_hiv/GTSP2715/analysis_data/GTSP2715-U5-3.R2.trim.primer.fastq.gz",
  "/home/kevin/projects/bar_hiv/GTSP2715/analysis_data/GTSP2715-U3-2.R2.trim.primer.fastq.gz",
  "/home/kevin/projects/bar_hiv/GTSP2716/analysis_data/GTSP2716-U5-2.R2.trim.primer.fastq.gz",
  "/home/kevin/projects/bar_hiv/GTSP2716/analysis_data/GTSP2716-U3-4.R2.trim.primer.fastq.gz",
  "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U5-4.R2.trim.primer.fastq.gz",
  "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U3-4.R2.trim.primer.fastq.gz")

R2 <- lapply(R2.files, function(file) {
  return(as.character(ShortRead::readFastq(file)@sread))
})

rls <- lapply(R2, function(reads) {
  get_common_read_lengths(reads)
})


R2 <- R2[nchar(R2) >= readLength]


consensus <- get_consensus(R2)

consensus_percentage <- get_consensus_percentage(R2, consensus)

names(consensus_percentage) <- strsplit(consensus, split="")[[1]]

endPos <- 0
for(i in seq(1, readLength-2)) {
  if(min(consensus_percentage[i:i+2]) < 0.9) {
    break
  }
  if(names(consensus_percentage)[i+1]=="C" & names(consensus_percentage)[i+2]=="A") {
    endPos <- i
  }
}

consensusLTR <- substr(consensus, 1, endPos)

print(args$LTR)
print(consensusLTR)

tmpBam <- Rsamtools::BamFile("/home/kevin/projects/bar_hiv/processing/A01.self.sorted.bam")
tmpAln <- as.data.frame(Rsamtools::scanBam(tmpBam)[[1]])
tmpAln$flag <- if_else(tmpAln$flag > 256, as.numeric(tmpAln$flag - 256), as.numeric(tmpAln$flag))

plotDNA(c(sort(tmpAln %>% filter(is.na(rname) & flag==141) %>% sample_n(100) %>% pull(seq)),
          sort(tmpAln %>% filter(!is.na(rname)) %>% sample_n(100) %>% pull(seq)),
          rep("AGATCCACAGATCAAGGATATCTTGTCTTCGTTGGGAGTGAATTAGCCCTTC",100),
          substr(A01.NFL, 8280, 8430), 100))


A01_blast_output_raw <- rbind(read.table("/home/kevin/projects/bar_hiv/processing/A01.R1.tbl",
                                         header=FALSE, skip=0, comment.char = "", sep = "\t") %>%
                                mutate(read="R1"),
                              read.table("/home/kevin/projects/bar_hiv/processing/A01.R2.tbl",
                                         header=FALSE, skip=0, comment.char = "", sep = "\t") %>%
                                mutate(read="R2"))
colnames(A01_blast_output_raw) <- c("qseqid",
                                "qseq",
                                "qlen",
                                "qstart",
                                "qend",
                                "sseqid",
                                "sseq",
                                "sstrand",
                                "slen",
                                "sstart",
                                "sned",
                                "evalue",
                                "bitscore",
                                "score",
                                "length",
                                "pident",
                                "nident",
                                "mismatch",
                                "btop",
                                "read")

A01_blast_output_paired <- A01_blast_output_raw %>% filter(read=="R1") %>%
  full_join(A01_blast_output_raw %>% filter(read=="R2"),
            by="qseqid", suffix=c(".R1",".R2")) %>%
  filter(sseqid.R1==sseqid.R2 | is.na(sseqid.R1) | is.na(sseqid.R2)) %>%
  mutate(sseqid=if_else(is.na(sseqid.R1), sseqid.R2, sseqid.R1),
         paired=if_else(is.na(sseqid.R1) | is.na(sseqid.R2), FALSE, TRUE)) %>%
  select(-sseqid.R1, -sseqid.R2)

A01_blast_output <- A01_blast_output_paired %>%
  group_by(qseqid) %>%
  mutate(max_bitscore.R1=max(bitscore.R1),
         max_bitscore.R2=max(bitscore.R2),
         best=if_else(paired,max_bitscore.R1==bitscore.R1 & max_bitscore.R2==bitscore.R2,
                      if_else(is.na(qseq.R1), max_bitscore.R2==bitscore.R2, max_bitscore.R1==bitscore.R1)),
         best_count=sum(best)) %>%
  ungroup()



q()



plotDNA(sample())