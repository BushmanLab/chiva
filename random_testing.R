library(dnar)
library(dnaplotr)
library(dplyr)
library(Biostrings)

HBX2_file <- "/home/kevin/dev/chiva/genomes/viral_genomes/HIV_HXB2.fasta"
HBX2 <- read.fa(HBX2_file)
HBX2$seq <- toupper(HBX2$seq)

sampleR1_file <- "/home/kevin/projects/bar_hiv/GTSP2718-U3-1/analysis_data/GTSP2718-U3-1.R1.fastq.gz"
sampleR2_file <- "/home/kevin/projects/bar_hiv/GTSP2718-U3-1/analysis_data/GTSP2718-U3-1.R2.fastq.gz"
sampleR1_file <- "/home/kevin/projects/bar_hiv/GTSP2718-U3-1/analysis_data/GTSP2718-U3-1.R1.trim.fastq.gz"
sampleR2_file <- "/home/kevin/projects/bar_hiv/GTSP2718-U3-1/analysis_data/GTSP2718-U3-1.R2.trim.primer.fastq.gz"
sampleR1 <- read.fastq(sampleR1_file)
sampleR2 <- read.fastq(sampleR2_file)

plotDNA(c(sort(sample(sampleR1$seq,100)),
          rep(substr(HBX2$seq,1,500),100),
          sort(sample(sampleR2$seq,100))
          ))

mostCommonR1 <- names(table(sampleR1$seq) %>% sort(decreasing=TRUE) %>% head(1))
mostCommonR2 <- names(table(sampleR2$seq) %>% sort(decreasing=TRUE) %>% head(1))

plotDNA(c(substr(HBX2$seq,326,650),
          revComp(mostCommonR2),
          mostCommonR2,
          revComp(substr(HBX2$seq,326,650))
))

plotDNA(c(rep(substr(HBX2$seq,326,650),100),
          sort(sample(sampleR1$seq,100)),
          sort(sample(sampleR1$seq,100)),
          rep(revComp(substr(HBX2$seq,326,650)),100)
))

pairwiseAlignment(mostCommonR1, HBX2$seq, type = "local")

pairwiseAlignment(c(sample(sampleR2$seq,100),mostCommonR1), HBX2$seq, scoreOnly=TRUE)
pairwiseAlignment(c(sample(sampleR2$seq,100),mostCommonR1), HBX2$seq, type="local", scoreOnly=TRUE)

random_genome_snippets <- mapply(function(x, y) {
    return(substring(HBX2$seq, x, y))
  },
  sample(1:50, 100, replace = T),
  sample(101:150, 100, replace = T)
  )

pairwiseAlignment(random_genome_snippets, HBX2$seq, type="local", scoreOnly = TRUE)
pairwiseAlignment(random_genome_snippets, HBX2$seq, type="local")








args <- list()
args$read1 <- "/data/sequencing/Illumina-archive/190108_M05588_0153_000000000-C7K9D/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz"
args$read2 <- "/data/sequencing/Illumina-archive/190108_M05588_0153_000000000-C7K9D/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz"
args$idx1 <- "/data/sequencing/Illumina-archive/190108_M05588_0153_000000000-C7K9D/Data/Intensities/BaseCalls/Undetermined_S0_L001_I1_001.fastq.gz"
args$outfolder <- "/home/kevin/projects/jones/demult_test/"
args$cores <- 32

args$read1 <- "/home/kevin/projects/bar_hiv/tmpr1.fastq.gz"
args$read2 <- "/home/kevin/projects/bar_hiv/tmpr2.fastq.gz"
args$idx1 <- "/home/kevin/projects/bar_hiv/tmpidx.fastq.gz"
args$outfolder <- "/home/kevin/projects/bar_hiv/demult_test/"
args$cores <- 1

args$idx2 <- "NA"
args$bc1 <- "R1"
args$bc2 <- "I1"
args$singleBarcode <- FALSE
args$bc1Mis <- 3
args$bc2Mis <- 1
args$manifest <- "/home/kevin/projects/jones/190108_Jones_HIV.sampleInfo.csv"
args$bc1Man <- "barcode1"
args$bc2Man <- "barcode2"
args$readNamePattern <- "[\\w\\:\\-\\+]+"
args$maxN <- 1
args$bc1Len <- 20
args$bc2Len <- 12
args$stat <- FALSE
args$chunkSize <- 100
args$poolreps <- FALSE
args$compress <- TRUE



demulti <- data.frame(
  "readType" = c("R1", "R2", "I1", "I2"),
  "path" = c(args$read1, args$read2, args$idx1, args$idx2)
)


demulti$bc1 <- grepl(args$bc1, demulti$readType)
demulti$bc2 <- grepl(args$bc2, demulti$readType)


if( demulti$readType[demulti$bc1] == demulti$readType[demulti$bc2] ){
  stop("Please select different read types for barcodes 1 and 2.\n")
}

if( demulti$readType[demulti$bc1] == "NA" ){
  stop("Barcode 1 is set to a read type that is not provided.\n")
}

if( demulti$readType[demulti$bc2] == "NA" ){
  stop("Barcode 2 is set to a read type that is not provided.\n")
}

if( args$singleBarcode ){
  demulti$bc2 <- FALSE
}

if( !is.null(args$maxMis) ){
  args$bc1Mis <- args$maxMis
  args$bc2Mis <- args$maxMis
}


input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply( 
    seq_along(args), 
    function(i){
      paste(args[[i]], collapse = ", ")
    }
  )
)

input_table <- input_table[
  match(
    c("manifest :", "idx1 :", "idx2 :", "read1 :", "read2 :", 
      "outfolder :", "stat :", "bc1 :", "bc2 :", "bc1Man :", "bc2Man :", 
      "bc1Len :", "bc2Len :", "bc1Mis :", "bc2Mis :", "cores :", "compress :", 
      "poolreps :", "singleBarcode :",  "readNamePattern :"
    ),
    input_table$Variables
  ),
  ]

cat("Demultiplex Inputs:\n")
print(
  x = data.frame(input_table, row.names = NULL), 
  right = FALSE, 
  row.names = FALSE
)

# Create output directory if not currently available ----
if( !file.exists(args$outfolder) ){
  
  attempt <- try(system(paste0("mkdir ", args$outfolder)))
  if(attempt == 1) stop("Cannot create output folder.\n")
  
}

# Check for required packages ----
required_packs <- c("stringr", "ShortRead", "Biostrings")
present_packs <- required_packs %in% row.names(installed.packages())

if( !all(present_packs) ){
  
  cat("Missing required r-packages:\n")
  print(
    data.frame(
      "R-Packages" = required_packs, 
      "Installed" = present_packs, 
      row.names = NULL
    ), right = FALSE, row.names = FALSE)
  
  stop("Check dependancies.\n")
  
}



file_ext <- unlist(strsplit(args$manifest, "\\."))
file_ext <- file_ext[length(file_ext)]

if( file_ext %in% c("yaml", "yml") ){
  
  if( !"yaml" %in% row.names(installed.packages()) ){
    stop("Package:yaml not loaded or installed.\n")
  }
  
  manifest <- yaml::yaml.load_file(args$manifest)
  
  if( args$singleBarcode ){
    
    samples_df <- data.frame(
      "sampleName" = names(manifest$samples),
      "bc1" = sapply( manifest$samples, function(x) x[args$bc1Man] ),
      row.names = NULL
    )
    
  }else{
    
    samples_df <- data.frame(
      "sampleName" = names(manifest$samples),
      "bc1" = sapply( manifest$samples, function(x) x[args$bc1Man] ),
      "bc2" = sapply( manifest$samples, function(x) x[args$bc2Man] ),
      row.names = NULL
    )
    
  }
  
}else{
  
  if( file_ext == "csv" ){
    manifest <- read.csv(args$manifest)
  }else if( file_ext == "tsv" ){
    manifest <- read.delim(args$manifest)
  }
  
  if( args$singleBarcode ){
    samples_df <- manifest[, c("sampleName", args$bc1Man)]
    names(samples_df) <- c("sampleName", "bc1")
  }else{
    samples_df <- manifest[, c("sampleName", args$bc1Man, args$bc2Man)]
    names(samples_df) <- c("sampleName", "bc1", "bc2")
  }
  
}


if( !args$singleBarcode ){
  
  unique_samples <- nrow(samples_df[,c("bc1", "bc2")]) == 
    nrow(unique(samples_df[,c("bc1", "bc2")]))
  
  if( !unique_samples ) stop("Ambiguous barcoding of samples. Please correct.\n")
  
}else{
  
  unique_samples <- length(samples_df[,c("bc1")]) == 
    length(unique(samples_df[,"bc1"]))
  
  if( !unique_samples ) stop("Ambiguous barcoding of samples. Please correct.\n")
  
}







i1f <- "/data/sequencing/Illumina-archive/190108_M05588_0153_000000000-C7K9D/Data/Intensities/BaseCalls/Undetermined_S0_L001_I1_001.fastq.gz"
i1 <- FastqStreamer(i1f, 10000)
r1f <- "/data/sequencing/Illumina-archive/190108_M05588_0153_000000000-C7K9D/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz"
r1 <- FastqStreamer(r1f, 10000)
r2f <- "/data/sequencing/Illumina-archive/190108_M05588_0153_000000000-C7K9D/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz"
r2 <- FastqStreamer(r2f, 10000)

while (length(i1.chunk <- yield(i1))) {
  bc1_reads <- yield(r1)
  
  all_indices <- stringr::str_extract(
    as.character(ShortRead::id(bc1_reads)), 
    args$readNamePattern
  )

  if( !all(table(all_indices) == 1) ){
    stop(
      "\n  Read names are not unique, check input sequence files or ",
      "adjust readNamePattern parameter.\n")
  }

  cat(paste("\nReads to demultiplex : ", length(bc1_reads), "\n"))
  
  if( args$cores > 1 ){
    
    bc1_proc_grps <- split(
      bc1_reads,
      ceiling( seq_along(bc1_reads) / (length(bc1_reads)/args$cores) )
    )
    
    split_indices <- split(
      all_indices,
      ceiling( seq_along(all_indices) / (length(bc1_reads)/args$cores) )
    )
    
    cluster <- parallel::makeCluster(min(c(parallel::detectCores(), args$cores)))

    BC1_parsed_list <-  parallel::clusterMap(
      cluster,
      function(reads, idx, parseIndexReads, samples_df, args){
        parseIndexReads(
          barcode.seqs = samples_df$bc1,
          reads = reads,
          indices = idx,
          barcode.length = args$bc1Len,
          max.mismatch = args$bc1Mis,
          max.N.count = args$maxN
        )
      },
      reads = bc1_proc_grps,
      idx = split_indices,
      MoreArgs = list(
        parseIndexReads = parseIndexReads, 
        samples_df = samples_df,
        args = args
      ),
      SIMPLIFY = FALSE
    )
    
    
    BC1_parsed <- lapply(
      names(BC1_parsed_list[[1]]), function(x){
        unlist(lapply(seq_along(BC1_parsed_list), function(i){
          BC1_parsed_list[[i]][[x]]
        }))
      }
    )
    
    names(BC1_parsed) <- names(BC1_parsed_list[[1]])
    rm(BC1_parsed_list, bc1_proc_grps)
    
    cat("\nbc1 breakdown:\n")
    print(
      data.frame(
        "bc1" = names(BC1_parsed),
        "Read Counts" = lengths(BC1_parsed)
      ),
      right = TRUE, 
      row.names = FALSE
    )
    
    if( !args$singleBarcode ){
      
      bc2_reads <- yield(r2)
      
      bc2_indices <- stringr::str_extract(
        as.character(ShortRead::id(bc2_reads)), 
        args$readNamePattern
      )
      
      if( !all(bc2_indices == all_indices) ){
        warning(
          "  Index reads are not in the same order. Sequencing files should ",
          "always be kept in order across read types.\n")
      }
      
      bc2_proc_grps <- split(
        bc2_reads,
        ceiling( seq_along(bc2_reads) / (length(bc2_reads)/args$cores) )
      )
      
      split_bc2_indices <- split(
        bc2_indices,
        ceiling( seq_along(bc2_reads) / (length(bc2_reads)/args$cores) )
      )
      
      BC2_parsed_list <- parallel::clusterMap(
        cluster,
        function(reads, idx, parseIndexReads, samples_df, args){
          parseIndexReads(
            barcode.seqs = samples_df$bc2,
            reads = reads,
            indices = idx,
            barcode.length = args$bc2Len,
            max.mismatch = args$bc2Mis,
            max.N.count = args$maxN
          )
        },
        reads = bc2_proc_grps,
        idx = split_bc2_indices,
        MoreArgs = list(
          parseIndexReads = parseIndexReads, 
          samples_df = samples_df,
          args = args
        ),
        SIMPLIFY = FALSE
      )
      
      BC2_parsed <- lapply(
        names(BC2_parsed_list[[1]]), function(x){
          unlist(lapply(seq_along(BC2_parsed_list), function(i){
            BC2_parsed_list[[i]][[x]]
          }))
        }
      )
      
      names(BC2_parsed) <- names(BC2_parsed_list[[1]])
      rm(BC2_parsed_list, bc2_proc_grps)
      
    }
    
    parallel::stopCluster(cluster)
    
  } else{
    
    BC1_parsed <-  parseIndexReads(
      barcode.seqs = samples_df$bc1, 
      reads = bc1_reads,
      indices = all_indices,
      barcode.length = args$bc1Len,
      max.mismatch = args$bc1Mis,
      max.N.count = args$maxN
    )
    
    cat("\nbc1 breakdown:\n")
    print(
      data.frame(
        "bc1" = names(BC1_parsed),
        "Read Counts" = lengths(BC1_parsed)
      ),
      right = TRUE, 
      row.names = FALSE
    )
    
    if( !args$singleBarcode ){
      
      bc2_reads <- r2.chunk
      
      bc2_indices <- stringr::str_extract(
        as.character(ShortRead::id(bc2_reads)), 
        args$readNamePattern
      )
      
      if( !all(bc2_indices == all_indices) ){
        warning(
          "  Index reads are not in the same order. Sequencing files should ",
          "always be kept in order across read types.\n")
      }
      
      BC2_parsed <- parseIndexReads(
        barcode.seqs = samples_df$bc2, 
        reads = bc2_reads,
        indices = bc2_indices,
        barcode.length = args$bc2Len,
        max.mismatch = args$bc2Mis,
        max.N.count = args$maxN
      )
      
    }
    
  }
  
  stop("end of iteration 1")
}
close(f)



  
  
 
  



demulti_func <- function(reads, read.type, cluster, args,
         multiplexed.data, writeDemultiplexedSequences){

  seqs <- reads@sread
  
  ids <- Biostrings::BStringSet(
    stringr::str_extract(
      as.character(reads@id), args$readNamePattern
    )
  )
  
  names(seqs) <- ids
  
  quals <- reads@quality@quality
  
  seqs <- split(
    seqs[match(multiplexed.data$index, as.character(ids))],
    multiplexed.data$sampleName
  )
  
  quals <- split(
    quals[match(multiplexed.data$index, as.character(ids))],
    multiplexed.data$sampleName
  )

  demultiplex <- parallel::clusterMap(
    cluster,
    writeDemultiplexedSequences,
    reads = seqs,
    quals = quals,
    samplename = names(seqs),
    MoreArgs = list(
      type = read.type,
      outfolder = args$outfolder,
      compress = args$compress
    )
  )
  
}

HXB2 <- as.character(ShortRead::sread(ShortRead::readFasta("/home/kevin/dev/chiva/genomes/viral_genomes/HIV_HXB2.fasta")))

r2Reads <- as.character(ShortRead::sread(ShortRead::readFastq("/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3035-4_u3.R2.fastq.gz")))
primers <- substr(r2Reads,1,8)


sampleGenome <- as.character(ShortRead::sread(ShortRead::readFasta("/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3035-4_u3.fasta")))
LTR <- substr(sampleGenome,1,700)

r1Reads <- as.character(ShortRead::sread(ShortRead::readFastq("/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3035-4_u3.R1.fastq.gz")))
r2Reads <- as.character(ShortRead::sread(ShortRead::readFastq("/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3035-4_u3.R2.fastq.gz")))
primerTrimmedr2Reads <- as.character(ShortRead::sread(ShortRead::readFastq("/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3035-4_u3.R2.trim.primer.fastq.gz")))

r1Sample <- sample(r1Reads, 1000)
r2Sample <- sample(r2Reads, 1000)
primerTrimmedr2Sample <- sample(primerTrimmedr2Reads, 1000)


plotDNA(c(r2Sample, rep(LTR, 1000)))
plotDNA(c(r2Sample, rep(revComp(LTR), 1000)))

plotDNA(c(sort(primerTrimmedr2Sample),
          rep(revComp(substr(LTR,1,54)), 1000),
          rep(revComp(substr(sampleGenome,9000,9139)), 1000)))

plotDNA(c(sort(r2Sample),
          rep(revComp(substr(LTR,1,62)), 1000),
          rep(revComp(substr(sampleGenome,9000,9147)), 1000)))

u5r2Reads <- as.character(ShortRead::sread(ShortRead::readFastq("/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3035-4_u5.R2.fastq.gz")))
u5r2Sample <- sample(u5r2Reads, 1000)
plotDNA(names(sort(table(u5r2Sample),decreasing = TRUE)))

old_U3_LTR <- read.csv("/home/kevin/projects/silicano/Silicano20190905.sampleInfo.csv", header = TRUE)[1,"LTR"]

old_U5_LTR <- read.csv("/home/kevin/projects/silicano/Silicano20190905.sampleInfo.csv", header = TRUE)[5,"LTR"]

plotDNA(c(sort(u5r2Sample),
          rep(revComp(old_U3_LTR), 1000),
          rep(substr(sampleGenome,9651,9700), 1000),
          rep(substr(LTR,566,700), 1000)))

trimmedu5r2Reads <- as.character(ShortRead::sread(ShortRead::readFastq("/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3035-4_u5.R2.trim.primer.fastq.gz")))
trimmedu5r2Sample <- sample(trimmedu5r2Reads, 1000)
plotDNA(names(sort(table(trimmedu5r2Sample),decreasing = TRUE)))

plotDNA(c(sort(trimmedu5r2Sample),
          rep(old_U5_LTR,1000),
          rep(substr(sampleGenome,9659,9700), 1000),
          rep(substr(LTR,574,700), 1000)))



get_consensus_sequence <- function(reads, len) {
  reads <- reads[nchar(reads) >= len]
  consensus <- paste0(sapply(seq(1,len), function(i) {
    return(names(sort(table(substr(reads,i,i)),decreasing = TRUE)[1]))
  }), collapse = "")
  return(consensus)
}

u3LTR <- get_consensus_sequence(primerTrimmedr2Sample, 57)
u5LTR <- get_consensus_sequence(trimmedu5r2Sample, 61)

plotDNA(c(sort(primerTrimmedr2Sample),
          rep(u3LTR,1000),
          rep(revComp(substr(LTR,1,54)), 1000),
          rep(revComp(substr(sampleGenome,9000,9139)), 1000)))

plotDNA(c(sort(trimmedu5r2Sample),
          rep(u5LTR,1000),
          rep(substr(sampleGenome,9659,9700), 1000),
          rep(substr(LTR,574,700), 1000)))





ltrOld <- sort(sample(as.character(ShortRead::sread(ShortRead::readFastq("/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3035-4_u3.R2.trim.ltr.fastq.gz"))),1000))
preLtrOld <- sort(sample(as.character(ShortRead::sread(ShortRead::readFastq("/home/kevin/projects/silicano/Silicano20190905/analysis_data/GTSP3035-4_u3.R2.trim.primer.fastq.gz"))),1000))
ltrNew <- sort(sample(as.character(ShortRead::sread(ShortRead::readFastq("/home/kevin/projects/silicano/Silicano20190905_New/analysis_data/GTSP3035-4_u3.R2.trim.ltr.fastq.gz"))),1000))
preLtr <- sort(sample(as.character(ShortRead::sread(ShortRead::readFastq("/home/kevin/projects/silicano/Silicano20190905_New/analysis_data/GTSP3035-4_u3.R2.trim.primer.fastq.gz"))),1000))
plotDNA(c(ltrOld,preLtrOld, ltrNew, preLtr))




