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

write.csv(result, file = args$output,
          quote = FALSE, row.names = FALSE, col.names = TRUE)

q()

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

readLength <- 
  

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
  mutate(max_nident.R1=max(nident.R1),
         max_nident.R2=max(nident.R2),
         best=if_else(paired,max_nident.R1==nident.R1 & max_nident.R2==nident.R2,
                      if_else(is.na(qseq.R1), max_nident.R2==nident.R2, max_nident.R1==nident.R1)),
         best_count=sum(best)) %>%
  ungroup()



q()



