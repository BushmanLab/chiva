# This script is meant to concatenate sequence files within a directory together
# into a single file. 
# 
# Usage Rscript concat_seq_files.R {input.directory} {output.fasta}
# 
# Files within the input directory should be in fasta or fasta.gz formats.
# 

library(tidyverse)


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
  "--fasta1", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--fasta2", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--blast1", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--blast2", nargs = 1, type = "character", help = desc$output
)

#parser$add_argument(
#  "--output", nargs = 1, type = "character", help = desc$output
#)

parser$add_argument(
  "--summary", nargs = 1, type = "character", help = desc$output
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

R1_reads <- ShortRead::readFasta(args$fasta1)
R1_reads_df <- data.frame(rawseq=R1_reads@sread,
                          seqid=R1_reads@id)

R2_reads <- ShortRead::readFasta(args$fasta2)
R2_reads_df <- data.frame(rawseq=R2_reads@sread,
                          seqid=R2_reads@id)


blast_output_raw_R1 <- read.table(args$blast1, header=FALSE, skip=0, comment.char = "", sep = "\t") %>%
  mutate(read="R1")
blast_output_raw_R2 <- read.table(args$blast2, header=FALSE, skip=0, comment.char = "", sep = "\t") %>%
  mutate(read="R2")
blast_colnames <-  c("qseqid",
                     "qseq",
                     "qlen",
                     "qstart",
                     "qend",
                     "sseqid",
                     "sseq",
                     "sstrand",
                     "slen",
                     "sstart",
                     "send",
                     "evalue",
                     "bitscore",
                     "score",
                     "length",
                     "pident",
                     "nident",
                     "mismatch",
                     "btop",
                     "read")
colnames(blast_output_raw_R1) <- blast_colnames
colnames(blast_output_raw_R2) <- blast_colnames

blast_output_raw_R1 <- blast_output_raw_R1 %>%
  group_by(qseqid) %>%
  mutate(max_bitscore=max(bitscore)) %>%
  ungroup()

blast_output_raw_R2 <- blast_output_raw_R2 %>%
  group_by(qseqid) %>%
  mutate(max_bitscore=max(bitscore)) %>%
  ungroup()

blast_output_paired <- blast_output_raw_R1 %>%
  full_join(blast_output_raw_R2,
            by=c("qseqid","sseqid"), suffix=c(".R1",".R2")) %>%
  group_by(qseqid) %>%
  # FIX THIS CODE
  # The max needs to be done to each R1 and R2 before the join
  # Otherwise, max can be set at NA
  mutate(best=max_bitscore.R1==bitscore.R1 & max_bitscore.R2==bitscore.R2) %>%
  filter(best) %>%
  mutate(best_count=n()) %>%
  ungroup() %>%
  select(-best) %>%
  left_join(R1_reads_df, by=c("qseqid"="seqid")) %>%
  rename(rawseq.R1=rawseq) %>%
  left_join(R2_reads_df, by=c("qseqid"="seqid")) %>%
  rename(rawseq.R2=rawseq)

good_paired_qseqids <- unique(blast_output_paired$qseqid)

blast_output_filtered_R1 <- blast_output_raw_R1 %>%
  filter(!qseqid %in% good_paired_qseqids) %>%
  mutate(best=max_bitscore==bitscore) %>%
  filter(best) %>%
  left_join(R1_reads_df, by=c("qseqid"="seqid"))

blast_output_filtered_R2 <- blast_output_raw_R2 %>%
  filter(!qseqid %in% good_paired_qseqids) %>%
  mutate(best=max_bitscore==bitscore) %>%
  filter(best)%>%
  left_join(R2_reads_df, by=c("qseqid"="seqid"))

blast_output_paired_unique_seqs <- blast_output_paired %>%
  group_by(sseqid, rawseq.R1, rawseq.R2) %>%
  summarize(count=n())

blast_output_R1_unique_seqs <- blast_output_filtered_R1 %>%
  group_by(sseqid, rawseq) %>%
  summarize(count=n())

blast_output_R2_unique_seqs <- blast_output_filtered_R2 %>%
  group_by(sseqid, rawseq) %>%
  summarize(count=n())



blast_summary_paired <- blast_output_paired_unique_seqs %>%
  group_by(sseqid) %>%
  summarize(uniq_seqs=n(), total_seqs=sum(count)) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(read="paired")

blast_summary_R1 <- blast_output_R1_unique_seqs %>%
  group_by(sseqid) %>%
  summarize(uniq_seqs=n(), total_seqs=sum(count)) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(read="R1")

blast_summary_R2 <- blast_output_R2_unique_seqs %>%
  group_by(sseqid) %>%
  summarize(uniq_seqs=n(), total_seqs=sum(count)) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(read="R2")

blast_summary <- rbind(blast_summary_paired, blast_summary_R1, blast_summary_R2)

paired_table <- gsub("summary.csv","paired.filtered.tbl",args$summary)
R1_table <- gsub("summary.csv","R1.filtered.tbl",args$summary)
R2_table <- gsub("summary.csv","R2.filtered.tbl",args$summary)
write.table(blast_output_paired, file=paired_table, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(blast_output_filtered_R1, file=R1_table, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(blast_output_filtered_R2, file=R2_table, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(blast_summary, file=args$summary, append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

q()

library(ggplot2)


blast_results %>%
  arrange(-uniq_seqs) %>%
  ggplot(aes(x=SampleID, y=sseqid, fill=uniq_seqs)) +
  geom_tile() +
  facet_grid(Patient ~ read, scales = "free") +
  theme (
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

blast_results %>%
  arrange(-total_seqs) %>%
  ggplot(aes(x=SampleID, y=sseqid, fill=total_seqs)) +
  geom_tile() +
  facet_grid(Patient ~ read, scales = "free") +
  theme (
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

blast_results %>%
  group_by(SampleID) %>%
  mutate(prop=uniq_seqs/sum(uniq_seqs)) %>%
  arrange(-prop) %>%
  ggplot(aes(x=SampleID, y=sseqid, fill=prop)) +
  geom_tile() +
  facet_grid(Patient ~ read, scales = "free") +
  theme (
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

blast_results %>%
  group_by(SampleID) %>%
  mutate(prop=total_seqs/sum(total_seqs)) %>%
  arrange(-prop) %>%
  ggplot(aes(x=SampleID, y=sseqid, fill=prop)) +
  geom_tile() +
  facet_grid(Patient ~ read, scales = "free") +
  theme (
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

blast_results %>%
  filter(Patient=="22") %>%
  arrange(-uniq_seqs) %>%
  ggplot(aes(x=SampleID, y=sseqid, fill=uniq_seqs)) +
  geom_tile() +
  facet_grid(Patient ~ read, scales = "free") +
  theme(
    axis.text.y = element_text(size=4)
  )


args$read1 <- "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U3-2.R1.blast.tbl"
args$read2 <- "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U3-2.R2.blast.tbl"


args$blast1 <- "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U5-2.R1.blast.tbl"
args$blast2 <- "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U5-2.R2.blast.tbl"
args$fasta1 <- "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U5-2.R1.blast.fasta"
args$fasta2 <- "/home/kevin/projects/bar_hiv/GTSP2722/analysis_data/GTSP2722-U5-2.R2.blast.fasta"


R1_reads <- ShortRead::readFasta(args$fasta1)
R1_reads_df <- data.frame(rawseq=R1_reads@sread,
                          seqid=R1_reads@id)

R2_reads <- ShortRead::readFasta(args$fasta2)
R2_reads_df <- data.frame(rawseq=R2_reads@sread,
                          seqid=R2_reads@id)


blast_output_raw_R1 <- read.table(args$blast1, header=FALSE, skip=0, comment.char = "", sep = "\t") %>%
  mutate(read="R1")
blast_output_raw_R2 <- read.table(args$blast2, header=FALSE, skip=0, comment.char = "", sep = "\t") %>%
  mutate(read="R2")
blast_colnames <-  c("qseqid",
                     "qseq",
                     "qlen",
                     "qstart",
                     "qend",
                     "sseqid",
                     "sseq",
                     "sstrand",
                     "slen",
                     "sstart",
                     "send",
                     "evalue",
                     "bitscore",
                     "score",
                     "length",
                     "pident",
                     "bitscore",
                     "mismatch",
                     "btop",
                     "read")
colnames(blast_output_raw_R1) <- blast_colnames
colnames(blast_output_raw_R2) <- blast_colnames

blast_output_paired_all <- blast_output_raw_R1 %>%
  full_join(blast_output_raw_R2,
            by=c("qseqid","sseqid"), suffix=c(".R1",".R2")) %>%
  group_by(qseqid) %>%
  mutate(max_bitscore.R1=max(bitscore.R1),
         max_bitscore.R2=max(bitscore.R2),
         best=max_bitscore.R1==bitscore.R1 & max_bitscore.R2==bitscore.R2,
         best_count=sum(best)) %>%
  ungroup()
blast_output_paired <- blast_output_paired_all %>%
  filter(best) %>%
  select(-best) %>%
  left_join(R1_reads_df, by=c("qseqid"="seqid")) %>%
  rename(rawseq.R1=rawseq) %>%
  left_join(R2_reads_df, by=c("qseqid"="seqid")) %>%
  rename(rawseq.R2=rawseq)

good_paired_qseqids <- unique(blast_output_paired$qseqid)

blast_output_filtered_R1 <- blast_output_raw_R1 %>%
  group_by(qseqid) %>%
  filter(!qseqid %in% good_paired_qseqids) %>%
  mutate(max_bitscore=max(bitscore),
         best=max_bitscore==bitscore,
         best_count=sum(best)) %>%
  filter(best) %>%
  left_join(R1_reads_df, by=c("qseqid"="seqid")) %>%
  ungroup()

blast_output_filtered_R2 <- blast_output_raw_R2 %>%
  group_by(qseqid) %>%
  filter(!qseqid %in% good_paired_qseqids) %>%
  mutate(max_bitscore=max(bitscore),
         best=max_bitscore==bitscore,
         best_count=sum(best)) %>%
  filter(best)%>%
  left_join(R2_reads_df, by=c("qseqid"="seqid")) %>%
  ungroup()

blast_output_paired_unique_seqs <- blast_output_paired %>%
  group_by(sseqid, rawseq.R1, rawseq.R2) %>%
  summarize(count=n())

blast_output_R1_unique_seqs <- blast_output_filtered_R1 %>%
  group_by(sseqid, rawseq) %>%
  summarize(count=n())

blast_output_R2_unique_seqs <- blast_output_filtered_R2 %>%
  group_by(sseqid, rawseq) %>%
  summarize(count=n())



blast_summary_paired <- blast_output_paired_unique_seqs %>%
  group_by(sseqid) %>%
  summarize(uniq_seqs=n(), total_seqs=sum(count)) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(read="paired")

blast_summary_R1 <- blast_output_R1_unique_seqs %>%
  group_by(sseqid) %>%
  summarize(uniq_seqs=n(), total_seqs=sum(count)) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(read="R1")

blast_summary_R2 <- blast_output_R2_unique_seqs %>%
  group_by(sseqid) %>%
  summarize(uniq_seqs=n(), total_seqs=sum(count)) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(read="R2")

blast_summary <- rbind(blast_summary_paired, blast_summary_R1, blast_summary_R2)







blast_output_filtered_R1 %>%
  group_by(best_count) %>%
  summarize(count=n()) %>%
  ungroup() %>%
  ggplot(aes(x=best_count, y=count)) +
  geom_bar(stat = "identity") +
  scale_y_log10()

blast_output_filtered_R2 %>%
  group_by(best_count) %>%
  summarize(count=n()) %>%
  ungroup() %>%
  ggplot(aes(x=best_count, y=count)) +
  geom_bar(stat = "identity") +
  scale_y_log10()

blast_output_paired %>%
  group_by(best_count) %>%
  summarize(count=n()) %>%
  ungroup() %>%
  ggplot(aes(x=best_count, y=count)) +
  geom_bar(stat = "identity") +
  scale_y_log10()






