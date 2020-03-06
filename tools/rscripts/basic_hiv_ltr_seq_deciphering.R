# This script is a short workflow to identify sequences for recent Silicano
# samples that will be used for filtering amplicons. The idea is to obtain the
# sequences from the provided data and if any gaps or incomplete sequeuences
# can be filled in with remaining HXB2 reference sequence.
# 
# Example by Christopher Nobles, Ph.D.


# Load options and packages ----
options(stringsAsFactors = FALSE, width = 120)
pkgs <- c("GenomicRanges", "Biostrings", "DECIPHER", "tidyverse")
loaded <- suppressMessages(sapply(pkgs, require, character.only = TRUE))
stopifnot(all(loaded))


# Read in relavent data ----
# This includes sequences provided, reference genonme (HXB2 for HIV-1, 
# subgroup B), and any additional informtion that would be needed to isolate 
# related sequences.

provided_seqs <- readDNAStringSet(
  "/home/nobles/Lab_Records/project_hiv_silicianoLab/viralSequences/LTR sequences for RB 101118.fasta"
)

hxb2 <- readDNAStringSet("/home/kevin/dev/chiva/genomes/viral_genomes/HIV_HXB2.fasta")

# I wrote this down a long time ago, from reference material. Useful.
hxb2_structure <- data.frame(
  feature = c(
    "5'LTR-U3", "5'LTR-R", "5'LTR-U5", "Lys_tRNA_PBS", "Packing_Loops",
    "Gag", "Pol", "Env", "3'LTR-U3", "3'LTR-R", "3'LTR-U5"),
  start = c(1, 456, 552, 635, 681, 790, 2085, 6225, 9086, 9541, 9637),
  end = c(455, 551, 634, 653, 789, 2292, 5096, 8795, 9540, 9636, 9719),
  frame = c(1, 1, 1, 1, 1, 1, 3, 3, 2, 2, 2)
)

hxb2_structure_gr <- GRanges(
  seqnames = "K03455|HIVHXB2CG",
  ranges = IRanges(start = hxb2_structure$start, end = hxb2_structure$end),
  strand = "+"
)

mcols(hxb2_structure_gr) <- hxb2_structure[,c("feature", "frame")]


# Isolate relative sequences and perform MSA ----
# Remove existing HXB2 references from given sequences. If you just take a look
# at the provided MSA, there are two different reference lengths. They likely 
# concatenated two MSA's, which is not helpful for our situation. We need one 
# LTR sequence, then we will pull the consensus of the last 50-60 bases from
# each side in a patient specific manner.
given_seqs <- provided_seqs[!grepl("hxb2", tolower(names(provided_seqs)))]

# Remove empty sequence indicators ("-") and return to a DNAStringSet
just_seqs <- data.frame(
    name = names(given_seqs), seqs = as.character(given_seqs)
  ) %>%
  mutate(seqs = str_remove_all(seqs, "-")) %>%
  deframe() %>%
  DNAStringSet()

# Add in the HXB2 LTR sequences, both 5' and 3'
ltr_gr <- hxb2_structure_gr[grepl("LTR", hxb2_structure_gr$feature)] %>%
  GenomicRanges::reduce()

ltr_seqs <- do.call(c, lapply(ltr_gr, function(x){
  DNAStringSet(hxb2, start = start(x), end = end(x))
}))

names(ltr_seqs) <- c("5'LTR-HXB2", "3'LTR-HXB2")

# Combine sequences and perform MSA using the DECIPHER package (or your 
# favorite)
msa_input <- c(ltr_seqs, just_seqs)
msa_seqs <- DECIPHER::AlignSeqs(msa_input)
  
# The results from the MSA seem a bit off, but it seems like it's only because
# of one of the sequences, "FSJJ384_T2_nefU3_memory_G12", but the sequences tend
# to agree with others identified. So why don't we hide that sequence for the
# msa. What I'm looking for here is to make sure the reference sequence isn't 
# broken in the primer to junction ranges (first and last 50 - 60 nts).
msa_input <- msa_input[names(msa_input) != "FSJJ384_T2_nefU3_memory_G12"]
msa_seqs <- DECIPHER::AlignSeqs(msa_input)

# Identify the ranges where the sequences match the LTRs and trim away all other
# sequences. These other sequences refer to internal vector sequence such as gag
# or nef.
ref_chars <- unlist(strsplit(as.character(msa_seqs[1]), ""))
first_nt <- match("T", ref_chars) #revcomp of CA would mean T comes first
last_nt <- length(ref_chars) - match("A", rev(ref_chars)) + 1
ltr_msa <- subseq(msa_seqs, start = first_nt, end = last_nt)

# Now with the clean MSA, we can select the regions we are amplifying, below are
# information on which primers were used (Hayley would know this too).
# Additionally, I always marked how far away the primers would be from the
# junctions given the HXB2 reference.

## PCR2 U3 Primer : CP07U3
u3_dist <- 59L

## PCR2 RU5 Primer : FSRU5
u5_dist <- 62L

# Just check to make sure you have clean reference alignment in the range of the
# primers
stopifnot(
  !grepl("-", subseq(ltr_msa[1], end = u3_dist)) &
  !grepl("-", subseq(ltr_msa[1], start = width(ltr_msa[1]) - u5_dist + 1))
)

# Isolate the regions of interest from the sequences
U3_seqs <- subseq(ltr_msa, end = u3_dist)
U5_seqs <- subseq(ltr_msa, start = unique(width(ltr_msa)) - u5_dist + 1)

# It looks like there is plenty of sequence to fill in the distance without the
# reference so let's remove the HXB2 reference sequences.
U3_seqs <- U3_seqs[!grepl("HXB2", names(U3_seqs))]
U5_seqs <- U5_seqs[!grepl("HXB2", names(U5_seqs))]

# !!!! It's important to remember that we've been working with U3_seqs in the 
# reverse complement orientation. So we should flip that around before we
# identify the final sequences to be used.
U3_seqs <- reverseComplement(U3_seqs)

# To me, it looks like the sequences only came from a single patient, "FSJJ384",
# so we should make a consensus sequence for each that captures all sequences 
# found. Here I set a really low threshold to get sequence variability.
U3_con_mat <- consensusMatrix(U3_seqs)[c("A", "C", "G", "T"),]
U5_con_mat <- consensusMatrix(U5_seqs)[c("A", "C", "G", "T"),]

U3_consensus <- names(IUPAC_CODE_MAP)[
  match(
    apply(U3_con_mat, 2, function(x){
      paste(c("A", "C", "G", "T")[x > 0], collapse = "")
    }),
    IUPAC_CODE_MAP
  )] %>%
  paste(collapse = "")

U5_consensus <- names(IUPAC_CODE_MAP)[
  match(
    apply(U5_con_mat, 2, function(x){
      paste(c("A", "C", "G", "T")[x > 0], collapse = "")
    }),
    IUPAC_CODE_MAP
  )] %>%
  paste(collapse = "")

# These sequences include the "CA" at the end, even if the consensus string
# mannaged to move in ambiguious bases. So for parameters into cHIVa, we would
# want to remove the terminal "CA" so that we can enforce it during processing.

U3_consensus <- substr(U3_consensus, start = 1, stop = nchar(U3_consensus) - 2)
U5_consensus <- substr(U5_consensus, start = 1, stop = nchar(U5_consensus) - 2)

# These last two sequences are what you will want to use in the LTR column of
# the sampleInfo.csv file.


