# Initial analysis of patient sequencing by High Throughput Method.
# commandline args : Rscript ....R inputfile output_dir install_dir
options(stringsAsFactors = FALSE, scipen = 99)

args <- commandArgs(trailingOnly = TRUE)
uniqueSitesFile <- args[1]
projectDir <- args[2]
chivaDir <- args[3]

packs <- c("stringr", "GenomicRanges", "igraph", "gintools", "dplyr", 
           "magrittr", "data.table", "Matrix", "BSgenome", "Biostrings",
           "reshape2")
packs_loaded <- suppressMessages(
  sapply(packs, require, character.only = TRUE))
if(!all(packs_loaded)){
  pander::pandoc.table(data.frame(
    "R-Packages" = names(packs_loaded), 
    "Loaded" = packs_loaded, 
    row.names = NULL))
  stop("Check dependancies.")
}

# Source additional functions (HARDCODED FOR NOW!!!)
source(file.path(chivaDir, "tools/rscripts/supporting_funcs.R"))

# Load unique sites data, by reads
uniq_reads <- fread(uniqueSitesFile, sep = ",") %>%
  as.data.frame() %>%
  mutate(
    chr = seqnames, 
    position = ifelse(strand == "+", start, end),
    breakpoint = ifelse(strand == "+", end, start),
    counts = rep(1, n())) %>%
  select(chr, position, breakpoint, strand, sampleName, counts) %>% 
  db_to_granges(split_sampleName = FALSE)

# Typically refinement of breakpoints should be done at the
# individual sample level, this helps increase the abundance range.
# For this study though, abundances are expected to all be very low.
# And therefore we can increase the chances of isolating crossovers
# if we refine the breakpoints of all samples together.

## By each sample independently.
#uniq_sites <- unlist(GRangesList(lapply(
#  split(uniq_reads, uniq_reads$samplename),
#  refine_breakpoints, counts = "counts")))

## All samples refined together
uniq_sites <- refine_breakpoints(uniq_reads, counts = "counts")
uniq_sites <- unique_granges(uniq_sites, counts.col = "counts")

uniq_sites <- standardize_sites(uniq_sites)
uniq_sites <- unique_granges(uniq_sites, counts.col = "counts")
uniq_sites$posid <- generate_posid(uniq_sites)

cond_sites <- condense_intsites(
  uniq_sites, grouping = "samplename", return.abundance = TRUE)

xofil_uniq_sites <- filter_crossovers(uniq_sites, "samplename", "counts")
xofil_cond_sites <- condense_intsites(
  xofil_uniq_sites, grouping = "samplename", return.abundance = TRUE
)

summary_tbl <- as.data.frame(xofil_cond_sites) %>% 
  mutate(sample = str_extract(samplename, "[\\w]+")) %>% 
  group_by(sample) %>% 
  summarise(
    reads = sum(counts), cells = sum(estAbund), sites = n_distinct(posid)
  )

# Summarize
saveRDS(
  list("unfil_uniq_sites" = uniq_sites, "fil_uniq_sites" = xofil_uniq_sites),
  file.path(projectDir, "standardized_uniq_sites.rds"))

write.csv(
  as.data.frame(cond_sites),
  file.path(projectDir, "condensed_sites.csv"),
  quote = FALSE, row.names = FALSE)

write.csv(
  as.data.frame(xofil_cond_sites),
  file.path(projectDir, "xofil_condensed_sites.csv"),
  quote = FALSE, row.names = FALSE)

read_site_matrix <- as.data.frame(uniq_sites) %>% 
  select(posid, samplename, counts) %>% 
  group_by(posid, samplename) %>% 
  summarise(counts = sum(counts)) %>% 
  melt(id.vars = c("posid", "samplename")) %>%
  acast(posid ~ samplename, value.var = "value", fill = 0)

write.csv(
  read_site_matrix, 
  file.path(projectDir, "read_site_matrix.csv"), 
  quote = FALSE, row.names = TRUE)

frag_site_matrix <- as.data.frame(cond_sites) %>% 
  select(posid, samplename, estAbund) %>% 
  group_by(posid, samplename) %>% 
  summarise(estAbund = sum(estAbund)) %>% 
  melt(id.vars = c("posid", "samplename")) %>%
  acast(posid ~ samplename, value.var = "value", fill = 0)

write.csv(
  frag_site_matrix, 
  file.path(projectDir, "/fragment_site_matrix.csv"), 
  quote = FALSE, row.names = TRUE)

write.csv(
  summary_tbl,
  file.path(projectDir, "/summary_table.csv"),
  quote = FALSE, row.names = FALSE)

q()
