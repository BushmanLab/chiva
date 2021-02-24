# Initial analysis of patient sequencing by High Throughput Method.
# commandline args : Rscript ....R inputfile output_dir install_dir
options(stringsAsFactors = FALSE, scipen = 99)

args <- commandArgs(trailingOnly = TRUE)

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
source(file.path(args[4], "tools/rscripts/supporting_funcs.R"))

# Load unique sites data, by reads
uniq_reads <- fread(args[1], sep = ",") %>%
  as.data.frame() %>%
  mutate(
    chr = seqnames, 
    position = ifelse(strand == "+", start, end),
    breakpoint = ifelse(strand == "+", end, start),
    counts = rep(1, n())) %>%
  select(chr, position, breakpoint, strand, sampleName, counts) %>% 
  db_to_granges(split_sampleName = FALSE)

# Refinement of breakpoints should be done at the
# individual sample level which helps increase the abundance range.

# uniq_sites <- unlist(GRangesList(lapply(
#  split(uniq_reads, uniq_reads$samplename),
#  refine_breakpoints, counts = "counts")))
# uniq_sites <- unique_granges(uniq_sites, sum.cols = "counts")

# Standardization of sites should be done at the patient level.

# supp_data <- read.table(args[5], sep = ",", header = TRUE)
# uniq_sites$GTSP <- substr(uniq_sites$samplename, 1, 8)
# uniq_sites$patient <- supp_data$patient[match(uniq_sites$GTSP, supp_data$specimen)]
# uniq_sites <- unlist(GRangesList(lapply(
#   split(uniq_sites, uniq_sites$patient),
#   standardize_sites, counts = "counts")))
# uniq_sites <- unique_granges(uniq_sites, sum.cols = "counts")
# uniq_sites$GTSP <- NULL
# uniq_sites$patient <- NULL


# For studies where abundances are expected to be low, we can increase
# the chances of isolating crossovers if we refine the breakpoints of all samples together.

uniq_sites <- refine_breakpoints(uniq_reads, counts = "counts")
uniq_sites <- unique_granges(uniq_sites, sum.cols = "counts")
uniq_sites <- standardize_sites(uniq_sites, counts = "counts")
uniq_sites <- unique_granges(uniq_sites, sum.cols = "counts")



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
  file.path(args[2], paste0("standardized_uniq_sites.",args[3],".rds")))

write.csv(
  as.data.frame(cond_sites),
  file.path(args[2], paste0("condensed_sites.",args[3],".csv")),
  quote = FALSE, row.names = FALSE)

write.csv(
  as.data.frame(xofil_cond_sites),
  file.path(args[2], paste0("xofil_condensed_sites.",args[3],".csv")),
  quote = FALSE, row.names = FALSE)

read_site_matrix <- as.data.frame(uniq_sites) %>% 
  select(posid, samplename, counts) %>% 
  group_by(posid, samplename) %>% 
  summarise(counts = sum(counts)) %>% 
  melt(id.vars = c("posid", "samplename")) %>%
  acast(posid ~ samplename, value.var = "value", fill = 0)

write.csv(
  read_site_matrix, 
  file.path(args[2], paste0("read_site_matrix.",args[3],".csv")),
  quote = FALSE, row.names = TRUE)

frag_site_matrix <- as.data.frame(cond_sites) %>% 
  select(posid, samplename, estAbund) %>% 
  group_by(posid, samplename) %>% 
  summarise(estAbund = sum(estAbund)) %>% 
  melt(id.vars = c("posid", "samplename")) %>%
  acast(posid ~ samplename, value.var = "value", fill = 0)

write.csv(
  frag_site_matrix, 
  file.path(args[2], paste0("fragment_site_matrix.",args[3],".csv")),
  quote = FALSE, row.names = TRUE)

write.csv(
  summary_tbl,
  file.path(args[2], paste0("summary_table.",args[3],".csv")),
  quote = FALSE, row.names = FALSE)

q()