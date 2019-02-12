# Notes
# Include a summary table at the beginning
# Barplot distributions, facet_wrap, split by patient, group celltype, order by
# timepoint. Top tables from each patient in sample set, includes sites sorted
# by abundance and read counts.


#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

options(stringsAsFactors = FALSE, scipen = 99, width = 999)
set.seed(1)


# Set up and gather command line arguments ----
parser <- argparse::ArgumentParser(
  description = "Generate a summary report for input run(s)."
)

parser$add_argument(
  "input", nargs = "+", type = "character",
  help = "Post-processing output .rds file(s). ie. standardized_uniq_sites.rds."
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", 
  help = "Output report file, extension not required."
)

parser$add_argument(
  "-c", "--config", nargs = "+", type = "character",
  help = "Run specific config file(s) in yaml format."
)

parser$add_argument(
  "-s", "--support", nargs = 1, type = "character",
  help = "Supplementary data input, csv or tsv format."
)

parser$add_argument(
  "-f", "--figures", action = "store_true",
  help = "Generate figures along with output report (pdf and png formats)."
)

parser$add_argument(
  "-d", "--data", action = "store_true",
  help = "Data to generate the report will be saved as an R image with output."
)

parser$add_argument(
  "-t", "--format", nargs = 1, type = "character", default = "html",
  help = "Output format for report. Either 'pdf' or 'html' (default)."
)

parser$add_argument(
  "--template", nargs = 1, type = "character", 
  default = "tools/rscripts/report_templates/HIV_report_template.Rmd",
  help = "File path to standard or custom report template."
)

parser$add_argument(
  "--install_path", nargs = 1, type = "character", default = "CHIVA_DIR",
  help = "Install directory path, do not change for normal applications."
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if( !dir.exists(args$install_path) ){
  root_dir <- Sys.getenv(args$install_path)
}else{
  root_dir <- args$install_path
}

if( !dir.exists(root_dir) ){
  stop(paste0("\n  Cannot find install path: ", args$install_path, ".\n"))
}else{
  args$install_path <- root_dir
}

if( length(args$input) != length(args$config) ){
  stop("\n  Must supply one config file for each input.")
}

report_formats <- c("html" = "html_document", "pdf" = "pdf_document")

if( !args$format %in% names(report_formats) ){
  stop("Please input either 'html' or 'pdf' for format.\n",
       "Other formats not supported.")
}

output_format <- report_formats[args$format]

## Resolve template file path.
if( file.exists(file.path(root_dir, args$template)) ){
  
  template_path <- normalizePath(file.path(root_dir, args$template))
  
}else if( file.exists(file.path(args$template)) ){
  
  template_path <- normalizePath(file.path(args$template))
  
}else{
  
  stop("\nCannot find template file: ", args$template, ".\n")
  
}


## Construct input table and print to terminal
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(seq_along(args), function(i){
    paste(args[[i]], collapse = ", ")
  })
)

input_table <- input_table[
  match(
    c("input :", "output :", "config :", "support :", 
      "figures :", "data :", "format :", "template :", 
      "install_path :"),
    input_table$Variables),
]

cat("\nReport Inputs:\n")

print(
  data.frame(input_table),
  right = FALSE, 
  row.names = FALSE
)


# Load dependancies ----
cat("\nLoading dependencies.\n")

add_packs <- c("hiAnnotator", "dplyr", "magrittr", "knitr")

add_packs_loaded <- suppressMessages(
  sapply(add_packs, require, character.only = TRUE)
)

if( !all(add_packs_loaded) ){
  
  print(
    data.frame(
      "R-Packages" = names(add_packs_loaded), 
      "Loaded" = add_packs_loaded
    ), 
    right = FALSE,
    row.names = FALSE
  )
  
  stop("Check dependancies.\n")
  
}

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

source(file.path(code_dir, "supporting_funcs.R"))

# Import metadata and consolidate into report objects ----
cat("Importing experimental data and configurations.\n\n")

## Load config files
configs <- lapply(args$config, function(x){
  if( file.exists(file.path(root_dir, x)) ){
    return(yaml::yaml.load_file(file.path(root_dir, x)))
  }else if( file.exists(x) ){
    return(yaml::yaml.load_file(x))
  }else{
    stop("\n  Cannot find config file: ", x, ".\n")
  }
})

names(configs) <- sapply(configs, "[[", "Run_Name")

## Load reference genome 
if( grepl(".fa", unique(sapply(configs, "[[", "Ref_Genome"))) ){
  
  if( !(
    file.exists(
      file.path(root_dir, unique(sapply(configs, "[[", "Ref_Genome")))
    ) | file.exists(unique(sapply(configs, "[[", "Ref_Genome")))
  ) ){
    stop("Specified reference genome file not found.")
  }
  
  ref_file_type <- ifelse(
    grepl(".fastq", unique(sapply(configs, "[[", "Ref_Genome"))), 
    "fastq", 
    "fasta"
  )
  
  if( file.exists(
    file.path(root_dir, unique(sapply(configs, "[[", "Ref_Genome"))) 
    ) ){

    ref_genome <- Biostrings::readDNAStringSet(
      filepath = file.path(
        root_dir, unique(sapply(configs, "[[", "Ref_Genome"))
      ),
      format = ref_file_type
    )
    
  }else{
    
    ref_genome <- Biostrings::readDNAStringSet(
      filepath = unique(sapply(configs, "[[", "Ref_Genome")), 
      format = ref_file_type
    )
  }
  
}else{
  
  ref_genome <- unique(sapply(configs, "[[", "Ref_Genome"))
  
  genome <- grep(
    pattern = ref_genome, 
    x = unique(BSgenome::installed.genomes()), 
    value = TRUE
  )
  
  if( length(genome) == 0 ){
    
    cat("\nInstalled genomes include:")
    print(unique(BSgenome::installed.genomes()))
    cat("\nSelected reference genome not in list.")
    stop("Error: Genome not available.")
    
  }else if( length(genome) > 1 ){
    
    cat("\nInstalled genomes include:")
    print(unique(BSgenome::installed.genomes()))
    cat(
      "\nPlease be more specific about reference genome.",
      "Multiple matches to input."
    )
    stop("Error: Multiple genomes requested.")
    
  }
  
  suppressMessages(library(genome, character.only = TRUE))
  
  ref_genome <- get(genome)
  
}


## Get versioning

soft_version <- as.character(read.delim(
  file = file.path(root_dir, ".version"), header = FALSE
))

build_version <- list.files(file.path(root_dir, "etc")) %>%
  grep(pattern = "build", x = ., value = TRUE) %>%
  stringr::str_extract(pattern = "v[0-9]+\\.[0-9]+.[0-9]+")

signature <- paste(
  unique(sort(unlist(lapply(configs, "[[", "signature")))), collapse = ", "
)
  
  
## Load reference files
ref_genes <- loadRefFiles(
  configs[[1]]$refGenes, 
  type = "GRanges", 
  freeze = configs[[1]]$Ref_Genome,
  root = root_dir
)

ref_exons <- extractExons(
  ref.genes = ref_genes, 
  exon.starts = ref_genes$exonStarts, 
  exon.ends = ref_genes$exonEnds, 
  gene.names = ref_genes$annot_sym
)

onco_genes <- loadRefFiles(
  configs[[1]]$oncoGeneList, 
  type = "gene.list", 
  freeze = configs[[1]]$Ref_Genome,
  root = root_dir
)

special_genes <- loadRefFiles(
  configs[[1]]$specialGeneList, 
  type = "gene.list", 
  freeze = config[[1]]$Ref_Genome,
  root = root_dir
)


## Combine sampleInfo files

sample_info <- dplyr::bind_rows(lapply(
      sapply(configs, "[[", "Sample_Info"), 
    function(x){
      
      if( file.exists(file.path(root_dir, x)) ){
        return(data.table::fread(file.path(root_dir, x), data.table = FALSE))
      }else if( file.exists(x) ){
        return(data.table::fread(x, data.table = FALSE))
      }else{
        stop("\n  Cannot find Sample_Info: ", x, ".\n")
      }

    }), .id = "run_set")

sample_name_col <- unique(sapply(configs, "[[", "Sample_Name_Column"))

if( length(sample_name_col) != 1 ){
  stop(
    "\n  SampleInfo files not in same format. Please revise for consistency.\n"
  )
}

sample_info$specimen <- stringr::str_extract(
  string = sample_info[,sample_name_col], 
  pattern = "[\\w]+"
)

specimen_levels <- unique(sample_info$specimen)

sample_info$specimen <- factor(sample_info$specimen, levels = specimen_levels)

required_sample_info_cols <- c(sample_name_col, "uniqueRegion")

if( !all(required_sample_info_cols %in% names(sample_info)) ){
  stop(
    "\n  Sample info files need to contain the following columns: ", 
    paste(required_sample_info_cols, collapse = ", "),
    ".\n"
  )
}

sample_info <- dplyr::rename(sample_info, "samplename" = sample_name_col)

## Load in supporting information ----
if( length(args$support) > 0 ){
  
  if( file.exists(file.path(root_dir, args$support)) ){
    support_path <- file.path(root_dir, args$support)
  }else if( file.exists(args$support) ){
    support_path <- args$support
  }else{
    stop("\n  Cannot find supporting data file: ", args$support, ".\n")
  }
  
  supp_data <- data.table::fread(support_path, data.table = FALSE)
  specimen_levels <- supp_data$specimen[supp_data$specimen %in% specimen_levels]
  patient_levels <- unique(supp_data$patient)
  
  supp_data <- dplyr::filter(supp_data, specimen %in% specimen_levels) %>%
    dplyr::mutate(
      specimen = factor(specimen, levels = specimen_levels),
      patient = factor(patient, levels = patient_levels)
    )
  
  sample_info <- dplyr::filter(sample_info, specimen %in% specimen_levels) %>%
    dplyr::mutate(
      specimen = factor(as.character(specimen), levels = specimen_levels)
    ) %>%
    dplyr::arrange(specimen)
  
}else{
  
  supp_data <- data.frame("specimen" = c(), "patient" = c())
  
}


## Consolidate supplementary data ----
if( length(args$support) > 0 ){
  
  cond_overview <- supp_data %>%
    dplyr::mutate(
      condition = vcollapse(
        d = dplyr::select(supp_data, -specimen, -patient), 
        sep = " - ", 
        fill = "NA"
      ),
      condition = factor(condition, levels = c(unique(condition), "Mock"))
    ) %>%
    dplyr::select(specimen, patient, condition)

}else{
  
  cond_overview <- sample_info %>%
    dplyr::mutate(
      condition = NA,
      patient = rep(NA, n())
    ) %>%
    dplyr::select(specimen, patient, condition)

}

## Read in experimental data and contatenate different sets ----
input_data <- lapply(args$input, function(x){
  if( file.exists(file.path(root_dir, x)) ){
    return(readRDS(file.path(root_dir, x)))
  }else if( file.exists(x) ){
    return(readRDS(x))
  }else{
    stop("\n  Cannot find input file: ", x, ".\n")
  }
})

names(input_data) <- names(configs)
data_type_names <- names(input_data[[1]])

input_data <- lapply(
  seq_along(input_data[[1]]), 
  function(i){
    dplyr::bind_rows(lapply(input_data, function(x){
        as.data.frame(x[[i]], row.names = NULL)
      }), 
      .id = "run.set"
    )
  }
)

names(input_data) <- data_type_names

input_data <- lapply(
  input_data, 
  function(x) x[x$samplename %in% sample_info$samplename,] 
)


# Analysis ----

cat("Starting analysis...\n")

# Determine dual detected sites ----
fil_sites <- input_data$fil_uniq_sites %>%
  dplyr::filter(seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))) %>%
  dplyr::left_join(
    dplyr::select(sample_info, samplename, specimen, uniqueRegion), 
    by = "samplename"
  ) %>%
  dplyr::mutate(
    specimen = stringr::str_extract(samplename, "[\\w]+"),
    specimen = factor(specimen, levels = specimen_levels),
    replicate = as.numeric(stringr::str_extract(samplename, "[0-9]+$"))
    #replicate = replicate - 4 * (ceiling(replicate/4) - 1) + 4 * run
  ) %>%
  dplyr::left_join(cond_overview, by = "specimen") %>%
  dplyr::mutate(
    patient = factor(patient, levels = unique(cond_overview$patient)),
    specimen = factor(specimen, levels = levels(cond_overview$specimen)),
    univ.id = generate_univID(posid, uniqueRegion),
    gene.id = assign_gene_id(
      seqnames = seqnames, 
      positions = as.numeric(stringr::str_extract(univ.id, "[0-9]+$")),
      reference = ref_genome,
      ref_genes = ref_genes, 
      onco_genes = onco_genes, 
      bad_actors = special_genes
      )
  )

cond_sites <- fil_sites %>%
  dplyr::group_by(
    patient, specimen, condition, univ.id, gene.id, replicate, uniqueRegion
  ) %>%
  dplyr::summarise(
    estAbund = n(), 
    reads = sum(counts)
  ) %>%
  dplyr::group_by(
    patient, specimen, condition, univ.id, gene.id, replicate
  ) %>%
  dplyr::summarise(
    estAbund = mean(estAbund),
    reads = sum(reads),
    detect = paste0(unique(uniqueRegion), collapse = ":")
  ) %>%
  dplyr::group_by(
    patient, specimen, condition, univ.id, gene.id
  ) %>%
  dplyr::summarise(
    estAbund = sum(estAbund), 
    reads = sum(reads),
    dual.detect = length(unique(unlist(strsplit(detect, ":")))) == 2,
    detect = paste(
      sort(unique(unlist(strsplit(detect, ":")))), collapse = ":"
    )
  ) %>%
  dplyr::ungroup()


# Specimen summary ----
# Summarize components and append to specimen table
summary_tbl <- cond_sites %>%
  dplyr::group_by(patient, specimen, condition) %>%
  dplyr::summarise(
    reads = sum(reads),
    infer.cells = round(sum(estAbund)),
    uniq.sites = dplyr::n_distinct(univ.id),
    dual.detect = sum(dual.detect),
    shannon = pop_calcs(estAbund, calc = "shannon"),
    gini = pop_calcs(estAbund, calc = "gini"),
    uc50 = pop_calcs(estAbund, calc = "uc50"),
    chao = calculateChao(as.character(S4Vectors::Rle(univ.id, estAbund)))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(specimen) %>%
  dplyr::select(-patient, -condition) %>%
  dplyr::full_join(cond_overview, ., by = "specimen", fill = 0)


# Additional analyses in development ----



cat("Analysis complete.\nStarting report generation...\n")

# Data passed to Rmd for report generation ----
set_names <- ifelse(
  length(configs) ==  1, 
  names(configs), 
  paste0(
    paste(names(configs)[seq_len(length(configs)-1)], collapse = ", "),
    ", and ", 
    names(configs)[length(configs)]
  )
)

# Normalize file output path
write(c(), file = args$output)
args$output <- normalizePath(args$output)
unlink(args$output)

output_path <- unlist(strsplit(args$output, "/"))
output_dir <- paste(output_path[seq_len(length(output_path)-1)], collapse = "/")
output_file <- output_path[length(output_path)]

if( args$format == "html" & !stringr::str_detect(output_file, ".html$") ){
  output_file <- paste0(output_file, ".html")
}

if( args$format == "pdf" & !stringr::str_detect(output_file, ".pdf$") ){
  output_file <- paste0(output_file, ".pdf")
}

if( args$data ){

  if( args$format == "html" ){
    
    save.image(file = file.path(
      output_dir, stringr::str_replace(output_file, ".html$", ".RData")
    )) 
    
  }else if( args$format == "pdf" ){
    
    save.image(file = file.path(
      output_dir, stringr::str_replace(output_file, ".pdf$", ".RData")
    )) 
    
  }
  
}

if( args$format == "html" ){
  
  css_path <- normalizePath(
    file.path(code_dir, "report_templates/template.css")
  )
  
  rmarkdown::render(
    input = template_path,
    output_format = output_format, 
    output_file = output_file,
    output_dir = output_dir,
    output_options = list("css" = css_path)
  )
  
}else{
  
  rmarkdown::render(
    input = template_path,
    output_format = output_format, 
    output_file = output_file,
    output_dir = output_dir
  )
  
}

q()
