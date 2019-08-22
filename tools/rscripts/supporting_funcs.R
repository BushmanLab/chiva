#' Assign gene ids given integration site locations and reference data.
#' @param seqnames character vector of chromosome names.
#' @param positions numeric vector of integration loci on chromosomes.
#' @param reference reference genome, i.e. hg38. Reference genome should contain
#' seqInfo data about the lengths of chromosomes and needs to be the same genome
#' that the integration sites were mapped and the `ref_genes` object is mapped 
#' as well.
#' @param ref_genes GRanges object with reference gene data.
#' @param onco_genes character vector of genes designated as cancer-related 
#' genes.
#' @param bad_actors character vector of genes designated as bad actors.
#' @param annotations logical. Should annotations be placed on gene ids? If
#' TRUE, then "*" will be placed for sites within genes, "~" for genes matching
#' to the `onco_genes` object, and "!" for genes matching to the `bad_actors` 
#' object.

assign_gene_id <- function(seqnames, positions, reference, ref_genes, 
                           onco_genes, bad_actors, annotations = TRUE){
  require(GenomicRanges)
  require(hiAnnotator)
  
  gr <- GRanges(
    seqnames = seqnames,
    ranges = IRanges(start = positions, width = 1),
    strand = "+",
    seqinfo = seqinfo(reference))
  
  # Annotate Sites with Gene names for within and nearest genes
  gr <- getSitesInFeature(
    gr, ref_genes, colnam = "in_gene", feature.colnam = "name2")
  gr <- getNearestFeature(
    gr, ref_genes, colnam = "nearest_gene", feature.colnam = "name2")
  
  ## Add gene marks ("*" for in_gene, "~" for onco gene, and "!" for bad_actors)
  gr$gene_id <- ifelse(
    gr$in_gene == "FALSE", gr$nearest_gene, gr$in_gene)
  gr$gene_id <- sapply(strsplit(gr$gene_id, ","), "[[", 1)
  
  gr$annot <- ifelse(
    gr$in_gene == "FALSE", rep("", length(gr)), rep("*", length(gr)))
  
  gr$annot <- ifelse(
    gr$gene_id %in% onco_genes, paste0(gr$annot, "~"), gr$annot)
  
  gr$annot <- ifelse(
    gr$gene_id %in% bad_actors, paste0(gr$annot, "!"), gr$annot)
  
  gr$annot <- ifelse(nchar(gr$annot) > 0, paste0(" ", gr$annot), gr$annot)
  
  if(annotations){
    return(paste0(gr$gene_id, gr$annot))
  }else{
    return(gr$gene_id)
  }
}

#' Compare sequences using multiple sequence alignment (DECIPHER) and converting
#' similar nucleotides to "." format. 
#' 
#' @usage compare_seqs(refSeq, testSeq, withRef = TRUE, gapOpening = -6)
#' 
#' @param refSeq character or XStringSet of length 1. The reference sequence to
#' compare against.
#' @param testSeq character vector or XStringSet. The sequence(s) to compare to 
#' to the reference.
#' @param withRef logical or character. Options: TRUE (default) will return the
#' reference sequence as the first sequence with all nucleotides. "only" will 
#' return only the reference sequence after the multiple sequence alignment. 
#' FALSE will return only the compared sequences. Returned sequences are always
#' in the same order as input.
#' @param gapOpening numeric or integer. Score for opening up a gap within the 
#' multiple sequence alignment.
#' 
#' @example 
#' ref <- "TCTAGCTTCTCAGCA"
#' dna <- c("TCTAGCTACTAAGCA", "TGTAGCTCCCAGCA", "TCTAGTTTCTCAACA")
#' compare_seqs(ref, dna) 
#' 
#' @author Christopher Nobles, Ph.D.

compare_seqs <- function(refSeq, testSeq, withRef = TRUE, 
                         gapOpening = -6, msa = TRUE){
  # Required packages
  loaded <- sapply(
    c("DECIPHER", "Biostrings", "stringr"), require, character.only = TRUE)
  stopifnot(all(loaded))
  
  refSeq <- as.character(refSeq)
  testSeq <- as.character(testSeq)
  
  alignedSeqs <- DNAStringSet(c(refSeq, testSeq))
  seqNames <- c("ref", paste0("test", 1:length(testSeq)))
  names(alignedSeqs) <- seqNames
  if(msa){
    alignedSeqs <- as.character(AlignSeqs(
      alignedSeqs, gapOpening = gapOpening, verbose = FALSE))
  }else{
    alignedSeqs <- as.character(alignedSeqs)
  }
  
  alignedMat <- as.matrix(
    sapply(alignedSeqs, function(x) unlist(strsplit(x, ""))))
  
  if(withRef == "only"){
    return(paste(alignedMat[,1], collapse = ""))
  }else if(withRef){
    return(c(
      paste(alignedMat[,1], collapse = ""),
      sapply(2:(length(testSeq) + 1), function(i){
        paste(sapply(1:nrow(alignedMat), function(j){
          ifelse(alignedMat[j,1] == alignedMat[j,i], ".", alignedMat[j,i])}),
          collapse = "")
      })))
  }else{
    return(sapply(2:(length(testSeq) + 1), function(i){
      paste(sapply(1:nrow(alignedMat), function(j){
        ifelse(alignedMat[j,1] == alignedMat[j,i], ".", alignedMat[j,i])}),
        collapse = "")
    }))
  }
}

#' Convert matches to GRanges objects
#' 
#' @param mindex mindex object, such as the output from 
#' Biostrings::vmatchPattern().
#' @param strand character, either "+" or "-", but indicating which strand was
#' matched.
#' @param ref GRanges object with seqInfo containing sequence lengths for 
#' sequences matched.

mindex_to_granges <- function(mindex, strand, ref){
  if(is.null(mindex@NAMES)) stop("NAMES column not found for seqnames.")
  ir <- unlist(mindex)
  unname(GRanges(
    seqnames = names(ir),
    ranges = ir,
    strand = strand,
    seqinfo = seqinfo(ref)))
}

#' Read psl.gz files into data.frame
#' 
#' @usage read_psl(pslfile)
#' @usage read_psl(pslfile, toNull = c("nCount", "qNumInsert", "qBaseInsert"))
#'  
#' @description Assuming psl.gz files don't have column header. Converts from 
#' 0-base coordinates to 1-base for all alignment information.
#' 
#' @param pslFile character vector of file name(s)
#' @param remove  character vector of column names to remove from final 
#' data.frame.
#' 
#' @return data.frame, psl table
#' @author Christopher Nobles, Ph.D.
#' 

read_psl <- function(pslFile, remove = NULL) {
  stopifnot(require("data.table"))
  stopifnot(all(file.exists(pslFile)))
  
  psl <- lapply(pslFile, function(f) {
    message("Reading ",f)
    data.table::fread( paste("zcat", f), sep="\t" )
  })
  psl <- data.table::rbindlist(psl)
  colnames(psl) <- psl_cols()
  
  # Convert to 1-base rather than 0-base (blat standard)
  psl$qStart <- psl$qStart + 1
  psl$tStart <- psl$tStart + 1
  
  if(!is.null(remove)) psl[, remove] <- NULL
  
  return(as.data.frame(psl))
}

#' Lists psl columns and classes
#' 
#' @usage psl_cols()
#' @usage psl_cols(with.class = TRUE)
#' 
#' @param with.class logical, include class of columns with names (FALSE by 
#' default)
#' 
#' @return character, list of two character vectors.
#' 
#' @author Christopher Nobles, Ph.D.

psl_cols <- function(with.class = FALSE){
  cols <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert",
            "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
            "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
            "blockCount", "blockSizes", "qStarts", "tStarts")
  cols.class <- c(rep("numeric",8), rep("character",2), rep("numeric",3),
                  "character", rep("numeric",4), rep("character",3))
  if(with.class){
    return(list(cols, cols.class))
  }else{
    return(cols)
  }
}

#' Quality filter psl read alignments
#' 
#' @usage quality_filter_psl(psl, qStartMax = NULL, globalIdentityMin = NULL, 
#' baseInsertMax = 5)
#' 
#' @description Using a data.frame of blat alignments (use read_psl() to convert
#' from alignment file to data.frame), filter out alignments that do not pass
#' certain thresholds.
#' 
#' @param psl Data.frame object of alignments (psl format) with alignment 
#' metrics of "matches", "repMatches", "qSize", "qStart", and "tBaseInsert".
#' @param qStartMax integer Maximum allowable nucleotide position to start read
#' alignment on the query sequence.
#' @param globalIdentityMin numeric Between 0 and 100 denoting the minimal 
#' percent of global identity for the alignment to pass the filter.
#' @param baseInsertMax integer Maximum number of allowable inserted nucleotides
#' within the alignment.
#' 
#' @return Subset of input data.frame passing filter criteria.
#' 
#' @author Christopher Nobles, Ph.D.
#' 

quality_filter_psl <- function(psl, qStartMax = NULL, 
                               globalIdentityMin = NULL, baseInsertMax = NULL){
  stopifnot(class(psl) == "data.frame")
  statsNeeded <- c("matches", "repMatches", "qSize", "qStart", "tBaseInsert")
  stopifnot(all(statsNeeded %in% names(psl)))
  
  gPctID <- 100 * (psl$matches + psl$repMatches) / psl$qSize
  
  if(length(qStartMax) > 0){
    psl <- subset(psl, psl$qStart <= qStartMax)}
  if(length(globalIdentityMin) > 0){
    psl <- subset(psl, gPctID >= globalIdentityMin)}
  if(length(baseInsertMax) > 0){
    psl <- subset(psl, psl$tBaseInsert <= baseInsertMax)}
  
  psl
}

#' Convert alignment information from psl format into a GRanges object
#' 
#' @usage psl_to_granges(psl, from, ranges = "t", remove = NULL)
#' 
#' @description From blat alignment information (psl), convert the data.frame
#' object to a GRanges object with the range information from either the target
#' alignments or query alignments. All range information should already be 
#' converted into 1-base coordinates from 0-base (use read_psl() for 
#' compatability). GRanges of query alignments will always have a strand 
#' designation of "+" as these coordinates refer to the alignment of input query
#' sequences on the target sequence. 
#' 
#' @param psl Data.frame in minimal psl format, including columns tName, tStart,
#' tEnd, strand, qName, qStart, qEnd. The former four are required for ranges = 
#' "t", while the latter three are required for ranges = "q".
#' @param from character Single name to designate alignment background / 
#' history.
#' @param ranges character Either "t" for target or "q" for query. Designates 
#' which alignment information will be transfered into the GRanges. Strand will 
#' always be positive for ranges = "q" as these coordinates only reference 
#' beginning and end of query sequence alignment.
#' @param remove character vector Column names to remove from metadata info.
#' 
#' @author Christopher Nobles, Ph.D.
#' 

psl_to_granges <- function(psl, from = NULL, ranges = "t", remove = NULL){
  stopifnot(ranges %in% c("q", "t")) # query or target
  psl$from <- from
  
  if(ranges == "q"){
    algns <- psl[,c("qName", "qStart", "qEnd")]
    algns$strand <- "+"
  }else{
    algns <- psl[,c("tName", "tStart", "tEnd", "strand")]
  }
  names(algns) <- c("names", "start", "end", "strand")
  
  algns.gr <- GRanges(
    seqnames = Rle(algns$names),
    ranges = IRanges(start = algns$start, end = algns$end),
    strand = Rle(algns$strand))
  
  if(!is.null(remove)) psl[, remove] <- NULL
  
  if(ncol(psl) > 1){
    mcols(algns.gr) <- psl
  }else if(ncols(psl) == 1){
    algns.gr$x <- psl[,1]
    names(mcols(algns.gr)) <- names(psl)
  }
  
  algns.gr
}

#' Filter and reassign reads associated with integration sites to suspected 
#' samples.
#' 
#' @usage filter_crossovers(
#'   sites, sample_col = "sampleName", 
#'   counts_col = "counts", cutoff_ratio = 100)
#'   
#' @description Filter and reassign reads (if single origin present) given a 
#' specific cutoff_ratio, or ratio between suspected origin and suspected 
#' crossover read counts. For reads to get c
#' 
#' @param sites GRanges object where each row is a single unique fragment within
#' a sample.
#' @param sample_col character name of column with sample information.
#' @param counts_col character name of column with read count information.
#' @param cutoff_ratio numeric minimum ratio between sample read counts for 
#' which to assign the reads to the numerator sample. Default is 100. i.e. If 
#' sample A had 200 reads assigned to a fragment and sample B had 5, then the 
#' ratio would be 200 / 5 or 40. Since 40 is less than the cutoff_ratio 100, the
#' reads from B would not be reassigned to A. If B instead had 1 read, then the 
#' ratio would change, 200 / 1 or 200. Since 200 is greater than the 
#' cutoff_ratio 100, the read from B is assigned to A. Reassignment occurs if 
#' ratio is equal to cutoff_ratio.
#' @param quiet logical Defaultly, the function will report the number of reads
#' and fragments reassigned. With TRUE, this behavior is silenced.
#' 
#' @author Christopher Nobles, Ph.D.
#' 

filter_crossovers <- function(sites, sample_col = NULL, counts_col = NULL, 
                              cutoff_ratio = 100, quiet = FALSE){
  # Require essential packages
  packs <- c("GenomicRanges", "igraph", "gintools")
  stopifnot(all(sapply(packs, require, character.only = TRUE)))
  
  # Check inputs
  stopifnot(!is.null(sample_col))
  stopifnot(!is.null(counts_col))
  sampleNames <- mcols(sites)[,sample_col]
  read_counts <- mcols(sites)[,counts_col]
  row_id <- paste0(
    seqnames(sites), ":", strand(sites), ":", start(sites), 
    ":", end(sites), ":", sampleNames)
  if(any(as.numeric(names(table(table(row_id)))) > 1)){
    stop("Each row needs to be a unique fragment (start and end) in a sample.")
  }
  
  # Identify the fragments in question for filtering / reassignment
  # Construct a matrix where rows are samples and columns are sites, that way
  # any colSum of a logical test for reads > 0 will tell how many samples a site
  # was identified in.
  assign_matrix <- table(data.frame(
      sampleNames = Rle(values = sampleNames, lengths = read_counts),
      pos_id = Rle(values = generate_posid(sites), lengths = read_counts)))

  crossover_pos_id <- colnames(assign_matrix)[
    which(colSums(assign_matrix > 0) > 0)]
  if(nrow(assign_matrix) > 1) {
    crossover_pos_id <- colnames(assign_matrix)[
      which(colSums(assign_matrix > 0) > 1)]
  }
  # If one sample, then just take all sites
  crossover_sites <- sites[generate_posid(sites) %in% crossover_pos_id]
  
  # For reads to be reassigned, they need to share the same site and breakpoint
  # The same strategy applied above is now applied to the breakpoints for each
  # site. This could likely be done in one large matrix, but here is separated
  # per integration site.
  crossover_sites <- split(crossover_sites, generate_posid(crossover_sites))

  filtered_frags <- do.call(rbind, lapply(crossover_sites, function(x){

    sams <- mcols(x)[,sample_col]
    reads <- mcols(x)[,counts_col]
    bp <- ifelse(strand(x) == "+", end(x), start(x))
    mat <- table(data.frame(
      samples = Rle(values = sams, lengths = reads),
      bp = Rle(values = bp, lengths = reads)))
    co_bp <- colnames(mat)[which(colSums(mat > 0) > 0)]
    if(nrow(mat)>1) {
      co_bp <- colnames(mat)[which(colSums(mat > 0) > 0)]
    }
    co_frags <- x[bp %in% co_bp]

    # With the isolated fragment that are crossing over, determine if their 
    # associated read count ratios are above the cutoff_ratio value.
    fil_frags <- data.frame(
        ori_sample = as.character(mcols(co_frags)[,sample_col]),
        pos_id = as.character(generate_posid(co_frags)),
        bp = ifelse(
          strand(co_frags) == "+", end(co_frags), start(co_frags)),
        counts = mcols(co_frags)[,counts_col],
        stringsAsFactors = FALSE
      ) %>%
      group_by(pos_id, bp) %>%
      mutate(
        adj_sample = ifelse(
          max(counts)/counts >= cutoff_ratio,
          ori_sample[which(counts == max(counts))],
          ori_sample)
      )
    if(nrow(mat)>1) {
      fil_frags <- fil_frags %>% filter(ori_sample != adj_sample)
    }

    return(fil_frags)
    
  }))

  # Transfer the reassignment data back into the dataset using the pos_id, bp,
  # and ori_sample data to identify the fragments. Then set their sample_col 
  # data to the adj_sample.
  adj_index <- sapply(seq_len(nrow(filtered_frags)), function(i){
    
    which(
      generate_posid(sites) == filtered_frags$pos_id[i] &
        ifelse(
          strand(sites) == "+", end(sites), 
          start(sites)) == filtered_frags$bp[i] &
        sampleNames == filtered_frags$ori_sample[i])
    
  })

  adj_sites <- sites
  
  mdata <- mcols(adj_sites)
  mdata[adj_index, sample_col] <- filtered_frags$adj_sample
  mcols(adj_sites) <- mdata
  
  # Message reassignments and produce output.
  filtered_frags <- split(
    filtered_frags, 
    paste0(filtered_frags$ori_sample, ":", filtered_frags$adj_sample)
  )
  
  if(!quiet){
    null <- sapply(filtered_frags, function(df){
      message(
        "Reassigned ", sum(df$counts), " read(s), representing ", 
        length(unique(df$pos_id)), " site(s), from ", unique(df$ori_sample), 
        " to ", unique(df$adj_sample))
    })
  }
  
  adj_sites <- unique_granges(adj_sites, sum.cols = counts_col)
  
  return(adj_sites)
  
}

#' Generate a universal position ID for vectors / viruses sequence from both
#' the U3 and U5 regions.
#' 
#' @usage generate_univID(pos.id, unique.region)
#' 
#' @description Similar to generate_posid from gintools, this function 
#' generates a string with positional information for the integration event,
#' independent of the sequencing direction. This function is based on the 
#' biology of lentiviral integration, with 5bp separation between the 
#' integration sites. The center base is selected as the integration position.
#' 
#' @param posID character position ID, as generated by gintools::generate_posid.
#' @param uniqueRegion character either "U3" or "U5" representing the targeted
#' sequencing region the data was generated from.
#' 
#' @author Christopher Nobles, Ph.D.
#' 

generate_univID <- function(pos.id, unique.region, delim = NULL){
  if(is.null(delim)){
    seqnames <- stringr::str_extract(pos.id, "[\\w]+")
    strand <- stringr::str_extract(pos.id, "[+-]")
    pos <- as.numeric(stringr::str_extract(pos.id, "[0-9]+$"))
    delim <- ""
  }else{
    parsed_id <- stringr::str_split(pos.id, pattern = delim, simplify = TRUE)
    seqnames <- parsed_id[,1]
    strand <- parsed_id[,2]
    pos <- as.numeric(parsed_id[,3])
  }
  
  stopifnot(all(strand %in% c("+", "-")))
  stopifnot(all(unique.region %in% c("U3", "U5")))
  
  # Change position to universal position
  pos <- ifelse(strand == "+", pos + 2, pos - 2)
  strand <- ifelse(
    unique.region == "U5", strand, ifelse(strand == "+", "-", "+"))
  
  # Generate univ id
  gintools::generate_posid(
    seqnames = seqnames, strand = strand, start = pos, end = pos, delim = delim)
}

#' Collapse row contents of a data.frame or matrix into single vector.
#'
#' \code{vcollapse} returns a single vector from input data.frame or matrix 
#' where row contents have been combined or collapsed, as with 
#' `paste(..., collaspe = "")`.
#'
#' @description Similar to python zip, `vzip` takes input vectors and merges
#' them together by their input order and index. A simple example is two numeric
#' vectors, A = c(1,1,1) and B = c(2,2,2). The output of vzip(A,B) would simply
#' be a single vector of c(1,2,1,2,1,2). Any number of vectors can be input, but
#' each input vector must be of the same length. Output vector class depends on
#' input vector consensus.
#'
#' @usage
#' vcollapse(d)
#' vcollapse(d, sep = "-", fill = "NA")
#'
#' @param d data.frame or matrix or object coercible to a matrix. Row contents
#' will be combined into a single output vector.
#' @param sep character used to separate collapsed contents.
#' @param fill character used to fill empty values within the coerced object.
#'
#' @examples
#' df <- data.frame(
#'   "A" = letters[1:5],
#'   "B" = 3:7,
#'   "C" = LETTERS[2:6])
#' vcollapse(df)
#' vcollapse(df, sep = "-")
#' vcollapse(df, sep = "-", fill = "z")
#'
#' @author Christopher Nobles, Ph.D.
#' @export
#'

vcollapse <- function(d, sep = "", fill = NULL){
  
  if( is.vector(d) ){
    stop("Function vcollapse() is not used on vectors, use paste(collapse = ...).")
  }
  
  if( any(sapply(seq_len(ncol(d)), function(i) class(d[,i])) == "factor") ){
    
    fct_idx <- which(
      sapply(seq_len(ncol(d)), function(i) class(d[,i])) == "factor"
    )
    
    mod_env <- new.env()
    mod_env$d <- d
    null <- lapply(fct_idx, function(i){ mod_env$d[,i] <- as.character(d[,i]) })
    d <- mod_env$d
    
  }
  
  if( class(d) != "matrix" ) d <- gsub(" ", "", as.matrix(d))
  mat <- d
  if( !is.null(fill) ) mat <- ifelse(is.na(mat), fill, mat)
  
  mat <- do.call(
    cbind, 
    lapply(seq_len(ncol(mat)), function(i){
      
      if( i < ncol(mat) ){
        cbind(mat[,i], rep(sep, nrow(mat)))
      }else{
        mat[,i]
      }
      
    })
  )
  
  mat <- cbind(mat, rep(">|<", nrow(mat)))
  
  if( any(is.na(mat)) ){
    stop("NA values detected in object, please use fill param for vcollapse.")
  }
  
  div_str <- stringr::str_c(t(mat), collapse = "")
  unlist(strsplit(div_str, split = ">\\|<"))
  
}

#' Calculation for descriptors of populations
#'
#' \code{pop_calcs} Calculations for describing features of populations.
#'
#' @description Given a population of integration sites, TCR/IGH sequences,
#' etc., a numerical value can be used to describe the diveristy (shannon index
#' or entropy), or the clonality (gini or clonality). These numerical values
#' allow for a comparison to be made between populations. Entropy and Clonality
#' are factors used by Adaptive Biotechnologies to describe their TCR / IGH
#' sequence data.
#'
#' @usage
#' pop_calcs(x, calc)
#'
#' @param x numeric a vector of abundances or proportions / frequencies.
#'
#' @param calc character one of the various population calcuations. Choices are:
#' "shannon", "gini", "entropy", "clonality", or "uc50".
#'
#' @examples
#' x <- sample(1:100, 50, replace = TRUE)
#'
#' pop_calcs(x, calc = "shannon")
#' pop_calcs(x, calc = "gini")
#' pop_calcs(x, calc = "entropy")
#' pop_calcs(x, calc = "clonality")
#' pop_calcs(x, calc = "uc50")
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @export

pop_calcs <- function(x, calc){
  if(calc == "shannon"){
    calc_shannon(x)
  }else if(calc == "gini"){
    calc_gini(x)
  }else if(calc == "entropy"){
    calc_entropy(x)
  }else if(calc == "clonality"){
    calc_clonality(x)
  }else if(calc == "uc50"){
    calc_uc50(x)
  }
}

#' @describeIn pop_calcs Calculate Shannon Diversity.
calc_shannon <- function(x, base = exp(1)){
  x <- x[!is.na(x)]
  x <- x/sum(x)
  shannon <- -x*log(x, base)
  shannon <- sum(shannon)
  shannon
}

#' @describeIn pop_calcs Calculate gini index.
calc_gini <- function(x, wt = NULL){
  if(is.null(wt)) wt <- rep(1, length(x))
  
  not_na <- which(!is.na(x))
  wt <- wt[not_na]
  x <- x[not_na]
  
  wt <- wt[order(x)]/sum(wt)  
  x <- x[order(x)]/sum(x)
  
  cum_wt <- cumsum(wt)
  cum_prod <- cumsum(x * wt)
  
  rel_prod <- cum_prod / cum_prod[length(cum_prod)]
  sum(rel_prod[-1] * cum_wt[-length(cum_wt)]) -
    sum(rel_prod[-length(rel_prod)] * cum_wt[-1])
}

#' @describeIn pop_calcs Calculate Entropy as defined by Adaptive Biotech.
calc_entropy <- function(x){
  x <- x[!is.na(x)]
  x <- x/sum(x)
  entropy <- -sum(x*log(x, base = 2))
  entropy
}

#' @describeIn pop_calcs Calcuate Clonality as defined by Adaptive Biotech.
calc_clonality <- function(x){
  x <- x[!is.na(x)]
  x <- x/sum(x)
  clonality <- 1+sum(x*log(x, base = length(x)))
  clonality
}

#' @describeIn pop_calcs Calculate UC50 or the number of must abundant clones
#' which make up half the population.
calc_uc50 <- function(x){
  stopifnot(is.vector(x) & is.numeric(x))
  x <- x[order(x)]
  accum <- sapply(1:length(x), function(i){sum(x[1:i])})
  length(accum[accum >= sum(x)/2])
}

#' Load reference files from config input format.
#'
loadRefFiles <- function(ref, type = "gene.list", freeze = NULL, root = NULL){
  
  stopifnot(type %in% c("gene.list", "GRanges", "data.frame"))
  
  if( !is.null(root) ) stopifnot(dir.exists(root))
  
  if( grepl(".rds$", ref$file) ){
    
    if( file.exists(file.path(root, ref$file)) ){
      ref_set <- readRDS(file.path(root, ref$file))
    }else if( file.exists(ref$file) ){
      ref_set <- readRDS(ref$file)
    }else{
      stop("\n  Cannot find file: ", ref$file, ".\n")
    }
    
  }else if( grepl(".RData$", ref$file) ){
    
    if( file.exists(file.path(root, ref$file)) ){
      rdata_path <- file.path(root, ref$file)
    }else if( file.exists(ref$file) ){
      rdata_path <- ref$file
    }else{
      stop("\n  Cannot find file: ", ref$file, ".\n")
    }
    
    ref_env <- new.env()
    load(rdata_path, envir = refs)
    ref_set <- ref_env[[ls(ref_env)]]
    
  }else if( 
    file.exists(file.path(root, ref$file)) | 
    file.exists(ref$file) | 
    grepl("^http", ref$file) 
  ){
    
    if( file.exists(file.path(root, ref$file)) ){
      ref_set <- data.table::fread(
        file.path(root, ref$file), data.table = FALSE
      )
    }else{
      ref_set <- data.table::fread(ref$file, data.table = FALSE)
    }
    
  }else{
    
    stopifnot( "hiAnnotator" %in% row.names(installed.packages()) )
    stopifnot( grepl(":", ref$file) )
    trackTable <- unlist(strsplit(ref$file, ":"))
    ucsc_session <- hiAnnotator::makeUCSCsession(freeze)
    
    stopifnot(
      trackTable[1] %in% rtracklayer::trackNames(ucsc_session) &
        trackTable[2] %in% rtracklayer::tableNames(
          rtracklayer::ucscTableQuery(ucsc_session, track = trackTable[1])
        )
    )
    
    ref_tbl <- rtracklayer::getTable(rtracklayer::ucscTableQuery(
      ucsc_session, track = trackTable[1], table = trackTable[2]
    ))
    
    ref_set <- rtracklayer::track(
      ucsc_session, name = trackTable[1], table = trackTable[2]
    )
    
    stopifnot( all(ref_tbl$name == ref_set$name) )
    ref_set <- GenomicRanges::granges(ref_set)
    GenomicRanges::mcols(ref_set) <- ref_tbl
    
  }
  
  if( !class(ref_set) %in% c("data.frame", "GRanges") ){
    stop("Import of reference data failed. Check input parameters.")
  }
  
  if( type == "GRanges" ){
    
    if( class(ref_set) == "GRanges" ){
      ref_set$annot_sym <- GenomicRanges::mcols(ref_set)[,ref$symbol]
      return(ref_set)
    }else{
      ref_set$annot_sym <- ref_set[,ref$symbol]
      return(GenomicRanges::makeGRangesFromDataFrame(
        ref_set, keep.extra.columns = TRUE
      ))
    }
    
  }else if( type == "gene.list" ){
    
    ref_set <- try(as.data.frame(ref_set), silent = TRUE)
    
    if( class(ref_set) == "data.frame" ){
      return(ref_set[,ref$symbol])
    }else{
      stop("Cannot coerce to geneList. Check input parameters.")
    }
    
  }else if( type == "data.frame" ){
    
    ref_set <- try(as.data.frame(ref_set), silent = TRUE)
    
    if( class(ref_set) == "data.frame" ){
      ref_set$annot_sym <- ref_set[,ref$symbol]
      return(ref_set)
      
    }else{
      
      stop("Cannot coerce to data.frame. Check input parameters.")
      
    }
    
  }
  
}

#' Functions for calculating Chao1 estimates of population size
jackIID <- function(ids, jrep = NULL, nrep = 10L){
  if(is.null(jrep)) jrep <- sample(rep(1:nrep,length=length(ids)))
  est0 <- vegan::estimateR(table(ids))
  jackrep <- jrep
  urepl <- unique(jrep)
  jackmat <- sapply(
    urepl, function(x) vegan::estimateR(table(ids[jackrep!=x])))
  pseudo <- length(urepl)*est0 - (length(urepl)-1)*jackmat
  rowMeans(pseudo)
}

#' jackknife biased or unbiased Chao estimator
#'
#' @param replicatedSites df with column posid
#' @return number population size estimate
calculateChao <- function(ids, biased=TRUE){
  if ( !biased ) { #regular Chao
    cluster.tab <- table(ids)
    return(round(estimateR(cluster.tab)["S.chao1"]))
  }
  round(jackIID(ids)["S.chao1"])
}

#' coverage plot functions
#' 
# Plots 
calc_coverage <- function(gr, resolution, diff.strands = TRUE){
  #Set up coverage gr
  strandless <- gr
  strand(strandless) <- "*"
  gr_ranges <- range(strandless)
  
  window_seqs <- lapply(gr_ranges, function(chr, res){
    seq(start(chr), end(chr), res)
  }, res = resolution)
  
  coverage_grl <- GRangesList(lapply(
    1:length(gr_ranges), function(i, gr_ranges, window_seqs){
      seqname <- seqnames(gr_ranges[i])
      window <- window_seqs[[i]]
      GRanges(
        seqnames = rep(seqname, length(window)),
        ranges = IRanges(
          start = window, width = rep(resolution, length(window))),
        strand = rep("*", length(window)))
    }, gr_ranges = gr_ranges, window_seqs = window_seqs))
  
  if(diff.strands){
    coverage_pos <- coverage_grl
    coverage_pos <- GRangesList(lapply(coverage_pos, function(x){
      strand(x) <- rep("+", length(x))
      x}))
    coverage_neg <- coverage_grl
    coverage_neg <- GRangesList(lapply(coverage_pos, function(x){
      strand(x) <- rep("-", length(x))
      x}))
  
    return(bind_rows(lapply(1:length(coverage_grl), function(i, gr){
      as.data.frame(coverage_grl[[i]], row.names = NULL) %>%
        select(seqnames, start, end, width) %>%
        mutate(
          readCountsPos = countOverlaps(coverage_pos[[i]], gr),
          readCountsNeg = countOverlaps(coverage_neg[[i]], gr)) %>%
        arrange(seqnames)
    }, gr = gr)))
  }else{
    return(bind_rows(lapply(1:length(coverage_grl), function(i, gr){
      as.data.frame(coverage_grl[[i]], row.names = NULL) %>%
        select(seqnames, start, end, width) %>%
        mutate(readCounts = countOverlaps(coverage_grl[[i]], gr)) %>%
        arrange(seqnames)
    }, gr = gr)))
  }
}

plot_coverage <- function(data, resolution = 10L){
  ggplot(data, aes(x = start)) + 
    geom_bar(
      aes(y = readCounts), 
      stat = "identity", width = resolution) +
    geom_abline(slope = 0, intercept = 0, color = "black") +
    theme_bw() +
    theme(
      legend.position = "None",
      strip.background = element_rect(fill = NA, linetype = 0),
      strip.text = element_text(size = 14, lineheight = 1.1),
      panel.border = element_rect(color = "white"),
      plot.background = element_rect(color = "white"),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_line(color = "white"),
      axis.title = element_text(color = "black", size = 14),
      axis.text = element_text(color = "black", size = 12),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"))
}

extractExons <- function(ref.genes, exon.starts, exon.ends, gene.names){
  
  stopifnot(
    length(ref.genes) == length(exon.starts) &
      length(exon.starts) == length(exon.ends) &
      length(exon.ends) == length(gene.names))
  
  exon.starts <- gsub(",$", "", exon.starts)
  exon.ends <- gsub(",$", "", exon.ends)
  
  exon.starts <- strsplit(exon.starts, ",")
  exon.ends <- strsplit(exon.ends, ",")
  
  stopifnot(all(sapply(
    seq_len(length(gene.names)),
    function(i) length(exon.starts[[i]]) == length(exon.ends[[i]])
  )))
  
  gr <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(
      values = as.character(GenomicRanges::seqnames(ref.genes)),
      lengths = lengths(exon.starts)
    ),
    ranges = IRanges::IRanges(
      start = as.numeric(unlist(exon.starts)),
      end = as.numeric(unlist(exon.ends))
    ),
    strand = S4Vectors::Rle(
      values = as.character(GenomicRanges::strand(ref.genes)),
      length = lengths(exon.starts)
    ),
    seqinfo = GenomeInfoDb::seqinfo(ref.genes)
  )
  
  gene_names <- as.character(S4Vectors::Rle(
    values = gene.names,
    lengths = lengths(exon.starts)
  ))
  
  exon_nums <- unlist(
    mapply(function(strand, exSt){
      if(strand == "+"){
        return(seq_along(exSt))
      }else{
        return(rev(seq_along(exSt)))
      }
    },
    strand = as.character(GenomicRanges::strand(ref.genes)),
    exSt = exon.starts
    )
  )
  
  exon_names <- paste0(gene_names, ":exon", sprintf("%03d", exon_nums))
  
  gr$exon.names <- exon_names
  gr
  
}

annotateFeatures <- function(df, reference, ref.genes, ref.exons, onco.genes,
                             special.genes, gene.name.col = "name2",
                             exon.name.col = "exon.names"){
  
  gr <- GenomicRanges::GRanges(
    seqnames = stringr::str_extract(df$Position, "[\\w]+"),
    ranges = IRanges::IRanges(
      start = as.numeric(stringr::str_extract(df$Position, "[0-9]+$")), 
      width = rep(1, nrow(df))
    ),
    strand = stringr::str_extract(df$Position, "[+-]"),
    seqinfo = GenomicRanges::seqinfo(reference)
  )
  
  if( length(gr) == 0 ){
    
    return(df)
    
  }else{
    
    # Annotate Sites with Gene names for within and nearest genes
    gr <- hiAnnotator::getNearestFeature(
      gr, ref.genes, colnam = "nearestGene", feature.colnam = gene.name.col
    )
    
    gr <- hiAnnotator::getSitesInFeature(
      gr, ref.genes, colnam = "inGene", feature.colnam = gene.name.col
    )
    
    gr <- hiAnnotator::getSitesInFeature(
      gr, ref.exons, colnam = "inExon", feature.colnam = exon.name.col
    )
    
    gr$exonIntron <- ifelse(
      gr$inGene != FALSE & gr$inExon != FALSE, 
      "exonic", 
      ifelse(
        gr$inGene != FALSE & gr$inExon == FALSE, 
        "intronic", 
        "intergenic"
      )
    )
    
    gr$isOnco <- gr$nearestGene %in% onco.genes
    gr$isSpecial <- gr$nearestGene %in% special.genes
    
    return(cbind(df, as.data.frame(GenomicRanges::mcols(gr))))
    
  }
  
}