
library(data.table)
library(rtracklayer)
library(parallel)
library(tictoc)

# Function that renames exon flanks bins in relation to splice sites
rename_bins <- function(gr) {
  
  gr$meta_bin <- as.numeric(NA)
  
  if (unique(strand(gr)) == "+" & unique(gr$type) == "upstream") {
    
    gr$meta_bin <- gr$bin - length(gr) - 1
    
    print(gr$meta_bin)
    
  } else if (unique(strand(gr)) == "+" & unique(gr$type) == "downstream") {
    
    gr$meta_bin <- gr$bin
    
  } else if (unique(strand(gr)) == "-" & unique(gr$type) == "upstream") {
    
    gr$meta_bin <- length(gr) - gr$bin + 1
    
  } else if (unique(strand(gr)) == "-" & unique(gr$type) == "downstream") {
    
    gr$meta_bin <- -gr$bin
  } else {
    warning("No matching condition met for strand and type.")
  }
  
  return(gr)
}

# Function to bin and rename exon flanks
prepare_exon_flanks <- function(exons.gr, window.size, bin.size, cores = 4) {
  
  # Number of bins in the window
  bin.n <- window.size / bin.size
  
  # Upstream to exons (genomic coordinates, ignoring strand)
  upstream.gr <- IRanges::shift(exons.gr, -(window.size))
  
  # Resetting the ranges for the upstream.gr; start of the exons.gr becomes the end; 3'ss or 5'ss
  ranges(upstream.gr) <- IRanges(start=start(upstream.gr), end=start(exons.gr) - 1) # exclude ss
  
  # Downstream to exons (genomic coordinates, ignoring strand)
  downstream.gr <- IRanges::shift(exons.gr, window.size)
  ranges(downstream.gr) <- IRanges(start=end(exons.gr) + 1, end=end(downstream.gr)) # exclude ss
  
  # Sanity check
  stopifnot(unique(width(upstream.gr)) == unique(width(downstream.gr)))
  stopifnot(unique(width(upstream.gr)) == window.size)
  
  # Binning upstream and downstream
  upstream.binned.gr <- unlist(tile(upstream.gr, n = bin.n))
  upstream.binned.gr$id <- rep(upstream.gr$id, each = bin.n)
  
  downstream.binned.gr <- unlist(tile(downstream.gr, n = bin.n))
  downstream.binned.gr$id <- rep(downstream.gr$id, each = bin.n)
  
  # Assigning bin numbers
  upstream.binned.gr$bin <- rep(1:bin.n, times = length(unique(upstream.gr$id)))
  downstream.binned.gr$bin <- rep(1:bin.n, times = length(unique(downstream.gr$id)))
  
  # Assign type so bins can be renamed accordingly
  upstream.binned.gr$type <- "upstream"
  downstream.binned.gr$type <- "downstream"
  
  # Split by exon ID to rename bins for each
  upstream.binned.grl <- split(upstream.binned.gr, upstream.binned.gr$id)
  downstream.binned.grl <- split(downstream.binned.gr, downstream.binned.gr$id)
  
  # Rename bins
  upstream.binned.renamed.gr <- unlist(GRangesList(mclapply(upstream.binned.grl, rename_bins, mc.cores = cores)))
  downstream.binned.renamed.gr <- unlist(GRangesList(mclapply(downstream.binned.grl, rename_bins, mc.cores = cores)))
  
  # Combine upstream and downstream
  exon_flanks.binned.renamed.gr <- c(upstream.binned.renamed.gr, downstream.binned.renamed.gr)
  
  # Clean up
  clean.exon_flanks.binned.renamed.gr <- exon_flanks.binned.renamed.gr
  clean.exon_flanks.binned.renamed.gr$type <- NULL
  clean.exon_flanks.binned.renamed.gr$bin <- clean.exon_flanks.binned.renamed.gr$meta_bin
  clean.exon_flanks.binned.renamed.gr$meta_bin <- NULL
  
  return(list(exon_flanks = clean.exon_flanks.binned.renamed.gr,
              debug_exon_flanks = exon_flanks.binned.renamed.gr))
}


## Alternative functions: either use the ones above or these ones; see usage below

# Function to process either upstream or downstream per exon
process_exon_flanks <- function(exon, direction) {
  id <- exon$id
  strand_val <- as.character(strand(exon))
  
  # Shift and define flank region
  if (direction == "upstream") {
    flank.gr <- IRanges::shift(exon, -window.size)
    ranges(flank.gr) <- IRanges(start = start(flank.gr), end = start(exon) - 1)
  } else if (direction == "downstream") {
    flank.gr <- IRanges::shift(exon, window.size)
    ranges(flank.gr) <- IRanges(start = end(exon) + 1, end = end(flank.gr))
  } else {
    stop("Invalid direction.")
  }
  
  # Tile into bins
  binned <- tile(flank.gr, n = bin.n)[[1]]
  binned$id <- rep(id, bin.n)
  binned$bin <- 1:bin.n
  binned$type <- direction
  
  # Rename bins
  binned <- rename_bins(binned)
  return(binned)

}

prepare_exon_flanks_simpler <- function(exons.gr, window.size, bin.size, cores = 4) {
  bin.n <- window.size / bin.size
  stopifnot(bin.n %% 1 == 0)  # ensure clean binning
  
  # Split input exons.gr by ID
  exons.grl <- split(exons.gr, exons.gr$id)
  
  # Process upstream and downstream in parallel
  upstream.binned <- mclapply(exons.grl, process_exon_flanks, direction = "upstream", mc.cores = cores)
  downstream.binned <- mclapply(exons.grl, process_exon_flanks, direction = "downstream", mc.cores = cores)
  
  # Combine and clean
  all_binned <- c(unlist(GRangesList(upstream.binned)), unlist(GRangesList(downstream.binned)))
  all_binned$bin <- all_binned$meta_bin
  all_binned$meta_bin <- all_binned$type <- NULL
  
  return(list(exon_flanks = all_binned,
              debug_exon_flanks = c(unlist(upstream.binned), unlist(downstream.binned))))
}




# test data prep
exons.df <- fread("~/Documents/projects/kennymatched.inverted.repeat_pair.alt.exons_Long.tsv")
ir_exons.df <- exons.df %>%
  dplyr::filter(type == "IRAlus-flanked alternative")

# Make a GRanges object out of the exons of interest
ir_exons.gr <- GRanges(seqnames = ir_exons.df$chrom,
                       ranges = IRanges(start = ir_exons.df$exon_start, end = ir_exons.df$exon_end),
                       strand = ir_exons.df$strand, id = ir_exons.df$id)
# we have exons!
#exons.gr <- ir_exons.gr[1:2]
exons.gr <- ir_exons.gr

# Set window and bin size of interest
window.nt <- 100 # flank size in nt
bin.nt <- 10    # bin size in nt
# The bin number is calculated by prepare_exon_flanks as window.nt/bin.nt



# Benchmarking speed 2 approaches: prepare_exon_flanks is faster

tic()
# Test mega-function 1
test1 <- prepare_exon_flanks(exons.gr, window.size = window.nt, bin.size = bin.nt)
toc()
# Access useful data
test1$exon_flanks


tic()
# Test mega-function 2
test2 <- prepare_exon_flanks_simpler(exons.gr, window.size = window.nt, bin.size = bin.nt)
toc()

# Access useful data
test2$exon_flanks




# Function for downstream analyses to find overlaps between data and bins
get_bin_overlaps <- function(windows.gr, data.gr, use_center = FALSE) {
  
  data.gr$bin <- as.numeric(NA)
  
  if (use_center) {
    
    # keep only the central 1 nt from your GRanges
    data.gr <- resize(data.gr, width = 1, fix = "center")
  }
  
  
  
  ol <- findOverlaps(data.gr, windows.gr)
  
  data.gr[queryHits(ol)]$bin <- windows.gr[subjectHits(ol)]$bin
  
  return(data.gr)
  
}


