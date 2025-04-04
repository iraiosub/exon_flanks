library(data.table)
library(rtracklayer)
library(parallel)

# Function that renames exon flanks bins in relation to splice sites
rename_bins <- function(gr) {
  
  gr$meta_bin <- as.numeric(NA)
  
  if (unique(strand(gr)) == "+" & unique(gr$type) == "upstream") {
    
    
    print(length(gr))
    print(gr$bin)
    
    gr$meta_bin <- gr$bin - length(gr) - 1
    
    print(gr$meta_bin)
    
  } else if (unique(strand(gr)) == "+" & unique(gr$type) == "downstream") {
    
    gr$meta_bin <- gr$bin
    
  } else if (unique(strand(gr)) == "-" & unique(gr$type) == "upstream") {
    
    gr$meta_bin <- length(gr) - gr$bin + 1
    
  } else if (unique(strand(gr)) == "-" & unique(gr$type) == "downstream") {
    
    gr$meta_bin <- -gr$bin
  } else {
    print ("no condition met")
  }
  
  return(gr)
}


# Function to bin and rename exon flanks
prepare_exon_flanks <- function(exons.gr, window.size, bin.size, cores = 4) {
  
  # Number of bins in the window
  bin.n <- window.size / bin.size
  
  # Upstream to CE
  upstream.gr <- IRanges::shift(exons.gr, -(window.size))
  
  # Resetting the ranges for the upstream.gr; start of the exons.gr becomes the end; 3'ss or 5'ss
  ranges(upstream.gr) <- IRanges(start=start(upstream.gr), end=start(exons.gr) - 1) # exclude actual ss
  
  # Downstream to CE
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
  
  # Assigning type
  upstream.binned.gr$type <- "upstream"
  downstream.binned.gr$type <- "downstream"
  
  # Split by exon ID to rename bins for each
  upstream.binned.grl <- split(upstream.binned.gr, upstream.binned.gr$id)
  downstream.binned.grl <- split(downstream.binned.gr, downstream.binned.gr$id)
  
  # Renaming bins
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

# Function for downstream analyses to find overlaps
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


# test data prep
exons.df <- fread("~/Documents/projects/kennymatched.inverted.repeat_pair.alt.exons_Long.tsv")
ir_exons.df <- exons.df %>%
  dplyr::filter(type == "IRAlus-flanked alternative")

# Make a GRanges object out of the exons of interest
ir_exons.gr <- GRanges(seqnames = ir_exons.df$chrom,
                       ranges = IRanges(start = ir_exons.df$exon_start, end = ir_exons.df$exon_end),
                       strand = ir_exons.df$strand, id = ir_exons.df$id)
# we have exons!
exons.gr <- ir_exons.gr[1:2]


# Set window and bin size of interest
window.nt <- 100 # flank size in nt
bin.nt <- 10    # bin size in nt
# The bin number is calculated by prepare_exon_flanks as window.nt/bin.nt

# Test mega-function
test <- prepare_exon_flanks(exons.gr, window.size = window.nt, bin.size = bin.nt)

# Access useful data
test$exon_flanks


#### Manual buiding of the function and testing

bin.n <- window.size/bin.size

#upstream to CE
#define the shift distance window
# shift_distance <- -window
upstream.gr <- GenomicRanges::shift(exons.gr, -window.size)

# Resetting the ranges for the upstream.gr; start of the exons.gr becomes the end; 3'ss or 5'ss
ranges(upstream.gr) <- IRanges(start=start(upstream.gr), end = start(exons.gr) - 1) # exclude actual ss

# downstream to CE
downstream.gr<- GenomicRanges::shift(exons.gr, window.size)
ranges(downstream.gr)<- IRanges(start=end(exons.gr) + 1, end =end(downstream.gr)) # exclude ss



# Now sanity check 
stopifnot(unique(width(upstream.gr)) == unique(width(downstream.gr)) )
stopifnot(unique(width(upstream.gr) == window.size))


# Now we can start binningm keeping exon ID
upstream.binned.gr <- unlist(GenomicRanges::tile(upstream.gr, n = bin.n))
upstream.binned.gr$id <- rep(upstream.gr$id, each = bin.n)

downstream.binned.gr <- unlist(GenomicRanges::tile(downstream.gr, n = bin.n))
downstream.binned.gr$id <- rep(downstream.gr$id, each = bin.n)

# Correctly adjust bins
# first initialise bins


upstream.binned.gr$bin <- rep(c(1:bin.n), 
                              times = length(unique(upstream.gr$id)))

downstream.binned.gr$bin <- rep(c(1:bin.n), 
                              times = length(unique(downstream.gr$id)))

# Convert to GRangesList to treat individual exons bins more easily/securely
upstream.binned.gr$type <- "upstream"
downstream.binned.gr$type <- "downstream"

# Split by exon ID to rename bins for each
upstream.binned.grl <- split(upstream.binned.gr, upstream.binned.gr$id)
downstream.binned.grl <- split(downstream.binned.gr, downstream.binned.gr$id)
# Apply renaming
upstream.binned.renamed.gr <- unlist(GRangesList(mclapply(upstream.binned.grl, rename_bins), mc.cores = 4))
downstream.binned.renamed.gr <- unlist(GRangesList(mclapply(downstream.binned.grl, rename_bins), mc.cores = 4))

# Now reunite upstream and downstream together and clean up
exon_flanks.binned.renamed.gr <- c(upstream.binned.renamed.gr, downstream.binned.renamed.gr)
clean.exon_flanks.binned.renamed.gr <- exon_flanks.binned.renamed.gr
clean.exon_flanks.binned.renamed.gr$type <- NULL
clean.exon_flanks.binned.renamed.gr$bin <- clean.exon_flanks.binned.renamed.gr$meta_bin
clean.exon_flanks.binned.renamed.gr$meta_bin <- NULL


# Now the final object for finding your overlaps is: clean.exon_flanks.binned.renamed.gr




