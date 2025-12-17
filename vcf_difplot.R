#!/usr/bin/env Rscript

# VCF Variant Position Plotting Script
# This script reads a tab-delimited file converted from VCF using GATK VariantsToTable
# and plots variant positions using ggplot2

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
})

# Constants
MAX_DISPLAY_ROWS <- 20  # Maximum number of positions to display in console output

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input tab-delimited file (required)", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="variant_plot.pdf",
              help="Output plot file [default=%default]", metavar="FILE"),
  make_option(c("-b", "--basename"), type="character", default=NULL,
              help="Baseline sample name", metavar="NAME"),
  make_option(c("-B", "--basecol"), type="integer", default=NULL,
              help="Baseline sample column position (1-based)", metavar="INT"),
  make_option(c("-c", "--copname"), type="character", default=NULL,
              help="Comparison sample name", metavar="NAME"),
  make_option(c("-C", "--copcol"), type="integer", default=NULL,
              help="Comparison sample column position (1-based)", metavar="INT"),
  make_option(c("-l", "--chrlength"), type="character", default=NULL,
              help="Chromosome length file (tab-delimited: CHROM LENGTH)", metavar="FILE"),
  make_option(c("-u", "--unit"), type="numeric", default=1e6,
              help="Chromosome length unit [default=%default]", metavar="NUM"),
  make_option(c("--baseHetcheck"), action="store_true", default=FALSE,
              help="Check if baseline sample is homozygous (e.g., A/A, G|G); ignore heterozygous positions"),
  make_option(c("--copHetcheck"), action="store_true", default=FALSE,
              help="Check if comparison sample is homozygous (e.g., A/A, G|G); ignore heterozygous positions"),
  make_option(c("--segmentColor"), type="character", default="red",
              help="Color for variant position segments [default=%default]", metavar="COLOR"),
  make_option(c("--segmentSize"), type="numeric", default=0.5,
              help="Thickness of variant position segments [default=%default]", metavar="NUM"),
  make_option(c("--chrBorderColor"), type="character", default="black",
              help="Color for chromosome borders [default=%default]", metavar="COLOR"),
  make_option(c("--chrBorderSize"), type="numeric", default=0.3,
              help="Thickness of chromosome borders [default=%default]", metavar="NUM")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="\nPlot variant positions from GATK VariantsToTable output\n\nExample usage:\n  Rscript vcf_difplot.R -i input.table -b sample1 -c sample2 -o output.pdf")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required (-i/--input)", call.=FALSE)
}

# Check if input file exists
if (!file.exists(opt$input)) {
  stop(paste("Input file does not exist:", opt$input), call.=FALSE)
}

# Validate baseline sample specification
if (is.null(opt$basename) && is.null(opt$basecol)) {
  stop("Either baseline sample name (-b/--basename) or column position (-B/--basecol) must be specified", call.=FALSE)
}

if (!is.null(opt$basename) && !is.null(opt$basecol)) {
  warning("Both baseline name and column specified. Using sample name.")
  opt$basecol <- NULL
}

# Validate comparison sample specification
if (is.null(opt$copname) && is.null(opt$copcol)) {
  stop("Either comparison sample name (-c/--copname) or column position (-C/--copcol) must be specified", call.=FALSE)
}

if (!is.null(opt$copname) && !is.null(opt$copcol)) {
  warning("Both comparison name and column specified. Using sample name.")
  opt$copcol <- NULL
}

# Read input file
cat("Reading input file:", opt$input, "\n")
data <- tryCatch({
  read.table(opt$input, header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="")
}, error = function(e) {
  stop(paste("Error reading input file:", e$message), call.=FALSE)
})

# Check required columns
required_cols <- c("CHROM", "POS")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse=", ")), call.=FALSE)
}

# Find GT columns (format: sampleID.GT)
gt_cols <- grep("\\.GT$", colnames(data), value=TRUE)
if (length(gt_cols) == 0) {
  stop("No GT columns found in input file. GT columns should be named as 'sampleID.GT'", call.=FALSE)
}

cat("Found", length(gt_cols), "samples with GT information\n")
sample_names <- sub("\\.GT$", "", gt_cols)
cat("Sample names:", paste(sample_names, collapse=", "), "\n")

# Determine baseline sample column
if (!is.null(opt$basename)) {
  base_col <- paste0(opt$basename, ".GT")
  if (!base_col %in% colnames(data)) {
    stop(paste("Baseline sample not found:", opt$basename, "\nAvailable samples:", paste(sample_names, collapse=", ")), call.=FALSE)
  }
  cat("Using baseline sample:", opt$basename, "\n")
} else {
  # Use column position
  if (opt$basecol < 1 || opt$basecol > length(gt_cols)) {
    stop(paste("Baseline column position out of range. Must be between 1 and", length(gt_cols)), call.=FALSE)
  }
  base_col <- gt_cols[opt$basecol]
  cat("Using baseline sample from column", opt$basecol, ":", sample_names[opt$basecol], "\n")
}

# Determine comparison sample column
if (!is.null(opt$copname)) {
  comp_col <- paste0(opt$copname, ".GT")
  if (!comp_col %in% colnames(data)) {
    stop(paste("Comparison sample not found:", opt$copname, "\nAvailable samples:", paste(sample_names, collapse=", ")), call.=FALSE)
  }
  cat("Using comparison sample:", opt$copname, "\n")
} else {
  # Use column position
  if (opt$copcol < 1 || opt$copcol > length(gt_cols)) {
    stop(paste("Comparison column position out of range. Must be between 1 and", length(gt_cols)), call.=FALSE)
  }
  comp_col <- gt_cols[opt$copcol]
  cat("Using comparison sample from column", opt$copcol, ":", sample_names[opt$copcol], "\n")
}

# Helper functions for genotype processing

# Parse genotype into alleles, handling missing data
# Input: genotype string like "A/T", "C|G", "./.", etc.
# Output: vector of alleles, or NULL if missing data
parse_genotype <- function(gt) {
  if (is.na(gt) || gt == "") {
    return(NULL)
  }
  
  # Replace | with / for consistent processing
  gt <- gsub("\\|", "/", gt)
  
  # Split by /
  alleles <- strsplit(gt, "/")[[1]]
  
  # Check for missing data
  if (any(alleles == ".")) {
    return(NULL)
  }
  
  # Validate: should have exactly 2 alleles for diploid genotypes
  if (length(alleles) != 2) {
    warning(paste("Unexpected number of alleles in genotype:", gt, "- expected 2, got", length(alleles)))
    return(NULL)
  }
  
  return(alleles)
}

# Normalize genotype: treat / and | as equivalent separators
# Input: genotype string like "A/T", "C|G", "./.", etc.
# Output: sorted alleles separated by "/" (e.g., "A/T" -> "A/T", "T|A" -> "A/T")
normalize_genotype <- function(gt) {
  alleles <- parse_genotype(gt)
  
  if (is.null(alleles)) {
    return(NA)
  }
  
  # Sort alleles to make "A/T" and "T/A" equivalent
  alleles_sorted <- sort(alleles)
  
  # Return normalized genotype
  return(paste(alleles_sorted, collapse="/"))
}

# Check if genotype is homozygous
# Input: genotype string like "A/A", "G|G", "A/T"
# Output: TRUE if homozygous, FALSE if heterozygous, NA if missing
is_homozygous <- function(gt) {
  alleles <- parse_genotype(gt)
  
  if (is.null(alleles)) {
    return(NA)
  }
  
  # Check if all alleles are the same
  return(length(unique(alleles)) == 1)
}

# Compare genotypes and mark variants
cat("Comparing genotypes...\n")

# Normalize genotypes for comparison
data$base_gt_norm <- sapply(data[[base_col]], normalize_genotype)
data$comp_gt_norm <- sapply(data[[comp_col]], normalize_genotype)

# Create initial filter: exclude positions with missing data (./.)
data$keep <- !is.na(data$base_gt_norm) & !is.na(data$comp_gt_norm)

cat("Positions after removing missing data (./.): ", sum(data$keep), "\n")

# Apply baseline homozygosity check if requested
if (opt$baseHetcheck) {
  cat("Applying baseline homozygosity check...\n")
  base_homo <- sapply(data[[base_col]], is_homozygous)
  data$keep <- data$keep & !is.na(base_homo) & base_homo
  cat("Positions after baseline homozygosity filter: ", sum(data$keep), "\n")
}

# Apply comparison homozygosity check if requested
if (opt$copHetcheck) {
  cat("Applying comparison homozygosity check...\n")
  comp_homo <- sapply(data[[comp_col]], is_homozygous)
  data$keep <- data$keep & !is.na(comp_homo) & comp_homo
  cat("Positions after comparison homozygosity filter: ", sum(data$keep), "\n")
}

# Filter data based on all criteria
data <- data[data$keep, ]

# Compare normalized genotypes
data$is_variant <- data$base_gt_norm != data$comp_gt_norm

cat("Total positions:", nrow(data), "\n")
cat("Variant positions:", sum(data$is_variant), "\n")
cat("Non-variant positions:", sum(!data$is_variant), "\n")

# Print first 20 positions that meet criteria
cat("\n=== First", MAX_DISPLAY_ROWS, "positions that meet filtering criteria ===\n")
if (nrow(data) > 0) {
  # Get first N rows (or all if less than N)
  n_display <- min(MAX_DISPLAY_ROWS, nrow(data))
  
  # Create display data frame with desired column names
  display_subset <- data.frame(
    CHROM = data$CHROM[1:n_display],
    POS = data$POS[1:n_display],
    Baseline_GT = data[[base_col]][1:n_display],
    Comparison_GT = data[[comp_col]][1:n_display],
    stringsAsFactors = FALSE
  )
  
  # Print as a formatted table
  print(display_subset, row.names=FALSE)
  
  if (nrow(data) > MAX_DISPLAY_ROWS) {
    cat(paste0("\n... (showing ", MAX_DISPLAY_ROWS, " of ", nrow(data), " total positions)\n"))
  }
} else {
  cat("No positions meet the filtering criteria.\n")
}
cat("========================================================\n\n")

# Get chromosome information
chromosomes <- unique(data$CHROM)
cat("Chromosomes found:", paste(chromosomes, collapse=", "), "\n")

# Determine chromosome lengths
if (!is.null(opt$chrlength)) {
  # Read chromosome length file
  if (!file.exists(opt$chrlength)) {
    stop(paste("Chromosome length file does not exist:", opt$chrlength), call.=FALSE)
  }
  
  cat("Reading chromosome length file:", opt$chrlength, "\n")
  
  # Intelligent separator detection
  # Read first few lines to detect separator (max 5 lines)
  first_lines <- readLines(opt$chrlength, n=5, warn=FALSE)
  
  # Count occurrences of different separators across all lines
  # gregexpr returns -1 when no match found, so we need to check for that
  count_separator <- function(lines, pattern) {
    sum(sapply(lines, function(x) {
      matches <- gregexpr(pattern, x)[[1]]
      if (matches[1] == -1) return(0)
      return(length(matches))
    }))
  }
  
  tab_count <- count_separator(first_lines, "\t")
  comma_count <- count_separator(first_lines, ",")
  semicolon_count <- count_separator(first_lines, ";")
  
  # Determine separator based on highest count with explicit priority
  # Priority when counts are tied: tab > comma > semicolon > whitespace
  detected_sep <- ""  # default to whitespace
  sep_name <- "whitespace"
  
  if (tab_count > 0 && tab_count >= comma_count && tab_count >= semicolon_count) {
    detected_sep <- "\t"
    sep_name <- "tab"
  } else if (comma_count > 0 && comma_count >= semicolon_count) {
    detected_sep <- ","
    sep_name <- "comma"
  } else if (semicolon_count > 0) {
    detected_sep <- ";"
    sep_name <- "semicolon"
  }
  
  cat("Detected separator in chromosome length file:", sep_name, "\n")
  
  chr_lengths <- tryCatch({
    read.table(opt$chrlength, header=FALSE, sep=detected_sep, stringsAsFactors=FALSE, col.names=c("CHROM", "LENGTH"))
  }, error = function(e) {
    stop(paste("Error reading chromosome length file:", e$message), call.=FALSE)
  })
  
  # Check if all chromosomes are in the length file
  missing_chrs <- setdiff(chromosomes, chr_lengths$CHROM)
  if (length(missing_chrs) > 0) {
    warning(paste("Some chromosomes not found in length file:", paste(missing_chrs, collapse=", "), 
                  "\nUsing max position for these chromosomes"))
  }
  
} else {
  # Use max position from data
  warning("No chromosome length file provided. Using maximum variant position for each chromosome.")
  chr_lengths <- aggregate(POS ~ CHROM, data=data, FUN=max)
  colnames(chr_lengths) <- c("CHROM", "LENGTH")
}

# Merge chromosome lengths with data
chr_info <- data.frame(CHROM = chromosomes, stringsAsFactors=FALSE)
chr_info <- merge(chr_info, chr_lengths, by="CHROM", all.x=TRUE)

# For chromosomes without length information, use max position from data
for (i in 1:nrow(chr_info)) {
  if (is.na(chr_info$LENGTH[i])) {
    chr_info$LENGTH[i] <- max(data$POS[data$CHROM == chr_info$CHROM[i]])
  }
}

# Sort chromosomes (numeric sorting for numeric chromosomes)
chr_info$chr_num <- suppressWarnings(as.numeric(chr_info$CHROM))
chr_info <- chr_info[order(chr_info$chr_num, chr_info$CHROM, na.last=TRUE), ]
chr_info$chr_order <- 1:nrow(chr_info)

# Scale lengths by unit
chr_info$LENGTH_scaled <- chr_info$LENGTH / opt$unit

cat("\nChromosome lengths (scaled by", opt$unit, "):\n")
print(chr_info[, c("CHROM", "LENGTH_scaled")])

# Prepare data for plotting
plot_data <- merge(data, chr_info[, c("CHROM", "chr_order", "LENGTH_scaled")], by="CHROM")
plot_data$POS_scaled <- plot_data$POS / opt$unit

# Filter for variants only
variant_data <- plot_data[plot_data$is_variant, ]

# Create plot
cat("\nGenerating plot...\n")

# Format unit label for better readability
unit_label <- if (opt$unit == 1e6) {
  "Mb"
} else if (opt$unit == 1e3) {
  "kb"
} else if (opt$unit == 1) {
  "bp"
} else {
  paste0(format(opt$unit, scientific = FALSE), " bp")
}

# Determine the correct parameter name for line width based on ggplot2 version
# ggplot2 >= 3.4.0 uses linewidth, older versions use size
ggplot2_version <- packageVersion("ggplot2")
use_linewidth <- ggplot2_version >= "3.4.0"

# Build the plot with version-appropriate parameters
if (use_linewidth) {
  p <- ggplot() +
    # Draw chromosome rectangles
    geom_rect(data=chr_info, 
              aes(xmin=0, xmax=LENGTH_scaled, ymin=chr_order-0.4, ymax=chr_order+0.4),
              fill="lightgray", color=opt$chrBorderColor, linewidth=opt$chrBorderSize) +
    # Draw variant positions as vertical segments
    geom_segment(data=variant_data,
                 aes(x=POS_scaled, xend=POS_scaled, y=chr_order-0.4, yend=chr_order+0.4),
                 color=opt$segmentColor, linewidth=opt$segmentSize, alpha=0.6) +
    scale_y_continuous(breaks=chr_info$chr_order, labels=chr_info$CHROM) +
    labs(x=paste0("Position (", unit_label, ")"),
         y="Chromosome",
         title="Variant Position Plot",
         subtitle=paste("Comparing", sub("\\.GT$", "", base_col), "vs", sub("\\.GT$", "", comp_col))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size=10))
} else {
  p <- ggplot() +
    # Draw chromosome rectangles
    geom_rect(data=chr_info, 
              aes(xmin=0, xmax=LENGTH_scaled, ymin=chr_order-0.4, ymax=chr_order+0.4),
              fill="lightgray", color=opt$chrBorderColor, size=opt$chrBorderSize) +
    # Draw variant positions as vertical segments
    geom_segment(data=variant_data,
                 aes(x=POS_scaled, xend=POS_scaled, y=chr_order-0.4, yend=chr_order+0.4),
                 color=opt$segmentColor, size=opt$segmentSize, alpha=0.6) +
    scale_y_continuous(breaks=chr_info$chr_order, labels=chr_info$CHROM) +
    labs(x=paste0("Position (", unit_label, ")"),
         y="Chromosome",
         title="Variant Position Plot",
         subtitle=paste("Comparing", sub("\\.GT$", "", base_col), "vs", sub("\\.GT$", "", comp_col))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size=10))
}

# Save plot
cat("Saving plot to:", opt$output, "\n")

# Determine output format based on file extension
output_ext <- tolower(tools::file_ext(opt$output))
if (output_ext == "pdf") {
  pdf(opt$output, width=12, height=max(4, nrow(chr_info) * 0.5))
} else if (output_ext == "png") {
  png(opt$output, width=1200, height=max(400, nrow(chr_info) * 50), res=100)
} else if (output_ext %in% c("jpg", "jpeg")) {
  jpeg(opt$output, width=1200, height=max(400, nrow(chr_info) * 50), res=100)
} else {
  # Default to PDF
  warning(paste("Unsupported output format:", output_ext, ". Using PDF instead."))
  opt$output <- sub(paste0("\\.", output_ext, "$"), ".pdf", opt$output)
  pdf(opt$output, width=12, height=max(4, nrow(chr_info) * 0.5))
}

print(p)
dev.off()

cat("\nDone! Plot saved to:", opt$output, "\n")
cat("Total variants plotted:", nrow(variant_data), "\n")
