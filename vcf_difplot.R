#!/usr/bin/env Rscript

# VCF Variant Position Plotting Script
# This script reads a tab-delimited file converted from VCF using GATK VariantsToTable
# and plots variant positions using ggplot2

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
  library(data.table)
})

# Constants
MAX_DISPLAY_ROWS <- 20   # Maximum number of positions to display in console output
MAX_PLOT_HEIGHT_IN <- 50 # Maximum plot height in inches to prevent rendering crashes

# Define command-line options
option_list <- list(
  make_option(c("-I", "--interactive"), action="store_true", default=FALSE,
              help="Enter interactive mode: prompts for each parameter with descriptions (all other flags are ignored)"),
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
              help="Thickness of chromosome borders [default=%default]", metavar="NUM"),
  make_option(c("-t", "--threads"), type="integer", default=1,
              help="[DEPRECATED] Ignored; vectorized processing is used instead", metavar="INT"),
  make_option(c("--output_table"), type="character", default=NULL,
              help="Optional output table file for variant positions used in plotting (tab-delimited: CHROM POS)", metavar="FILE"),
  make_option(c("--CMplot"), action="store_true", default=FALSE,
              help="Use CMplot R package for density plotting instead of ggplot2"),
  make_option(c("--CMplot_bin_size"), type="numeric", default=1e6,
              help="Bin size (bp) for CMplot density calculation [default=%default]", metavar="NUM"),
  make_option(c("--CMplot_col"), type="character", default="darkgreen,yellow,red",
              help="Comma-separated colors for CMplot density gradient [default=%default]", metavar="COLORS"),
  make_option(c("--CMplot_dpi"), type="integer", default=300,
              help="DPI for CMplot raster output [default=%default]", metavar="INT"),
  make_option(c("--CMplot_width"), type="numeric", default=9,
              help="Width in inches for CMplot output [default=%default]", metavar="NUM"),
  make_option(c("--CMplot_height"), type="numeric", default=6,
              help="Height in inches for CMplot output [default=%default]", metavar="NUM"),
  make_option(c("--CMplot_main"), type="character", default="Variant Density Plot",
              help="Title for CMplot density plot [default=%default]", metavar="TITLE")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="\nPlot variant positions from GATK VariantsToTable output\n\nExample usage:\n  Rscript vcf_difplot.R -i input.table -b sample1 -c sample2 -o output.pdf")
opt <- parse_args(opt_parser)

# ---------------------------------------------------------------------------
# Interactive mode
# ---------------------------------------------------------------------------
# Prompts the user for each parameter with a description.
# Required items loop until valid input is provided.
# Optional items accept Enter to keep the default.
run_interactive_mode <- function() {
  cat("============================================================\n")
  cat("  VCF Difplot -- Interactive Parameter Setup\n")
  cat("  Press Enter to accept the default value shown in [brackets].\n")
  cat("  Required parameters (*) must receive a non-empty value.\n")
  cat("============================================================\n\n")

  # Helper: read a line from the user's keyboard using bash's "read -e", which
  # enables GNU readline editing including Tab file-path completion, arrow-key
  # history navigation, and Ctrl+C / Ctrl+D exit handling.
  #
  # The multi-line label/description is printed via cat() first; only the
  # trailing "> " prompt is handled by bash so that readline knows where to
  # draw the cursor.  The user's input is captured through a temp file so that
  # R's stdout pipe is not disturbed.
  read_line <- function(prompt) {
    # Split prompt: print everything up to (not including) the final "> "
    # via cat, then let bash's "read -e" display "    > " itself.
    prompt_display <- sub(" > $", "", prompt)
    cat(prompt_display)
    flush(stdout())

    tmp <- tempfile(fileext = ".readline")
    on.exit(unlink(tmp), add = TRUE)

    ret <- system(
      sprintf("bash -c %s",
        shQuote(paste0(
          "read -e -r -p '    > ' _rl_input 2>/dev/tty",
          " && printf '%s\\n' \"$_rl_input\" > ", shQuote(tmp)
        ))
      ),
      wait = TRUE
    )

    if (ret != 0L) {   # non-zero exit = Ctrl+C or bash not available
      cat("\n[!] Interrupted by user. Exiting.\n")
      quit(save = "no", status = 1L, runLast = FALSE)
    }
    if (!file.exists(tmp)) return("")
    lines <- readLines(tmp, warn = FALSE)
    if (length(lines) == 0L) return("")
    trimws(lines[[1L]])
  }

  # Helper: prompt for a required character/path value; loops until non-empty
  prompt_required <- function(label, description) {
    repeat {
      val <- read_line(paste0("(*) ", label, "\n    ", description, "\n    > "))
      if (nchar(val) > 0) return(val)
      cat("    [!] This field is required. Please enter a value.\n")
    }
  }

  # Helper: prompt for an optional value; returns default when empty
  prompt_optional <- function(label, description, default_val) {
    disp <- if (is.null(default_val)) "none" else as.character(default_val)
    val <- read_line(paste0("( ) ", label, "\n    ", description,
                            "\n    [default: ", disp, "] > "))
    if (nchar(val) == 0) return(default_val)
    val
  }

  # Helper: prompt for an optional yes/no (logical) value; default is FALSE
  prompt_yesno <- function(label, description, default_val = FALSE) {
    disp <- if (default_val) "yes" else "no"
    repeat {
      val <- tolower(read_line(paste0("( ) ", label, "\n    ", description,
                                      "\n    [default: ", disp, "] (yes/no) > ")))
      if (nchar(val) == 0) return(default_val)
      if (val %in% c("y", "yes", "true",  "1")) return(TRUE)
      if (val %in% c("n", "no",  "false", "0")) return(FALSE)
      cat("    [!] Please enter yes or no.\n")
    }
  }

  # Helper: parse and validate a numeric value from a string
  parse_numeric_input <- function(s) {
    v <- suppressWarnings(as.numeric(s))
    if (is.na(v)) return(NULL)
    v
  }

  result <- list(interactive = TRUE)

  # --- Required: input file ---
  available_samples <- character(0)   # sample names read from file header
  repeat {
    result$input <- prompt_required(
      "Input file (-i / --input)  [REQUIRED]",
      "Tab-delimited file produced by GATK VariantsToTable (must exist)."
    )
    if (!file.exists(result$input)) {
      cat("    [!] File not found:", result$input, "-- please try again.\n")
      next
    }
    # Read header line only to extract sample names; support plain-text and gzip
    hdr <- tryCatch({
      is_gz <- grepl("\\.gz$", result$input, ignore.case = TRUE)
      con   <- if (is_gz) gzfile(result$input, "rt") else file(result$input, "rt")
      on.exit(try(close(con), silent = TRUE), add = TRUE)
      h <- readLines(con, n = 1L, warn = FALSE)
      close(con)
      if (length(h) == 0L) character(0L) else strsplit(h, "\t")[[1L]]
    }, error = function(e) character(0L))
    gt_hdr <- grep("\\.GT$", hdr, value = TRUE)
    if (length(gt_hdr) == 0L) {
      cat("    [!] No '.GT' columns found in that file.\n")
      cat("    [!] The input must be a GATK VariantsToTable output (tab-delimited),\n")
      cat("    [!] not a raw VCF file. Please try again.\n")
      next
    }
    available_samples <- sub("\\.GT$", "", gt_hdr)
    cat("    [i] Found", length(available_samples), "sample(s):",
        paste(available_samples, collapse = ", "), "\n")
    break
  }

  # --- Optional: output file ---
  result$output <- prompt_optional(
    "Output plot file (-o / --output)",
    "Path for the output PDF/PNG/JPG plot.",
    "variant_plot.pdf"
  )

  # --- Baseline sample (name OR column) ---
  sample_hint <- if (length(available_samples) > 0L)
    paste0("Available samples: ", paste(available_samples, collapse = ", "), ".")
  else
    "Name of the baseline sample as it appears in the column header (e.g. sample1)."

  cat("\n  Baseline sample identification -- provide EITHER a name OR a column number.\n")
  result$basename <- prompt_optional(
    "Baseline sample name (-b / --basename)",
    sample_hint,
    NULL
  )
  if (is.null(result$basename) || nchar(result$basename) == 0) {
    result$basename <- NULL
    repeat {
      raw <- prompt_optional(
        "Baseline sample column position (-B / --basecol)",
        "1-based index of the baseline GT column among all GT columns.",
        NULL
      )
      if (is.null(raw)) {
        cat("    [!] You must specify either a baseline sample name or column position.\n")
        next
      }
      v <- suppressWarnings(as.integer(raw))
      if (!is.na(v) && v >= 1) { result$basecol <- v; break }
      cat("    [!] Please enter a positive integer.\n")
    }
  } else {
    result$basecol <- NULL
  }

  # --- Comparison sample (name OR column) ---
  cat("\n  Comparison sample identification -- provide EITHER a name OR a column number.\n")
  result$copname <- prompt_optional(
    "Comparison sample name (-c / --copname)",
    sample_hint,
    NULL
  )
  if (is.null(result$copname) || nchar(result$copname) == 0) {
    result$copname <- NULL
    repeat {
      raw <- prompt_optional(
        "Comparison sample column position (-C / --copcol)",
        "1-based index of the comparison GT column among all GT columns.",
        NULL
      )
      if (is.null(raw)) {
        cat("    [!] You must specify either a comparison sample name or column position.\n")
        next
      }
      v <- suppressWarnings(as.integer(raw))
      if (!is.na(v) && v >= 1) { result$copcol <- v; break }
      cat("    [!] Please enter a positive integer.\n")
    }
  } else {
    result$copcol <- NULL
  }

  # --- Optional: chromosome length file ---
  cat("\n")
  result$chrlength <- prompt_optional(
    "Chromosome length file (-l / --chrlength)",
    "Two-column file (CHROM LENGTH). If omitted, max variant position is used per chromosome.",
    NULL
  )
  if (!is.null(result$chrlength) && nchar(result$chrlength) == 0) result$chrlength <- NULL

  # --- Optional: unit ---
  repeat {
    raw <- prompt_optional(
      "Chromosome length unit (-u / --unit)",
      "Numeric divisor for position axis (e.g. 1e6 = Mb, 1e3 = kb, 1 = bp).",
      "1e6"
    )
    v <- parse_numeric_input(if (is.null(raw)) "1e6" else raw)
    if (!is.null(v) && v > 0) { result$unit <- v; break }
    cat("    [!] Please enter a positive number.\n")
  }

  # --- Optional: homozygosity filters ---
  cat("\n")
  result$baseHetcheck <- prompt_yesno(
    "Baseline homozygosity filter (--baseHetcheck)",
    "Ignore positions where the baseline genotype is heterozygous.",
    FALSE
  )
  result$copHetcheck <- prompt_yesno(
    "Comparison homozygosity filter (--copHetcheck)",
    "Ignore positions where the comparison genotype is heterozygous.",
    FALSE
  )

  # --- Optional: plot aesthetics ---
  cat("\n  Plot aesthetics (press Enter to keep defaults):\n")
  result$segmentColor <- prompt_optional(
    "Variant segment color (--segmentColor)",
    "Any R color name or hex code for variant position marks.",
    "red"
  )

  repeat {
    raw <- prompt_optional(
      "Variant segment thickness (--segmentSize)",
      "Line width for variant marks (e.g. 0.5).",
      "0.5"
    )
    v <- parse_numeric_input(if (is.null(raw)) "0.5" else raw)
    if (!is.null(v) && v > 0) { result$segmentSize <- v; break }
    cat("    [!] Please enter a positive number.\n")
  }

  result$chrBorderColor <- prompt_optional(
    "Chromosome border color (--chrBorderColor)",
    "Any R color name or hex code for chromosome outlines.",
    "black"
  )

  repeat {
    raw <- prompt_optional(
      "Chromosome border thickness (--chrBorderSize)",
      "Line width for chromosome outlines (e.g. 0.3).",
      "0.3"
    )
    v <- parse_numeric_input(if (is.null(raw)) "0.3" else raw)
    if (!is.null(v) && v > 0) { result$chrBorderSize <- v; break }
    cat("    [!] Please enter a positive number.\n")
  }

  # --- Optional: output table ---
  cat("\n")
  result$output_table <- prompt_optional(
    "Variant position output table (--output_table)",
    "Optional path to write a tab-delimited CHROM/POS table of plotted variants.",
    NULL
  )
  if (!is.null(result$output_table) && nchar(result$output_table) == 0) result$output_table <- NULL

  # --- Optional: CMplot ---
  cat("\n")
  result$CMplot <- prompt_yesno(
    "Use CMplot for density plot (--CMplot)",
    "Use the CMplot R package to draw a density plot instead of ggplot2.",
    FALSE
  )
  if (isTRUE(result$CMplot)) {
    repeat {
      raw <- prompt_optional(
        "CMplot bin size (--CMplot_bin_size)",
        "Genomic window size in bp for density calculation (e.g. 1e6).",
        "1e6"
      )
      v <- parse_numeric_input(if (is.null(raw)) "1e6" else raw)
      if (!is.null(v) && v > 0) { result$CMplot_bin_size <- v; break }
      cat("    [!] Please enter a positive number.\n")
    }
    result$CMplot_col <- prompt_optional(
      "CMplot density colors (--CMplot_col)",
      "Comma-separated color names for the density gradient (low to high).",
      "darkgreen,yellow,red"
    )
    repeat {
      raw <- prompt_optional(
        "CMplot DPI (--CMplot_dpi)",
        "Resolution for raster output files (e.g. 300).",
        "300"
      )
      v <- suppressWarnings(as.integer(raw))
      if (!is.na(v) && v > 0) { result$CMplot_dpi <- v; break }
      cat("    [!] Please enter a positive integer.\n")
    }
    repeat {
      raw <- prompt_optional("CMplot width (--CMplot_width)", "Plot width in inches.", "9")
      v <- parse_numeric_input(if (is.null(raw)) "9" else raw)
      if (!is.null(v) && v > 0) { result$CMplot_width <- v; break }
      cat("    [!] Please enter a positive number.\n")
    }
    repeat {
      raw <- prompt_optional("CMplot height (--CMplot_height)", "Plot height in inches.", "6")
      v <- parse_numeric_input(if (is.null(raw)) "6" else raw)
      if (!is.null(v) && v > 0) { result$CMplot_height <- v; break }
      cat("    [!] Please enter a positive number.\n")
    }
    result$CMplot_main <- prompt_optional(
      "CMplot title (--CMplot_main)",
      "Title displayed on the density plot.",
      "Variant Density Plot"
    )
  } else {
    result$CMplot_bin_size <- 1e6
    result$CMplot_col      <- "darkgreen,yellow,red"
    result$CMplot_dpi      <- 300L
    result$CMplot_width    <- 9
    result$CMplot_height   <- 6
    result$CMplot_main     <- "Variant Density Plot"
  }

  # Deprecated field kept for compatibility
  result$threads <- 1L

  # -----------------------------------------------------------------------
  # Build and display the equivalent non-interactive command so the user
  # can reproduce this run directly from the command line next time.
  # -----------------------------------------------------------------------
  script_path <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0L) sub("^--file=", "", file_arg[1L]) else "vcf_difplot.R"
  }, error = function(e) "vcf_difplot.R")

  cmd <- c("Rscript", shQuote(script_path))
  cmd <- c(cmd, "-i", shQuote(result$input))
  cmd <- c(cmd, "-o", shQuote(result$output))
  if (!is.null(result$basename))  cmd <- c(cmd, "-b", shQuote(result$basename))
  if (!is.null(result$basecol))   cmd <- c(cmd, "-B", result$basecol)
  if (!is.null(result$copname))   cmd <- c(cmd, "-c", shQuote(result$copname))
  if (!is.null(result$copcol))    cmd <- c(cmd, "-C", result$copcol)
  if (!is.null(result$chrlength)) cmd <- c(cmd, "-l", shQuote(result$chrlength))
  if (!is.null(result$unit))      cmd <- c(cmd, "-u", result$unit)
  if (isTRUE(result$baseHetcheck)) cmd <- c(cmd, "--baseHetcheck")
  if (isTRUE(result$copHetcheck))  cmd <- c(cmd, "--copHetcheck")
  cmd <- c(cmd, "--segmentColor",   shQuote(result$segmentColor))
  cmd <- c(cmd, "--segmentSize",    result$segmentSize)
  cmd <- c(cmd, "--chrBorderColor", shQuote(result$chrBorderColor))
  cmd <- c(cmd, "--chrBorderSize",  result$chrBorderSize)
  if (!is.null(result$output_table)) cmd <- c(cmd, "--output_table", shQuote(result$output_table))
  if (isTRUE(result$CMplot)) {
    cmd <- c(cmd, "--CMplot")
    cmd <- c(cmd, "--CMplot_bin_size", result$CMplot_bin_size)
    cmd <- c(cmd, "--CMplot_col",      shQuote(result$CMplot_col))
    cmd <- c(cmd, "--CMplot_dpi",      result$CMplot_dpi)
    cmd <- c(cmd, "--CMplot_width",    result$CMplot_width)
    cmd <- c(cmd, "--CMplot_height",   result$CMplot_height)
    cmd <- c(cmd, "--CMplot_main",     shQuote(result$CMplot_main))
  }

  cat("\n============================================================\n")
  cat("  Parameters confirmed. Starting analysis...\n")
  cat("  Equivalent command:\n\n")
  cat("    ", paste(cmd, collapse = " "), "\n")
  cat("============================================================\n\n")

  result
}

# If --interactive (-I) was requested, replace opt with interactively gathered values
if (isTRUE(opt$interactive)) {
  opt <- run_interactive_mode()
}

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
message("Reading input file: ", opt$input)
data <- tryCatch({
  as.data.frame(fread(opt$input, sep="\t", header=TRUE))
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

message("Found ", length(gt_cols), " samples with GT information")
sample_names <- sub("\\.GT$", "", gt_cols)
message("Sample names: ", paste(sample_names, collapse=", "))

# Determine baseline sample column
if (!is.null(opt$basename)) {
  base_col <- paste0(opt$basename, ".GT")
  if (!base_col %in% colnames(data)) {
    stop(paste("Baseline sample not found:", opt$basename, "\nAvailable samples:", paste(sample_names, collapse=", ")), call.=FALSE)
  }
  message("Using baseline sample: ", opt$basename)
} else {
  # Use column position
  if (opt$basecol < 1 || opt$basecol > length(gt_cols)) {
    stop(paste("Baseline column position out of range. Must be between 1 and", length(gt_cols)), call.=FALSE)
  }
  base_col <- gt_cols[opt$basecol]
  message("Using baseline sample from column ", opt$basecol, ": ", sample_names[opt$basecol])
}

# Determine comparison sample column
if (!is.null(opt$copname)) {
  comp_col <- paste0(opt$copname, ".GT")
  if (!comp_col %in% colnames(data)) {
    stop(paste("Comparison sample not found:", opt$copname, "\nAvailable samples:", paste(sample_names, collapse=", ")), call.=FALSE)
  }
  message("Using comparison sample: ", opt$copname)
} else {
  # Use column position
  if (opt$copcol < 1 || opt$copcol > length(gt_cols)) {
    stop(paste("Comparison column position out of range. Must be between 1 and", length(gt_cols)), call.=FALSE)
  }
  comp_col <- gt_cols[opt$copcol]
  message("Using comparison sample from column ", opt$copcol, ": ", sample_names[opt$copcol])
}

# Vectorized genotype helpers

# Return a logical vector: TRUE where a diploid GT string is missing or contains
# wildcard alleles and should be excluded from analysis.
gt_is_missing <- function(gt) {
  gt_norm <- gsub("\\|", "/", gt)
  allele1 <- sub("/.*$", "", gt_norm)
  allele2 <- sub("^[^/]*/", "", gt_norm)
  is.na(gt) | gt == "" |
    allele1 == "." | allele2 == "." |
    allele1 == "*" | allele2 == "*" |
    !grepl("/", gt_norm, fixed=TRUE)
}

# Normalize genotype: treat / and | as equivalent separators, sort alleles.
# Input: character vector of genotype strings like "A/T", "C|G", "./.", etc.
# Output: character vector of sorted alleles separated by "/" (NA for missing/wildcard)
normalize_gt_vec <- function(gt) {
  gt_norm <- gsub("\\|", "/", gt)
  allele1 <- sub("/.*$", "", gt_norm)
  allele2 <- sub("^[^/]*/", "", gt_norm)
  is_missing <- gt_is_missing(gt)
  swap <- !is_missing & allele2 < allele1
  a1 <- allele1
  a2 <- allele2
  a1[swap] <- allele2[swap]
  a2[swap] <- allele1[swap]
  result <- paste(a1, a2, sep="/")
  result[is_missing] <- NA
  result
}

# Check if genotype is homozygous (vectorized).
# Input: character vector of genotype strings
# Output: logical vector (NA for missing/wildcard genotypes)
is_homo_vec <- function(gt) {
  gt_norm <- gsub("\\|", "/", gt)
  allele1 <- sub("/.*$", "", gt_norm)
  allele2 <- sub("^[^/]*/", "", gt_norm)
  is_missing <- gt_is_missing(gt)
  result <- allele1 == allele2
  result[is_missing] <- NA
  result
}

# Compare genotypes and mark variants
message("Comparing genotypes...")

data$base_gt_norm <- normalize_gt_vec(data[[base_col]])
data$comp_gt_norm <- normalize_gt_vec(data[[comp_col]])

# Exclude positions with missing data (./.) or wildcards (*/*)
data$keep <- !is.na(data$base_gt_norm) & !is.na(data$comp_gt_norm)
message("Positions after removing missing data (./.) and wildcards (*/*): ", sum(data$keep))

# Apply baseline homozygosity check if requested
if (opt$baseHetcheck) {
  message("Applying baseline homozygosity check...")
  base_homo <- is_homo_vec(data[[base_col]])
  data$keep <- data$keep & !is.na(base_homo) & base_homo
  message("Positions after baseline homozygosity filter: ", sum(data$keep))
}

# Apply comparison homozygosity check if requested
if (opt$copHetcheck) {
  message("Applying comparison homozygosity check...")
  comp_homo <- is_homo_vec(data[[comp_col]])
  data$keep <- data$keep & !is.na(comp_homo) & comp_homo
  message("Positions after comparison homozygosity filter: ", sum(data$keep))
}

# Filter data based on all criteria
data <- data[data$keep, ]
gc()

# Compare normalized genotypes
data$is_variant <- data$base_gt_norm != data$comp_gt_norm

message("Total positions: ", nrow(data))
message("Variant positions: ", sum(data$is_variant))
message("Non-variant positions: ", sum(!data$is_variant))

# Print first 20 variant positions (where base and cop differ)
variant_data <- data[data$is_variant, ]
cat("\n=== First", MAX_DISPLAY_ROWS, "variant positions (base != cop) ===\n")
if (nrow(variant_data) > 0) {
  # Get first N rows (or all if less than N)
  n_display <- min(MAX_DISPLAY_ROWS, nrow(variant_data))
  
  # Create display data frame with desired column names
  display_subset <- head(variant_data[, c("CHROM", "POS", base_col, comp_col)], n_display)
  colnames(display_subset) <- c("CHROM", "POS", "Baseline_GT", "Comparison_GT")
  
  # Print as a formatted table
  print(display_subset, row.names=FALSE)
  
  if (nrow(variant_data) > MAX_DISPLAY_ROWS) {
    cat("\n... (showing", MAX_DISPLAY_ROWS, "of", nrow(variant_data), "total variant positions)\n")
  }
} else {
  cat("No variant positions found (all positions have matching genotypes).\n")
}
cat("========================================================\n\n")

# Get chromosome information
chromosomes <- unique(data$CHROM)
message("Chromosomes found: ", paste(chromosomes, collapse=", "))

# Determine chromosome lengths
if (!is.null(opt$chrlength)) {
  # Read chromosome length file
  if (!file.exists(opt$chrlength)) {
    stop(paste("Chromosome length file does not exist:", opt$chrlength), call.=FALSE)
  }
  
  message("Reading chromosome length file: ", opt$chrlength)
  chr_lengths <- tryCatch({
    as.data.frame(fread(opt$chrlength, header=FALSE, col.names=c("CHROM", "LENGTH")))
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

# Sort chromosomes using a key that handles numeric names, X, Y, MT and other
# non-numeric names in conventional genomic order.
chr_sort_key <- function(chrom) {
  # Strip leading "chr"/"Chr"/"CHR" prefix for comparison
  stripped <- sub("^[Cc][Hh][Rr]", "", chrom)
  num <- suppressWarnings(as.numeric(stripped))
  ifelse(!is.na(num), num,
    ifelse(toupper(stripped) == "X",  10000L,
    ifelse(toupper(stripped) == "Y",  10001L,
    ifelse(toupper(stripped) %in% c("MT", "M"), 10002L,
    10003L  # everything else: secondary sort by chrom name (alphabetical) via order()
  ))))
}
chr_info <- chr_info[order(chr_sort_key(chr_info$CHROM), chr_info$CHROM), ]
chr_info$chr_order <- seq_len(nrow(chr_info))

# Scale lengths by unit
chr_info$LENGTH_scaled <- chr_info$LENGTH / opt$unit

cat("\nChromosome lengths (scaled by", opt$unit, "):\n")
print(chr_info[, c("CHROM", "LENGTH_scaled")])

# Prepare data for plotting
plot_data <- merge(data, chr_info[, c("CHROM", "chr_order", "LENGTH_scaled")], by="CHROM")
plot_data$POS_scaled <- plot_data$POS / opt$unit

# Filter for variants only
variant_data <- plot_data[plot_data$is_variant, ]

# Write output table if requested
if (!is.null(opt$output_table)) {
  message("Writing variant position table to: ", opt$output_table)
  write.table(variant_data[, c("CHROM", "POS")], file=opt$output_table,
              sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  message("Variant position table saved: ", nrow(variant_data), " positions written")
}

# Create plot
message("\nGenerating plot...")

# Guard: if no variant positions exist, generate an empty chromosome frame and warn
if (nrow(variant_data) == 0) {
  message("Warning: No variant positions found. Generating empty chromosome frame plot.")
}

if (isTRUE(opt$CMplot)) {
  # ---------------------------------------------------------------------------
  # CMplot density plot branch
  # ---------------------------------------------------------------------------
  suppressPackageStartupMessages(library(CMplot))

  # Build a data frame in CMplot format: SNP, Chromosome, Position
  cmplot_data <- data.frame(
    SNP        = paste(variant_data$CHROM, variant_data$POS, sep = "_"),
    Chromosome = variant_data$CHROM,
    Position   = variant_data$POS,
    stringsAsFactors = FALSE
  )

  # Parse comma-separated color string into a character vector
  chr_den_col <- trimws(strsplit(opt$CMplot_col, ",")[[1]])

  # Determine output format from file extension
  output_ext <- tolower(tools::file_ext(opt$output))
  if (!output_ext %in% c("pdf", "png", "jpg", "jpeg")) {
    warning(paste("Unsupported output format:", output_ext, ". Using pdf instead."))
    output_ext <- "pdf"
    opt$output <- paste0(tools::file_path_sans_ext(opt$output), ".pdf")
  }
  # CMplot uses "jpg" not "jpeg"
  if (output_ext == "jpeg") output_ext <- "jpg"

  # CMplot saves to the working directory; temporarily switch to the output dir
  output_dir  <- dirname(opt$output)
  output_base <- tools::file_path_sans_ext(basename(opt$output))
  if (!nzchar(output_dir)) output_dir <- "."

  old_wd <- getwd()
  tryCatch(setwd(output_dir), error = function(e)
    stop(paste("Cannot change to output directory:", output_dir), call. = FALSE))
  on.exit(setwd(old_wd), add = TRUE)

  message("Saving CMplot density plot to directory: ", normalizePath(output_dir))
  message("Output file name prefix: ", output_base)

  CMplot(cmplot_data,
         plot.type    = "d",
         bin.size     = opt$CMplot_bin_size,
         chr.den.col  = chr_den_col,
         file         = output_ext,
         file.name    = output_base,
         dpi          = opt$CMplot_dpi,
         main         = opt$CMplot_main,
         file.output  = TRUE,
         verbose      = TRUE,
         width        = opt$CMplot_width,
         height       = opt$CMplot_height)

  setwd(old_wd)
  message("\nDone! CMplot output saved.")
  message("Total variants used for density plot: ", nrow(variant_data))

} else {
  # ---------------------------------------------------------------------------
  # Default ggplot2 plotting branch
  # ---------------------------------------------------------------------------

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
  message("Saving plot to: ", opt$output)

  # Determine output format based on file extension
  output_ext <- tolower(tools::file_ext(opt$output))
  plot_height_inches <- min(MAX_PLOT_HEIGHT_IN, max(4, nrow(chr_info) * 0.5))
  plot_height_pixels <- min(MAX_PLOT_HEIGHT_IN * 100L, max(400L, nrow(chr_info) * 50L))
  if (output_ext == "pdf") {
    pdf(opt$output, width=12, height=plot_height_inches)
  } else if (output_ext == "png") {
    png(opt$output, width=1200, height=plot_height_pixels, res=100)
  } else if (output_ext %in% c("jpg", "jpeg")) {
    jpeg(opt$output, width=1200, height=plot_height_pixels, res=100)
  } else {
    # Default to PDF
    warning(paste("Unsupported output format:", output_ext, ". Using PDF instead."))
    opt$output <- sub(paste0("\\.", output_ext, "$"), ".pdf", opt$output)
    pdf(opt$output, width=12, height=plot_height_inches)
  }

  print(p)
  if (dev.cur() > 1) dev.off()

  message("\nDone! Plot saved to: ", opt$output)
  message("Total variants plotted: ", nrow(variant_data))
}
