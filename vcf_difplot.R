#!/usr/bin/env Rscript

# VCF Variant Position Plotting Script
# Reads a tab-delimited file from GATK VariantsToTable and plots variant
# positions where two samples' genotypes differ.

# --- Load libraries and early version checks --------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
  library(data.table)
})

if (packageVersion("ggplot2") < "3.4.0") {
  stop("This script requires ggplot2 >= 3.4.0. Update: install.packages('ggplot2')", call. = FALSE)
}

# --- Constants ---------------------------------------------------------------
MAX_DISPLAY_ROWS <- 20

# --- Command-line options ----------------------------------------------------
option_list <- list(
  make_option(c("-I", "--interactive"), action="store_true", default=FALSE,
              help="Enter interactive mode: prompts for each parameter"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input tab-delimited file (required)", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="variant_plot.pdf",
              help="Output plot file (pdf/png/jpg/svg) [default=%default]", metavar="FILE"),
  make_option(c("-b", "--basename"), type="character", default=NULL,
              help="Baseline sample name", metavar="NAME"),
  make_option(c("-B", "--basecol"), type="integer", default=NULL,
              help="Baseline GT column index (1-based, among GT columns only)", metavar="INT"),
  make_option(c("-c", "--copname"), type="character", default=NULL,
              help="Comparison sample name", metavar="NAME"),
  make_option(c("-C", "--copcol"), type="integer", default=NULL,
              help="Comparison GT column index (1-based, among GT columns only)", metavar="INT"),
  make_option(c("-l", "--chrlength"), type="character", default=NULL,
              help="Chromosome length file (2+ cols; first two used: CHROM LENGTH). Supports .fai", metavar="FILE"),
  make_option(c("-u", "--unit"), type="numeric", default=1e6,
              help="Position unit divisor (1e6=Mb, 1e3=kb, 1=bp) [default=%default]", metavar="NUM"),
  make_option(c("--baseHetcheck"), action="store_true", default=FALSE,
              help="Only include positions where baseline is homozygous"),
  make_option(c("--copHetcheck"), action="store_true", default=FALSE,
              help="Only include positions where comparison is homozygous"),
  make_option(c("--segmentColor"), type="character", default="red",
              help="Color for variant segments [default=%default]", metavar="COLOR"),
  make_option(c("--segmentSize"), type="numeric", default=0.5,
              help="Thickness of variant segments [default=%default]", metavar="NUM"),
  make_option(c("--chrBorderColor"), type="character", default="black",
              help="Color for chromosome borders [default=%default]", metavar="COLOR"),
  make_option(c("--chrBorderSize"), type="numeric", default=0.3,
              help="Thickness of chromosome borders [default=%default]", metavar="NUM"),
  make_option(c("--segmentAlpha"), type="numeric", default=0.6,
              help="Alpha transparency for variant segments (0-1) [default=%default]", metavar="NUM"),
  make_option(c("--output_table"), type="character", default=NULL,
              help="Write variant positions (CHROM POS) to this file", metavar="FILE"),
  # CMplot SNP-density plot options
  make_option(c("--cmplot"), action="store_true", default=FALSE,
              help="Additionally generate a CMplot SNP-density plot (requires CMplot package)"),
  make_option(c("--cmplot_bin_size"), type="numeric", default=1e6,
              help="CMplot density window size in bp [default=%default]", metavar="NUM"),
  make_option(c("--cmplot_den_col"), type="character", default="darkgreen,yellow,red",
              help="CMplot density gradient colors, comma-separated [default=%default]", metavar="COLORS"),
  make_option(c("--cmplot_dpi"), type="integer", default=300,
              help="CMplot output DPI [default=%default]", metavar="INT"),
  make_option(c("--cmplot_width"), type="numeric", default=NULL,
              help="CMplot figure width in inches (default: auto)", metavar="NUM"),
  make_option(c("--cmplot_height"), type="numeric", default=NULL,
              help="CMplot figure height in inches (default: auto)", metavar="NUM"),
  make_option(c("--cmplot_main"), type="character", default=NULL,
              help="CMplot figure title (default: auto-generated)", metavar="TEXT")
)

opt_parser <- OptionParser(option_list=option_list,
  description="\nPlot variant positions from GATK VariantsToTable output\n\nExample:\n  Rscript vcf_difplot.R -i input.table -b sample1 -c sample2 -o output.pdf")
opt <- parse_args(opt_parser)

# Conditionally load CMplot if requested
if (isTRUE(opt$cmplot)) {
  if (!requireNamespace("CMplot", quietly = TRUE)) {
    stop("--cmplot requires the CMplot package. Install: install.packages('CMplot')", call. = FALSE)
  }
  suppressPackageStartupMessages(library(CMplot))
}

# --- Helper functions --------------------------------------------------------

# Read column names from a tab-delimited file (plain or gzipped).
# fread auto-detects compression (requires R.utils for .gz).
read_header_cols <- function(filepath) {
  names(fread(filepath, nrows = 0L, sep = "\t"))
}

# Resolve a sample specification (name OR column index) to a GT column name.
resolve_sample <- function(name, col, gt_cols, sample_names, label) {
  if (!is.null(name)) {
    if (!is.null(col)) warning("Both ", label, " name and column specified. Using name.")
    idx <- match(name, sample_names)
    if (is.na(idx)) {
      stop(label, " sample not found: '", name, "'\nAvailable: ",
           paste(sample_names, collapse = ", "), call. = FALSE)
    }
  } else if (!is.null(col)) {
    if (col < 1 || col > length(gt_cols)) {
      stop(label, " GT column index out of range (1-", length(gt_cols), ")", call. = FALSE)
    }
    idx <- col
  } else {
    stop("Specify ", label, " sample name or GT column index", call. = FALSE)
  }
  message("Using ", label, " sample: ", sample_names[idx])
  gt_cols[idx]
}

# Unified genotype parser. Processes only unique GT values for speed, then maps
# back via match(). Returns list(norm, homo, ploidy, missing) aligned to input.
parse_gt <- function(gt) {
  gt <- as.character(gt)
  u <- unique(gt)
  s <- gsub("|", "/", u, fixed = TRUE)
  al <- strsplit(s, "/", fixed = TRUE)

  bad <- vapply(al, function(a) {
    length(a) == 0L || any(a == "") || any(a == ".") || any(a == "*")
  }, logical(1L), USE.NAMES = FALSE)
  bad <- bad | is.na(u) | u == ""

  norm <- vapply(al, function(a) {
    if (length(a) == 1L) a <- c(a, a)
    paste(sort(a), collapse = "/")
  }, character(1L), USE.NAMES = FALSE)
  norm[bad] <- NA_character_

  homo <- vapply(al, function(a) length(unique(a)) == 1L, logical(1L), USE.NAMES = FALSE)
  homo[bad] <- NA

  ploidy <- lengths(al)
  ploidy[bad] <- NA_integer_

  i <- match(gt, u)
  list(norm = norm[i], homo = homo[i], ploidy = ploidy[i], missing = bad[i])
}

# Chromosome sort key: numeric first, then X, Y, MT/M, then others.
chr_sort_key <- function(chrom) {
  s <- toupper(sub("^[Cc][Hh][Rr]", "", chrom))
  num <- suppressWarnings(as.numeric(s))
  special <- c(X = 1e4, Y = 1e4 + 1, M = 1e4 + 2, MT = 1e4 + 2)
  ifelse(!is.na(num), num, ifelse(s %in% names(special), special[s], 1e4 + 3))
}

# Validate a color string (case-insensitive, accepts R names and hex).
validate_color <- function(color_str, param_name) {
  if (is.null(color_str)) return(invisible(NULL))
  ok <- tryCatch({ col2rgb(color_str); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) {
    stop("Invalid color for ", param_name, ": '", color_str,
         "'. Use an R color name or hex code (e.g. '#FF5733').", call. = FALSE)
  }
}

# Validate output path/extension early. Returns list(path, ext).
validate_output <- function(output_path) {
  ext <- tolower(tools::file_ext(output_path))
  supported <- c("pdf", "png", "jpg", "jpeg", "svg")
  if (ext == "") {
    warning("Output filename has no extension. Appending '.pdf'.")
    return(list(path = paste0(output_path, ".pdf"), ext = "pdf"))
  } else if (!ext %in% supported) {
    warning("Unsupported output format: '", ext, "'. Using PDF instead.")
    return(list(path = sub("\\.[^.]*$", ".pdf", output_path), ext = "pdf"))
  }
  list(path = output_path, ext = ext)
}

# --- Interactive mode --------------------------------------------------------
run_interactive_mode <- function() {
  cat("============================================================\n")
  cat("  VCF Difplot -- Interactive Parameter Setup\n")
  cat("  Press Enter to accept the default value shown in [brackets].\n")
  cat("  Required parameters (*) must receive a non-empty value.\n")
  cat("============================================================\n\n")

  # On Unix with bash+TTY: bash readline for Tab completion / history.
  # Otherwise: readLines(file("stdin")) works reliably under Rscript.
  has_bash <- Sys.which("bash") != "" && .Platform$OS.type != "windows"
  has_tty <- if (has_bash) system("test -t 0", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0L else FALSE
  use_bash <- has_bash && has_tty

  read_line <- function(prompt) {
    cat(sub("\\s*> $", "", prompt))
    flush(stdout())

    if (use_bash) {
      tmp <- tempfile(fileext = ".readline")
      on.exit(unlink(tmp), add = TRUE)
      ret <- system(sprintf("bash -c %s", shQuote(paste0(
        "read -e -r -p '    > ' _rl_input 2>/dev/tty",
        " && printf '%s\\n' \"$_rl_input\" > ", shQuote(tmp)))), wait = TRUE)
      if (ret != 0L) { cat("\n[!] Interrupted. Exiting.\n"); quit(save = "no", status = 1L) }
      if (!file.exists(tmp)) return("")
      lines <- readLines(tmp, warn = FALSE)
      if (length(lines) == 0L) return("")
      return(trimws(lines[[1L]]))
    }

    # file("stdin") is correct under Rscript; stdin() refers to the script itself.
    val <- tryCatch(readLines(file("stdin"), n = 1L, warn = FALSE), error = function(e) character(0L))
    if (length(val) == 0L) { cat("\n[!] EOF. Exiting.\n"); quit(save = "no", status = 1L) }
    trimws(val[[1L]])
  }

  prompt_required <- function(label, desc) {
    repeat {
      val <- read_line(paste0("(*) ", label, "\n    ", desc, "\n    > "))
      if (nchar(val) > 0) return(val)
      cat("    [!] Required. Please enter a value.\n")
    }
  }
  prompt_optional <- function(label, desc, def) {
    disp <- if (is.null(def)) "none" else as.character(def)
    val <- read_line(paste0("( ) ", label, "\n    ", desc, "\n    [default: ", disp, "] > "))
    if (nchar(val) == 0) return(def)
    val
  }
  prompt_yesno <- function(label, desc, def = FALSE) {
    repeat {
      val <- tolower(read_line(paste0("( ) ", label, "\n    ", desc,
                                      "\n    [default: ", if (def) "yes" else "no", "] (y/n) > ")))
      if (nchar(val) == 0) return(def)
      if (val %in% c("y", "yes", "1")) return(TRUE)
      if (val %in% c("n", "no", "0")) return(FALSE)
      cat("    [!] Enter yes or no.\n")
    }
  }
  prompt_positive <- function(label, desc, def) {
    repeat {
      v <- suppressWarnings(as.numeric(prompt_optional(label, desc, as.character(def))))
      if (!is.na(v) && v > 0) return(v)
      cat("    [!] Enter a positive number.\n")
    }
  }
  prompt_sample <- function(label, flag_n, flag_c, hint) {
    cat("\n  ", label, " sample -- name OR column number.\n")
    nm <- prompt_optional(paste0(label, " name (", flag_n, ")"), hint, NULL)
    if (!is.null(nm) && nchar(nm) > 0) return(list(name = nm, col = NULL))
    repeat {
      raw <- prompt_optional(paste0(label, " GT col index (", flag_c, ")"),
                             "1-based index among GT columns.", NULL)
      if (is.null(raw)) { cat("    [!] Must specify name or index.\n"); next }
      v <- suppressWarnings(as.numeric(raw))
      if (!is.na(v) && v >= 1 && v == floor(v)) return(list(name = NULL, col = as.integer(v)))
      cat("    [!] Enter a positive integer.\n")
    }
  }

  result <- list(interactive = TRUE)
  samples <- character(0)

  repeat {
    result$input <- prompt_required("Input file (-i)  [REQUIRED]",
      "GATK VariantsToTable output (must exist).")
    if (!file.exists(result$input)) { cat("    [!] Not found:", result$input, "\n"); next }
    hdr <- tryCatch(read_header_cols(result$input), error = function(e) character(0L))
    gt_hdr <- grep("\\.GT$", hdr, value = TRUE)
    if (length(gt_hdr) == 0L) { cat("    [!] No .GT columns. Not a VariantsToTable file?\n"); next }
    samples <- sub("\\.GT$", "", gt_hdr)
    cat("    [i]", length(samples), "sample(s):", paste(samples, collapse = ", "), "\n")
    break
  }

  result$output <- prompt_optional("Output file (-o)", "pdf/png/jpg/svg.", "variant_plot.pdf")
  hint <- paste0("Available: ", paste(samples, collapse = ", "), ".")
  bs <- prompt_sample("Baseline", "-b", "-B", hint)
  result$basename <- bs$name; result$basecol <- bs$col
  cs <- prompt_sample("Comparison", "-c", "-C", hint)
  result$copname <- cs$name; result$copcol <- cs$col

  cat("\n")
  result$chrlength <- prompt_optional("Chr length file (-l)", "2+ cols (CHROM LENGTH). Supports .fai.", NULL)
  result$unit <- prompt_positive("Unit (-u)", "1e6=Mb, 1e3=kb, 1=bp.", 1e6)
  cat("\n")
  result$baseHetcheck <- prompt_yesno("Baseline homozygosity (--baseHetcheck)", "Skip het in baseline.", FALSE)
  result$copHetcheck <- prompt_yesno("Comparison homozygosity (--copHetcheck)", "Skip het in comparison.", FALSE)
  cat("\n  Aesthetics:\n")
  result$segmentColor <- prompt_optional("Segment color", "R name or hex.", "red")
  result$segmentSize <- prompt_positive("Segment size", "Line width.", 0.5)
  result$segmentAlpha <- prompt_positive("Segment alpha", "0-1.", 0.6)
  result$chrBorderColor <- prompt_optional("Border color", "R name or hex.", "black")
  result$chrBorderSize <- prompt_positive("Border size", "Line width.", 0.3)
  cat("\n")
  result$output_table <- prompt_optional("Output table (--output_table)", "CHROM/POS file path.", NULL)

  # CMplot density plot
  cat("\n")
  result$cmplot <- prompt_yesno("CMplot SNP-density plot (--cmplot)",
    "Additionally generate a CMplot density plot of variant positions.", FALSE)
  if (isTRUE(result$cmplot)) {
    result$cmplot_bin_size <- prompt_positive("  Density bin size (--cmplot_bin_size)", "Window in bp.", 1e6)
    result$cmplot_den_col <- prompt_optional("  Density colors (--cmplot_den_col)", "Comma-separated.", "darkgreen,yellow,red")
    result$cmplot_dpi <- prompt_positive("  DPI (--cmplot_dpi)", "Output resolution.", 300)
    result$cmplot_main <- prompt_optional("  Title (--cmplot_main)", "Plot title.", NULL)
    result$cmplot_width <- NULL
    result$cmplot_height <- NULL
  }

  # Equivalent command
  sp <- tryCatch({ a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
    if (length(f)) sub("^--file=", "", f[1L]) else "vcf_difplot.R" }, error = function(e) "vcf_difplot.R")
  cmd <- c("Rscript", shQuote(sp), "-i", shQuote(result$input), "-o", shQuote(result$output))
  if (!is.null(result$basename))  cmd <- c(cmd, "-b", shQuote(result$basename))
  if (!is.null(result$basecol))   cmd <- c(cmd, "-B", result$basecol)
  if (!is.null(result$copname))   cmd <- c(cmd, "-c", shQuote(result$copname))
  if (!is.null(result$copcol))    cmd <- c(cmd, "-C", result$copcol)
  if (!is.null(result$chrlength)) cmd <- c(cmd, "-l", shQuote(result$chrlength))
  cmd <- c(cmd, "-u", result$unit)
  if (isTRUE(result$baseHetcheck)) cmd <- c(cmd, "--baseHetcheck")
  if (isTRUE(result$copHetcheck))  cmd <- c(cmd, "--copHetcheck")
  cmd <- c(cmd, "--segmentColor", shQuote(result$segmentColor), "--segmentSize", result$segmentSize,
           "--segmentAlpha", result$segmentAlpha, "--chrBorderColor", shQuote(result$chrBorderColor),
           "--chrBorderSize", result$chrBorderSize)
  if (!is.null(result$output_table)) cmd <- c(cmd, "--output_table", shQuote(result$output_table))
  if (isTRUE(result$cmplot)) {
    cmd <- c(cmd, "--cmplot", "--cmplot_bin_size", result$cmplot_bin_size,
             "--cmplot_den_col", shQuote(result$cmplot_den_col), "--cmplot_dpi", result$cmplot_dpi)
    if (!is.null(result$cmplot_main)) cmd <- c(cmd, "--cmplot_main", shQuote(result$cmplot_main))
  }
  cat("\n============================================================\n")
  cat("  Equivalent command:\n\n    ", paste(cmd, collapse = " "), "\n")
  cat("============================================================\n\n")
  result
}

if (isTRUE(opt$interactive)) opt <- run_interactive_mode()

# --- Parameter validation (shared by CLI and interactive) --------------------
if (is.null(opt$input)) { print_help(opt_parser); stop("Input file required (-i)", call. = FALSE) }
if (!file.exists(opt$input)) stop("Input file not found: ", opt$input, call. = FALSE)
if (is.null(opt$basename) && is.null(opt$basecol))
  stop("Specify baseline sample (-b name or -B index)", call. = FALSE)
if (is.null(opt$copname) && is.null(opt$copcol))
  stop("Specify comparison sample (-c name or -C index)", call. = FALSE)
if (!is.numeric(opt$unit) || opt$unit <= 0) stop("--unit must be positive", call. = FALSE)
if (!is.numeric(opt$segmentSize) || opt$segmentSize <= 0) stop("--segmentSize must be positive", call. = FALSE)
if (!is.numeric(opt$chrBorderSize) || opt$chrBorderSize <= 0) stop("--chrBorderSize must be positive", call. = FALSE)
if (!is.numeric(opt$segmentAlpha) || opt$segmentAlpha < 0 || opt$segmentAlpha > 1)
  stop("--segmentAlpha must be 0-1", call. = FALSE)
validate_color(opt$segmentColor, "--segmentColor")
validate_color(opt$chrBorderColor, "--chrBorderColor")

# Validate CMplot parameters if enabled
if (isTRUE(opt$cmplot)) {
  if (!is.numeric(opt$cmplot_bin_size) || opt$cmplot_bin_size <= 0)
    stop("--cmplot_bin_size must be a positive number", call. = FALSE)
  if (!is.null(opt$cmplot_dpi) && (!is.numeric(opt$cmplot_dpi) || opt$cmplot_dpi < 1))
    stop("--cmplot_dpi must be a positive integer", call. = FALSE)
  # Parse density colors
  cmplot_den_cols <- trimws(strsplit(opt$cmplot_den_col, ",")[[1]])
  for (cc in cmplot_den_cols) validate_color(cc, "--cmplot_den_col")
}

# Validate output early (before expensive processing)
out_info <- validate_output(opt$output)
opt$output <- out_info$path
output_ext <- out_info$ext

# --- Read header, resolve samples --------------------------------------------
message("Reading header: ", opt$input)
all_cols <- tryCatch(read_header_cols(opt$input), error = function(e)
  stop("Error reading input: ", e$message, call. = FALSE))

gt_cols <- grep("\\.GT$", all_cols, value = TRUE)
if (length(gt_cols) == 0)
  stop("No GT columns found. Need GATK VariantsToTable output with 'sample.GT' columns.", call. = FALSE)
sample_names <- sub("\\.GT$", "", gt_cols)
message("Found ", length(gt_cols), " sample(s): ", paste(sample_names, collapse = ", "))

base_col <- resolve_sample(opt$basename, opt$basecol, gt_cols, sample_names, "baseline")
comp_col <- resolve_sample(opt$copname, opt$copcol, gt_cols, sample_names, "comparison")
if (base_col == comp_col)
  stop("Baseline and comparison are the SAME sample (", base_col, ").", call. = FALSE)

base_name <- sub("\\.GT$", "", base_col)
comp_name <- sub("\\.GT$", "", comp_col)

# --- Read data (only needed columns) -----------------------------------------
message("Reading data...")
data <- tryCatch(
  as.data.frame(fread(opt$input, sep = "\t", header = TRUE,
                      select = unique(c("CHROM", "POS", base_col, comp_col)),
                      colClasses = list(character = "CHROM"))),
  error = function(e) stop("Error reading input: ", e$message, call. = FALSE))

if (!all(c("CHROM", "POS") %in% colnames(data)))
  stop("Missing CHROM or POS column.", call. = FALSE)
message("Read ", nrow(data), " positions")

# --- Parse and compare genotypes ---------------------------------------------
message("Comparing genotypes...")
b <- parse_gt(data[[base_col]])
cmp <- parse_gt(data[[comp_col]])

valid_ploidy <- unique(c(b$ploidy[!is.na(b$ploidy)], cmp$ploidy[!is.na(cmp$ploidy)]))
if (any(valid_ploidy == 1L)) {
  n <- sum(b$ploidy == 1L | cmp$ploidy == 1L, na.rm = TRUE)
  warning(n, " position(s) with haploid genotype(s). Expanded to diploid (A -> A/A).")
}
if (any(valid_ploidy > 2L)) {
  n <- sum(b$ploidy > 2L | cmp$ploidy > 2L, na.rm = TRUE)
  warning(n, " position(s) with polyploid genotype(s). Compared as sorted allele set.")
}

keep <- !b$missing & !cmp$missing
message("After removing missing/malformed: ", sum(keep))
if (opt$baseHetcheck) { keep <- keep & !is.na(b$homo) & b$homo; message("After base hom filter: ", sum(keep)) }
if (opt$copHetcheck)  { keep <- keep & !is.na(cmp$homo) & cmp$homo; message("After comp hom filter: ", sum(keep)) }

data <- data[keep, ]
is_variant <- b$norm[keep] != cmp$norm[keep]

message("Total: ", nrow(data), " | Variants: ", sum(is_variant), " | Non-variants: ", sum(!is_variant))

# Preview
cat("\n=== First", MAX_DISPLAY_ROWS, "variants (", base_name, "!=", comp_name, ") ===\n")
if (sum(is_variant) > 0) {
  vi <- which(is_variant)[1:min(MAX_DISPLAY_ROWS, sum(is_variant))]
  print(data.frame(CHROM = data$CHROM[vi], POS = data$POS[vi],
                   Baseline_GT = data[[base_col]][vi], Comparison_GT = data[[comp_col]][vi]),
        row.names = FALSE)
  if (sum(is_variant) > MAX_DISPLAY_ROWS) cat("... (", sum(is_variant), "total )\n")
} else cat("No variants found.\n")
cat("========================================================\n\n")

# --- Chromosome lengths ------------------------------------------------------
chromosomes <- unique(data$CHROM)

if (!is.null(opt$chrlength)) {
  if (!file.exists(opt$chrlength)) stop("Length file not found: ", opt$chrlength, call. = FALSE)
  message("Reading length file: ", opt$chrlength)
  # select=1:2 handles .fai (5 cols) and any multi-column format
  chr_lengths <- tryCatch(
    as.data.frame(fread(opt$chrlength, header = FALSE, select = 1:2,
                        col.names = c("CHROM", "LENGTH"), colClasses = list(character = 1))),
    error = function(e) stop("Error reading length file: ", e$message, call. = FALSE))

  chr_lengths$LENGTH <- suppressWarnings(as.numeric(chr_lengths$LENGTH))
  chr_lengths <- chr_lengths[!is.na(chr_lengths$LENGTH) & chr_lengths$LENGTH > 0, ]
  dups <- duplicated(chr_lengths$CHROM)
  if (any(dups)) { warning("Removing ", sum(dups), " duplicate CHROM from length file."); chr_lengths <- chr_lengths[!dups, ] }
  if (nrow(chr_lengths) == 0) stop("No valid entries in length file.", call. = FALSE)

  # Unified: length-file chromosomes UNION data chromosomes
  all_chrs <- unique(c(chr_lengths$CHROM, chromosomes))
  miss <- setdiff(chromosomes, chr_lengths$CHROM)
  if (length(miss) > 0) warning("Not in length file (using max POS): ", paste(miss, collapse = ", "))
  chr_info <- merge(data.frame(CHROM = all_chrs, stringsAsFactors = FALSE), chr_lengths, by = "CHROM", all.x = TRUE)
} else {
  if (nrow(data) == 0) stop("All positions filtered out; provide -l to plot empty frame.", call. = FALSE)
  warning("No length file. Using max observed position per chromosome.")
  chr_info <- data.frame(CHROM = chromosomes, LENGTH = NA_real_, stringsAsFactors = FALSE)
}

# Fill NA lengths from data (vectorized)
if (nrow(data) > 0 && any(is.na(chr_info$LENGTH))) {
  obs_max <- tapply(data$POS, data$CHROM, max)
  na_i <- is.na(chr_info$LENGTH)
  chr_info$LENGTH[na_i] <- as.numeric(obs_max[chr_info$CHROM[na_i]])
}
chr_info$LENGTH[is.na(chr_info$LENGTH)] <- 0

# Sort genomically
chr_info <- chr_info[order(chr_sort_key(chr_info$CHROM), chr_info$CHROM), ]
chr_info$chr_order <- seq_len(nrow(chr_info))
chr_info$LENGTH_scaled <- chr_info$LENGTH / opt$unit

cat("\nChromosomes (", nrow(chr_info), "):\n")
print(chr_info[, c("CHROM", "LENGTH_scaled")], row.names = FALSE)

# --- Plot data ---------------------------------------------------------------
if (nrow(data) > 0) {
  data$is_variant <- is_variant
  plot_data <- merge(data, chr_info[, c("CHROM", "chr_order", "LENGTH_scaled")], by = "CHROM")
  plot_data$POS_scaled <- plot_data$POS / opt$unit
  plot_data <- plot_data[order(plot_data$chr_order, plot_data$POS), ]
  # POS overflow check
  over <- plot_data$POS > plot_data$LENGTH_scaled * opt$unit
  if (any(over)) warning(sum(over), " variant(s) beyond chromosome length. Check assembly version.")
  plot_variants <- plot_data[plot_data$is_variant, ]
} else {
  plot_variants <- data.frame(CHROM = character(0), POS = numeric(0), POS_scaled = numeric(0),
                              chr_order = numeric(0), is_variant = logical(0), stringsAsFactors = FALSE)
}

# Output table (genomic order)
if (!is.null(opt$output_table)) {
  message("Writing table: ", opt$output_table)
  fwrite(plot_variants[, c("CHROM", "POS")], opt$output_table, sep = "\t")
  message("  ", nrow(plot_variants), " positions written")
}

# --- Plot --------------------------------------------------------------------
message("\nGenerating plot...")
if (nrow(plot_variants) == 0) message("  (empty frame - no variants)")

unit_label <- if (opt$unit == 1e6) "Mb" else if (opt$unit == 1e3) "kb" else if (opt$unit == 1) "bp" else
  paste0("\u00d7", format(opt$unit, scientific = FALSE), " bp")

p <- ggplot() +
  geom_rect(data = chr_info,
            aes(xmin = 0, xmax = LENGTH_scaled, ymin = chr_order - 0.4, ymax = chr_order + 0.4),
            fill = "lightgray", color = opt$chrBorderColor, linewidth = opt$chrBorderSize) +
  geom_segment(data = plot_variants,
               aes(x = POS_scaled, xend = POS_scaled, y = chr_order - 0.4, yend = chr_order + 0.4),
               color = opt$segmentColor, linewidth = opt$segmentSize, alpha = opt$segmentAlpha) +
  scale_y_reverse(breaks = chr_info$chr_order, labels = chr_info$CHROM) +
  labs(x = paste0("Position (", unit_label, ")"), y = "Chromosome",
       title = "Variant Position Plot",
       subtitle = paste("Comparing", base_name, "vs", comp_name)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 10))

# ggsave: auto-detects device from extension, handles open/close
plot_h <- min(50, max(4, nrow(chr_info) * 0.5))
message("Saving: ", opt$output)
ggsave(opt$output, p, width = 12, height = plot_h, dpi = 100, limitsize = FALSE,
       device = if (output_ext == "svg") grDevices::svg else NULL)

message("\nDone! ", nrow(plot_variants), " variants plotted -> ", opt$output)

# --- CMplot SNP-density plot (optional) --------------------------------------
if (isTRUE(opt$cmplot)) {
  if (nrow(plot_variants) == 0) {
    warning("--cmplot skipped: no variant positions to plot density for.")
  } else {
    message("\nGenerating CMplot SNP-density plot...")

    # Construct CMplot input: 3 columns (SNP name, Chromosome, Position)
    cmplot_input <- data.frame(
      SNP = paste0(plot_variants$CHROM, "_", plot_variants$POS),
      Chr = plot_variants$CHROM,
      Pos = plot_variants$POS,
      stringsAsFactors = FALSE
    )

    # Derive output file name: same directory, add _density suffix
    out_dir <- normalizePath(dirname(opt$output), mustWork = TRUE)
    out_base <- sub("\\.[^.]*$", "", basename(opt$output))
    # CMplot constructs filename as: Marker_Density.<file.name>.<ext>
    # So we pass only the name part and temporarily setwd to the output directory.
    cmplot_file_ext <- if (output_ext %in% c("jpg", "jpeg")) "jpg" else output_ext
    cmplot_name <- paste0(out_base, "_density")

    # Title
    cmplot_main <- if (!is.null(opt$cmplot_main)) opt$cmplot_main else
      paste("SNP Density:", base_name, "vs", comp_name)

    # Density colors
    cmplot_den_cols <- trimws(strsplit(opt$cmplot_den_col, ",")[[1]])

    # Call CMplot (temporarily switch to output directory)
    old_wd <- setwd(out_dir)
    on.exit(setwd(old_wd), add = TRUE)

    CMplot(cmplot_input,
      plot.type = "d",
      bin.size = opt$cmplot_bin_size,
      chr.den.col = cmplot_den_cols,
      chr.labels = chr_info$CHROM,
      chr.pos.max = TRUE,
      main = cmplot_main,
      main.cex = 1.2,
      file = cmplot_file_ext,
      file.name = cmplot_name,
      file.output = TRUE,
      dpi = opt$cmplot_dpi,
      width = opt$cmplot_width,
      height = opt$cmplot_height,
      verbose = FALSE
    )

    setwd(old_wd)
    # CMplot outputs: Marker_Density.<file.name>.<ext>
    cmplot_out <- file.path(out_dir, paste0("Marker_Density.", cmplot_name, ".", cmplot_file_ext))
    message("CMplot density plot saved: ", cmplot_out)
  }
}
