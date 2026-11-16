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

# Built-in colour presets. Each entry is a list with:
#   segment  - variant segment colour
#   border   - chromosome rectangle border colour
#   density  - 2-3 colour gradient for the density plot (low -> high)
#   palettes - discrete colours used for the chromosome bands / density legend
# All values are plain R colour names or hex codes, so no extra package is
# needed. Names are chosen to match common journal figure styles.
COLOR_PRESETS <- list(
  # Nature/Cell style: muted blue-gray body, vivid accent for hits
  "nature" = list(
    label   = "Nature (muted blue-gray, orange accents)",
    segment = "#D1495B",
    border  = "#2B2D42",
    density = c("#F1F3F5", "#8D99AE", "#2B2D42"),
    palettes = c("#8D99AE", "#2B2D42")
  ),
  # Science/AAAS style: teal-to-crimson, high contrast
  "science" = list(
    label   = "Science (teal to crimson)",
    segment = "#C8102E",
    border  = "#1D3557",
    density = c("#F7F7F7", "#2A9D8F", "#C8102E"),
    palettes = c("#2A9D8F", "#C8102E")
  ),
  # Lancet / NEJM style: classic navy + brick red
  "lancet" = list(
    label   = "Lancet (navy and brick red)",
    segment = "#B22222",
    border  = "#003049",
    density = c("#FFF8F0", "#4A7C9B", "#003049"),
    palettes = c("#4A7C9B", "#003049")
  ),
  # Colorblind-safe (Okabe-Ito derived), good for print and screen
  "okabe_ito" = list(
    label   = "Okabe-Ito (colour-blind safe)",
    segment = "#D55E00",
    border  = "#0072B2",
    density = c("#F0F0F0", "#56B4E9", "#0072B2"),
    palettes = c("#56B4E9", "#0072B2")
  ),
  # Viridis-like sequential, perceptually uniform
  "viridis" = list(
    label   = "Viridis (perceptually uniform)",
    segment = "#FDE725",
    border  = "#440154",
    density = c("#440154", "#21918C", "#FDE725"),
    palettes = c("#440154", "#21918C")
  ),
  # Grayscale, for journals that forbid colour
  "grayscale" = list(
    label   = "Grayscale (print-safe)",
    segment = "#B2182B",
    border  = "#333333",
    density = c("#F7F7F7", "#999999", "#252525"),
    palettes = c("#999999", "#252525")
  ),
  # Legacy: the original red-on-lightgray look
  "classic" = list(
    label   = "Classic (red on light gray)",
    segment = "red",
    border  = "black",
    density = c("#F5F5F5", "#FDBB84", "#B30000"),
    palettes = c("#FDBB84", "#B30000")
  )
)

DEFAULT_PRESET <- "nature"

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
  make_option(c("--preset"), type="character", default=DEFAULT_PRESET,
              help=paste0("Colour preset: ", paste(names(COLOR_PRESETS), collapse=", "),
                          " [default=%default]"), metavar="NAME"),
  make_option(c("--segmentColor"), type="character", default=NULL,
              help="Color for variant segments (overrides --preset)", metavar="COLOR"),
  make_option(c("--segmentSize"), type="numeric", default=0.05,
              help="Thickness of variant segments [default=%default]", metavar="NUM"),
  make_option(c("--chrBorderColor"), type="character", default=NULL,
              help="Color for chromosome borders (overrides --preset)", metavar="COLOR"),
  make_option(c("--chrBorderSize"), type="numeric", default=0.3,
              help="Thickness of chromosome borders [default=%default]", metavar="NUM"),
  make_option(c("--segmentAlpha"), type="numeric", default=0.6,
              help="Alpha transparency for variant segments (0-1) [default=%default]", metavar="NUM"),
  make_option(c("--output_table"), type="character", default=NULL,
              help="Write variant positions + genotypes to this file", metavar="FILE"),
  make_option(c("--table_min_len_diff"), type="integer", default=NULL,
              help=paste("Only write --output_table rows whose LEN_DIFF is at least",
                         "this value. 1 drops length-neutral substitutions (SNPs and",
                         "equal-length indels). Default: keep every row"), metavar="INT"),
  # SNP-density plot options (drawn with ggplot2, no CMplot dependency)
  make_option(c("--density"), action="store_true", default=FALSE,
              help="Additionally generate a SNP-density plot"),
  make_option(c("--density_bin_size"), type="numeric", default=1e6,
              help="Density window size in bp [default=%default]", metavar="NUM"),
  make_option(c("--density_den_col"), type="character", default=NULL,
              help="Density gradient colors, comma-separated (overrides --preset)", metavar="COLORS"),
  make_option(c("--density_dpi"), type="integer", default=300,
              help="Density plot output DPI [default=%default]", metavar="INT"),
  make_option(c("--density_width"), type="numeric", default=NULL,
              help="Density plot width in inches (default: auto)", metavar="NUM"),
  make_option(c("--density_height"), type="numeric", default=NULL,
              help="Density plot height in inches (default: auto)", metavar="NUM"),
  make_option(c("--density_main"), type="character", default=NULL,
              help="Density plot title (default: auto-generated)", metavar="TEXT"),
  make_option(c("--list_presets"), action="store_true", default=FALSE,
              help="List available colour presets and exit"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE,
              help="Print version and exit")
)

SCRIPT_VERSION <- "2.1.0"

opt_parser <- OptionParser(option_list=option_list,
  description="\nPlot variant positions from GATK VariantsToTable output\n\nExample:\n  Rscript vcf_difplot.R -i input.table -b sample1 -c sample2 -o output.pdf")
opt <- parse_args(opt_parser)

if (isTRUE(opt$version)) {
  cat("vcf_difplot.R", SCRIPT_VERSION, "\n")
  quit(save = "no", status = 0L)
}

if (isTRUE(opt$list_presets)) {
  cat("Available colour presets:\n\n")
  for (nm in names(COLOR_PRESETS)) {
    p <- COLOR_PRESETS[[nm]]
    cat(sprintf("  %-12s %s\n", nm, p$label))
    cat(sprintf("  %-12s segment=%s  border=%s  density=%s\n\n",
                "", p$segment, p$border, paste(p$density, collapse=" -> ")))
  }
  cat("Use --preset NAME to select one, or override individual colours with\n")
  cat("--segmentColor / --chrBorderColor / --density_den_col.\n")
  quit(save = "no", status = 0L)
}

# Resolve the colour preset. Individual colour flags override preset values.
if (is.null(COLOR_PRESETS[[opt$preset]])) {
  stop("Unknown --preset '", opt$preset, "'. Available: ",
       paste(names(COLOR_PRESETS), collapse = ", "), call. = FALSE)
}
preset <- COLOR_PRESETS[[opt$preset]]

if (is.null(opt$segmentColor))    opt$segmentColor    <- preset$segment
if (is.null(opt$chrBorderColor))  opt$chrBorderColor  <- preset$border
if (is.null(opt$density_den_col)) opt$density_den_col <- paste(preset$density, collapse = ",")

# Parse the density colors once; reused by validation and the density call.
density_den_cols <- trimws(strsplit(opt$density_den_col, ",")[[1]])

# --- Helper functions --------------------------------------------------------

# Read column names from a tab-delimited file.
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
#
# Malformed genotypes (empty allele, '.', '*', consecutive separators) are
# marked missing. A trailing/leading separator such as "A/" or "/A" is
# malformed too: strsplit() would silently turn "A/" into a single allele and
# the value would be mistaken for a haploid genotype.
parse_gt <- function(gt) {
  gt <- as.character(gt)
  u <- unique(gt)
  s <- gsub("|", "/", u, fixed = TRUE)

  # Detect a leading/trailing separator before splitting, so "A/" is not
  # reinterpreted as haploid "A".
  edge_sep <- grepl("^/|/$", s)
  al <- strsplit(s, "/", fixed = TRUE)

  bad <- vapply(al, function(a) {
    length(a) == 0L || any(a == "") || any(a == ".") || any(a == "*")
  }, logical(1L), USE.NAMES = FALSE)
  bad <- bad | edge_sep | is.na(u) | u == ""

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

# Collapse a homozygous genotype to a single allele: "A/A" -> "A", "G|G" -> "G".
# Heterozygous / polyploid genotypes are returned unchanged so no information is
# lost when the homozygosity filters are off.
collapse_homo_gt <- function(gt) {
  gt <- as.character(gt)
  u <- unique(gt)
  out <- vapply(u, function(g) {
    if (is.na(g) || !nzchar(g)) return(NA_character_)
    al <- strsplit(gsub("|", "/", g, fixed = TRUE), "/", fixed = TRUE)[[1L]]
    if (length(al) >= 1L && all(al == al[1L])) al[1L] else g
  }, character(1L), USE.NAMES = FALSE)
  out[match(gt, u)]
}

# Length of the allele sequence carried by a genotype, used for the
# "how many bp does this differ by" column. Two identical alleles count once.
gt_allele_len <- function(gt) {
  gt <- as.character(gt)
  u <- unique(gt)
  out <- vapply(u, function(g) {
    if (is.na(g) || !nzchar(g)) return(NA_integer_)
    al <- strsplit(gsub("|", "/", g, fixed = TRUE), "/", fixed = TRUE)[[1L]]
    al <- unique(al)
    if (length(al) == 0L) return(NA_integer_)
    as.integer(nchar(al[1L]))
  }, integer(1L), USE.NAMES = FALSE)
  out[match(gt, u)]
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
    # file_path_sans_ext() only strips the basename's extension, so a dotted
    # directory component (e.g. "my.dir/plot.xyz") is preserved.
    return(list(path = paste0(tools::file_path_sans_ext(output_path), ".pdf"), ext = "pdf"))
  }
  list(path = output_path, ext = ext)
}

# Create the parent directory of a path if it does not exist yet.
ensure_dir <- function(path, what = "output") {
  d <- dirname(path)
  if (!nzchar(d) || dir.exists(d)) return(invisible(TRUE))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d))
    stop("Cannot create ", what, " directory: ", d, call. = FALSE)
  message("Created ", what, " directory: ", d)
  invisible(TRUE)
}

# --- SNP density plot (self-contained, no CMplot) ----------------------------

# Count variants per fixed-width window along each chromosome.
# Returns a data.frame with CHROM, start, end, mid, count, chr_order.
compute_density <- function(variants, chr_info, bin_size) {
  if (nrow(variants) == 0L) {
    return(data.frame(CHROM = character(0), start = numeric(0), end = numeric(0),
                      mid = numeric(0), count = numeric(0), chr_order = numeric(0),
                      stringsAsFactors = FALSE))
  }
  # Window index within each chromosome; the last (partial) window is kept.
  win <- floor(variants$POS / bin_size)
  key <- paste(variants$CHROM, win, sep = "\r")
  agg <- data.frame(
    CHROM = tapply(variants$CHROM, key, `[`, 1L),
    win   = tapply(win, key, `[`, 1L),
    count = as.integer(table(key)),
    stringsAsFactors = FALSE
  )
  agg$start <- agg$win * bin_size
  # Clamp the window end to the chromosome length so the last window does not
  # spill past the chromosome body.
  chr_len <- chr_info$LENGTH[match(agg$CHROM, chr_info$CHROM)]
  agg$end <- pmin(agg$start + bin_size, chr_len)
  agg$mid <- (agg$start + agg$end) / 2
  agg$chr_order <- chr_info$chr_order[match(agg$CHROM, chr_info$CHROM)]
  agg <- agg[order(agg$chr_order, agg$start), ]
  rownames(agg) <- NULL
  agg
}

# Draw the density plot: one horizontal band per chromosome, filled with a
# colour gradient proportional to the local variant count.
plot_density <- function(dens, chr_info, pos_unit, unit_label, den_cols,
                         title, subtitle, border_color) {
  den_cols <- if (length(den_cols) == 1L) c(den_cols, den_cols) else den_cols
  max_count <- max(dens$count, na.rm = TRUE)
  if (!is.finite(max_count) || max_count <= 0) max_count <- 1

  # Chromosome body (light track) so empty regions are still visible.
  chr_info$LENGTH_scaled <- chr_info$LENGTH / pos_unit

  ggplot() +
    # 1) light track underneath, so regions with no variants still show
    geom_rect(data = chr_info,
              aes(xmin = 0, xmax = LENGTH_scaled,
                  ymin = chr_order - 0.4, ymax = chr_order + 0.4),
              fill = "grey96", colour = NA) +
    # 2) density windows on top of the track
    geom_rect(data = dens,
              aes(xmin = start / pos_unit, xmax = end / pos_unit,
                  ymin = chr_order - 0.4, ymax = chr_order + 0.4,
                  fill = count),
              colour = NA) +
    # 3) outline drawn LAST so the density fill cannot cover it. Drawing the
    #    border on the underlying track would leave gaps wherever a window
    #    covers it.
    geom_rect(data = chr_info,
              aes(xmin = 0, xmax = LENGTH_scaled,
                  ymin = chr_order - 0.4, ymax = chr_order + 0.4),
              fill = NA, colour = border_color, linewidth = 0.3) +
    scale_fill_gradientn(colours = den_cols, limits = c(0, max_count),
                         name = "Variants\nper window") +
    scale_y_reverse(breaks = chr_info$chr_order, labels = chr_info$CHROM) +
    labs(x = paste0("Position (", unit_label, ")"), y = NULL,
         title = title, subtitle = subtitle) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 9),
          legend.key.height = grid::unit(1.1, "lines"))
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
    if (length(val) == 0L) {
      cat("\n[!] No more input available on stdin. Exiting.\n")
      cat("    Interactive mode needs a terminal (TTY). When piping answers in,\n")
      cat("    non-interactive use is simpler:\n")
      cat("      Rscript vcf_difplot.R -i FILE -b BASE -c COMP -o OUT\n")
      quit(save = "no", status = 1L)
    }
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
  prompt_unit_interval <- function(label, desc, def) {
    repeat {
      v <- suppressWarnings(as.numeric(prompt_optional(label, desc, as.character(def))))
      if (!is.na(v) && v >= 0 && v <= 1) return(v)
      cat("    [!] Enter a number between 0 and 1.\n")
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
  preset_names <- paste(names(COLOR_PRESETS), collapse = ", ")
  result$preset <- prompt_optional("Colour preset (--preset)",
    paste0("One of: ", preset_names, ". Leave empty for default."), DEFAULT_PRESET)
  if (!result$preset %in% names(COLOR_PRESETS)) {
    cat("    [!] Unknown preset, using", DEFAULT_PRESET, "\n")
    result$preset <- DEFAULT_PRESET
  }
  preset_sel <- COLOR_PRESETS[[result$preset]]
  result$segmentColor <- prompt_optional("Segment color", "R name or hex. Empty = preset.",
                                         preset_sel$segment)
  result$segmentSize <- prompt_positive("Segment size", "Line width.", 0.05)
  result$segmentAlpha <- prompt_unit_interval("Segment alpha", "0-1.", 0.6)
  result$chrBorderColor <- prompt_optional("Border color", "R name or hex. Empty = preset.",
                                           preset_sel$border)
  result$chrBorderSize <- prompt_positive("Border size", "Line width.", 0.3)
  cat("\n")
  result$output_table <- prompt_optional("Output table (--output_table)", "CHROM/POS/GT file path.", NULL)
  if (!is.null(result$output_table)) {
    result$table_min_len_diff <- prompt_optional(
      "  Min LEN_DIFF (--table_min_len_diff)",
      "Drop rows below this length difference. 1 drops SNPs/equal-length indels. Empty = keep all.",
      NULL)
    if (!is.null(result$table_min_len_diff)) {
      v <- suppressWarnings(as.integer(result$table_min_len_diff))
      result$table_min_len_diff <- if (is.na(v) || v < 0) NULL else v
    }
  }

  # SNP-density plot
  cat("\n")
  result$density <- prompt_yesno("SNP-density plot (--density)",
    "Additionally generate a variant-density plot.", FALSE)
  if (isTRUE(result$density)) {
    result$density_bin_size <- prompt_positive("  Density bin size (--density_bin_size)", "Window in bp.", 1e6)
    result$density_den_col <- prompt_optional("  Density colors (--density_den_col)",
      "Comma-separated. Empty = preset.", paste(preset_sel$density, collapse = ","))
    result$density_dpi <- prompt_positive("  DPI (--density_dpi)", "Output resolution.", 300)
    result$density_main <- prompt_optional("  Title (--density_main)", "Plot title.", NULL)
    result$density_width <- NULL
    result$density_height <- NULL
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
  cmd <- c(cmd, "--preset", shQuote(result$preset),
           "--segmentColor", shQuote(result$segmentColor), "--segmentSize", result$segmentSize,
           "--segmentAlpha", result$segmentAlpha, "--chrBorderColor", shQuote(result$chrBorderColor),
           "--chrBorderSize", result$chrBorderSize)
  if (!is.null(result$output_table)) cmd <- c(cmd, "--output_table", shQuote(result$output_table))
  if (!is.null(result$table_min_len_diff))
    cmd <- c(cmd, "--table_min_len_diff", result$table_min_len_diff)
  if (isTRUE(result$density)) {
    cmd <- c(cmd, "--density", "--density_bin_size", result$density_bin_size,
             "--density_den_col", shQuote(result$density_den_col), "--density_dpi", result$density_dpi)
    if (!is.null(result$density_main)) cmd <- c(cmd, "--density_main", shQuote(result$density_main))
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

if (!is.null(opt$table_min_len_diff)) {
  if (!is.numeric(opt$table_min_len_diff) || opt$table_min_len_diff < 0)
    stop("--table_min_len_diff must be a non-negative integer", call. = FALSE)
  if (is.null(opt$output_table))
    warning("--table_min_len_diff has no effect without --output_table.")
}

# Validate density-plot parameters if enabled
if (isTRUE(opt$density)) {
  if (!is.numeric(opt$density_bin_size) || opt$density_bin_size <= 0)
    stop("--density_bin_size must be a positive number", call. = FALSE)
  if (!is.null(opt$density_dpi) && (!is.numeric(opt$density_dpi) || opt$density_dpi < 1))
    stop("--density_dpi must be a positive integer", call. = FALSE)
  for (cc in density_den_cols) validate_color(cc, "--density_den_col")
}

# Validate output early (before expensive processing)
out_info <- validate_output(opt$output)
opt$output <- out_info$path
output_ext <- out_info$ext
ensure_dir(opt$output, "output")
if (!is.null(opt$output_table)) ensure_dir(opt$output_table, "output table")

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

# POS must be a non-negative integer. Validate before any arithmetic so a bad
# value fails here with a clear message instead of deep inside plotting.
data$POS <- suppressWarnings(as.numeric(data$POS))
if (anyNA(data$POS)) {
  warning(sum(is.na(data$POS)), " row(s) with non-numeric POS removed.")
  data <- data[!is.na(data$POS), , drop = FALSE]
}
if (any(data$POS < 0)) {
  warning(sum(data$POS < 0), " row(s) with negative POS removed.")
  data <- data[data$POS >= 0, , drop = FALSE]
}
frac <- data$POS != floor(data$POS)
if (any(frac)) {
  warning(sum(frac), " row(s) with non-integer POS rounded down.")
  data$POS[frac] <- floor(data$POS[frac])
}
message("Read ", nrow(data), " positions")

# Remember the chromosomes present in the raw input so an all-filtered run can
# still draw an empty frame instead of aborting.
raw_chromosomes <- unique(data$CHROM)

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

# Collapse duplicate CHROM/POS rows. The same position can appear on several
# rows (e.g. multi-allelic records split by VariantsToTable); without this the
# plot draws identical segments on top of each other and --output_table gets
# duplicate rows.
n_before <- nrow(data)
if (n_before > 0L) {
  dup <- duplicated(paste(data$CHROM, data$POS, sep = "\r"))
  if (any(dup)) {
    message("Collapsing ", sum(dup), " duplicate CHROM/POS row(s) (",
            n_before, " -> ", sum(!dup), ").")
    data <- data[!dup, , drop = FALSE]
    is_variant <- is_variant[!dup]
  }
}

n_total <- nrow(data)
n_var <- sum(is_variant)
message("Total: ", n_total, " | Variants: ", n_var, " | Non-variants: ", n_total - n_var)

# Preview
cat("\n=== First", MAX_DISPLAY_ROWS, "variants (", base_name, "!=", comp_name, ") ===\n")
if (n_var > 0) {
  vi <- which(is_variant)[seq_len(min(MAX_DISPLAY_ROWS, n_var))]
  print(data.frame(CHROM = data$CHROM[vi], POS = data$POS[vi],
                   Baseline_GT = data[[base_col]][vi], Comparison_GT = data[[comp_col]][vi]),
        row.names = FALSE)
  if (n_var > MAX_DISPLAY_ROWS) cat("... (", n_var, "total )\n")
} else cat("No variants found.\n")
cat("========================================================\n\n")

# --- Chromosome lengths ------------------------------------------------------
chromosomes <- unique(data$CHROM)

if (!is.null(opt$chrlength)) {  if (!file.exists(opt$chrlength)) stop("Length file not found: ", opt$chrlength, call. = FALSE)
  message("Reading length file: ", opt$chrlength)
  # select=1:2 handles .fai (5 cols) and any multi-column format
  chr_lengths <- tryCatch(
    as.data.frame(fread(opt$chrlength, header = FALSE, select = 1:2,
                        col.names = c("CHROM", "LENGTH"), colClasses = list(character = 1))),
    error = function(e) stop("Error reading length file: ", e$message, call. = FALSE))

  chr_lengths$LENGTH <- suppressWarnings(as.numeric(chr_lengths$LENGTH))
  bad_len <- is.na(chr_lengths$LENGTH) | chr_lengths$LENGTH <= 0
  if (any(bad_len)) {
    warning("Ignoring ", sum(bad_len), " non-positive/unparseable LENGTH entr",
            if (sum(bad_len) == 1L) "y" else "ies", " in length file: ",
            paste(chr_lengths$CHROM[bad_len], collapse = ", "))
    chr_lengths <- chr_lengths[!bad_len, ]
  }
  dups <- duplicated(chr_lengths$CHROM)
  if (any(dups)) { warning("Removing ", sum(dups), " duplicate CHROM from length file."); chr_lengths <- chr_lengths[!dups, ] }
  if (nrow(chr_lengths) == 0) stop("No valid entries in length file.", call. = FALSE)

  # Length-file chromosomes are the authoritative set; data chromosomes missing
  # from it fall back to their max observed position.
  all_chrs <- unique(c(chr_lengths$CHROM, chromosomes))
  miss <- setdiff(chromosomes, chr_lengths$CHROM)
  if (length(miss) > 0)
    warning("Chromosome(s) absent from length file, using max observed POS: ",
            paste(miss, collapse = ", "))
  chr_info <- merge(data.frame(CHROM = all_chrs, stringsAsFactors = FALSE), chr_lengths, by = "CHROM", all.x = TRUE)
} else {
  if (nrow(data) == 0) {
    # Nothing survived filtering. Draw an empty frame using the chromosomes
    # seen in the raw input (or the length file) instead of aborting, so the
    # user still gets a plot showing the genomic axes.
    warning("No positions survived filtering; drawing an empty frame. ",
            "Supply -l to control the chromosome set.")
    chromosomes <- if (length(raw_chromosomes) > 0) raw_chromosomes else character(0)
    if (length(chromosomes) == 0)
      stop("No usable positions and no chromosome information available; ",
           "supply -l to draw an empty frame.", call. = FALSE)
  } else {
    warning("No length file. Using max observed position per chromosome.")
  }
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
  # match() instead of merge(): merge() re-sorts rows and is ~17x slower on
  # multi-million-row inputs (measured 17.6s vs 1.0s at 8M rows).
  idx <- match(data$CHROM, chr_info$CHROM)
  plot_data <- data
  plot_data$is_variant <- is_variant
  plot_data$chr_order <- chr_info$chr_order[idx]
  plot_data$LENGTH <- chr_info$LENGTH[idx]
  plot_data$POS_scaled <- plot_data$POS / opt$unit
  plot_data <- plot_data[order(plot_data$chr_order, plot_data$POS), ]

  # Overflow check against the raw LENGTH, avoiding the divide-then-multiply
  # round trip through LENGTH_scaled.
  over <- plot_data$POS > plot_data$LENGTH
  if (any(over)) warning(sum(over), " variant(s) beyond chromosome length. Check assembly version.")

  # Keep the genotype columns alongside the plot coordinates so --output_table
  # can report what each sample actually carried at every variant position.
  plot_variants <- plot_data[plot_data$is_variant,
                             c("CHROM", "POS", "POS_scaled", "chr_order",
                               base_col, comp_col)]
} else {
  plot_variants <- data.frame(CHROM = character(0), POS = numeric(0), POS_scaled = numeric(0),
                              chr_order = numeric(0), stringsAsFactors = FALSE)
  plot_variants[[base_col]] <- character(0)
  plot_variants[[comp_col]] <- character(0)
}
n_pv <- nrow(plot_variants)

# --- Merge overlapping segments ----------------------------------------------
# Every variant is drawn as a vertical segment of width `segmentSize` at
# (chr_order, POS_scaled). Each segment costs a full path in the PDF content
# stream, so a plot with millions of variants produces a huge file that takes
# ages to open. Segments that overlap horizontally, however, cover the same
# pixels as one wider rectangle: three 0.05-wide segments spanning
# (1.00,1.05) (1.02,1.07) (1.06,1.11) look the same as a single (1.00,1.11)
# rectangle. Overlapping neighbours are therefore coalesced into one rectangle
# per contiguous run, which cuts the path count dramatically while keeping the
# painted area identical.
#
# The segments are drawn with alpha, so where they overlap the colour stacks up
# and gets darker. That is accounted for by compositing the n stacked copies
# over the chromosome body and painting the merged rectangle with the resulting
# opaque colour:
#   C_out = C_seg * (1 - (1-a)^n) + C_bg * (1-a)^n
# A run of length 1 reproduces the original single segment exactly.
#
# plot_variants itself is left untouched, so the --output_table export still
# lists every variant position.
PLOT_WIDTH_IN   <- 12
PLOT_DPI        <- 100
# Fraction of the figure width occupied by the panel: the rest is the left axis
# labels, tick marks and margins.
PANEL_FRACTION  <- 0.80

# Colour of the chromosome body, i.e. what the segments are painted onto.
CHR_FILL <- "lightgray"

plot_segments <- plot_variants
if (n_pv > 1L) {
  # Segment half-width in x-axis units. `segmentSize` is a ggplot linewidth in
  # millimetres; the panel spans PANEL_FRACTION of PLOT_WIDTH_IN inches and
  # covers max_len x-units, so the conversion is a plain ratio.
  max_len <- max(chr_info$LENGTH_scaled, na.rm = TRUE)
  if (!is.finite(max_len) || max_len <= 0) max_len <- 1
  panel_in <- PLOT_WIDTH_IN * PANEL_FRACTION
  half_w <- (opt$segmentSize / 25.4) / 2 / panel_in * max_len

  # Walk the segments in genomic order and start a new run whenever the next
  # segment does not touch the current run (gap larger than the line width).
  # The data is already sorted by (chr_order, POS_scaled).
  same_chr <- c(FALSE, diff(plot_variants$chr_order) == 0L)
  prev_pos <- c(NA_real_, plot_variants$POS_scaled[-n_pv])
  touches  <- same_chr & ((plot_variants$POS_scaled - half_w) <= (prev_pos + half_w))
  touches[is.na(touches)] <- FALSE

  run_id <- cumsum(!touches)
  n_runs <- run_id[n_pv]

  if (n_runs < n_pv) {
    # A run is one contiguous band of touching segments. Painting a long run as
    # a single flat rectangle would erase the density information: on a
    # chromosome densely packed with variants the whole strip would become one
    # uniform colour, whereas the original shows darker and lighter stretches.
    # Each run is therefore cut into slices about one pixel wide, and every
    # slice gets its own colour, so the shading survives while the path count
    # drops from one-per-variant to a few hundred per chromosome.
    max_len_slice <- max_len
    x_per_px <- max_len_slice / (PLOT_WIDTH_IN * PLOT_DPI * PANEL_FRACTION)

    first <- !duplicated(run_id)
    last  <- !duplicated(run_id, fromLast = TRUE)
    run_xmin <- plot_variants$POS_scaled[first] - half_w
    run_xmax <- plot_variants$POS_scaled[last] + half_w

    # Slice index inside each run: 0, 1, 2, ... for successive pixel-wide bins.
    run_start <- run_xmin[run_id]
    slice_of <- floor((plot_variants$POS_scaled - run_start) / x_per_px)
    slice_key <- paste(run_id, slice_of, sep = "\r")

    s_first <- !duplicated(slice_key)
    s_last  <- !duplicated(slice_key, fromLast = TRUE)
    plot_segments <- data.frame(
      chr_order = plot_variants$chr_order[s_first],
      # Slice edges are the bin boundaries, clamped to the run's real extent so
      # the band never grows beyond the segments it replaces. A slice never
      # spans two runs, because run_id is part of slice_key.
      xmin = run_start[s_first] + slice_of[s_first] * x_per_px,
      xmax = pmin(run_xmax[run_id][s_last],
                  run_start[s_first] + (slice_of[s_last] + 1) * x_per_px),
      stringsAsFactors = FALSE
    )

    # Average number of segments covering a point of each slice. Group the
    # segments by slice once (split() is a single pass), then sweep each
    # slice's edge events, so the whole step stays linear in the input size.
    seg_x <- plot_variants$POS_scaled
    slice_f <- factor(slice_key, levels = unique(slice_key))
    groups <- split(seg_x, slice_f)
    plot_segments$n_eff <- vapply(groups, function(gx) {
      lo <- gx - half_w
      hi <- gx + half_w
      ev_x <- c(lo, hi)
      ev_d <- c(rep(1, length(lo)), rep(-1, length(hi)))
      ord <- order(ev_x)
      ev_x <- ev_x[ord]; ev_d <- ev_d[ord]
      depth <- cumsum(ev_d)
      widths <- diff(ev_x)
      depth <- depth[-length(depth)]
      keep <- widths > 0
      if (!any(keep)) return(1)
      sum(widths[keep] * depth[keep]) / sum(widths[keep])
    }, numeric(1L), USE.NAMES = FALSE)

    seg_rgb <- col2rgb(opt$segmentColor)[, 1]
    bg_rgb  <- col2rgb(CHR_FILL)[, 1]
    a <- opt$segmentAlpha
    # (1-a)^n in log space, so a large n cannot underflow prematurely.
    stay <- if (a >= 1) rep(0, nrow(plot_segments)) else
      exp(plot_segments$n_eff * log1p(-a))
    mixed <- outer(1 - stay, seg_rgb) + outer(stay, bg_rgb)
    plot_segments$fill <- rgb(mixed[, 1], mixed[, 2], mixed[, 3], maxColorValue = 255)

    message("  Merged overlapping segments (", n_pv, " -> ",
            nrow(plot_segments), " drawn).")
  } else {
    # Nothing overlaps: every segment stays its own rectangle, exactly as wide
    # as the original line. The drawing layer below always works on xmin/xmax,
    # so the plain segment coordinates are converted here.
    plot_segments <- data.frame(
      chr_order = plot_variants$chr_order,
      xmin = plot_variants$POS_scaled - half_w,
      xmax = plot_variants$POS_scaled + half_w,
      stringsAsFactors = FALSE
    )
  }
} else if (n_pv == 1L) {
  # A single variant: no merging to do, but the rectangle columns are still
  # needed by the drawing layer.
  max_len <- max(chr_info$LENGTH_scaled, na.rm = TRUE)
  if (!is.finite(max_len) || max_len <= 0) max_len <- 1
  half_w <- (opt$segmentSize / 25.4) / 2 / (PLOT_WIDTH_IN * PANEL_FRACTION) * max_len
  plot_segments <- data.frame(
    chr_order = plot_variants$chr_order,
    xmin = plot_variants$POS_scaled - half_w,
    xmax = plot_variants$POS_scaled + half_w,
    stringsAsFactors = FALSE
  )
} else {
  # No variants at all: keep an empty frame with the columns the layer expects.
  plot_segments <- data.frame(
    chr_order = numeric(0), xmin = numeric(0), xmax = numeric(0),
    stringsAsFactors = FALSE
  )
}

# --- Plot --------------------------------------------------------------------
message("\nGenerating plot...")
if (n_pv == 0) message("  (empty frame - no variants)")

unit_label <- if (opt$unit == 1e6) "Mb" else if (opt$unit == 1e3) "kb" else if (opt$unit == 1) "bp" else
  paste0("\u00d7", format(opt$unit, scientific = FALSE), " bp")

# Segment appearance: merged runs carry a pre-composited opaque colour in
# `fill`; otherwise the plain segment colour and alpha are used, which is the
# original single-segment appearance.
if (!is.null(plot_segments$fill)) {
  seg_colour <- plot_segments$fill
  seg_alpha  <- 1
} else {
  seg_colour <- opt$segmentColor
  seg_alpha  <- opt$segmentAlpha
}

p <- ggplot() +
  geom_rect(data = chr_info,
            aes(xmin = 0, xmax = LENGTH_scaled, ymin = chr_order - 0.4, ymax = chr_order + 0.4),
            fill = CHR_FILL, color = opt$chrBorderColor, linewidth = opt$chrBorderSize) +
  # One rectangle per run of touching segments. A run of length 1 is drawn as a
  # rectangle exactly as wide as the original segment, so the picture is
  # unchanged; longer runs replace a stack of overlapping segments.
  geom_rect(data = plot_segments,
            aes(xmin = xmin, xmax = xmax, ymin = chr_order - 0.4, ymax = chr_order + 0.4),
            colour = seg_colour, fill = seg_colour, linewidth = 0, alpha = seg_alpha) +
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
ggsave(opt$output, p, width = PLOT_WIDTH_IN, height = plot_h, dpi = PLOT_DPI, limitsize = FALSE,
       device = if (output_ext == "svg") grDevices::svg else NULL)

message("\nDone! ", n_pv, " variants plotted -> ", opt$output)

# Output table (genomic order). Written only after the plot succeeded so a
# failed run does not leave a half-finished deliverable behind.
if (!is.null(opt$output_table)) {
  message("Writing table: ", opt$output_table)

  base_gt <- plot_variants[[base_col]]
  comp_gt <- plot_variants[[comp_col]]

  # Homozygous genotypes are reported as a single allele ("A/A" -> "A"), which
  # is what the homozygosity-filtered workflow always produces.
  out_tbl <- data.frame(
    CHROM = plot_variants$CHROM,
    POS = format(plot_variants$POS, scientific = FALSE, trim = TRUE),
    base = collapse_homo_gt(base_gt),
    comp = collapse_homo_gt(comp_gt),
    stringsAsFactors = FALSE
  )
  names(out_tbl)[3:4] <- c(base_col, comp_col)

  # Absolute allele-length difference between the two samples.
  len_base <- gt_allele_len(base_gt)
  len_comp <- gt_allele_len(comp_gt)
  out_tbl$LEN_DIFF <- abs(len_base - len_comp)

  # Optional length-difference filter. LEN_DIFF == 0 covers SNPs and
  # equal-length substitutions, which carry no length change at all.
  n_before_filter <- nrow(out_tbl)
  if (!is.null(opt$table_min_len_diff)) {
    out_tbl <- out_tbl[!is.na(out_tbl$LEN_DIFF) &
                       out_tbl$LEN_DIFF >= opt$table_min_len_diff, , drop = FALSE]
    message("  LEN_DIFF >= ", opt$table_min_len_diff, ": ",
            nrow(out_tbl), " of ", n_before_filter, " positions kept (",
            n_before_filter - nrow(out_tbl), " removed)")
    if (n_before_filter > 0L && nrow(out_tbl) == 0L)
      warning("--table_min_len_diff ", opt$table_min_len_diff,
              " removed every position. SNPs and equal-length indels all have ",
              "LEN_DIFF = 0; use 0 to keep them.", call. = FALSE)
  }

  fwrite(out_tbl, opt$output_table, sep = "\t", quote = FALSE, na = "NA")
  message("  ", nrow(out_tbl), " positions written (columns: CHROM, POS, ",
          base_col, ", ", comp_col, ", LEN_DIFF)")
}

# --- SNP-density plot (optional, self-contained) -----------------------------
if (isTRUE(opt$density)) {
  if (n_pv == 0) {
    warning("--density skipped: no variant positions to plot density for.")
  } else {
    message("\nGenerating SNP-density plot...")

    dens <- compute_density(plot_variants, chr_info, opt$density_bin_size)
    dens <- dens[!is.na(dens$chr_order), ]

    # Same directory as the main plot, "_density" suffix, same format.
    out_dir  <- dirname(opt$output)
    out_base <- tools::file_path_sans_ext(basename(opt$output))
    dens_ext <- output_ext
    dens_out <- file.path(out_dir, paste0(out_base, "_density.", dens_ext))

    dens_title <- if (!is.null(opt$density_main)) opt$density_main else
      paste("SNP Density:", base_name, "vs", comp_name)
    dens_sub <- paste0("Window: ",
                       format(opt$density_bin_size, big.mark = ",", scientific = FALSE),
                       " bp  |  max ", max(dens$count), " variants/window")

    # Only label chromosomes that actually carry variants, otherwise the y axis
    # shows length-file chromosomes that have no data.
    chr_sub <- chr_info[chr_info$CHROM %in% unique(dens$CHROM), ]
    chr_sub$chr_order <- seq_len(nrow(chr_sub))
    dens$chr_order <- chr_sub$chr_order[match(dens$CHROM, chr_sub$CHROM)]

    pd <- plot_density(dens, chr_sub, opt$unit, unit_label, density_den_cols,
                       dens_title, dens_sub, opt$chrBorderColor)

    dens_h <- if (!is.null(opt$density_height)) opt$density_height else
      min(50, max(3, nrow(chr_sub) * 0.45))
    dens_w <- if (!is.null(opt$density_width)) opt$density_width else 12

    message("Saving: ", dens_out)
    ggsave(dens_out, pd, width = dens_w, height = dens_h,
           dpi = opt$density_dpi, limitsize = FALSE,
           device = if (dens_ext == "svg") grDevices::svg else NULL)
    message("SNP-density plot saved: ", dens_out)
  }
}
