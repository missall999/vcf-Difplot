# VCF Variant Position Difference Plot

An R script project for plotting variant positions from VCF files, comparing genotypes between two samples.

## Overview

This tool reads a tab-delimited file converted from VCF format and creates visualizations showing variant positions across chromosomes. It compares genotypes between a baseline sample and a comparison sample, highlighting positions where they differ.

## Prerequisites

- R (version 3.6 or higher)
- R packages:
  - `ggplot2` (>= 3.4.0)
  - `optparse`
  - `data.table`
  - `R.utils` (optional, required only for reading gzipped input files)
- GATK (for converting VCF to tab-delimited format)

## Installation

Install required R packages:

```r
install.packages(c("ggplot2", "optparse", "data.table"))
# Optional: for gzipped input support
install.packages("R.utils")
```

## Workflow

### Step 1: Convert VCF to Tab-Delimited Format

Use GATK's VariantsToTable to convert your VCF file to a tab-delimited format:

```bash
gatk VariantsToTable \
   -V input.vcf \
   -F CHROM -F POS -GF GT \
   -O output.table
```

This command extracts:
- `CHROM`: Chromosome name
- `POS`: Position
- `GT`: Genotype for each sample (creates columns like `sampleID.GT`)

### Step 2: Plot Variant Positions

Run the R script to generate the plot:

```bash
Rscript vcf_difplot.R -i output.table -b baseline_sample -c comparison_sample -o variant_plot.pdf
```

## Usage

```bash
Rscript vcf_difplot.R [options]
```

### Required Arguments

- `-i, --input FILE`: Input tab-delimited file (required)

### Sample Selection (choose one method for each)

**Baseline Sample:**
- `-b, --basename NAME`: Baseline sample name
- `-B, --basecol INT`: Baseline GT column index (1-based, among GT columns only)

**Comparison Sample:**
- `-c, --copname NAME`: Comparison sample name
- `-C, --copcol INT`: Comparison GT column index (1-based, among GT columns only)

> **Note:** The column index refers to the position among GT columns only, not all columns in the file. For example, if the file has columns `CHROM POS sample1.GT sample2.GT`, then `-B 1` selects `sample1.GT` and `-B 2` selects `sample2.GT`.

### Optional Arguments

- `-o, --output FILE`: Output plot file (default: `variant_plot.pdf`). Supports PDF, PNG, JPEG, SVG.
- `-l, --chrlength FILE`: Chromosome length file (first two columns used: CHROM, LENGTH)
  - Supports `.fai` format (samtools faidx output, 5 columns) — only the first two are read
  - Separator is automatically detected (tab, comma, semicolon, or whitespace)
  - If not provided, uses maximum observed position per chromosome (warning issued)
- `-u, --unit NUM`: Position unit divisor (default: 1e6 for Mb; use 1e3 for kb, 1 for bp)
- `--baseHetcheck`: Only include positions where baseline is homozygous
- `--copHetcheck`: Only include positions where comparison is homozygous
- `--output_table FILE`: Write variant positions (CHROM, POS) to a tab-delimited file, sorted by genomic coordinate
- `-I, --interactive`: Interactive mode — prompts for each parameter with descriptions

### Visualization Customization

- `--segmentColor COLOR`: Color for variant segments (default: `red`)
- `--segmentSize NUM`: Thickness of variant segments (default: `0.5`)
- `--segmentAlpha NUM`: Transparency of variant segments, 0-1 (default: `0.6`)
- `--chrBorderColor COLOR`: Color for chromosome borders (default: `black`)
- `--chrBorderSize NUM`: Thickness of chromosome borders (default: `0.3`)

All color parameters accept R color names (case-insensitive) or hex codes (e.g., `"#FF5733"`).

### Genotype Handling

The script properly handles GATK VariantsToTable genotype formats:
- Supports both `/` and `|` as separators (phased and unphased)
- Treats `A/T` and `T|A` as equivalent (normalizes by sorting alleles)
- Automatically filters out positions with missing data (`./.`), wildcards (`*/*`), or malformed genotypes
- **Haploid support**: Single-allele genotypes (e.g., male chrX/chrY `A`) are expanded to homozygous diploid (`A/A`) for correct comparison
- **Polyploid support**: Genotypes with 3+ alleles are sorted and compared as allele sets
- Can optionally filter for homozygous positions only
- **Smart chromosome sorting**: Natural genomic order (Chr1, Chr2, ..., Chr10, X, Y, MT)
- **Y-axis orientation**: Chr1 at top (conventional genome browser layout)

## Example

```bash
Rscript vcf_difplot.R \
  -i variants.table \
  -b sample1 \
  -c sample2 \
  -o comparison.pdf \
  -l chr_lengths.txt \
  -u 1000000 \
  --baseHetcheck \
  --copHetcheck \
  --segmentColor red \
  --segmentSize 0.5 \
  --segmentAlpha 0.6 \
  --chrBorderColor black \
  --chrBorderSize 0.3 \
  --output_table positions.tsv
```

### Chromosome Length File Format

The first two columns are used (CHROM, LENGTH). Extra columns are ignored, so `.fai` files work directly:

**Tab-delimited:**
```
chr1	248956422
chr2	242193529
chr3	198295559
```

**samtools .fai (5 columns — only first two used):**
```
chr1	248956422	112	80	81
chr2	242193529	252092603	80	81
```

**Comma-delimited:**
```
chr1,248956422
chr2,242193529
```

The script automatically detects the separator.

## Output

The script generates a plot where:
- Each chromosome is a horizontal rectangle (light gray), with chr1 at the top
- Variant positions (where genotypes differ) are vertical lines (red by default)
- The x-axis shows position (scaled by the specified unit)
- The y-axis lists chromosomes in genomic order

Console output includes:
- Summary statistics (total positions, variants, non-variants)
- First 20 variant positions with genotypes
- Chromosome length information
- Warnings for data issues (positions beyond chromosome length, haploid/polyploid genotypes, etc.)

![image](example/output.png)

## Features

- **Automatic Sample Detection**: Reads GT column names to identify available samples
- **Flexible Sample Selection**: Specify samples by name or GT column index
- **Chromosome Length Handling**: Supports custom length files (including .fai) or auto-detects from data
- **Multiple Output Formats**: PDF, PNG, JPEG, SVG (via ggsave)
- **Robust Error Handling**: Early validation of all parameters before expensive processing
- **Memory Efficient**: Only reads needed columns from input file
- **Fast Genotype Parsing**: Processes unique GT values only, then maps back (10-50x faster on large datasets)
- **Interactive Mode**: Guided parameter setup with equivalent command output

## Error Handling

The script validates early and fails fast:
- ggplot2 version check at startup
- Output format validation before data processing
- Color and numeric parameter validation
- Input file existence and column structure
- Sample name/index resolution with clear error messages
- Same-sample detection (stops with error)
- Chromosome length file validation (duplicates, non-numeric, header detection)
- POS overflow detection (warns if variants exceed chromosome length)

## Notes

- GT columns must be named in the format `sampleID.GT`
- Missing and malformed genotypes are automatically excluded
- Chromosomes are sorted in genomic order (numeric, then X, Y, MT)
- Plot height auto-adjusts based on chromosome count
- Gzipped input (.gz, .bgz) requires the `R.utils` package

## License

This project is open source and available for use and modification.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.
