# VCF Variant Position Difference Plot

An R script project for plotting variant positions from VCF files, comparing genotypes between two samples.

## Overview

This tool reads a tab-delimited file converted from VCF format and creates visualizations showing variant positions across chromosomes. It compares genotypes between a baseline sample and a comparison sample, highlighting positions where they differ.

## Prerequisites

- R (version 3.6 or higher)
- R packages:
  - `ggplot2`
  - `optparse`
- GATK (for converting VCF to tab-delimited format)

## Installation

Install required R packages:

```r
install.packages(c("ggplot2", "optparse"))
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
- `-B, --basecol INT`: Baseline sample column position (1-based)

**Comparison Sample:**
- `-c, --copname NAME`: Comparison sample name
- `-C, --copcol INT`: Comparison sample column position (1-based)

### Optional Arguments

- `-o, --output FILE`: Output plot file (default: `variant_plot.pdf`)
- `-l, --chrlength FILE`: Chromosome length file (CHROM LENGTH)
  - If not provided, uses maximum variant position (warning will be issued)
  - Separator is automatically detected (supports tab, comma, semicolon, or whitespace)
- `-u, --unit NUM`: Chromosome length unit (default: 1e6 for Mb)
- `--baseHetcheck`: Check if baseline sample is homozygous; ignore heterozygous positions
  - Only positions where baseline is homozygous (e.g., A/A, G|G) will be included
- `--copHetcheck`: Check if comparison sample is homozygous; ignore heterozygous positions
  - Only positions where comparison is homozygous (e.g., A/A, G|G) will be included

### Visualization Customization

- `--segmentColor COLOR`: Color for variant position segments (default: `red`)
  - Accepts any valid R color name or hex code (e.g., "blue", "#FF5733")
- `--segmentSize NUM`: Thickness of variant position segments (default: `0.5`)
- `--chrBorderColor COLOR`: Color for chromosome borders (default: `black`)
  - Accepts any valid R color name or hex code
- `--chrBorderSize NUM`: Thickness of chromosome borders (default: `0.3`)

### Genotype Handling

The script properly handles GATK VariantsToTable genotype formats:
- Supports both `/` and `|` as separators (phased and unphased)
- Treats `A/T` and `T|A` as equivalent (normalizes for comparison)
- Automatically filters out positions with missing data (`./.`)
- Can optionally filter for homozygous positions only

## Examples

### Example 1: Using Sample Names

```bash
Rscript vcf_difplot.R \
  -i variants.table \
  -b sample1 \
  -c sample2 \
  -o comparison.pdf
```

### Example 2: Using Column Positions

```bash
Rscript vcf_difplot.R \
  -i variants.table \
  -B 1 \
  -C 2 \
  -o comparison.pdf
```

### Example 3: With Chromosome Length File and Custom Unit (kb)

```bash
Rscript vcf_difplot.R \
  -i variants.table \
  -b sample1 \
  -c sample2 \
  -l chr_lengths.txt \
  -u 1000 \
  -o comparison.pdf
```

This example uses kilobase (kb) units instead of the default megabase (Mb) units.

### Example 4: Filter for Homozygous Baseline Positions

```bash
Rscript vcf_difplot.R \
  -i variants.table \
  -b sample1 \
  -c sample2 \
  --baseHetcheck \
  -o homozygous_baseline.pdf
```

This example only includes positions where the baseline sample is homozygous (e.g., A/A, G|G).

### Example 5: Filter for Both Homozygous Positions

```bash
Rscript vcf_difplot.R \
  -i variants.table \
  -b sample1 \
  -c sample2 \
  --baseHetcheck \
  --copHetcheck \
  -o both_homozygous.pdf
```

This example only includes positions where both samples are homozygous.

### Example 6: Custom Colors and Line Thickness

```bash
Rscript vcf_difplot.R \
  -i variants.table \
  -b sample1 \
  -c sample2 \
  --segmentColor blue \
  --segmentSize 1.0 \
  --chrBorderColor darkgray \
  --chrBorderSize 0.5 \
  -o custom_colors.pdf
```

This example customizes the appearance with blue variant segments (thicker) and dark gray chromosome borders.

### Example 7: Using Hex Color Codes

```bash
Rscript vcf_difplot.R \
  -i variants.table \
  -b sample1 \
  -c sample2 \
  --segmentColor "#FF5733" \
  --segmentSize 0.8 \
  -o hex_colors.pdf
```

This example uses a hex color code for the variant segments.

### Chromosome Length File Format

The chromosome length file should have two columns (no header). The separator is automatically detected (tab, comma, semicolon, or whitespace):

**Tab-delimited:**
```
chr1	248956422
chr2	242193529
chr3	198295559
...
```

**Comma-delimited:**
```
chr1,248956422
chr2,242193529
chr3,198295559
```

**Space-delimited:**
```
chr1 248956422
chr2 242193529
chr3 198295559
```

The script will automatically detect and use the appropriate separator.

## Output

The script generates a plot where:
- Each chromosome is represented as a horizontal rectangle (light gray by default)
- Variant positions (where genotypes differ) are shown as vertical lines (red by default)
- The x-axis shows position (scaled by the specified unit)
- The y-axis lists chromosomes

## Features

- **Automatic Sample Detection**: Reads GT column names to identify available samples
- **Flexible Sample Selection**: Specify samples by name or column position
- **Chromosome Length Handling**: Supports custom length files or auto-detects from data
- **Multiple Output Formats**: Supports PDF, PNG, and JPEG
- **Robust Error Handling**: Validates input files, parameters, and data structure
- **Informative Messages**: Provides detailed progress and summary information

## Error Handling

The script includes comprehensive error checking for:
- Missing or inaccessible input files
- Missing required columns (CHROM, POS, GT)
- Invalid sample names or column positions
- Missing chromosome length data
- Invalid file formats

## Notes

- GT columns must be named in the format `sampleID.GT`
- Missing genotypes are automatically excluded from comparison
- Chromosomes are sorted numerically when possible, alphabetically otherwise
- The plot height automatically adjusts based on the number of chromosomes

## License

This project is open source and available for use and modification.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.