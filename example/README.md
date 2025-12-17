# Example Data

This directory contains example data files for testing the vcf_difplot.R script.

## Files

### example_data.table

A sample tab-delimited file that simulates the output from GATK VariantsToTable. This file contains:

- **CHROM**: Chromosome identifier (chr1-chr5)
- **POS**: Position on the chromosome
- **sample1.GT**: Genotype for sample1 (in GATK format: A/T, C|G, ./., etc.)
- **sample2.GT**: Genotype for sample2 (in GATK format: A/T, C|G, ./., etc.)

The data includes 17 positions across 5 chromosomes with varying genotypes, including:
- Homozygous genotypes (e.g., A/A, G|G)
- Heterozygous genotypes (e.g., A/T, C|G)
- Missing data (./.)
- Both phased (|) and unphased (/) separators

### chr_lengths.txt

A tab-delimited file containing chromosome lengths for human chromosomes 1-5. Format:

```
chromosome_name<TAB>length_in_bp
```

These lengths are based on the GRCh38 human reference genome.

### chr_lengths_comma.txt

The same chromosome length data but comma-delimited, to demonstrate automatic separator detection:

```
chromosome_name,length_in_bp
```

## Running Examples

### Example 1: Basic Plot Using Sample Names

```bash
Rscript ../vcf_difplot.R \
    -i example_data.table \
    -b sample1 \
    -c sample2 \
    -o basic_plot.pdf
```

This will:
- Read the example data file
- Compare sample1 (baseline) vs sample2 (comparison)
- Generate a PDF plot showing variant positions
- Use auto-detected chromosome lengths (with warning)

### Example 2: Plot with Chromosome Length File

```bash
Rscript ../vcf_difplot.R \
    -i example_data.table \
    -b sample1 \
    -c sample2 \
    -l chr_lengths.txt \
    -o plot_with_lengths.pdf
```

This will:
- Use the provided chromosome length file
- Generate more accurate visualization with proper chromosome scales

### Example 3: Using Column Positions

```bash
Rscript ../vcf_difplot.R \
    -i example_data.table \
    -B 1 \
    -C 2 \
    -l chr_lengths.txt \
    -o plot_by_column.pdf
```

This will:
- Select baseline sample using column position (1st GT column = sample1)
- Select comparison sample using column position (2nd GT column = sample2)

### Example 4: PNG Output

```bash
Rscript ../vcf_difplot.R \
    -i example_data.table \
    -b sample1 \
    -c sample2 \
    -l chr_lengths.txt \
    -o plot_output.png
```

This generates a PNG image instead of PDF.

### Example 5: Custom Unit (kb instead of Mb)

```bash
Rscript ../vcf_difplot.R \
    -i example_data.table \
    -b sample1 \
    -c sample2 \
    -l chr_lengths.txt \
    -u 1000 \
    -o plot_kb_units.pdf
```

This displays positions in kilobases (kb) instead of megabases (Mb).

### Example 6: Filter for Homozygous Baseline

```bash
Rscript ../vcf_difplot.R \
    -i example_data.table \
    -b sample1 \
    -c sample2 \
    -l chr_lengths.txt \
    --baseHetcheck \
    -o plot_homozygous_base.pdf
```

This includes only positions where sample1 (baseline) is homozygous.

### Example 7: Filter for Both Homozygous

```bash
Rscript ../vcf_difplot.R \
    -i example_data.table \
    -b sample1 \
    -c sample2 \
    -l chr_lengths.txt \
    --baseHetcheck \
    --copHetcheck \
    -o plot_both_homozygous.pdf
```

This includes only positions where both samples are homozygous.

### Example 8: Custom Colors and Line Thickness

```bash
Rscript ../vcf_difplot.R \
    -i example_data.table \
    -b sample1 \
    -c sample2 \
    -l chr_lengths.txt \
    --segmentColor blue \
    --segmentSize 1.0 \
    --chrBorderColor darkgray \
    --chrBorderSize 0.5 \
    -o plot_custom_style.pdf
```

This demonstrates custom visualization with blue variant segments and dark gray borders.

### Example 9: Using Comma-Delimited Chromosome Length File

```bash
Rscript ../vcf_difplot.R \
    -i example_data.table \
    -b sample1 \
    -c sample2 \
    -l chr_lengths_comma.txt \
    -o plot_with_comma_sep.pdf
```

This demonstrates automatic separator detection using the comma-delimited chromosome length file.

## Expected Output

The generated plot will show:
- 5 horizontal gray bars (one for each chromosome)
- Red vertical lines indicating positions where sample1 and sample2 have different genotypes
- X-axis showing position (in Mb by default)
- Y-axis showing chromosome names

The console output will display:
- Summary statistics (total positions, variant/non-variant counts)
- **Table of first 20 positions** that meet the filtering criteria, showing:
  - CHROM: Chromosome name
  - POS: Position on chromosome
  - Baseline_GT: Genotype of baseline sample
  - Comparison_GT: Genotype of comparison sample
- Processing information and warnings

For the example data, you should see:

**Without filters (basic run):**
- Total positions: 15 (after removing 2 positions with ./.)
- Variant positions vary based on genotype differences
- Positions with ./. in either column are automatically excluded

**With --baseHetcheck:**
- Only positions where sample1 is homozygous (A/A, G/G, T|T, etc.)
- Excludes heterozygous positions (A/T, C/G, etc.)

**With both --baseHetcheck and --copHetcheck:**
- Only positions where both samples are homozygous
- Smallest subset of positions

Note: The script normalizes genotypes, so A/T and T|A are treated as equivalent.

## Creating Your Own Test Data

To create your own test data from a real VCF file:

```bash
gatk VariantsToTable \
   -V your_file.vcf \
   -F CHROM -F POS -GF GT \
   -O your_data.table
```

Then run the plotting script on your data:

```bash
Rscript ../vcf_difplot.R \
    -i your_data.table \
    -b your_baseline_sample \
    -c your_comparison_sample \
    -o your_plot.pdf
```
