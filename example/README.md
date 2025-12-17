# Example Data

This directory contains example data files for testing the vcf_difplot.R script.

## Files

### example_data.table

A sample tab-delimited file that simulates the output from GATK VariantsToTable. This file contains:

- **CHROM**: Chromosome identifier (chr1-chr5)
- **POS**: Position on the chromosome
- **sample1.GT**: Genotype for sample1
- **sample2.GT**: Genotype for sample2

The data includes 15 variant positions across 5 chromosomes with varying genotypes.

### chr_lengths.txt

A tab-delimited file containing chromosome lengths for human chromosomes 1-5. Format:

```
chromosome_name<TAB>length_in_bp
```

These lengths are based on the GRCh38 human reference genome.

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

## Expected Output

The generated plot will show:
- 5 horizontal gray bars (one for each chromosome)
- Red vertical lines indicating positions where sample1 and sample2 have different genotypes
- X-axis showing position (in Mb by default)
- Y-axis showing chromosome names

For the example data, you should see:
- chr1: 3 variants (at positions 100kb, 500kb, 1000kb)
- chr2: 2 variants (at positions 150kb, 600kb)
- chr3: 2 variants (at positions 200kb, 800kb)
- chr4: 2 variants (at positions 100kb, 500kb)
- chr5: 1 variant (at position 250kb)

Total: 10 variant positions out of 15 total positions in the file.

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
