#!/bin/bash

# Test script for vcf_difplot.R
# This script runs various test cases to verify the functionality

echo "=== VCF Difplot Test Suite ==="
echo ""

# Check if R is installed
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R to run tests."
    exit 1
fi

# Check if required R packages are installed
echo "Checking R package dependencies..."
Rscript -e "if (!require('ggplot2', quietly=TRUE)) stop('ggplot2 package not installed')" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: ggplot2 package not installed. Install with: install.packages('ggplot2')"
    exit 1
fi

Rscript -e "if (!require('optparse', quietly=TRUE)) stop('optparse package not installed')" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: optparse package not installed. Install with: install.packages('optparse')"
    exit 1
fi

echo "✓ All required packages installed"
echo ""

# Test 1: Help message
echo "Test 1: Help message"
Rscript vcf_difplot.R --help
echo ""

# Test 2: Missing input file (should fail)
echo "Test 2: Missing input file (should fail)"
Rscript vcf_difplot.R -b sample1 -c sample2 2>&1 | grep -q "Input file is required"
if [ $? -eq 0 ]; then
    echo "✓ Test passed: Correctly reports missing input file"
else
    echo "✗ Test failed: Should report missing input file"
fi
echo ""

# Test 3: Non-existent input file (should fail)
echo "Test 3: Non-existent input file (should fail)"
Rscript vcf_difplot.R -i nonexistent.table -b sample1 -c sample2 2>&1 | grep -q "does not exist"
if [ $? -eq 0 ]; then
    echo "✓ Test passed: Correctly reports non-existent file"
else
    echo "✗ Test failed: Should report non-existent file"
fi
echo ""

# Test 4: Missing baseline sample parameter (should fail)
echo "Test 4: Missing baseline sample parameter (should fail)"
Rscript vcf_difplot.R -i example/example_data.table -c sample2 2>&1 | grep -q "baseline sample"
if [ $? -eq 0 ]; then
    echo "✓ Test passed: Correctly reports missing baseline sample"
else
    echo "✗ Test failed: Should report missing baseline sample"
fi
echo ""

# Test 5: Missing comparison sample parameter (should fail)
echo "Test 5: Missing comparison sample parameter (should fail)"
Rscript vcf_difplot.R -i example/example_data.table -b sample1 2>&1 | grep -q "comparison sample"
if [ $? -eq 0 ]; then
    echo "✓ Test passed: Correctly reports missing comparison sample"
else
    echo "✗ Test failed: Should report missing comparison sample"
fi
echo ""

# Test 6: Basic usage with sample names
echo "Test 6: Basic usage with sample names"
Rscript vcf_difplot.R \
    -i example/example_data.table \
    -b sample1 \
    -c sample2 \
    -o test_output1.pdf
if [ -f test_output1.pdf ]; then
    echo "✓ Test passed: Plot generated successfully"
    rm test_output1.pdf
else
    echo "✗ Test failed: Plot not generated"
fi
echo ""

# Test 7: Usage with column positions
echo "Test 7: Usage with column positions"
Rscript vcf_difplot.R \
    -i example/example_data.table \
    -B 1 \
    -C 2 \
    -o test_output2.pdf
if [ -f test_output2.pdf ]; then
    echo "✓ Test passed: Plot generated with column positions"
    rm test_output2.pdf
else
    echo "✗ Test failed: Plot not generated with column positions"
fi
echo ""

# Test 8: Usage with chromosome length file
echo "Test 8: Usage with chromosome length file"
Rscript vcf_difplot.R \
    -i example/example_data.table \
    -b sample1 \
    -c sample2 \
    -l example/chr_lengths.txt \
    -o test_output3.pdf
if [ -f test_output3.pdf ]; then
    echo "✓ Test passed: Plot generated with chromosome length file"
    rm test_output3.pdf
else
    echo "✗ Test failed: Plot not generated with chromosome length file"
fi
echo ""

# Test 9: PNG output format
echo "Test 9: PNG output format"
Rscript vcf_difplot.R \
    -i example/example_data.table \
    -b sample1 \
    -c sample2 \
    -o test_output4.png
if [ -f test_output4.png ]; then
    echo "✓ Test passed: PNG plot generated"
    rm test_output4.png
else
    echo "✗ Test failed: PNG plot not generated"
fi
echo ""

# Test 10: Invalid sample name (should fail)
echo "Test 10: Invalid sample name (should fail)"
Rscript vcf_difplot.R \
    -i example/example_data.table \
    -b invalid_sample \
    -c sample2 2>&1 | grep -q "not found"
if [ $? -eq 0 ]; then
    echo "✓ Test passed: Correctly reports invalid sample name"
else
    echo "✗ Test failed: Should report invalid sample name"
fi
echo ""

echo "=== Test Suite Complete ==="
