# Differential Gene Expression and Gene Annotation Analysis

## Repository Name

`dge-gene-annotation-analysis`

## Description

This repository contains Python implementations for differential gene expression (DGE) analysis, gene annotation enrichment analysis, and motif sequence finding using Position-Specific Scoring Matrices (PSSMs). The project covers statistical analysis of gene expression data, functional annotation of differentially expressed genes, and visualization of sequence motifs.

## Features

- **Differential Gene Expression Analysis**: Compute arithmetic means, variances, log fold changes, p-values, and perform multiple testing correction using Bonferroni adjustment.
- **Gene Functional Annotation Enrichment**: Associate Gene Ontology (GO) terms with differentially expressed genes, compute GO term frequencies, and perform hypergeometric tests for enrichment analysis.
- **Position-Specific Scoring Matrix (PSSM)**: Compute and visualize PSSMs to identify and analyze motifs in biological sequences.

## Files and Functions

### Differential Gene Expression Analysis

1. **Reading Data**:
   - `read_expression_data(filepath)`: Reads the expression data from a tab-separated file into a data structure.

2. **Calculations**:
   - `calculate_means_variances()`: Computes arithmetic means and variances for expression values across conditions and time points.
   - `get_LFC()`: Calculates log fold change (logFC) between conditions.
   - `get_pval()`: Computes p-values using a gene-specific t-test.
   - `adjust_pval()`: Applies Bonferroni correction to p-values.

3. **Filtering and Reporting**:
   - `filter_probes()`: Identifies top probes with the lowest adjusted p-values.
   - Report includes probe identifiers, gene symbols, log fold changes, and p-values.

### Gene Functional Annotation Enrichment Analysis

1. **GO Annotation**:
   - `associate_go_terms()`: Associates GO terms with differentially expressed genes.

2. **Annotation Analysis**:
   - `common_least_common_go_terms()`: Returns the most and least common GO identifiers.
   - `hypergeometric_test()`: Computes hypergeometric p-values for GO term enrichment.
   - `adjust_pval()`: Applies Bonferroni correction to hypergeometric p-values.

3. **Reporting**:
   - List of enriched biological processes with the lowest adjusted p-values.

### Position-Specific Scoring Matrix (PSSM)

1. **PSSM Computation**:
   - `get_FrequencyMatrix()`: Computes the frequency matrix from aligned sequences.
   - `get_CorrectedFrequencyMatrix()`: Computes the corrected frequency matrix with pseudo-weight adjustment.
   - `get_ScoringMatrix()`: Computes the scoring matrix from the corrected frequency matrix.

2. **Visualization**:
   - `plot_SequenceLogo()`: Plots a sequence logo for the frequency or scoring matrix.

## Getting Started

### Prerequisites

Make sure you have Python and the following libraries installed:

- `numpy`
- `scipy`
- `matplotlib`

Install these dependencies using pip:

```bash
pip install numpy scipy matplotlib
