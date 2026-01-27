# Auditing RNA-seq Data Quality Using Unsupervised Learning
A Case Study on GTEx Brain-Cortex Samples

## Project Overview

High-throughput RNA-seq datasets are widely reused for biological discovery, yet they often contain hidden technical and biological anomalies that can distort downstream analysis. This project treats RNA-seq analysis not as prediction, but as data auditing.

>**Core Question:**
>*Do RNA-seq samples behave as expected, given known biology and sequencing principles?*

Using GTEx Brain Cortex RNA-seq data, this project combines standard RNA-seq preprocessing, exploratory data analysis, and unsupervised machine learning to identify samples that deviate from expected behavior and to evaluate how these anomalies affect biological conclusions.

## Objectives

- Perform rigorous preprocessing and normalization of RNA-seq data  
- Explore global and local structure using PCA and UMAP  
- Detect anomalous samples using unsupervised machine learning  
- Validate anomalies using biological metadata and marker genes  
- Assess the impact of anomalies on downstream differential expression

## Dataset

Source: GTEx v8 (Genotype-Tissue Expression Project)

Tissue: Brain – Cortex

Data Types Used:
- Gene-level RNA-seq read counts  
- Sample attributes  
- Subject phenotypes

Subset used in this study:
- 70 Brain Cortex RNA-seq samples  
- Balanced across sex  
- Mostly unique donors (69 donors)

## Workflow Summary

1. Metadata filtering:
   - Select Brain Cortex samples  
   - Merge sample-level and subject-level metadata  
   - Remove samples with missing key attributes

2. Gene count processing:
   - Load GTEx gene-level counts  
   - Subset to selected samples

3. Gene Filtering:
   - Remove lowly expressed genes to reduce noise
   
4. Normalization:
   - Apply variance-stabilizing or log-based normalization

5. Exploratory analysis:
   - PCA for global structure  
   - UMAP for local neighborhood structure

6. Anomaly detection:
   - Isolation Forest (global anomalies)  
   - Local Outlier Factor (local inconsistencies)  
   - Consensus anomaly scoring

7. Biological validation:
   - Examine known marker genes  
   - Inspect donor-level metadata  
   - Interpret anomalies biologically

8. Impact assessment:
   - Differential expression before vs after anomaly removal

## Why Unsupervised Learning?

- No ground-truth labels exist for “bad” RNA-seq samples  
- Anomalies are rare and heterogeneous  
- Unsupervised methods:
  - scale well to high-dimensional gene expression
  - do not require labels
  - complement classical QC metrics

Machine learning is used as an assistant, not a decision-maker, since all flagged samples are biologically interpreted.

## Methods & Tools Used

This project combines established RNA-seq processing practices with interpretable, unsupervised machine learning methods. All tools were selected to prioritize robustness, transparency, and biological interpretability.

**1. RNA-seq Processing & Statistical Analysis**
- R
- tidyverse
- edgeR
- DESeq2

**2. Exploratory Data Analysis**
- Principal Component Analysis (PCA)
- UMAP

**3. Unsupervised Anomaly Detection**
- Isolation Forest
- Local Outlier Factor (LOF)
- Consensus Anomaly scoring

**4. Data Source**
- GTEx v8

## Key Takeaway

This project demonstrates how unsupervised learning can be applied responsibly to audit RNA-seq datasets, improving data reliability and downstream biological interpretation.
