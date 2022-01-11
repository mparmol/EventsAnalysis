
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **README**

## **R code for: Parras-Moltó & , 2022**

This document describes the use of
Antibiotic_resistance_genes_analyzer.R. The evolutionary process of
Antibiotic resistance genes (ARGs) could be studied by the flow of
genetic information between different phyla bacteria, what we call as
events. Events are described as node-points within the tree with two
branches that belong to at least two different phyla. Within this script
we will use a tree based on the sequence aligment of ARGs of each class
and their related taxonomy to define where events happens.

    R version 4.0.3 (2020-10-10)

### **Dependencies**

#### **Conda**

A Conda environment is required to run the script.

Miniconda3 installation guide:

``` r
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
```

Environment installation:

``` r
conda env create -f Antibiotic_analyze.yml
```

#### **R Packages**

· `adephylo` v1.1.11. This package is devoted to exploratory analysis of
phylogenetic comparative data. It re-implements and extends phylogenetic
procedures from the ade4 package (which are now deprecated).

· `ape` v5.5. ape provides functions for reading, writing, manipulating,
analysing, and simulating phylogenetic trees and DNA sequences,
computing DNA distances, translating into AA sequences, estimating trees
with distance-based methods, and a range of methods for comparative
analyses and analysis of diversification. Functionalities are also
provided for programming new phylogenetic methods.

· `cowplot` v1.1.1. The cowplot package is a simple add-on to ggplot. It
provides various features that help with creating publication-quality
figures, such as a set of themes, functions to align plots and arrange
them into complex compound figures, and functions that make it easy to
annotate plots and or mix plots with images.

· `data.table` v1.14.0. data.table inherits from data.frame. It offers
fast and memory efficient: file reader and writer, aggregations,
updates, equi, non-equi, rolling, range and interval joins, in a short
and flexible syntax, for faster development.

· `genbankr` v1.18.0. Provides utilities to parse GenBank and GenPept
files into data structures which integrate with the Bioconductor
ecosystem..

· `ggnewscale` v0.4.5. Use multiple fill and colour scales in ’ggplot2.

· `ggplot2` v3.3.5. A system for ‘declaratively’ creating graphics,
based on “The Grammar of Graphics”. You provide the data, tell ‘ggplot2’
how to map variables to aesthetics, what graphical primitives to use,
and it takes care of the details.

· `ggtree` v2.4.2. Visualizing phylogenetic tree and heterogenous
associated data based on grammar of graphics ggtree provides functions
for visualizing phylogenetic tree and its associated data in R.

· `gtools` v3.9.2. The gtools R package provides functions to assist in
R programming.

· `phangorn` v2.7.1. Phylogenetic analysis in R (Estimation of
phylogenetic trees and networks using Maximum Likelihood, Maximum
Parsimony, Distance methods & Hadamard conjugation).

· `randomcoloR` v1.1.0.1. Simple methods to generate attractive random
colors.

· `reshape2` v1.4.4. Flexibly restructure and aggregate data using just
two functions: melt and ‘dcast’ (or ‘acast’).

· `stringr` v1.4.0. A consistent, simple and easy to use set of wrappers
around the fantastic ‘stringi’ package. All function and argument names
(and positions) are consistent, all functions deal with “NA”’s and zero
length vectors in the same way, and the output from one function is easy
to feed into the input of another.

· `taxonomizr` v0.8.0. Functions for assigning taxonomy to NCBI
accession numbers and taxon IDs based on NCBI’s accession2taxid and
taxdump files. This package allows the user to downloads NCBI data dumps
and create a local database for fast and local taxonomic assignment.

#### **Databases creation**

``` r
### 1 - Taxonomy database

/storage/parras/Taxa_ncbi/accessionTaxa.sql

### 2 - Get accession number and names for each analyzed sequence

/home/parras/RESULTADOS_FINALES_NO_BORRAR/total_ids.txt

### 3 - Get an original dataframe from all sequence taxonomy

/home/parras/RESULTADOS_FINALES_NO_BORRAR/complete_taxa_file.txt

### 4 - Download all information related to each sequence, like isolation_source

/storage/parras/databaseR/Tablas_taxa/Total_isolation.txt

### 5 - Download ResFinder database file (nucleotide format db)

/storage/parras/ResFinder_db/AntibioticDatabase
```

#### **Auxiliary files**

1- isolation_criteria.txt. Each isolation cluster is defined in this
file.

``` r
/home/parras/RESULTADOS_FINALES_NO_BORRAR/isolation_criteria.txt
```

2- isolation_exclusion_criteria.txt. In this list we have all the terms
that are excluded from any group to avoid ambiguity.

``` r
/home/parras/RESULTADOS_FINALES_NO_BORRAR/isolation_exclusion_criteria.txt
```
