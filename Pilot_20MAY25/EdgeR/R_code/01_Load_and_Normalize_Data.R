# This is the initial data cleaning and normalization script for differential expression for the pilot experiment from Konrad's mouse organoids
# Written by KM on 20MAY25

# First load packages ----
library(edgeR) # For differential expression calculation
library(dplyr) # For general data manipulation
library(AnnotationDbi) # To annotate with Entrez IDs
library(org.Hs.eg.db) # To annotate with Entrez IDs
library(GO.db) # For Gene Ontology analysis 
library(ghibli) # For pretty colors
library(ggplot2) # For plots (e.g. Volcano plot)
pacman::p_load(here,  
               tidyverse, 
               janitor, # Cleaning column names  
               scales, # Transform axis scales   
               ggrepel) # Optimize plot label separation 
library(GseaVis) # For nicer GSEA plot

# Read in the raw data and metadata ----
## Raw data
raw.data <- read.delim("../results/star_salmon/salmon.merged.gene_counts.tsv", check.names=FALSE, stringsAsFactors=FALSE,row.names="gene_id")
## Metadata
metadata <- read.delim("input/organoid_pilot_metadata.txt",check.names=FALSE, stringsAsFactors=FALSE)
### View the metadata
metadata
### Make "Condition" (wildtype or mutant) a factor
metadata$Condition <- factor(metadata$Condition, levels=c("WT","Mut"))
### View the metadata
str(metadata)
table(metadata)
#### CAUTION 
#### edgeR metadata is assigned by the order of the rows in the raw gene expression table
#### This is fine for us because we want WT to be the reference, and WT was listed first
#### We need to be cautious that we are setting the appropriate reference

# Create edgeR Object ----
# Create edgeR list object (called DGEList)
y <- DGEList(counts=raw.data,group=metadata$Condition)
## View our edgeR object
str(y)

# Filter out lowly expressed genes ----
keep <- filterByExpr(y)
## This is a logical vector
## We can see the number of genes kept by counting how many are `TRUE`
## And how many dropped by counting `FALSE`
summary(keep)

# Re-calculate library size after filtering
y <- y[keep,,keep.lib.sizes=FALSE]
## Confirm that the genes we wanted to keep are kept
## You can do this by comparing the length of the `gene_name` variable matches the `TRUE` count from above
str(y)

# Normalization ----
## Calculates a scaling factor (trimmed mean of M-values, TMM)
## Multiplies this normalization factor by the original library size to get the effective (normalized) library size used for downstream analysis
## norm.factor < 1 = small number of high-count genes are monopolizing the library, library size is scaled down/counts are scaled up
## norm.factor > 1 = library size scaled up/counts scaled down 
y <- normLibSizes(y)
## You can see the norm.factors below:
y$samples

# Save the data ----
## Save the normalized data in RDS format
saveRDS(y, file = "input/clean_data.rds")