# This is the initial data cleaning and normalization script for differential expression for the pilot experiment from Konrad's mouse organoids
# Written by KM on 20MAY25

# First load packages ----
# These are all from Bioconductor
library(edgeR) # For differential expression calculation. We will create an edgeR object and use the package to normalize our data.
library(AnnotationDbi) # To annotate with Entrez IDs
library(org.Mm.eg.db) # To annotate with Mouse Entrez IDs

# Read in the raw data and metadata ----
## Raw data
raw.data <- read.delim("../results/star_salmon/salmon.merged.gene_counts.tsv", check.names=FALSE, stringsAsFactors=FALSE)
### View the raw data structure
str(raw.data)

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

# Create edgeR DGEList Object ----
y <- DGEList(counts=raw.data,group=metadata$Condition)
## View our edgeR object
str(y)

# Annotating with ENTREZ IDs to do GO and pathway enrichment ----
# This is kind of complex, but to do pathway enrichment we need to add other gene IDs (called Entrez IDs)

# Need ENTREZ ID to be the row name, but cannot have duplicate row names (some genes have unique Ensembl IDs but not gene names/Entrez IDs)
# Finding the Ensembl IDs currently in the database 
idfound<-y$genes$gene_id %in% mappedRkeys(org.Mm.egENSEMBL)
y<-y[idfound,]
dim(y)

# Then adding ENTREZ ID
egENSEMBL<-toTable(org.Mm.egENSEMBL)
head(egENSEMBL)

m<-match(y$genes$gene_id, egENSEMBL$ensembl_id)
y$genes$EntrezGene<-egENSEMBL$gene_id[m]
head(y$genes)

# Dropping duplcated genes (both gene names and Entrez IDs)
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$gene_name)
y.unique <- y[!d,]
e <- duplicated(y.unique$genes$EntrezGene)
y.unique <- y.unique[!e,]
nrow(y.unique)

# Making the Entrez ID the row name (and dropping duplicative data from EntrezGene column)
rownames(y.unique$counts) <- rownames(y.unique$genes) <- y.unique$genes$EntrezGene
y.unique$genes$EntrezGene <- NULL

# Filter out lowly expressed genes ----
keep <- filterByExpr(y)
## This is a logical vector
## We can see the number of genes kept by counting how many are `TRUE`
## And how many dropped by counting `FALSE`
summary(keep)

# To filter out some noise, we are going to set a stricter filtering threshold
keep2 <- rowSums( cpm(y) > 0.5 ) >=2
summary(keep2)

# We will take the `keep2` filtering threshold going forward


# Re-calculate library size after filtering
y <- y[keep2,,keep.lib.sizes=FALSE]
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
saveRDS(y, file = "input/normalized_data.rds")




