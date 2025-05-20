# This is the initial data cleaning and normalization script for differential expression for the pilot experiment from Konrad's mouse organoids
# Written by KM on 20MAY25

# First load packages ----
library(edgeR) # For differential expression calculation. We will create an edgeR object and use the package to normalize our data.

# Read in the raw data and metadata ----
## Raw data
raw.data <- read.delim("../results/star_salmon/salmon.merged.gene_counts.tsv", check.names=FALSE, stringsAsFactors=FALSE)
### View the raw data
view(raw.data)
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


# Setting model design, defining condition as the covariate of interest
## Defining the grouping variable (condition)
group <- factor(c("WT","Mut"))

design <- model.matrix(~group)
design.reduced <- matrix(1,2,1)
dgList <- estimateGLMCommonDisp(data, design.reduced,
                                method="deviance", robust = TRUE, subset=NULL)




# Annotating with ENTREZ IDs to do GO and pathway enrichment

# Need ENTREZ ID to be the row name, but cannot have duplicate row names (some genes have unique Ensembl IDs but not gene names/Entrez IDs)
# Finding the Ensembl IDs currently in the database 
idfound<-dgList$genes$gene_id %in% mappedRkeys(org.Mm.egENSEMBL)
dgList<-dgList[idfound,]
dim(dgList)

# Then adding ENTREZ ID
egENSEMBL<-toTable(org.Mm.egENSEMBL)
head(egENSEMBL)

m<-match(dgList$genes$gene_id, egENSEMBL$ensembl_id)
dgList$genes$EntrezGene<-egENSEMBL$gene_id[m]
head(dgList$genes)

# Dropping duplcated genes (both gene names and Entrez IDs)
o <- order(rowSums(dgList$counts), decreasing=TRUE)
dgList <- dgList[o,]
d <- duplicated(dgList$genes$gene_name)
dgList.unique <- dgList[!d,]
e <- duplicated(dgList.unique$genes$EntrezGene)
dgList.unique <- dgList.unique[!e,]
nrow(dgList.unique)

# Making the Entrez ID the row name (and dropping duplicative data from EntrezGene column)
rownames(dgList.unique$counts) <- rownames(dgList.unique$genes) <- dgList.unique$genes$EntrezGene
dgList.unique$genes$EntrezGene <- NULL

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
saveRDS(y, file = "input/normalized_data.rds")
