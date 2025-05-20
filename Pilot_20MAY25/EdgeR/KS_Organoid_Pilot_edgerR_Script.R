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

# This is the initial data cleaning and normalization script for differential expression for the pilot experiment from Konrad's mouse organoids

# Read in the raw data and metadata ----
## Raw data
raw.data <- read.delim("../results/star_salmon/salmon.merged.gene_counts.tsv", check.names=FALSE, stringsAsFactors=FALSE,row.names="gene_id")
## Metadata
metadata <- read.delim("organoid_pilot_metadata.txt",check.names=FALSE, stringsAsFactors=FALSE)
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


# Create edgeR list object (called DGEList)
y <- DGEList(counts=raw.data,group=metadata$Condition)
## View our edgeR object
str(y)

# Filter out lowly expressed genes
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

# Normalization 
## Calculates a scaling factor (trimmed mean of M-values, TMM)
## Multiplies this normalization factor by the original library size to get the effective (normalized) library size used for downstream analysis
## norm.factor < 1 = small number of high-count genes are monopolizing the library, library size is scaled down/counts are scaled up
## norm.factor > 1 = library size scaled up/counts scaled down 
y <- normLibSizes(y)
## You can see the norm.factors below:
y$samples

# Setting model design, defining treatment as the covariate of interest 
design <- model.matrix(~metadata$Treatment)

# Standard analysis would be to fit a glm quasi-likelihood model like so...
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)

# ...But this experiment has no replicates since it's a pilot analysis so we are using the Likelihood Ratio Test as suggested by the vignette
# We first need to estimate the dispersion. Setting robust=T identifies outliers to the mean dispersion trend.
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

# And then we fit the LRT model
fit <- glmFit(y, design)
lrt.2vs1 <- glmLRT(fit, coef=2)
topTags(lrt.2vs1)


# Plotting multidimensional scaling plot (MDS) -- 2D scatterplot to illustrate expression differences between samples
# Using top 500 genes to calculate pairwise distances 
plotMDS(y, top=500)

# Transforming to log counts per million (undefined values are avoided and poorly-defined log FCs for low counts are shrunk to zero)
logcpm <- cpm(y, log=TRUE)

# Plotting log-fold change against log counts per million, with DEGs highlighted 
summary(decideTests(lrt.2vs1))
plotMD(lrt.2vs1)
abline(h=c(-1, 1), col="blue")

####

# Tumor (treated vs untreated) analysis ----
# Using option 3 from section 2.12 of the vignette (What to do if you have no replicates)
# Reloading data to make Ensembl ID the first column instead of the row name
raw.data <- read.delim("salmon.merged.gene_counts.tsv", check.names=FALSE, stringsAsFactors=FALSE)
# Metadata
## Just doing separate sample sheets for Tumor and Normal
metadata <- read.delim("lenvatinib_metadata_tumor.txt",check.names=FALSE, stringsAsFactors=FALSE)
metadata$Condition <- factor(metadata$Condition, levels=c("T","T_Tx")) # Actual (1v1)x2 experiment

raw.tum <- raw.data[c("gene_id","gene_name","T.C","T.100")]

group <- factor(c("T","T_Tx"))

dgList <- DGEList(counts=raw.tum,group=metadata$Condition)
## Default filtering
keep <- filterByExpr(dgList)
summary(keep)
dgList <- dgList[keep,,keep.lib.sizes=FALSE]
## Stricter filtering?
keep2 <- rowSums( cpm(dgList) > 0.5 ) >=2
summary(keep2)
dgList <- dgList[keep2,,keep.lib.sizes=FALSE]
dgList <- normLibSizes(dgList)
dgList$samples

design <- model.matrix(~group)
design.reduced <- matrix(1,2,1)
dgList <- estimateGLMCommonDisp(dgList, design.reduced,
                                method="deviance", robust = TRUE, subset=NULL)
fit <- glmFit(dgList, design)
## Normal untreated should be the ref, and then set coef 2:4 to correspond to N_tx, T, and T_tx, respectively
## Doing tumor only, tumor vs treated tumor
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

# Annotating with ENTREZ IDs to do GO and pathway enrichment

# Need ENTREZ ID to be the row name, but cannot have duplicate row names (some genes have unique Ensembl IDs but not gene names/Entrez IDs)
# Finding the Ensembl IDs currently in the database 
idfound<-dgList$genes$gene_id %in% mappedRkeys(org.Hs.egENSEMBL)
dgList<-dgList[idfound,]
dim(dgList)

# Then adding ENTREZ ID
egENSEMBL<-toTable(org.Hs.egENSEMBL)
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

# Running differential expression
group <- factor(c("T","T_Tx"))

dgList.unique <- DGEList(counts=raw.tum,group=metadata$Condition)
keep <- filterByExpr(dgList.unique)
summary(keep)
dgList.unique <- dgList.unique[keep,,keep.lib.sizes=FALSE]
## Stricter filtering
keep2 <- rowSums( cpm(dgList.unique) > 0.1 ) >=2
summary(keep2)
dgList.unique <- dgList.unique[keep2,,keep.lib.sizes=FALSE]
dgList.unique <- normLibSizes(dgList.unique)
dgList.unique$samples
#barplot(dgList.unique$samples$lib.size*1e-6, names=1:4, ylab="Library size (millions)")
## We got quite poor coverage on the normals; they are topping off the sequencing now

design <- model.matrix(~group)
design.reduced <- matrix(1,2,1)
dgList.unique <- estimateGLMCommonDisp(dgList.unique, design.reduced,
                                method="deviance", robust = TRUE, subset=NULL)
fit <- glmFit(dgList.unique, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)



# Plotting log-fold change against log counts per million, with DEGs highlighted 
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")


# Gene ontology (GO) and pathway analysis
go <- goana(lrt, species = "Hs")
topGO(go)
keg <- kegga(lrt, species = "Hs")
topKEGG(keg)


# IFN stuff
## Dmitriy had said that we might expect to see upregulation of Type I IFN
## Here is the Type-I IFN gene set list (GO:0034340)
my_string <- "10010,103,10410,10581,11140,11274,11277,117854,149095,201292,23170,253943,25939,25942,26058,26065,29072,29110,29883,3093,338376,340061,3429,3434,3439,3440,3441,3442,3443,3444,3445,3446,3447,3448,3449,3451,3452,3454,3455,3456,3467,3654,3661,3665,3716,389856,401397,406991,4321,4599,4615,4938,4939,4940,51441,54502,55337,56339,56832,57506,57663,5770,5771,5777,5781,58509,64135,64343,6472,6609,6672,6772,6773,7187,7297,7474,79594,81030,81844,84166,8519,8772,9023,90933,91445,9470,9636,9641,9755,9759"
ifn_entrez <- strsplit(my_string, ",")[[1]]
print(ifn_entrez)

IFNlist <- rownames(dgList.unique)%in%ifn_entrez
summary(IFNlist)

fry(dgList.unique, index=IFNlist, design=design)
barcodeplot(lrt2$table$logFC, index=IFNlist)

camera(dgList2, index=IFNlist, design=design2)

## Volcano plots ----
# Results for DEGs -- defining up/down regulated genes
res.sig <- topTags(lrt, n=nrow(lrt$table), adjust.method = "fdr")$table
res.sig$topDE <- "NA"
res.sig$topDE <- as.character(res.sig$topDE)

res.sig$topDE[res.sig$logFC > 1 & res.sig$FDR < 0.05] <- "Up"
res.sig$topDE[res.sig$logFC < -1 & res.sig$FDR < 0.05] <- "Down"
res.sig$topDE[res.sig$topDE!="Up" & res.sig$topDE!="Down"] <- "ns"
res.sig$topDE <- as.factor(res.sig$topDE)

res.sig %>% distinct(topDE) %>% pull()

write.csv(res.sig, "TumorT_TumorTx_sig_filtered.csv")

# Make nice colors
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")
ghibli_colors
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[7], ghibli_colors[4])

cols <- c("Up" = ghibli_colors[4], "Down" = ghibli_colors[3], "ns" = "grey") 
sizes <- c("Up" = 2, "Down" = 2, "ns" = 1) 
alphas <- c("Up" = 1, "Down" = 1, "ns" = 0.5)

ggplot(data=res.sig, aes(x=logFC, y=-log10(FDR), fill=topDE, size=topDE, alpha=topDE)) + 
  geom_point() +
  theme_minimal() +
  geom_hline(yintercept=1.3, linetype="dashed", color="gray") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed",
             color="gray") +
  #xlim(-10,10) +
  geom_point(shape=21,
             color="black") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  labs(title = "Gene expression changes comparing Lenvatinib-treated and untreated tumor tissue",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

head(res.sig)

# Good tutorial for customizing volcano: https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/ 

# Highlighting specific genes on volcano

# Cytokine-cytokine receptor interaction pathway genes
cytogenes <- c("ACVRL1","CXCL6","CSF3","CXCL8","CCL3L3","IL24","CXCR4","TNFRSF11B","CSF2RA","CXCL5","IL1RL1","CCL3","TNFSF11","CCL2","CCL18","IL10","IL11","IL33","TGFB2","CCL21","INHBA","IL17RB","BMP4","IL1A","CXCL10","CXCL11","IL6","CXCL12","IL1B","IL3RA","ACKR3")

up_cyto_genes <- res.sig %>% filter(gene_name %in% cytogenes  & res.sig$topDE=="Up")
down_cyto_genes <- res.sig %>% filter(gene_name %in% cytogenes  & res.sig$topDE=="Down")

# Up/down cytokine genes
updown_genes <- rbind(up_cyto_genes, down_cyto_genes)

ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_cyto_genes,
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + # For upregulated genes
  geom_point(data = down_cyto_genes,
             shape = 21,
             size = 2, 
             fill = "steelblue", 
             colour = "black") + # For downregulated genes
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = updown_genes,   # Combine both up and down gene lists
                   aes(label = gene_name), 
                   force = 2,
                   nudge_y = 1) + # Add labels last to appear as the top layer
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 5)),
  #                    limits = c(-8, 8)) +
 # scale_y_continuous(limits = c(0, 11)) +
  labs(title = "Gene expression changes comparing Lenvatinib-treated and untreated tumor tissue",
       subtitle = "Significantly differentially expressed genes in KEGG 'Cytokine-cytokine receptor interaction' pathway (adjusted pathway enrichment p-value=1.09E-09)",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin treated tumor") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

## Resized for figure
ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_cyto_genes,
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + # For upregulated genes
  geom_point(data = down_cyto_genes,
             shape = 21,
             size = 2, 
             fill = "steelblue", 
             colour = "black") + # For downregulated genes
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = updown_genes,   # Combine both up and down gene lists
                   aes(label = gene_name, size = 30 ), 
                   force = 4,
                   nudge_y = 2,
                   show.legend=F) + # Add labels last to appear as the top layer
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 5)),
  #                    limits = c(-8, 8)) +
  # scale_y_continuous(limits = c(0, 11)) +
  labs(title = "DEGs comparing Lenvatinib-treated and untreated tumor",
       subtitle = "KEGG 'Cytokine-cytokine receptor interaction' pathway (adj. p-value=1.09E-09)",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin treated tumor") +
  theme_bw(base_size = 25) + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
## 1500 x 900

## IL-17 pathway genes
IL17genes <- c("CXCL6","CSF3","CXCL8","MMP1","MMP3","PTGS2","MUC5B","MMP9","CXCL5","IL17RB","MUC5AC","CXCL10","MAPK11","IL6","IL1B","CCL2","S100A9","S100A8")

up_IL17_genes <- res.sig %>% filter(gene_name %in% IL17genes  & res.sig$topDE=="Up")
down_IL17_genes <- res.sig %>% filter(gene_name %in% IL17genes  & res.sig$topDE=="Down")

# Up/down cytokine genes
updown_genes <- rbind(up_IL17_genes, down_IL17_genes)

ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_IL17_genes, # For upregulated genes
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + 
  geom_point(data = down_IL17_genes, # For downregulated genes
             shape = 21,
             size = 2, 
             fill = "steelblue", 
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = updown_genes,   # Combine both up and down gene lists
                   aes(label = gene_name), 
                   force = 2,
                   nudge_y = 1) + # Add labels last to appear as the top layer
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 5)),
  #                    limits = c(-8, 8)) +
  # scale_y_continuous(limits = c(0, 11)) +
  labs(title = "Gene expression changes comparing Lenvatinib-treated and untreated tumor tissue",
       subtitle = "Significantly differentially expressed genes in KEGG 'IL-17 signaling pathway' (adjusted pathway enrichment p-value=1.09E-09)",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin treated tumor") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


## Pathways in cancer genes
cancergenes <- c("CXCL8","LAMC3","NOTCH4","FLT4","LAMA4","CXCR4","LAMC2","PTGS2","CSF2RA","DLL4","EDNRB","BDKRB2","BDKRB1","NQO1","PDGFRA","TGFB2","MMP1","WNT5A","FN1","WNT7A","IGF1","AXIN2","GNG11","MMP9","PGF","BMP4","IL6","CXCL12","COL4A1","IL3RA","GSTA1")

up_cancer_genes <- res.sig %>% filter(gene_name %in% cancergenes  & res.sig$topDE=="Up")
down_cancer_genes <- res.sig %>% filter(gene_name %in% cancergenes  & res.sig$topDE=="Down")

# Up/down cytokine genes
updown_genes <- rbind(up_cancer_genes, down_cancer_genes)

ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_cancer_genes, # For upregulated genes
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + 
  geom_point(data = down_cancer_genes, # For downregulated genes
             shape = 21,
             size = 2, 
             fill = "steelblue", 
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = updown_genes,   # Combine both up and down gene lists
                   aes(label = gene_name, size = 30 ), 
                   force = 4,
                   nudge_y = 2,
                   show.legend=F) + # Add labels last to appear as the top layer
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 5)),
  #                    limits = c(-8, 8)) +
  # scale_y_continuous(limits = c(0, 11)) +
  labs(title = "DEGs comparing Lenvatinib-treated and untreated tumor",
       subtitle = "KEGG 'Pathways in cancer' (adjusted pathway enrichment p-value=6.16E-04)",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin treated tumor") +
  theme_bw(base_size = 25) + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

## GSEA Results ----
# Ran GSEA on all the tumor up/downregulated genes using the classic GSEA software (downloaded locally) 
# Making nicer plots of the GSEA results
# Documentation here: https://junjunlab.github.io/gseavis-manual/parse-output-results-from-gsea-software.html

# Loading in the output from the GSAE done on all genes 
integrated <- readGseaFile(filePath = "GSEA/output/Lenvatinib_RNA_Tumor_All_Genes.GseaPreranked.1738084316379/")
str(integrated)

# check enrichment data
head(integrated$meta[1:3,1:5])
# check gene lists
head(integrated$glist)
# check terms
head(integrated$gset[1:3,])

# Pathways of interest
setid <- c("KEGG_PATHWAYS_IN_CANCER",
           "KEGG_CHEMOKINE_SIGNALING_PATHWAY")


# plot
gseaNb(filePath = "GSEA/output/Lenvatinib_RNA_Tumor_All_Genes.GseaPreranked.1738084316379/",
       geneSetID = setid[2])

lapply(1:2, function(x){
  gseaNb(filePath = integrated,
         geneSetID = setid[x],
         addPval = T,
         pvalX = 0.6,
         pvalY = 0.6)
}) -> plist
# combine
cowplot::plot_grid(plotlist = plist,nrow = 1,align = 'hv')

# new style plot
gseaNb(filePath = integrated,
       geneSetID = setid[1],
       newGsea = T,
       addGene = T,
       markTopgene = T)

gseaNb(filePath=integrated,
       geneSetID = setid[1],
       addPval=T,
       pvalX=0.8,
       pvalY=0.4)



# Normal (treated vs untreated) analysis ----
# Using option 3 from section 2.12 of the vignette (What to do if you have no replicates)
# Reloading data to make Ensembl ID the first column instead of the row name
raw.data <- read.delim("salmon.merged.gene_counts.tsv", check.names=FALSE, stringsAsFactors=FALSE)
# Metadata
## Just doing separate sample sheets for Tumor and Normal
metadata <- read.delim("lenvatinib_metadata_normal.txt",check.names=FALSE, stringsAsFactors=FALSE)
metadata$Condition <- factor(metadata$Condition, levels=c("N","N_Tx")) # Actual (1v1)x2 experiment

raw.norm <- raw.data[c("gene_id","gene_name","N.C","N.100")]

group <- factor(c("N","N_Tx"))

dgList <- DGEList(counts=raw.norm,group=metadata$Condition)
keep <- filterByExpr(dgList)
dgList <- dgList[keep,,keep.lib.sizes=FALSE]
## Stricter filtering?
keep2 <- rowSums( cpm(dgList) > 0.5 ) >=2
dgList <- dgList[keep2,,keep.lib.sizes=FALSE]
dgList <- normLibSizes(dgList)
dgList$samples

design <- model.matrix(~group)
design.reduced <- matrix(1,2,1)
dgList <- estimateGLMCommonDisp(dgList, design.reduced,
                                method="deviance", robust = TRUE, subset=NULL)
fit <- glmFit(dgList, design)
## Normal only, treated vs untreated
lrt <- glmLRT(fit, coef=2)
topTags(lrt)


# Annotating with ENTREZ IDs to do GO and pathway enrichment

# Need ENTREZ ID to be the row name, but cannot have duplicate row names (some genes have unique Ensembl IDs but not gene names/Entrez IDs)
# Finding the Ensembl IDs currently in the database 
idfound<-dgList$genes$gene_id %in% mappedRkeys(org.Hs.egENSEMBL)
dgList<-dgList[idfound,]
dim(dgList)

# Then adding ENTREZ ID
egENSEMBL<-toTable(org.Hs.egENSEMBL)
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

# Running differential expression
group <- factor(c("N","N_Tx"))

dgList.unique <- DGEList(counts=raw.norm,group=metadata$Condition)
keep <- filterByExpr(dgList.unique)
summary(keep)
dgList.unique <- dgList.unique[keep,,keep.lib.sizes=FALSE]
## Stricter filtering
keep2 <- rowSums( cpm(dgList.unique) > 0.1 ) >=2
summary(keep2)
dgList.unique <- dgList.unique[keep2,,keep.lib.sizes=FALSE]
dgList.unique <- normLibSizes(dgList.unique)
dgList.unique$samples

design <- model.matrix(~group)
design.reduced <- matrix(1,2,1)
dgList.unique <- estimateGLMCommonDisp(dgList.unique, design.reduced,
                                       method="deviance", robust = TRUE, subset=NULL)
fit <- glmFit(dgList.unique, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)


# Plotting log-fold change against log counts per million, with DEGs highlighted 
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")


# Gene ontology (GO) and pathway analysis
go <- goana(lrt, species = "Hs")
topGO(go)
keg <- kegga(lrt, species = "Hs")
topKEGG(keg)


## Volcano plots ----
# Results for DEGs -- defining up/down regulated genes
res.sig <- topTags(lrt, n=nrow(lrt$table), adjust.method = "fdr")$table
res.sig$topDE <- "NA"
res.sig$topDE <- as.character(res.sig$topDE)

res.sig$topDE[res.sig$logFC > 1 & res.sig$FDR < 0.05] <- "Up"
res.sig$topDE[res.sig$logFC < -1 & res.sig$FDR < 0.05] <- "Down"
res.sig$topDE[res.sig$topDE!="Up" & res.sig$topDE!="Down"] <- "ns"
res.sig$topDE <- as.factor(res.sig$topDE)

res.sig %>% distinct(topDE) %>% pull()

write.csv(res.sig, "NormalN_NormalTx_sig_filtered.csv")

# Make nice colors
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")
ghibli_colors
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[7], ghibli_colors[4])

cols <- c("Up" = ghibli_colors[4], "Down" = ghibli_colors[3], "ns" = "grey") 
sizes <- c("Up" = 2, "Down" = 2, "ns" = 1) 
alphas <- c("Up" = 1, "Down" = 1, "ns" = 0.5)

ggplot(data=res.sig, aes(x=logFC, y=-log10(FDR), fill=topDE, size=topDE, alpha=topDE)) + 
  geom_point() +
  theme_minimal() +
  geom_hline(yintercept=1.3, linetype="dashed", color="gray") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed",
             color="gray") +
  #xlim(-10,10) +
  geom_point(shape=21,
             color="black") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  labs(title = "Gene expression changes comparing Lenvatinib-treated and untreated normal tissue",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

head(res.sig)

# Highlighting specific genes on volcano

# Cytokine-cytokine receptor interaction pathway genes
cytogenes <- c("GDF10","IL11","IL33","CSF3","CXCL8","TNFRSF12A","TNFSF15","IL24","LIF","CXCR4","CXCL1","PRL","INHBA","CX3CL1","CXCL5","IL1RL1","BMP2","IL6","CCL2","IL7R")

up_cyto_genes <- res.sig %>% filter(gene_name %in% cytogenes  & res.sig$topDE=="Up")
down_cyto_genes <- res.sig %>% filter(gene_name %in% cytogenes  & res.sig$topDE=="Down")

# Up/down cytokine genes
updown_genes <- rbind(up_cyto_genes, down_cyto_genes)

ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_cyto_genes,
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + # For upregulated genes
  geom_point(data = down_cyto_genes,
             shape = 21,
             size = 2, 
             fill = "steelblue", 
             colour = "black") + # For downregulated genes
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = updown_genes,   # Combine both up and down gene lists
                   aes(label = gene_name), 
                   force = 2,
                   nudge_y = 1) + # Add labels last to appear as the top layer
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 5)),
  #                    limits = c(-8, 8)) +
  # scale_y_continuous(limits = c(0, 11)) +
  labs(title = "Gene expression changes comparing Lenvatinib-treated and untreated normal tissue",
       subtitle = "Significantly differentially expressed genes in KEGG 'Cytokine-cytokine receptor interaction' pathway (adjusted pathway enrichment p-value=2.26E-03)",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin treated normal") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

## Bigger
ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_cyto_genes,
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + # For upregulated genes
  geom_point(data = down_cyto_genes,
             shape = 21,
             size = 2, 
             fill = "steelblue", 
             colour = "black") + # For downregulated genes
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = updown_genes,   # Combine both up and down gene lists
                   aes(label = gene_name, size=30), 
                   force = 4,
                   nudge_y = 1,
                   show.legend=F) + # Add labels last to appear as the top layer
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 5)),
  #                    limits = c(-8, 8)) +
  # scale_y_continuous(limits = c(0, 11)) +
  labs(title = "DEGs comparing Lenvatinib-treated and untreated normal tissue",
       subtitle = "KEGG 'Cytokine-cytokine receptor interaction' pathway (adj. p-value=2.26E-03)",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin treated normal") +
  theme_bw(base_size=25) + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
## 1500 x 900

## IL-17 pathway genes
IL17genes <- c("FOSL1","CSF3","IL6","CXCL8","MMP1","CCL2","CXCL1","PTGS2","MMP9","CXCL5","MAPK13")

up_IL17_genes <- res.sig %>% filter(gene_name %in% IL17genes  & res.sig$topDE=="Up")
down_IL17_genes <- res.sig %>% filter(gene_name %in% IL17genes  & res.sig$topDE=="Down")

# Up/down cytokine genes
updown_genes <- rbind(up_IL17_genes, down_IL17_genes)

ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_IL17_genes, # For upregulated genes
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + 
  geom_point(data = down_IL17_genes, # For downregulated genes
             shape = 21,
             size = 2, 
             fill = "steelblue", 
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = updown_genes,   # Combine both up and down gene lists
                   aes(label = gene_name), 
                   force = 2,
                   nudge_y = 1) + # Add labels last to appear as the top layer
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 5)),
  #                    limits = c(-8, 8)) +
  # scale_y_continuous(limits = c(0, 11)) +
  labs(title = "Gene expression changes comparing Lenvatinib-treated and untreated normal tissue",
       subtitle = "Significantly differentially expressed genes in KEGG 'IL-17 signaling pathway' (adjusted pathway enrichment p-value=1.49E-03)",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin treated normal") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

