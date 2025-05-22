# Load packages 

# Bioconductor packages
library(edgeR) # For differential expression calculation
library(GO.db) # For Gene Ontology analysis 
library(GseaVis) # For nicer GSEA plot

# Other packages 
library(dplyr) # For general data manipulation
library(ghibli) # For pretty colors
library(ggplot2) # For plots (e.g. Volcano plot)
pacman::p_load(here,  
               tidyverse, 
               janitor, # Cleaning column names  
               scales, # Transform axis scales   
               ggrepel) # Optimize plot label separation 


# This is the data analysis script for differential expression for the pilot experiment from Konrad's mouse organoids

# Read in the normalized data and metadata ----
data <- readRDS(file = "input/normalized_data.rds")
metadata <- read.delim("input/organoid_pilot_metadata.txt",check.names=FALSE, stringsAsFactors=FALSE)

# Running differential expression ----

## Model design ----
# Because this experiment has no replicates since it's a pilot analysis, we are using the Likelihood Ratio Test as suggested by the vignette
# Using option 3 from section 2.12 of the vignette (What to do if you have no replicates)

# Define our grouping variable (WT vs Mut)
group <- factor(c("WT","Mut"))

design <- model.matrix(~group)
design.reduced <- matrix(1,2,1)
data <- estimateGLMCommonDisp(data, design.reduced,
                                       method="deviance", robust = TRUE, subset=NULL)

# Fit our differential expression model ----
fit <- glmFit(data, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)


# Plotting log-fold change against log counts per million, with DEGs highlighted 
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")


# Gene ontology (GO) and pathway analysis
go <- goana(lrt, species = "Mm")
topGO(go)
keg <- kegga(lrt, species = "Mm")
topKEGG(keg)


## Volcano plots ----
### Results for DEGs -- defining up/down regulated genes ----
res.sig <- topTags(lrt, n=nrow(lrt$table), adjust.method = "fdr")$table
res.sig$topDE <- "NA"
res.sig$topDE <- as.character(res.sig$topDE)

res.sig$topDE[res.sig$logFC > 1 & res.sig$FDR < 0.05] <- "Up"
res.sig$topDE[res.sig$logFC < -1 & res.sig$FDR < 0.05] <- "Down"
res.sig$topDE[res.sig$topDE!="Up" & res.sig$topDE!="Down"] <- "ns"
res.sig$topDE <- as.factor(res.sig$topDE)

res.sig %>% distinct(topDE) %>% pull()

write.csv(res.sig, "output/Sig_DEGs_NoCre_vs_PlusCre_P53mut.csv")

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
  labs(title = "Gene expression changes comparing \nP53 Wild-type (pre-Cre) and P53 mutant (post-Cre) uterine horn organoids",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

head(res.sig)

# Good tutorial for customizing volcano: https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/ 

### Highlighting specific genes on volcano ----
#### Top Genes ----
# Top 15 differentially expressed genes (taken manually from results .csv file)
topgenes <- c("Cldn2","Trf","Gpx3","Mmp7","Alb","Havcr1","Ccn2","Elapor1","Cpe","Dmbt1","Kctd12",
              "Obscn","Armcx3","Kcnj16","Inpp5d")

up_top_genes <- res.sig %>% filter(gene_name %in% topgenes  & res.sig$topDE=="Up")
down_top_genes <- res.sig %>% filter(gene_name %in% topgenes  & res.sig$topDE=="Down")

# Up/down Rap1 Signaling genes
updown_top_genes <- rbind(up_top_genes, down_top_genes)

ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_top_genes,
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + # For upregulated genes
  geom_point(data = down_top_genes,
             shape = 21,
             size = 2, 
             fill = "steelblue", 
             colour = "black") + # For downregulated genes
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = updown_top_genes,   # Combine both up and down gene lists
                   aes(label = gene_name), 
                   force = 2,
                   nudge_y = 1) + # Add labels last to appear as the top layer
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 5)),
  #                    limits = c(-8, 8)) +
  # scale_y_continuous(limits = c(0, 11)) +
  labs(title = "Gene expression changes comparing \nP53 Wild-type (pre-Cre) and P53 mutant (post-Cre) uterine horn organoids",
       subtitle = "Top 15 Differentially Expressed Genes",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin WT Organoid") +
  theme_bw(base_size = 12) + # Select theme with a white background, scale up font size slightly   
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


#### Pathways ----
# I put the top DEGs (adjusted p-value <= 0.05) into Enrichr 
# I am showing the genes in a few of the top pathways

# Rap1 Signaling Pathway
rap1path <- c("Vav3","Itgam","Itgb2","Pdgfb","Vegfb","Pard6g","Itgal","Thbs1","Gnai1","Adcy9","Plcb4",
             "Akt3","Id1","Kit","Tln2","Rac2","Fgfr3","Fgfr2","Rapgef4")

up_rap1_genes <- res.sig %>% filter(gene_name %in% rap1path  & res.sig$topDE=="Up")
down_rap1_genes <- res.sig %>% filter(gene_name %in% rap1path  & res.sig$topDE=="Down")

# Up/down Rap1 Signaling genes
updown_rap1_genes <- rbind(up_rap1_genes, down_rap1_genes)

ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_rap1_genes,
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + # For upregulated genes
  geom_point(data = down_rap1_genes,
             shape = 21,
             size = 2, 
             fill = "steelblue", 
             colour = "black") + # For downregulated genes
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = updown_rap1_genes,   # Combine both up and down gene lists
                   aes(label = gene_name), 
                   force = 2,
                   nudge_y = 1) + # Add labels last to appear as the top layer
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 5)),
  #                    limits = c(-8, 8)) +
  # scale_y_continuous(limits = c(0, 11)) +
  labs(title = "Gene expression changes comparing \nP53 Wild-type (pre-Cre) and P53 mutant (post-Cre) uterine horn organoids",
       subtitle = "Significantly differentially expressed genes in KEGG 'Rap1 Signaling' pathway \n(adjusted pathway enrichment p-value=0.01)",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin WT Organoid") +
  theme_bw(base_size = 12) + # Select theme with a white background, scale up font size slightly   
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


# Chemokine Signaling Pathway
chemokines <- c("Vav3","Ncf1","Ppbp","Cxcl2","Cxcl15","Gnai1","Prex1","Cxcl10","Pak1","Adcy9","Plcb4","Ccl6","Ccl5","Akt3","Cxcr2","Rac2","Ccl2")

up_chemo_genes <- res.sig %>% filter(gene_name %in% chemokines  & res.sig$topDE=="Up")
down_chemo_genes <- res.sig %>% filter(gene_name %in% chemokines  & res.sig$topDE=="Down")

# Up/down Chemokine Signaling genes
updown_chemo_genes <- rbind(up_chemo_genes, down_chemo_genes)

ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_chemo_genes,
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + # For upregulated genes
  geom_point(data = down_chemo_genes,
             shape = 21,
             size = 2, 
             fill = "steelblue", 
             colour = "black") + # For downregulated genes
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = updown_chemo_genes,   # Combine both up and down gene lists
                   aes(label = gene_name), 
                   force = 2,
                   nudge_y = 1) + # Add labels last to appear as the top layer
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 5)),
  #                    limits = c(-8, 8)) +
  # scale_y_continuous(limits = c(0, 11)) +
  labs(title = "Gene expression changes comparing \nP53 Wild-type (pre-Cre) and P53 mutant (post-Cre) uterine horn organoids",
       subtitle = "Significantly differentially expressed genes in WikiPathways 'Chemokine Signaling' pathway \n(adjusted pathway enrichment p-value=2.88E-03)",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin WT Organoid") +
  theme_bw(base_size = 12) + # Select theme with a white background, scale up font size slightly   
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


# Regulation of IGF Transport and Uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)
igfbp <- c("Rcn1","Sdc2","Alb","Tnc","Fn1","Igf2","Apoe","Prss23","Kng1")

up_igfbp_genes <- res.sig %>% filter(gene_name %in% igfbp  & res.sig$topDE=="Up")
down_igfbp_genes <- res.sig %>% filter(gene_name %in% igfbp  & res.sig$topDE=="Down")

# Up/down Chemokine Signaling genes
updown_igfbp_genes <- rbind(up_chemo_genes, down_chemo_genes)

ggplot(data = res.sig,
       aes(x = logFC,
           y = -log10(FDR))) + 
  geom_point(aes(colour = topDE), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = up_chemo_genes,
             shape = 21,
             size = 2, 
             fill = "firebrick", 
             colour = "black") + # For upregulated genes
  geom_point(data = down_chemo_genes,
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
  labs(title = "Gene expression changes comparing \nP53 Wild-type (pre-Cre) and P53 mutant (post-Cre) uterine horn organoids",
       subtitle = "Significantly differentially expressed genes in Reactome 'Regulation of IGF Transport and Uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)' pathway \n(adjusted pathway enrichment p-value=2.18E-04)",
       x = "log(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression change\nin WT Organoid") +
  theme_bw(base_size = 12) + # Select theme with a white background, scale up font size slightly   
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


