# Pilot Project Results
Results from initial pilot project, with RNA-seq alignment and quantification performed on `20 May 2025`.

## Experiment Details
- **Samples**
  - **Species:** Mouse (Mus musculus)
  - **Malignancy:** Uterine horn (endometrial)
  - `Sample 1` Wild-type (pre-Cre) for P53 mutant
  - `Sample 2` Mutant (plus-Cre) P53
- **Timeframe**
  -   RNA isolated and sent for sequencing on `10 April 2025`
  -   Sequencing returned `29 April 2025` (KM received download link `01 May 2025`)
  -   Nextflow [NF-Core RNA-seq](https://github.com/nf-core/rnaseq) pipeline run `20 May 2025`
  -   *Next up*: KM to complete EdgeR analysis in R

## Contents
- **Pipeline QC:** Download [MultiQC](multiqc_report.html) html file to view the aggregated quality control results for the pipeline
- **Quantified Reads:**
  - Use the results from the [star_salmon](results/star_salmon) folder: *[salmon.merged.gene_counts.tsv](results/star_salmon/salmon.merged.gene_counts.tsv)*
  - This is the input for differential expression analysis
- **Differential expression analysis:** In the [edgeR](EdgeR) folder.
- *Everything else is semi-optional, worry about these first.*
