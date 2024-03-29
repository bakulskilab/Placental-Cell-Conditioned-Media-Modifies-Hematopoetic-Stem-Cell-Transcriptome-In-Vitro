# Placental cell conditioned media modifies hematopoietic stem cell transcriptome in vitro

## Citation Information

Harris SM, Su AL, Dou JF, Loch-Caruso R, Elkin ER, Bakulski KM. 2022. Placental cell conditioned media modifies hematopoietic stem cell transcriptome in vitro. In review

This Github repository contains the data management and analytic scripts to produce the following manuscript: Placental cell conditioned media modifies hematopoietic stem cell transcriptome in vitro

## Abstract

Background: Factors secreted from the placenta into maternal and fetal circulation are important for maternal and fetal development during pregnancy. Maternal hematopoietic stem cells maintain maternal blood supply. The placenta is an early site of fetal hematopoiesis, and placental hematopoietic stem cells are involved in the early stages of fetal blood cell differentiation. Cross-talk between placental cells and hematopoietic stem cells represents an important but understudied phenomenon of pregnancy. The impact of toxicant exposure on placental-immune cell communication is poorly understood. The goals of this study were to 1) determine if factors secreted from placental cells alter transcriptomic responses in hematopoietic stem cells in vitro and 2) if monoethylhexyl phthalate (MEHP), the major metabolite of the pollutant diethylhexyl phthalate, modifies these effects.

Methods: Using in vitro cell line models of hematopoietic stem cells (K-562) and placental syncytiotrophoblasts (differentiated BeWo), we treated K-562 cells for 24 hours with media conditioned by incubation with BeWo cells, media conditioned with BeWo cells treated with 10 µM MEHP, or unconditioned controls (n = 4 replicates for each group). We extracted K-562 cell RNA and performed RNA sequencing. We then conducted differential gene expression and pathway analysis by treatment group and accounted for multiple comparisons using false discovery rate (FDR). 

Results: K-562 cells treated with BeWo cell conditioned media differentially expressed 173 genes (FDR<0.05 and fold-change>2.0), relative to vehicle control. Cells treated with BeWo media upregulated TPM4 2.4 fold (FDR=1.8x10-53) and S1PR3 3.3 fold (FDR=1.6x10-40) compared to controls. Upregulated genes were enriched for biological processes including stem cell maintenance (“somatic stem cell population maintenance”), cell proliferation (“positive regulation of endothelial cell proliferation”) and immune processes (“myeloid leukocyte cytokine production”). Downregulated genes were enriched for protein translation (“mitochondrial translation”) and transcriptional regulation (“RNA processing”). Cells treated with media from BeWo cells treated with MEHP upregulated (FDR<0.05) eight genes, including genes involved in fat/lipid metabolism (PLIN2, fold-change: 1.4; CPT1A, fold-change: 1.4) and iron uptake (TFRC, fold-change: 1.3), relative to the BeWo conditioned media treatment.

Conclusion: Hematopoietic stem cells are responsive to media that has been conditioned by placental cells, potentially impacting processes related to stem cell maintenance and proliferation, which may represent placental-immune communication important for development. The metabolite MEHP only had a modest impact on these responses at the single concentration tested. Future directions will investigate components of placental cell media (hormones, microvesicles, and proteins) contributing to hematopoietic cell signals. 


## Script Files

1_fastQC.pbs: using fastQC and multiQC for raw RNA sequencing data

2_STAR.pbs: aligning raw data to human genome

3_QoRTs.pbs: quality control of the aligned data

4_featureCounts.pbs: generating gene count data

5_file_prep.R: data preparation; checking gene IDs and mappling them to gene symbol

6_diff_exp_analysis.Rmd: differential gene expression analysis
