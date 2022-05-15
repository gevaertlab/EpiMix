# EpiMix
## An intergrative tool for epigenomic subtyping using DNA methylation

<hr>
### Introduction
EpiMix is a comprehensive tool for the integrative analysis of DNA methylation data and gene expression data (**Figure 1**). EpiMix enables automated data downloading (from TCGA or GEO), pre-processing, methylation modeling, interactive visualization and functional annotation. To identify hypo- or hypermethylated genes, EpiMix uses a beta mixture model to identify the methylation states of each CpG site and compares the DNA methylation in an experimental group to a control group. The output from EpiMix is the functional changes in DNA methylation that is associated with gene expression. 

EpiMix incorporates specialized algorithms to identify functional DNA methylation at various genetic elements, including proximal cis-regulatory elements of protein-coding genes, distal enhancers, and genes encoding microRNAs (miRNAs) or lncRNAs. There are four alternative analytic modes for modeling DNA methylation at different genetic elements:

* **Regular**: cis-regulatory elements within or immediately surrounding protein-coding genes.
* **Enhancer**: distal enhancers.
* **miRNA**: miRNA-coding genes.
* **lncRNA**: lncRNA-coding genes.

![**Figure 1**. Overview of EpiMix's workflow. EpiMix incorporates four functional modules: downloading, preprocessing, methylation modeling and functional analysis. The methylation modeling module enables four alternative analytic modes, including "Regular", "Enhancer", "miRNA" and "lncRNA". These analytic modes target DNA methylation analysis on different genetic elements.](vignettes/figures/Workflow.png)
<hr>

### Interactive web tool

EpiMix is also available as an interactive web tool: [https://epimix.stanford.edu](https://epimix.stanford.edu)

<hr>

### Authors
<hr>

Yuanning Zheng, John Jun, Kevin Brennan and Olivier Gevaert<br>
Stanford Center for Biomedical Informatics Research (BIMR)<br>
Department of Medicine<br>
1265 Welch Road<br>
Stanford CA, 94305-5479

<hr>
