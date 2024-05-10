---
title: Analysis of epigenetic signals in cell-free DNA
author: Michael Noe
date: 2024-04-26
---

Analyses of epigenetic signals in cell-free DNA
===============================================


The following sections provide details for the analyses performed in:

*DNA methylation and gene expression as potential determinants of genome-wide cell-free DNA fragmentation*

# Citation
Michaël Noë, Dimitrios Mathios, Akshaya V. Annapragada, Shashikant Koul, Zacharia H. Foda, Jamie Medina, Stephen Cristiano, Christopher Cherry, Daniel C. Bruhm, Noushin Niknafs, Vilmos Adleff, Leonardo Ferreira, Hari Easwaran, Stephen Baylin, Jillian Phallen, Robert B. Scharpf, and Victor E. Velculescu.
DNA methylation and gene expression as potential determinants of genome-wide cell-free DNA fragmentation


# Abstract
Circulating cell-free DNA (cfDNA) is emerging as a diagnostic avenue for cancer detection, but the characteristics and origins of cfDNA fragmentation in the blood are poorly understood.  We evaluated the effect of DNA methylation and gene expression on naturally occurring genome-wide cfDNA fragmentation through analysis of plasma from 969 individuals, including 182 with cancer.  cfDNA fragment ends occurring at preferred locations genome-wide more frequently contained CCs or CGs, and fragments ending with CGs or CCGs were enriched or depleted, respectively, at methylated CpG positions, consistent with structural models showing increased interaction of methylated CG fragment ends with nucleosomes.  Higher levels and larger sizes of cfDNA fragments were independently associated with regions of CpG methylation and reduced gene expression, and reflected differences in cfDNA fragmentation in tissue-specific pathways.  The effects of methylation and expression on cfDNA coverage were validated by analyses of human cfDNA in mice implanted with isogenic tumors with or without the mutant IDH1 chromatin modifier.  Tumor-related hypomethylation and increased gene expression were associated with global decrease in cfDNA fragment size that may explain the overall smaller cfDNA fragments observed in human cancers.  Cancer-specific methylation at CpGs of pancreatic cancer patients was associated with genome-wide changes in cfDNA fragment ends in patients with cancers.  These results provide a connection between epigenetic changes and cfDNA fragmentation that may have implications for disease detection.

Table of contents
-----------------

* [Citation](#citation)
* [Abstract](#abstract)
* [Introduction](#introduction)
* [Available data used in analysis](#available-data-used-in-analysis)
* [Pre-processing](#pre-processing)
* [Figure 1](#figure-1)
* [Figure 2](#figure-2)
* [Figure 3](#figure-3)
* [Figure 4](#figure-4)
* [Supplementary Figure 1](#supplementary-figure-1)
* [Supplementary Figure 2](#supplementary-figure-2)
* [Supplementary Figure 3](#supplementary-figure-3)
* [Supplementary Figure 4](#supplementary-figure-4)
* [Supplementary Figure 5](#supplementary-figure-5)
* [Supplementary Figure 6](#supplementary-figure-6)
* [Supplementary Figure 7](#supplementary-figure-7)
* [Supplementary Figure 8](#supplementary-figure-8)
* [Supplementary Figure 9](#supplementary-figure-9)
* [Supplementary Figure 10](#supplementary-figure-10)
* [Supplementary Figure 12](#supplementary-figure-12)
* [Supplementary Figure 13](#supplementary-figure-13)
* [Supplementary Figure 14](#supplementary-figure-14)
* [Supplementary Figure 15](#supplementary-figure-15)

Introduction
------------

Some of the analyses described in this repository relay on data from previous studies (Christiano et al., Nature, 2019 and Mathios et al., Nature Communications, 2021) that have been deposited at the database of Genotypes and Phenotypes (dbGaP) and the European Genome-Phenome Archive (EGA). The GitHub-repositories for these supporting studies provide additional details regarding how the data was processed ([pre-processing](https://github.com/cancer-genomics/reproduce_lucas_wflow/tree/master/code/preprocessing) in the '[reproduce_lucas_wflow](https://github.com/cancer-genomics/reproduce_lucas_wflow)'-GitHub). These objects were saved as '.rds' file and referred to as analytic data.  Some of the analytic data is provided as `GRanges` objects.  Unlike BED-files (not used here), GRanges objects contain the start-position of the fragment, while BED-files use the position before the start-position as the first base in the fragment.

For the analyses desribed in this README, we start from analytic data saved in '.rds'-files.  To construct the Figures in the paper, we often generate temporary files that are too big or numerous to upload to the GitHub repository. These files are stored in a top-level directory called 'cfepigenetics_data'. The temporary files stored there are used as input for scripts that generate a final summary-file that are small enough to store in the [data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)-folder of this repository.

Available data used in analysis
-------------------------------

Datasets used in this study that have been made publicly available from previous publications:

* Cristiano et al., Nature, 2019: the database of Genotypesand Phenotypes (dbGaP, study ID 34536).
* Mathios et al., Nature Communications, 2021: the European Genome-Phenome Archive (EGA, EGAS00001005340).

Methylation datasets of cell-free DNA:

* Moss et al., Nature Communications, 2018: the NCBI Gene Expression Omnibus (GEO, GSE1221261) database repository.
* Loyfer et al., Nature, 2023: the NCBI Gene Expression Omnibus (GEO, GSE186458) database repository.

Gene expression datasets of white-blood cells (myeloid):

* Uhlen et al., Science, 2015: data is available at [The Human Protein Atlas](https://v13.proteinatlas.org/download/rna.csv.zip)
* [Genotype-Tissue Expression Project (GTEx)](https://gtexportal.org/home)

Gene expression and DNA methylation datasets of other tumors and healthy white blood cells ([Supplementary Figure 15](#supplementary-figure-15)) were downloaded from [The Cancer Genome Atlas (TCGA)](https://portal.gdc.cancer.gov/).


Pre-processing
--------------

Processing of raw sequencing data (fastq files) to generate GRanges-objects containing information about cell-free DNA fragment alignments (as defined by 'chromosome', 'start' and 'end'), we refer to the [GitHub-repository](https://github.com/cancer-genomics/reproduce_lucas_wflow) from the paper: Mathios et al., Nature Communications, 2021. In the [code-folder](https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/master/code), there is a [pre-processing-folder](https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/master/code/pre-processing), containing the scripts for pre-processing.

* [fastp.sh](https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/master/code/preprocessing/fastp.sh): We used '[fastp](https://github.com/OpenGene/fastp)' specifically to trim adapters, when the cell-free DNA fragments were shorter than the amount of read-out cycles, which would have led to reading into the adapter on the other side of the DNA-strand.
* [align.sh](https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/master/code/preprocessing/align.sh): We used '[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)' to align the reads of the fastq-files to HG19.
* [post_alignment.sh](https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/master/code/preprocessing/post_alignment.sh): We used '[Sambamba](https://lomereiter.github.io/sambamba/)' to flag duplicates and to extract cell-free DNA fragment information, like positions (as defined by 'chromosome', 'start' and 'end') and MAPQ-values, which are written to a BED-file.
* [bed_to_granges.sh](https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/master/code/preprocessing/bed_to_granges.sh): We used a custom R-script ([01-bed_to_granges.r](https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/master/code/preprocessing/01-bed_to_granges.r)), to transform the BED-file into a GRanges-object, saved as a '.rds'-file. During this process we also filtered cell-free DNA fragments with MAPQ-value <= 30 to exclude low-quality alignments.

Figure 1
--------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Fig_1.jpg" width = "400">

* [Pre_Figure1.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure1.rmd): contains a step-by-step guide with scripts to process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and to summarize these files (all samples).  The summary-file has been uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Figure1.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Figure1.rmd): processes the provided summarized file (analytic data) to generate Figure 1.
* Figure 1D contains a 3D rendering of a nucleosome, as captured from the Protein Data Bank (structure [7COW](https://www.rcsb.org/structure/7cow)).


Figure 2
--------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Fig_2.jpg" width = "400">

* [Pre_Figure2.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure2.rmd): contains a step-by-step to process the raw data (GRanges-objects; per sample) to intermediary files (per sample) as well as a summarization step.  The summary-file is provided with this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Figure2.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Figure2.rmd): processes the provided summarized file (analytic data) to generate Figure 2.

Figure 3
--------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Fig_3AB.jpg" width = "400">
<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Fig_3CD.jpg" width = "400">
<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Fig_3E.png" width = "400">
<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Fig_3F1.png" width = "400">
<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Fig_3F2.png" width = "400">
<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Fig_3G.jpg" width = "400">
<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Fig_3H.jpg" width = "400">

* [Pre_Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure3.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Figure3.rmd): processes the provided summarized file (analytic data) to generate Figure 3 (except Figure 3E and 3F).
* [Figure3EF.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Figure3EF.rmd): processes the provided summarized file (analytic data) to generate Figure 3E and 3F.

Figure 4
--------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Fig_4.jpg" width = "400">

* [Pre_Figure4.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure4.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Pre_Figure4_Model_Ensemble.Rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure4_Model_Ensemble.Rmd): code that uses the summarized file (all samples) to generate the prediction model and generate a summarized file of this model. The code uses a 10-fold cross validation method, which needs to be manually changed for every fold (indicated in the code).
* [Figure4.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Figure4.rmd): processes the provided summarized file (analytic data) to generate Figure 4.

Supplementary Figure 1
----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_1.jpg" width = "400">

* [Pre_Supplementary_Figure1.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure1.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository. This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)). ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)).
* [Supplementary_Figure1.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure1.rmd): processes the provided summarized file file (analytic data) to generate Supplementary Figure 1.

Supplementary Figure 2
----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_2.jpg" width = "400">

This figure is a variation on Figure 2A, using different beta-value cut-offs to define 'methylated' and 'unmethylated'.
* [Pre_Figure2.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure2.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Supplementary_Figure2.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure2.rmd): processes the provided summarized file file (analytic data) to generate Supplementary Figure 2.

Supplementary Figure 3
----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_3.jpg" width = "400">

* [Pre_Supplementary_Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Supplementary_Figure3.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Supplementary_Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure3.rmd): processes the provided summarized file file (analytic data) to generate Supplementary Figure 3.

Supplementary Figure 4
----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_4.jpg" width = "400">

* [Pre_Supplementary_Figure4.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Supplementary_Figure4.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Supplementary_Figure4.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure4.rmd): processes the provided summarized file file (analytic data) to generate Supplementary Figure 4.

Supplementary Figure 5
----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_5.jpg" width = "400">

* [Pre_Supplementary_Figure5.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Supplementary_Figure5.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Supplementary_Figure5.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure5.rmd): processes the provided summarized file file (analytic data) to generate Supplementary Figure 5.

Supplementary Figure 6
----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_6.jpg" width = "400">

* [Pre_Supplementary_Figure4.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Supplementary_Figure4.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Supplementary_Figure6.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure6.rmd): processes the provided summarized file file (analytic data) to generate Supplementary Figure 6.

Supplementary Figure 7
----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_7.jpg" width = "400">

* [Pre_Figure2.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure2.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Supplementary_Figure7.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure7.rmd): processes the provided summarized file file (analytic data) to generate Supplementary Figure 7.

Supplementary Figure 8
----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_8.jpg" width = "400">

* [Pre_Figure2.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure2.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Supplementary_Figure8.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure8.rmd): processes the provided summarized file file (analytic data) to generate Supplementary Figure 8.

Supplementary Figure 9
----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_9A.jpg" width = "400">
<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_9B.jpg" width = "400">
<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_9C.jpg" width = "400">

* [Pre_Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure3.rmd): generates a matrix, connecting CpG-islands (and beta-values) to the nearest transcription start site (TSS) (and gene-expression values), which is uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Figure3.rmd): populates the previously generated matrix with information about coverage, fragment-size and nucleosome positioning ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)).
* [Supplementary_Figure9.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure9.rmd): processes the provided summarized file file (analytic data) to generate Supplementary Figure 9.


Supplementary Figure 10
-----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_10.jpg" width = "400">

* [Pre_Supplementary_Figure10.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Supplementary_Figure10.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Supplementary_Figure10.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure10.rmd): process provided summarized file (analytic data) to generate Supplementary Figure 10.

Supplementary Figure 11
-----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_11.jpg" width = "400">

* [Supplementary_Figure11.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure11.rmd): process provided summarized file (analytic data) to generate Supplementary Figure 11.

Supplementary Figure 12
----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_12.jpg" width = "400">

* [Pre_Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure3.rmd): generates a matrix, connecting CpG-islands (and beta-values) to the nearest transcription start site (TSS) (and gene-expression values), which is uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)).
* [Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Figure3.rmd): populates the previously generated matrix with information about coverage, fragment-size and nucleosome positioning ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)).
* [Supplementary_Figure12.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure12.rmd): process provided summarized file (analytic data) to generate Supplementary Figure 12.

Supplementary Figure 13
-----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_13.jpg" width = "400">

* [Pre_Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure3.rmd): contains a step-by-step guide which scripts will process the raw data (GRanges-objects; per sample) to intermediary files (per sample) and summarize them (all samples) into a summary-file, uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Supplementary_Figure13.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure13.rmd): process provided summarized file (analytic data) to generate Supplementary Figure 13.

Supplementary Figure 14
-----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_14A.jpg" width = "400">
<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_14B.jpg" width = "400">

* [Pre_Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Pre_Figure3.rmd): generates a matrix, connecting CpG-islands (and beta-values) to the nearest transcription start site (TSS) (and gene-expression values), which is uploaded to this repository ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)). This script requires the raw data (after [pre-processing](#pre-processing) the data from [Cristiano et al. and Mathios et al.](#required-data)).
* [Figure3.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Figure3.rmd): populates the previously generated matrix with information about coverage, fragment-size and nucleosome positioning ([data](https://github.com/cancer-genomics/cfepigenetics/blob/main/data)).
* [Supplementary_Figure14.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure14.rmd): process provided summarized file (analytic data) to generate Supplementary Figure 14.

Supplementary Figure 15
-----------------------

<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_15A.jpg" width = "400">
<img src="https://github.com/cancer-genomics/cfepigenetics/blob/main/output/Supplementary_Fig_15B.jpg" width = "400">

 [Supplementary_Figure15.rmd](https://github.com/cancer-genomics/cfepigenetics/blob/main/analysis/Supplementary_Figure15.rmd): process provided summarized file (analytic data) to generate Supplementary Figure 15.


Session information
-------------------


```r
pander::pander(sessionInfo())
```

**R version 4.3.2 (2023-10-31)**

**Platform:** aarch64-apple-darwin20 (64-bit)

**locale:**
en_US.UTF-8||en_US.UTF-8||en_US.UTF-8||C||en_US.UTF-8||en_US.UTF-8

**attached base packages:**
_stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**loaded via a namespace (and not attached):**
_digest(v.0.6.33)_, _R6(v.2.5.1)_, _fastmap(v.1.1.1)_, _xfun(v.0.41)_, _cachem(v.1.0.8)_, _knitr(v.1.45)_, _htmltools(v.0.5.7)_, _rmarkdown(v.2.25)_, _lifecycle(v.1.0.4)_, _cli(v.3.6.1)_, _pander(v.0.6.5)_, _sass(v.0.4.7)_, _jquerylib(v.0.1.4)_, _compiler(v.4.3.2)_, _rstudioapi(v.0.15.0)_, _tools(v.4.3.2)_, _evaluate(v.0.23)_, _bslib(v.0.6.0)_, _Rcpp(v.1.0.11)_, _yaml(v.2.3.7)_, _rlang(v.1.1.2)_ and _jsonlite(v.1.8.7)_
