# Cell Line specific gene Identification Pipeline (CLIP)

## CLIP Rationale
CLIP enables systematic meta-analysis and integration of multi-modal omics datasets collected for multiple research laboratories. CLIP is a “bottom-up” approach  to meta-analysis and data integration. To boost the statistical power toward finding robust and reproducible signals, the CLIP framework accounts for the substantial variability in the consistency of the various types of modalities between laboratories. 

First, we define the notion of a “Cancer Cell line Specific” (CCS) gene; as a gene that exhibits a molecular attribute unique to a particular cancer cell line. In other words, a CCS gene is unique to a given cell line in reference to all the other cell lines in the particular dataset, i.e. having a context-specific property, and therefore possibly related to a specific cancer subtype or cellular function. Statistically, a gene that has the tendency to be located towards the extremes of the distribution in any given dataset is considered as a CCS gene. For instance, the expression of ERBB2 gene is much higher in HER2 driven breast cancer cell lines, compared to cell lines from other tissue types. Thus, in all HER2+ cell lines, ERBB2 is a CCS gene in the gene expression modality. 

<p align="center"> 
<img src="images/figure1.png" style="width: 5%; height: 5%"/>​
</p>

## CLIP Flowchart

<p align="center"> 
<img src="images/figure2.png" style="width: 10%; height: 20%"/>​
</p>

## CLIP pipeline execution

### Step 1

```r
# Code for CLIP exection
source("clip-meta/R/clip/clip.R")

###########################################################################
### Step 1a: Quantile Normalization and Scaling of continuous datasets  ###
### Outlier Evidence Score (OES) scores                                 ###
###########################################################################
# The first step involves loading all processed data sets and Quantile normalization and Outlier Evidence Score (OES) estimates

source("clip-meta/R/clip/clip_functions.R")
# source functions utilized in CLIP pipeline 

#load processed datasets for each modality 
load("clip-meta/data/CNV.RData")
load("clip-meta/data/GEXP.RData")
load("clip-meta/data/FUNC.RData")
load("clip-meta/data/METH.RData")
load("clip-meta/data/MUT.RData")
load("clip-meta/data/PEXP.RData")
load("clip-meta/data/PHOS.RData")
load("clip-meta/dataTAS.RData")
..

######################################################################
### Step 1b: Proportion scores in binary datasets                  ###
######################################################################
# For binary datasets estimate proportion scores
..
# Save all processed files 
save(file="clip-meta/processed_files/clip/Normalized_Modalities.RData",compress="bzip2",
     FUNC.QUANT.LIST, METH.QUANT.LIST, GEXP.QUANT.LIST, PEXP.QUANT.LIST, 
     PHOS.QUANT.LIST, TAS.QUANT.LIST, CNV.DEL.LIST, CNV.AMP.LIST, MUT.LIST)

``` 

### Step 2
```r
######################################################################
### Step 2: CLIP - meta analysis of breast cancer cell lines       ###
######################################################################
# Meta-analysis and integration of OES scores from multiple studies for each breast cancer cell lines

#load normalized and scaled data in Step 1
load("clip-meta/processed_files/clip/Normalized_Modalities.RData")

..

# Save all processed files for next step
save(file = "clip-meta/processed_files/clip//CLIP_Integrated_Breast.RData",
     brca.mut.meta, brca.cnv_del.meta, brca.cnv_amp.meta, brca.meth.meta, brca.gexp.meta,
     brca.pexp.meta, brca.phos.meta, brca.func.meta, brca.tas.meta, compress="bzip2")


``` 

### Step 3
```r
#########################################################################
### Step 3: CLIP - Rank Product Integration                          ####
### Identification of rCCS genes                                     ####
#########################################################################
# Integration of Cell line Specific genes (CCS) genes from multiple modalites to identify robust CCS (rCCS) for each breast cancer cell lines

#load data from Step 2
load("clip-meta/processed_files/clip//CLIP_Integrated_Breast.RData")

..
``` 

### Benchmarking CLIP 

```r
#Source the following R file to run CLIP benchmarking
source("clip-meta/R/clip/clip_benchmark.R")
..
``` 

### Reproduciblity analyses

```r
#Source the following R files to run correlation analyses between datasets
source("clip-meta/R/reproducibility_analysis/reproducibility_analysis_functions.R")
source("clip-meta/R/reproducibility_analysis/reproducibility_analysis.R")
..
``` 

## Citation

    @article{Jaiswal2020,
        author = {Jaiswal, A. and Gautam, P. and Pietilä, A. and Timmonen, S. and Zenz, T. and Nordström, N. and Akimov, Y. and Sipari, N. and Tanoli, Z. and Lehti, K. and Wennerberg, K. and Aittokallio, T.},
        title = {Multi-modal meta-analysis of cancer cell line omics profiles identifies ECHDC1 as a novel breast tumor suppressor},
        year = {2020},
        doi = {10.1101/2020.01.31.929372},
        publisher = {Cold Spring Harbor Laboratory},
        URL = {https://www.biorxiv.org/content/10.1101/2020.01.31.929372v1.abstract},
        journal = {bioRxiv}
    }



