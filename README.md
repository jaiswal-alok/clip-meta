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


# clip-meta
