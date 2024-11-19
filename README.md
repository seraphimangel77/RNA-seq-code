# Scripts and files used in "Cancer cells hijack orphan ESCRT subunit to fuel oncogenic growth"
## Data source

1. RNA-seq data of ExtendedDataFig.7f, have been deposited in the Genome Sequence Archive under the accession number PRJCA020571.

2. The Pan Cancer RNA-seq and gene copy number variations of HRS, STAM1, and STAM2 were obtained from the University of California Santa Cruz (UCSC) Xena browser (http://xenabrowser.net/).

3. The MS proteomics data have been deposited to the ProteomeXchange Consortium (http://proteomecentral.proteomexchange.org) via the iProX partner repository under the dataset identifier PXD045939.

## Multivariate Cox analysis
The Fig.3f_code script uses the forestploter package to perform multivariate Cox analysis on 54 patient-derived HCC samples (group_cox.tsv).
## Quantitative MS analysis
The Fig.2e_Fig. 5h_Fig. 7a_Fig. 8a_code script uses limma and ggVolcano to perform Quantitative MS analysis on JHH-7 cells.
## Correlation analysis
The Fig.3b_code script performs correlation analysis on HCC tumor data (Input Fig 3b.xlsx).
## **HRS**, **STAM1**, and **STAM2** expression levels
The ExtendedDataFig.1a_code script analyzes the expression levels of **HRS**, **STAM1**, and **STAM2** in the TCGA Pan-Cancer dataset (HGS_STAM_STAM2_expression.tsv).
## The frequency of gene copy number variation
The ExtendedDataFig.1b_code script analyzes the gene copy number variation of HRS (HGS_CNV.tsv), STAM1 (STAM_CNV.tsv), and STAM2 (STAM2_CNV.tsv) in the TCGA Pan-Cancer dataset.
## Heatmap and GSEA
The ExtendedDataFig.7f_code script generates a heatmap and performs GSEA analysis on the RNA-seq data (PRJCA020571).
