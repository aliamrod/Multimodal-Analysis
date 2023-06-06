# Multimodal-Analysis of transcriptional regulatory landscape of developing mouse cortex (snRNA + snATAC Pipeline)

![upload_gif](https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/e0ed5460-3481-487c-9670-c58466d627fa)

## Description of Scripts and Analysis

**snRNA-seq Analysis**

**01_snRNA_CellRanger**

**02_snRNA_CellBender**

**03_snRNA_DoubletFinder**

**04_snRNA_QC** - Quality control and standard pre-processing of snRNA data using Seurat v3.1.5. 

**05_snRNA_UMAP** - Louvain clustering + UMAP visualization of snRNA data following transformation.

**06_snRNA_Labelling** - LabelTransfer, SingleR, Manual labelling methodologies.

**07_snRNA_DE - Assignment of cluster identities based on DE genes (based on MAST, Finak et al., 2015)** 

DESeq2 detects count outliers using Cook's distance and removes these genes from analysis. It also automatically removes genes whose mean of normalized counts is below a threshold determined by an optimization procedure. Removing these genes with low counts improves the detection power by making the multiple testing adjustment of the p-values less severe.


**08_snRNA_monocle3** - Pseudotime analysis utilizing Monocle3. Analyses include model fitting to identify gene expression changes as a function of pseudotime.



<p align="center">
<img width="821" alt="Screen Shot 2023-05-27 at 8 29 24 PM" src="https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/7b343b31-995c-46e5-a501-539df1cb6fed">

<img width="878" alt="Screen Shot 2023-06-06 at 11 38 14 AM" src="https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/9554b6a0-064b-4e17-9a71-ac726b5f45b6">

  
  <img width="2029" alt="Screen Shot 2023-06-06 at 11 39 26 AM" src="https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/f2b93562-33a7-4c24-bf8c-d6d7ca94f22c">




![IMG-5504](https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/a296e82a-8041-40bb-a7f7-076f9fcddd25)

<p align="center">
  <img width="568" alt="Screen Shot 2023-03-03 at 2 52 15 PM" src="https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/fbd6e40c-9f8f-4c64-9889-0cdfe8e42881">


------------------------------------------------------------------------------------------------


**snATAC-seq Analysis**

**01_snATAC_QC** - Quality control of snATAC data based on number of fragments and TSS enrichment present.

**02_snATAC_aggregateBin** - Final QC and code for generating count matrices based on genomic bins (for initial clustering), gene bodies and promoters.

**03_snATAC** - 
  Includes: Normalization => Initial clustering based on fragments from fixed-size genome wide bins.

            Peak Normalization => Final peak calling based on initial clusters to generate high-quality peak set; used for final clustering and visualization.

            Identification of chromVARMotifs => Computing motif accessibility deviations using chromVAR (Schep et al., 2017) implemented in Signac.

            Compatation of Gene Scores => Computing gene activity scores utilizing Cicero; used for subsequent integration analyses. 

            Locating unique peaks => Identification of cluster specific peaks.
  
<img width="825" alt="Screen Shot 2023-05-27 at 8 30 00 PM" src="https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/00919d18-f45d-40c9-9a51-78a620a4fe8c">

  
  
<img width="562" alt="Screen Shot 2023-06-06 at 11 40 34 AM" src="https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/763510e3-5ae1-4db8-a8e1-fc3351c6124e">

  <p align="center"><img width="877" alt="Screen Shot 2023-06-06 at 11 42 37 AM" src="https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/f971f8af-d3b7-4d7c-9be0-da62079676b4">

  
  
Last Update: 01 May, 2023.
