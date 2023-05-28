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





![IMG-5504](https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/a296e82a-8041-40bb-a7f7-076f9fcddd25)



------------------------------------------------------------------------------------------------


**snATAC-seq Analysis**

01_snATAC_QC - Quality control of snATAC data based on number of fragments and TSS enrichment present.

02_snATAC_aggregateBin - Final QC and code for generating count matrices based on genomic bins (for initial clustering), gene bodies and promoters.

03_snATAC_Normalization - Initial clustering based on fragments from fixed-size genome wide bins.

04_snATAC_peakNormalization - Final peak calling based on initial clusters to generate high-quality peak set; used for final clustering and visualization.

05_snATAC_chromVARMotifs - Computing motif accessibility deviations using chromVAR (Schep et al., 2017) implemented in Signac.

06_snATAC_Compute_Gene_Scores - Computing gene activity scores utilizing Cicero; used for subsequent integration analyses. 

07_snATAC_Cluster_Unique_Peaks - Identification of cluster specific peaks.



<img width="825" alt="Screen Shot 2023-05-27 at 8 30 00 PM" src="https://github.com/aliamrod/Multimodal-Analysis/assets/62684338/65bed62b-fcd0-4a08-a95b-e3baaa229202">

Last Update: 01 May, 2023.
