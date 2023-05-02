# Load required packages onto BioHPC.

# Open terminal and enter the following command-line statements:
# $ module purge && module load shared slurm python/3.7.x-anaconda
# $ module load R/4.1.1-gccmkl
# $ module load hdf5_18/1.8.17
# $ module load gcc/8.3.0
# $ module load htslib
# $ module   load gsl/2.4
# $ module load macs/2.1.2
# $ module load rstudio-desktop/1.1.456

# Open R/4.1.1-gccmkl module.
# $ rstudio
library(devtools)
library(BiocManager)
library(dplyr)
library(Seurat)
library(ggplot2)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(tidyverse)  
rm(list = ls())



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########## SAMPLE AO2 ##########
# Doublet Finder
rm(list=ls())
library(devtools)
library(BiocManager)
library(dplyr)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(tidyverse)  
# Doublet Finder
rm(list=ls())
library(devtools)
library(BiocManager)
library(dplyr)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(tidyverse)  


setwd("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO2")


list.files() 

# Create counts matrix
#cts_sample_ao2 <- ReadMtx(mtx = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/matrix.mtx", 
#        features = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/features.tsv",
#        cells = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/barcodes.tsv")
cts_sample_ao2<- Read10X_h5("AO2_cellbender_filtered.h5")
cts_sample_ao2[1:10, 1:10] # inspect count sparse matrix

# Create Seurat object. ---------------------------------------------------------------------------------------
ao2.seurat<- CreateSeuratObject(counts = cts_sample_ao2)

# QC and Filtering. ---------------------------------------------------------------------------------------
ao2.seurat$mitoPercent<- PercentageFeatureSet(ao2.seurat, pattern = '^mt-')

# Removing cells which report low number of molecules + low number of genes + high mitochondrial count.
#ao2.seurat.filtered<-  
#  subset(ao2.seurat, subset = nCount_RNA > 800 &
#           nFeature_RNA > 500 & 
#           mitoPercent < 10)

# Pre-process standard workflow. ---------------------------------------------------------------------------------------
ao2.seurat.filtered<- NormalizeData(object = ao2.seurat.filtered)
ao2.seurat.filtered<- FindVariableFeatures(object = ao2.seurat.filtered)
ao2.seurat.filtered<- ScaleData(object = ao2.seurat.filtered)
ao2.seurat.filtered<- RunPCA(object = ao2.seurat.filtered)
ElbowPlot(ao2.seurat.filtered) # captures variation data in accordance to dimensionality
ao2.seurat.filtered<- FindNeighbors(object = ao2.seurat.filtered, dims = 1:20)
ao2.seurat.filtered<- FindClusters(object = ao2.seurat.filtered)
ao2.seurat.filtered<- RunUMAP(object = ao2.seurat.filtered, dims = 1:20)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
# Find optimal pK value (bimodality coefficient and the highest value = optimal pK value). 
sweep.res.list_ao2 <- paramSweep_v3(ao2.seurat.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_ao2 <- summarizeSweep(sweep.res.list_ao2, GT = FALSE)
bcmvn_ao2 <- find.pK(sweep.stats_ao2)

# BCmetric on y-axis; pK values on x-axis. pK value = optimal pK value)
ggplot(bcmvn_ao2, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_ao2 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- ao2.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)   
nExp_poi <- round(0.076*nrow(ao2.seurat.filtered@meta.data))  ## Assuming 3.9% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# Run DoubletFinder. 
ao2.seurat.filtered <- doubletFinder_v3(ao2.seurat.filtered, 
                                        PCs = 1:20, 
                                        pN = 0.25, 
                                        pK = pK, 
                                        nExp = nExp_poi.adj,
                                        reuse.pANN = FALSE, sct = FALSE)
# Look at predictions. Singlet vs. Doublet. 
View(ao2.seurat.filtered@meta.data)

# Visualize doublets.
unfiltered.dim<- DimPlot(ao2.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_450")

# number of singlets and doublets
table(ao2.seurat.filtered@meta.data$DF.classifications_0.25_0.3_450)

vlnplot<- VlnPlot(ao2.seurat.filtered, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.3_450")

# Remove all predicted doublets from data. 

## data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
data.filt.ao2 = ao2.seurat.filtered[, ao2.seurat.filtered@meta.data[, "DF.classifications_0.25_0.3_450"] == "Singlet"]
dim(data.filt.ao2)


filtered.dim<- DimPlot(data.filt.ao, reduction='umap', group.by = "DF.classifications_0.25_0.3_450")

# Plots
(unfiltered.dim / filtered.dim) | vlnplot
# Save data
# Save QC-filtered data for further analysis. Create output directory results and save data to that folder.

# dir.create("data/results", showWarnings = F)
saveRDS(data.filt.ao6, "/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ao6_NEW_doub_removed.rds")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########## SAMPLE AO5 ##########
# Doublet Finder
rm(list=ls())
library(devtools)
library(BiocManager)
library(dplyr)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(tidyverse)  


setwd("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO5")


list.files() 

# Create counts matrix
#cts_sample_ao2 <- ReadMtx(mtx = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/matrix.mtx", 
#        features = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/features.tsv",
#        cells = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/barcodes.tsv")
cts_sample_ao5<- Read10X_h5("AO5_cellbender_filtered.h5")
cts_sample_ao5[1:10, 1:10] # inspect count sparse matrix


# Create Seurat object. ---------------------------------------------------------------------------------------
ao5.seurat<- CreateSeuratObject(counts = cts_sample_ao5)

# QC and Filtering. ---------------------------------------------------------------------------------------
ao5.seurat$mitoPercent<- PercentageFeatureSet(ao5.seurat, pattern = '^mt-')

# Removing cells which report low number of molecules + low number of genes + high mitochondrial count.
ao5.seurat.filtered<-  
  subset(ao5.seurat, subset = nCount_RNA > 800 &
           nFeature_RNA > 500 & 
           mitoPercent < 10)

# Pre-process standard workflow. ---------------------------------------------------------------------------------------
ao5.seurat.filtered<- NormalizeData(object = ao5.seurat.filtered)
ao5.seurat.filtered<- FindVariableFeatures(object = ao5.seurat.filtered)
ao5.seurat.filtered<- ScaleData(object = ao5.seurat.filtered)
ao5.seurat.filtered<- RunPCA(object = ao5.seurat.filtered)
ElbowPlot(ao5.seurat.filtered) # captures variation data in accordance to dimensionality
ao5.seurat.filtered<- FindNeighbors(object = ao5.seurat.filtered, dims = 1:20)
ao5.seurat.filtered<- FindClusters(object = ao5.seurat.filtered)
ao5.seurat.filtered<- RunUMAP(object = ao5.seurat.filtered, dims = 1:20)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------

# Find optimal pK value (bimodality coefficient and the highest value = optimal pK value). 
sweep.res.list_ao5 <- paramSweep_v3(ao5.seurat.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_ao5 <- summarizeSweep(sweep.res.list_ao5, GT = FALSE)
bcmvn_ao5 <- find.pK(sweep.stats_ao5)

# BCmetric on y-axis; pK values on x-axis. pK value = optimal pK value)
ggplot(bcmvn_ao5, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_ao5 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- ao5.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)   
nExp_poi <- round(0.039*nrow(ao5.seurat.filtered@meta.data))  ## Assuming 3.9% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# Run DoubletFinder. 
ao5.seurat.filtered <- doubletFinder_v3(ao5.seurat.filtered, 
                                        PCs = 1:20, 
                                        pN = 0.25, 
                                        pK = pK, 
                                        nExp = nExp_poi.adj,
                                        reuse.pANN = FALSE, sct = FALSE)
# Look at predictions. Singlet vs. Doublet. 
View(ao5.seurat.filtered@meta.data)

# Visualize doublets.
unfiltered.dim<- DimPlot(ao5.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_499")

# number of singlets and doublets
table(ao5.seurat.filtered@meta.data$DF.classifications_0.25_0.3_499)

vlnplot<- VlnPlot(ao5.seurat.filtered, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.3_499")

# Remove all predicted doublets from data. 

## data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
data.filt.ao5 = ao5.seurat.filtered[, ao5.seurat.filtered@meta.data[, "DF.classifications_0.25_0.3_499"] == "Singlet"]
dim(data.filt.ao5)


filtered.dim<- DimPlot(data.filt.ao5, reduction='umap', group.by = "DF.classifications_0.25_0.3_499")

# Plots
(unfiltered.dim / filtered.dim) | vlnplot
# Save data
# Save QC-filtered data for further analysis. Create output directory results and save data to that folder.

# dir.create("data/results", showWarnings = F)
saveRDS(data.filt.ao5, "/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ao5_doub_removed.rds")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########## SAMPLE AO6 ##########
# Doublet Finder
rm(list=ls())
library(devtools)
library(BiocManager)
library(dplyr)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(tidyverse)  


setwd("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO6")


list.files() 

# Create counts matrix
#cts_sample_ao2 <- ReadMtx(mtx = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/matrix.mtx", 
#        features = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/features.tsv",
#        cells = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/barcodes.tsv")
cts_sample_ao6<- Read10X_h5("AO6_cellbender_filtered.h5")
cts_sample_ao6[1:10, 1:10] # inspect count sparse matrix


# Create Seurat object. ---------------------------------------------------------------------------------------
ao6.seurat<- CreateSeuratObject(counts = cts_sample_ao6)

# QC and Filtering. ---------------------------------------------------------------------------------------
ao6.seurat$mitoPercent<- PercentageFeatureSet(ao6.seurat, pattern = '^mt-')

# Removing cells which report low number of molecules + low number of genes + high mitochondrial count.
ao6.seurat.filtered<-  
  subset(ao6.seurat, subset = nCount_RNA > 800 &
           nFeature_RNA > 500 & 
           mitoPercent < 10)

# Pre-process standard workflow. ---------------------------------------------------------------------------------------
ao6.seurat.filtered<- NormalizeData(object = ao6.seurat.filtered)
ao6.seurat.filtered<- FindVariableFeatures(object = ao6.seurat.filtered)
ao6.seurat.filtered<- ScaleData(object = ao6.seurat.filtered)
ao6.seurat.filtered<- RunPCA(object = ao6.seurat.filtered)
ElbowPlot(ao6.seurat.filtered) # captures variation data in accordance to dimensionality
ao6.seurat.filtered<- FindNeighbors(object = ao6.seurat.filtered, dims = 1:20)
ao6.seurat.filtered<- FindClusters(object = ao6.seurat.filtered)
ao6.seurat.filtered<- RunUMAP(object = ao6.seurat.filtered, dims = 1:20)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
# Find optimal pK value (bimodality coefficient and the highest value = optimal pK value). 
sweep.res.list_ao6 <- paramSweep_v3(ao6.seurat.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_ao6 <- summarizeSweep(sweep.res.list_ao6, GT = FALSE)
bcmvn_ao6 <- find.pK(sweep.stats_ao6)

# BCmetric on y-axis; pK values on x-axis. pK value = optimal pK value)
ggplot(bcmvn_ao6, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_ao6 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- ao6.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)   
nExp_poi <- round(0.039*nrow(ao6.seurat.filtered@meta.data))  ## Assuming 3.9% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# Run DoubletFinder. 
ao6.seurat.filtered <- doubletFinder_v3(ao6.seurat.filtered, 
                                        PCs = 1:20, 
                                        pN = 0.25, 
                                        pK = pK, 
                                        nExp = nExp_poi.adj,
                                        reuse.pANN = FALSE, sct = FALSE)
# Look at predictions. Singlet vs. Doublet. 
View(ao6.seurat.filtered@meta.data)

# Visualize doublets.
unfiltered.dim<- DimPlot(ao6.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_692")

# number of singlets and doublets
table(ao6.seurat.filtered@meta.data$DF.classifications_0.25_0.005_692)

vlnplot<- VlnPlot(ao6.seurat.filtered, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.005_692")

# Remove all predicted doublets from data. 

## data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
data.filt.ao6 = ao6.seurat.filtered[, ao6.seurat.filtered@meta.data[, "DF.classifications_0.25_0.005_692"] == "Singlet"]
dim(data.filt.ao6)


filtered.dim<- DimPlot(data.filt.ao6, reduction='umap', group.by = "DF.classifications_0.25_0.005_692")

# Plots
(unfiltered.dim / filtered.dim) | vlnplot
# Save data
# Save QC-filtered data for further analysis. Create output directory results and save data to that folder.

# dir.create("data/results", showWarnings = F)
saveRDS(data.filt.ao6, "/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ao6_doub_removed.rds")


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########## SAMPLE AO7 ##########
# Doublet Finder
rm(list=ls())
library(devtools)
library(BiocManager)
library(dplyr)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(tidyverse)  


setwd("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO7")


list.files() 

# Create counts matrix
#cts_sample_ao2 <- ReadMtx(mtx = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/matrix.mtx", 
#        features = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/features.tsv",
#        cells = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/barcodes.tsv")
cts_sample_ao7<- Read10X_h5("AO7_cellbender_filtered.h5")
cts_sample_ao7[1:10, 1:10] # inspect count sparse matrix


# Create Seurat object. ---------------------------------------------------------------------------------------
ao7.seurat<- CreateSeuratObject(counts = cts_sample_ao7)

# QC and Filtering. ---------------------------------------------------------------------------------------
ao7.seurat$mitoPercent<- PercentageFeatureSet(ao7.seurat, pattern = '^mt-')

# Removing cells which report low number of molecules + low number of genes + high mitochondrial count.
ao7.seurat.filtered<-  
  subset(ao7.seurat, subset = nCount_RNA > 800 &
           nFeature_RNA > 500 & 
           mitoPercent < 10)

# Pre-process standard workflow. ---------------------------------------------------------------------------------------
ao7.seurat.filtered<- NormalizeData(object = ao7.seurat.filtered)
ao7.seurat.filtered<- FindVariableFeatures(object = ao7.seurat.filtered)
ao7.seurat.filtered<- ScaleData(object = ao7.seurat.filtered)
ao7.seurat.filtered<- RunPCA(object = ao7.seurat.filtered)
ElbowPlot(ao7.seurat.filtered) # captures variation data in accordance to dimensionality
ao7.seurat.filtered<- FindNeighbors(object = ao7.seurat.filtered, dims = 1:20)
ao7.seurat.filtered<- FindClusters(object = ao7.seurat.filtered)
ao7.seurat.filtered<- RunUMAP(object = ao7.seurat.filtered, dims = 1:20)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
# Find optimal pK value (bimodality coefficient and the highest value = optimal pK value). 
sweep.res.list_ao7 <- paramSweep_v3(ao7.seurat.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_ao7 <- summarizeSweep(sweep.res.list_ao7, GT = FALSE)
bcmvn_ao7 <- find.pK(sweep.stats_ao7)

# BCmetric on y-axis; pK values on x-axis. pK value = optimal pK value)
ggplot(bcmvn_ao7, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_ao7 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- ao7.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)   
nExp_poi <- round(0.039*nrow(ao7.seurat.filtered@meta.data))  ## Assuming 3.9% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# Run DoubletFinder. 
ao7.seurat.filtered <- doubletFinder_v3(ao7.seurat.filtered, 
                                        PCs = 1:20, 
                                        pN = 0.25, 
                                        pK = pK, 
                                        nExp = nExp_poi.adj,
                                        reuse.pANN = FALSE, sct = FALSE)
# Look at predictions. Singlet vs. Doublet. 
View(ao7.seurat.filtered@meta.data)

# Visualize doublets.
unfiltered.dim<- DimPlot(ao7.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.06_452")

# number of singlets and doublets
table(ao7.seurat.filtered@meta.data$DF.classifications_0.25_0.06_452)

vlnplot<- VlnPlot(ao7.seurat.filtered, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.06_452")

# Remove all predicted doublets from data. 

## data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
data.filt.ao7 = ao7.seurat.filtered[, ao7.seurat.filtered@meta.data[, "DF.classifications_0.25_0.06_452"] == "Singlet"]
dim(data.filt.ao7)


filtered.dim<- DimPlot(data.filt.ao7, reduction='umap', group.by = "DF.classifications_0.25_0.06_452")

# Plots
(unfiltered.dim / filtered.dim) | vlnplot
# Save data
# Save QC-filtered data for further analysis. Create output directory results and save data to that folder.

# dir.create("data/results", showWarnings = F)
saveRDS(data.filt.ao7, "/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ao7_doub_removed.rds")



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########## SAMPLE AO8 ##########
# Doublet Finder
rm(list=ls())
library(devtools)
library(BiocManager)
library(dplyr)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(tidyverse)  


setwd("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO6")


list.files() 

# Create counts matrix
#cts_sample_ao2 <- ReadMtx(mtx = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/matrix.mtx", 
#        features = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/features.tsv",
#        cells = "/work/Neuroinformatics_Core/s204365/snRNA/snRNA108/SAMPLE_AO24_FLX/outs/raw_feature_bc_matrix/barcodes.tsv")
cts_sample_ao8<- Read10X_h5("AO8_cellbender_filtered.h5")
cts_sample_ao8[1:10, 1:10] # inspect count sparse matrix


# Create Seurat object. ---------------------------------------------------------------------------------------
ao8.seurat<- CreateSeuratObject(counts = cts_sample_ao8)

# QC and Filtering. ---------------------------------------------------------------------------------------
ao8.seurat$mitoPercent<- PercentageFeatureSet(ao8.seurat, pattern = '^mt-')

# Removing cells which report low number of molecules + low number of genes + high mitochondrial count.
ao8.seurat.filtered<-  
  subset(ao8.seurat, subset = nCount_RNA > 800 &
           nFeature_RNA > 500 & 
           mitoPercent < 10)

# Pre-process standard workflow. ---------------------------------------------------------------------------------------
ao8.seurat.filtered<- NormalizeData(object = ao8.seurat.filtered)
ao8.seurat.filtered<- FindVariableFeatures(object = ao8.seurat.filtered)
ao8.seurat.filtered<- ScaleData(object = ao8.seurat.filtered)
ao8.seurat.filtered<- RunPCA(object = ao8.seurat.filtered)
ElbowPlot(ao8.seurat.filtered) # captures variation data in accordance to dimensionality
ao8.seurat.filtered<- FindNeighbors(object = ao8.seurat.filtered, dims = 1:20)
ao8.seurat.filtered<- FindClusters(object = ao8.seurat.filtered)
ao8.seurat.filtered<- RunUMAP(object = ao8.seurat.filtered, dims = 1:20)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
# Find optimal pK value (bimodality coefficient and the highest value = optimal pK value). 
sweep.res.list_ao8 <- paramSweep_v3(ao8.seurat.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_ao8 <- summarizeSweep(sweep.res.list_ao8, GT = FALSE)
bcmvn_ao8 <- find.pK(sweep.stats_ao8)

# BCmetric on y-axis; pK values on x-axis. pK value = optimal pK value)
ggplot(bcmvn_ao8, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_ao8 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- ao8.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)   
nExp_poi <- round(0.039*nrow(ao8.seurat.filtered@meta.data))  ## Assuming 3.9% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# Run DoubletFinder. 
ao8.seurat.filtered <- doubletFinder_v3(ao8.seurat.filtered, 
                                        PCs = 1:20, 
                                        pN = 0.25, 
                                        pK = pK, 
                                        nExp = nExp_poi.adj,
                                        reuse.pANN = FALSE, sct = FALSE)
# Look at predictions. Singlet vs. Doublet. 
View(ao8.seurat.filtered@meta.data)

# Visualize doublets.
unfiltered.dim<- DimPlot(ao8.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.21_307")

# number of singlets and doublets
table(ao8.seurat.filtered@meta.data$DF.classifications_0.25_0.21_307)

vlnplot<- VlnPlot(ao8.seurat.filtered, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.21_307")

# Remove all predicted doublets from data. 

## data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
data.filt.ao8 = ao8.seurat.filtered[, ao8.seurat.filtered@meta.data[, "DF.classifications_0.25_0.21_307"] == "Singlet"]
dim(data.filt.ao8)


filtered.dim<- DimPlot(data.filt.ao8, reduction='umap', group.by = "DF.classifications_0.25_0.21_307")

# Plots
(unfiltered.dim / filtered.dim) | vlnplot
# Save data
# Save QC-filtered data for further analysis. Create output directory results and save data to that folder.

# dir.create("data/results", showWarnings = F)
saveRDS(data.filt.ao8, "/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ao8_doub_removed.rds")
