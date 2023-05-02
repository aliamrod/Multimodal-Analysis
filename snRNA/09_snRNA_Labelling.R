# snRNA-seq Labelling
rm(list=ls())
setwd('/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/input_dir')
library(devtools)
library(tidyverse)
library(knitr)
library(patchwork)
library(Seurat)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(ggpubr)
library(stringr)
library(knitr)

umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)


#Load reference set.
# Read in reference dataset (Arlotta et al)
ref<-Read10X_h5("~/Downloads/GSM4635081_P1_S2_filtered_gene_bc_matrices_h5.h5")
ref<-CreateSeuratObject(ref)

# Load query set.
setwd('/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir')
query<-readRDS("./query_0404.rds")

# Renormalize.
query<- NormalizeData(query, normalization.method="LogNormalize", scale.factor=10000)

# Scale data.
query<- FindVariableFeatures(query, 
                             selection.method="vst",
                             nfeatures=2000)
genes<-rownames(query)
query<-ScaleData(query, features=genes)
query<-RunPCA(query,
              features=VariableFeatures(object=query),
              npcs=60)

# Manual Labelling.
canonical_markers <- list(
  'Radial glia' = c("Pax6", "Mki67", "Emx2"),
  'Astrocytes' = c("Apoe", "Aldh1l1", "Slc1a3"),
  'oRG' = c("Hopx", "Tnx", "Lgals2"),
  'Caudal ganglionic eminence (CGE)-derived interneurons' = c("Pax6", "Sp8", "Cxcl14", "Htr3a", "Nr2f2"),
  'Medial ganglionic eminence (MGE)-derived interneurons' = c("Sst", "Npy", "Lhx6", "Nxph1", "Nxph2", "Nkx2-1"),
  'CPN' = c("Satb2", "Inhba", "Frmd4b"),
  'CFuPN' = c("Bcl11b", "Crym", "Tle4"),
  'Cajal-Retzius Cells' = c("Nrg1", "Car10", "Edil3", "Flrt3", "Cdh4"),
  'Migrating Neurons' = c("Flna", "Dcx"),
  'Immature Neurons' = c("Sox2", "Nestin", "Cr", "Tuc-4", "Cb", "Dcx", "Pcna", "Ki-67"),
  'Stellate, Layer 4' = c("Rorb"),
  'CPN, Layers 5&6' = c("Rorb", "Fam19a2", "Cdh13", "Igsf21", "Gnb4"),
  'Ependymocytes' = c("Foxj1", "Pifo", "Dynlrb2"),
  'Red Blood Cells' = c("Car2", "Hemgn", "Hba1", "Hba2"),
  'SCPN' = c("Thsd7a"),
  'CPN, Layers 2&3' = c("Lamp5", "Cux2", "Glra3"),
  'CThPN, Layer 6' = c("Ntng2"),
  'Cycling Glial Cells' = c("Mki67", "Top2a", "Birc5"),
  'Apical Progenitors' = c("Cspg4", "Pdgfrb"),
  'Vascular and Leptomeningeal Cells (VLMC)' = c("Col1a1", "Vtn", "Lgals1"),
  'Oligodendrocyte Precursor Cells' = c("Olig1", "Olig2", "Pdgfra"),
  'Endothelial Cells' = c("Cldn5", "Mcam"),
  'Pericytes' = c("Cspg4", "Pdgfrb"),
  'Intermediate progenitors' = c("Eomes", "Neurog2", "Btg2", "Ppp1r17", "Tmem158"),
  'Interneurons(Inhib)' = c("Gad2","Gad1","Sp9","Npy","Sst","Cck", "Dlx1", "Dlx2", "Slc32a1"),
  'Neurons' = c("Satb2", "Mef2c", "Rbfox3", "Dcx"),
  'UL neurons' = c("Cux1","Cux2", "Satb2", "Rorb"),
  'DL neurons' = c("Bcl11b", "Foxp2","Tle4","Tbr1"),
  'L1' = c('Reln'),
  'Astrocytes' = c("Aqp4", "Aldh1l1"),
  'Microglia' = c("Trem2", "Mpeg1", "P2ry12", "Tmem119", "Aif1"),
  'Oligodendrocytes' = c("Olig2", "Pdgfra", "Sox10")
)

# Heatmap
library(viridis)
DefaultAssay(query)<-"RNA"
# DoHeatmap(query, group.by ="seurat_clusters", features=as.character(unlist(canonical_markers))) 

# Dot Plot
library(scCustomize)
DefaultAssay(query)<-'RNA'
Clustered_DotPlot(seurat_object = query, features = canonical_markers)

# Manual Labelling
new_ident <- setNames(c("Interneurons", "Neurons",
                        "CThPN, Layer 6",
                        "CThPN, Layer 6",
                        "CThPN, Layer 6",
                        "CThPN, Layer 6",
                        "CfuPN",
                        "Neurons", "CfuPN", "Neurons", 
                        "CfuPN",
                        "Interneurons",
                        "Cycling Radial Glial Cells", "Interneurons", 
                        "CfuPN", "CPN, Layers 5&6", "Oligodendrocytes",
                        "Interneurons",
                        "ScPN",
                        "CfuPN",
                        "CPN, Layers 5&6",
                        "CThPN, Layer 6",
                        "Cycling Radial Glial Cells", "Interneurons" ,"Astrocytes",
                        "CfuPN",
                        "Astrocytes",
                        "Interneurons",
                        "Unknown",
                        "Cycling Radial Glial Cells",
                        "CPN, Layers 5&6",
                        "Interneurons",
                        "Stellate, Layer 4",
                        "Cajal-Retzius Cells",
                        "Microglia",
                        "Ependymocytes",
                        "Ependymocytes", 
                        "Cycling Radial Glial Cells"),
                      levels(query))
query_02 <- RenameIdents(query, new_ident)
DimPlot(query_02, reduction = "umap", label = TRUE, label.box=T, label.size=2)+NoLegend()+ggtitle("Clusters in UMAP Space")

# Manual annotations
cluster_annotations <- list(
  '0' = "Interneurons",
  '1' = 'Neurons',
  '2' = 'CThPN, Layer 6',
  '3' = "CThPN, Layer 6",
  '4' = "CThPN, Layer 6",
  '5' = "CThPN, Layer 6",
  '6' = "CfuPN",
  '7' = "Neurons", '8' = "CfuPN", '9' = "Neurons", 
  '10' = "CfuPN",
  '11' = "Interneurons",
  '12' =  "Cycling Radial Glial Cells",'13'= "Interneurons", 
  '14' =  "CfuPN", '15' = "CPN, Layers 5&6", '16'="Oligodendrocytes",
  '17' =  "Interneurons",
  '18'= "ScPN",
  '19' = "CfuPN",
  '20' ="CPN",
  '21' =  "CThPN, Layer 6",
  '22'="Cycling Radial Glial Cells",'23'= "Interneurons" ,'24'="Astrocytes",
  '25'= "CfuPN",
  '26'="Astrocytes",
  '27'="Interneurons",
  '28' = "Cajal-Retzius Cells",
  '29' = "Cycling Radial Glial Cells",
  '30'="CPN, Layers 5&6",
  '31'= "Interneurons",
  '32' = "Stellate, Layer 4",
  '33'= "Cajal-Retzius Cells",
  '34' =  "Microglia",
  '35' = "Ependymocytes",
  '36' =  "Ependymocytes", 
  '37' =  "Cycling Radial Glial Cells"
)


# Manual annotations
cluster_annotations_cellClass <- list(
  '0' = "Inhibitory Neurons",
  '1' = 'Neurons',
  '2' = 'Excitatory Neurons',
  '3' = "Excitatory Neurons",
  '4' = "Excitatory Neurons",
  '5' = "Excitatory Neurons",
  '6' = "Excitatory Neurons",
  '7' = "Neurons", '8' = "Excitatory Neurons", '9' = "Neurons", 
  '10' = "Excitatory Neurons",
  '11' = "Inhibitory Neurons",
  '12' =  "Cycling Radial Glial Cells",'13'= "Inhibitory Neurons", 
  '14' =  "Excitatory Neurons", '15' = "Excitatory Neurons", '16'="Oligodendrocytes",
  '17' =  "Inhibitory Neurons",
  '18'= "Excitatory Neurons",
  '19' = "Excitatory Neurons",
  '20' ="Excitatory Neurons",
  '21' =  "Excitatory Neurons",
  '22'="Cycling Radial Glial Cells",'23'= "Inhibitory Neurons" ,'24'="Astrocytes",
  '25'= "Excitatory Neurons",
  '26'="Astrocytes",
  '27'="Inhibitory Neurons",
  '28' = "Unknown",
  '29' = "Cycling Radial Glial Cells",
  '30'="Excitatory Neurons",
  '31'= "Inhibitory Neurons",
  '32' = "Excitatory Neurons",
  '33'= "Cajal-Retzius Cells",
  '34' =  "Microglia",
  '35' = "Ependymocytes",
  '36' =  "Ependymocytes", 
  '37' =  "Cycling Radial Glial Cells"
)

# add CellType to seurat metadata
query_02$cellClass <- unlist(cluster_annotations_cellClass[query_02$seurat_clusters])
DimPlot(query_02, reduction = "umap", label = TRUE, label.box=T, label.size=2, group.by="cellClass")+NoLegend()+ggtitle("Clusters in UMAP Space")

saveRDS(query_02, "./labelled_query_02_0410.rds")

# Process Seurat file.
DefaultAssay(ref)<-"RNA"

ref <- FindVariableFeatures(
  ref,
  selection.method = "vst",
  nfeatures = 4000
)

all.genes <- rownames(ref)
ref <- ScaleData(ref, features = all.genes)

ref <- RunPCA(
  ref,
  features = VariableFeatures(object = ref),
  npcs=50
)

ref <- RunUMAP(ref, reduction = "pca", dims = 1:30)
ref <- FindNeighbors(ref, reduction = "pca", dims = 1:30)
ref <- FindClusters(ref,resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))

# Determine percent of variation associated with each PC
pct<- seur.ref[["pca"]]@stdev/sum(ref[["pca"]]@stdev)*100

# Calculate cumulative percents for each PC.
cumu<-cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and %variation associated with the PC as less than 5.
co1<-which(cumu>90&pct<5)[1]

# Determine the difference between variation of PC and subsequent PC.
# co2 = last point where change of % of variation is more than 0.1%.
co2<- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# Minimum of co1 versus co2.
pcs<-min(co1,co2)
numPCs=pcs

print("---Make tSNE.---")
ref<-RunTSNE(object=ref,dims.use=1:numPCs, check_duplicates=F)
write.table(ref[["pca"]]@feature.loadings,'./VST_PCA.txt',quote=F,row.names=T,col.names=T,sep='\t')

#Read reference meta.data
ref.path<-"/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/input_dir/metaData_scDevSC.txt"
ref.tab<- read.delim(file = ref.path, header=T, sep="\t")
ref.subset1<- ref.tab[grep("P1_S2", ref.tab$NAME), ]
ref.subset2<- ref.subset1[c(1,17)]

# Remove string from rows
#ref.subset2$NAME<- gsub("P1_S1_", "", as.character(ref.subset2$NAME))
ref.subset2$NAME<- gsub("P1_S2_", "", as.character(ref.subset2$NAME))

colnames(ref.subset2)<- c("CellID", "CellType")
ref.subset2$CellID<- sub("$", "-1", as.character(ref.subset2$CellID))

ref.subset3<- as.data.frame(ref.subset2)
ref.df<- ref.subset3

new_frame<- merge(as.data.frame(ref@meta.data), ref.df,
                  by.x = "row.names", by.y = "CellID")

row.names(new_frame)<- new_frame$Row.names
new_frame$Row.names <- NULL

# Create a vector.
newcelltype<- new_frame$CellType
names(newcelltype)<- row.names(new_frame)

# View.
table(newcelltype)

# Integrate 'newcelltype' vector into P1 meta.data. 
ref$newCellType<-newcelltype
table(ref@meta.data$newCellType)

Idents(ref)<- ref$newCellType

# LabelTransfer: transfer label onto new frame.
reference<-ref
query<-query

DefaultAssay(reference)<-"RNA"
DefaultAssay(query)<-"RNA"

anchors<-FindTransferAnchors(reference=reference, query=query, dims=1:30, npcs=30)
predictions<-TransferData(anchorset=anchors, refdata=reference$newCellType, dims=1:30)
query$newCellType_transfer<-predictions$predicted.id

# LabelTransfer final UMAP.
plot1 <- UMAPPlot(query, label=T)+NoLegend()
plot2 <- UMAPPlot(query, group.by="newCellType_transfer", label.box=T, label.size=3,label=T,repel=T)+ggtitle("CellType Transfer (Arlotta et al)")
plot1 | plot2


# Save object.
saveRDS(reference, "/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/input_dir/ARLOTTA_REF.rds")
saveRDS(query, "/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/input_dir/label_transfer_02.rds")
