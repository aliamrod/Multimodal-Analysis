# snRNA-seq QC
rm(list=ls())
input_dir<- '/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/input_dir/'
setwd(input_dir)
library(devtools)
#install.packages("stringfish")
#remotes::install_cran("qs", type = "source", configure.args = "--with-simd=AVX2")
library(tidyverse)
library(knitr)
library(patchwork)
library(Seurat)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(ggpubr)
#library(qs)

#CellRanger vs. CellBender vs. DoubletFinder data.

#AO1
ao1_cellranger<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO1/outs/filtered_feature_bc_matrix.h5")
ao1_cellbender<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO1/AO1_CELLBENDER2_filtered.h5")
ao1<- readRDS("./ao1_doub_removed.rds")

#AO2
ao2_cellranger<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO2/outs/filtered_feature_bc_matrix.h5")
ao2_cellbender<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO6/AO6_cellbender_filtered.h5")
ao2<- readRDS("./ao2_doub_removed.rds")

#AO3
ao3_cellranger<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO3/outs/filtered_feature_bc_matrix.h5")
ao3_cellbender<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO3/AO3_CELLBENDER2_filtered.h5")
ao3<- readRDS("./ao3_doub_removed.rds")

#AO4
ao4_cellranger<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO4/outs/filtered_feature_bc_matrix.h5")
ao4_cellbender<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO4/AO4_CELLBENDER2_filtered.h5")
ao4<- readRDS("./ao4_doub_removed.rds")

#AO5
ao5_cellranger<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO5/outs/filtered_feature_bc_matrix.h5")
ao5_cellbender<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO5/AO5_cellbender_filtered.h5")
ao5<- readRDS("./ao5_doub_removed.rds")

#AO6
ao6_cellranger<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO6/outs/filtered_feature_bc_matrix.h5")
ao6_cellbender<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO6/AO6_cellbender_filtered.h5")
ao6<- readRDS("./ao6_doub_removed.rds")

#AO7
ao7_cellranger<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO7/outs/filtered_feature_bc_matrix.h5")
ao7_cellbender<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO7/AO7_cellbendr_filtered.h5")
ao7<- readRDS("./ao7_doub_removed.rds")

#AO8
ao8_cellranger<- Read10X_h5("RNA/AO8/outs/filtered_feature_bc_matrix.h5")
ao8_cellbender<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO8/AO8_cellbender_filtered.h5")
ao8<- readRDS("./ao8_doub_removed.rds")

# Common nuclei.
length(intersect(colnames(ao1_cellranger), colnames(ao1_cellbender)))
length(intersect(colnames(ao2_cellranger), colnames(ao2_cellbender)))
length(intersect(colnames(ao3_cellranger), colnames(ao3_cellbender)))
length(intersect(colnames(ao4_cellranger), colnames(ao4_cellbender)))
length(intersect(colnames(ao5_cellranger), colnames(ao5_cellbender)))
length(intersect(colnames(ao6_cellranger), colnames(ao6_cellbender)))
length(intersect(colnames(ao7_cellranger), colnames(ao7_cellbender)))
length(intersect(colnames(ao8_cellranger), colnames(ao8_cellbender)))

library(ggplot2)
graph_plot_QC <- read.table(
  header=TRUE, text='Category        SampleID NucleiCount
1   a_CellRanger       AO1      24868
2   a_CellRanger       AO2      21370
3   a_CellRanger       AO3      22801
4   a_CellRanger       AO4      15102
5   a_CellRanger       AO5       15794
6   a_CellRanger       AO6       20934
7   a_CellRanger       AO7      11995
8   a_CellRanger        AO8      9240
9   b_CellBender      AO1       15097
10  b_CellBender      AO2      10826
11  b_CellBender      AO3      22860
12  b_CellBender      AO4      14875
13  b_CellBender      AO5      14991
14  b_CellBender      AO6      19647
15  b_CellBender      AO7      13443
16  b_CellBender      AO8      8554
17 c_DoubletFinder    AO1      14058
18 c_DoubletFinder    AO2      10067
19 c_DoubletFinder    AO3      21234
20 c_DoubletFinder    AO4      13806
21 c_DoubletFinder    AO5      14441
22 c_DoubletFinder    AO6      18230
23 c_DoubletFinder    AO7      12493
24 c_DoubletFinder    AO8      7950 
  ')

ggplot(graph_plot_QC, aes(SampleID, NucleiCount, fill = Category)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1")+ggtitle("Detected Nuclei Count")+DarkTheme()

ao1$Genotype<-"FLOX"
ao2$Genotype<-"CKO"
ao3$Genotype<-"CKO"
ao4$Genotype<-"FLOX"
ao5$Genotype<-"FLOX"
ao6$Genotype<-"FLOX"
ao7$Genotype<-"CKO"
ao8$Genotype<-"CKO"

ao1$SampleID<-"AO1"
ao2$SampleID<-"AO2"
ao3$SampleID<-"AO3"
ao4$SampleID<-"AO4"
ao5$SampleID<-"AO5"
ao6$SampleID<-"AO6"
ao7$SampleID<-"AO7"
ao8$SampleID<-"AO8"

ao1$Sample<-"AO1_FLOX"
ao2$Sample<-"AO2_CKO"
ao3$Sample<-"AO3_CKO"
ao4$Sample<-"AO4_FLOX"
ao5$Sample<-"AO5_FLOX"
ao6$Sample<-"AO6_FLOX"
ao7$Sample<-"AO7_CKO"
ao8$Sample<-"AO8_CKO"

ao1$library.batch<-"1"
ao2$library.batch<-"1"
ao3$library.batch<-"2"
ao4$library.batch<-"2"
ao5$library.batch<-"3"
ao6$library.batch<-"3"
ao7$library.batch<-"3"
ao8$library.batch<-"3"

# Merge by genotype + harmonize.
library(Seurat)
library(cowplot)
library(harmony)

cko.merged<-merge(ao2, y=c(ao3,ao7,ao8),
                  add.cell.ids=c("AO2", "AO3", "AO7", "AO8"))

flox.merged<-merge(ao1, y=c(ao4,ao5,ao6),
                   add.cell.ids=c("AO1", "AO4", "AO5", "AO6"))

cko.merged<- NormalizeData(cko.merged, normalization.method="LogNormalize", scale.factor=10000)
flox.merged<- NormalizeData(flox.merged, normalization.method="LogNormalize", scale.factor=10000)

# Scale data.
print("---Identify Variable Features---")
cko.merged<- FindVariableFeatures(cko.merged, 
                                  selection.method="vst",
                                  nfeatures=2000)
flox.merged<- FindVariableFeatures(flox.merged,
                                   selection.method="vst",
                                   nfeatures=2000)
print("---Regression---")
cko.genes<-rownames(cko.merged)
flox.genes<-rownames(flox.merged)

print("---Scale Data---")
cko.merged<-ScaleData(cko.merged, features=cko.genes)
flox.merged<-ScaleData(flox.merged, features=flox.genes)

print("---Run PCA---")
cko.merged<-RunPCA(cko.merged,
                  features=VariableFeatures(object=cko.merge),
                  npcs=60)
flox.merge<-RunPCA(flox.merged,
                   features=VariableFeatures(object=flox.merge),
                   npcs=60)

# Harmonize.
options(repr.plot.height = 2.5, repr.plot.width = 6)
cko.merged <- cko.merge %>% 
  RunHarmony("library.batch", plot_convergence = TRUE)

flox.merged <- flox.merge %>% 
  RunHarmony("library.batch", plot_convergence = TRUE)

# Integration.
library(Seurat)
library(patchwork)

merged<- merge(cko.merged, y=c(flox.merged), add.cell.ids = c("CKO", "FLOX"))

# Split the dataset into a list of two seurat objects (FLOX and CKO)
list <- SplitObject(merged, split.by = "Genotype")

# Normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)

anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
cortex.merged <- IntegrateData(anchorset = anchors)

#-------------

cortex.merged<- Add_Cell_Complexity_Seurat(seurat_object = cortex.merged)
cortex.merged@meta.data$group[cortex.merged@meta.data$SampleID =="AO1" |
                                cortex.merged@meta.data$SampleID == "AO4" |
                                cortex.merged@meta.data$SampleID == "AO5" |
                                cortex.merged@meta.data$SampleID == "AO6"]<- "FLOX"

cortex.merged@meta.data$group[cortex.merged@meta.data$SampleID == "AO2" |
                                cortex.merged@meta.data$SampleID == "AO3" |
                                cortex.merged@meta.data$SampleID == "AO7" |
                                cortex.merged@meta.data$SampleID == "AO8"]<- "cKO"

###### Correlation plots between replicates and samples #####


# Grouped by: cKO vs. FLOX groups.
qc_1<- QC_Plots_Genes(seurat_object = cortex.merged, low_cutoff = 300, high_cutoff = 8000, group.by = "group", pt.size = 0)
qc_2<- QC_Plots_UMIs(seurat_object = cortex.merged, low_cutoff = 1000, high_cutoff = 55000, group.by = "group", pt.size = 0)
ggarrange(qc_1,qc_2,
          labels = c("A","B"))


#Grouped by: CKO vs. FLX groups.
p1<- QC_Plots_Genes(seurat_object = cortex.merged, low_cutoff = 0,high_cutoff = 6000, group.by="group", pt.size=0)
p2<- QC_Plots_UMIs(seurat_object = cortex.merged, low_cutoff = 1000, high_cutoff = 55000,group.by="group",pt.size=0)
p3<- QC_Plots_Mito(seurat_object = cortex.merged, high_cutoff = 7.5, group.by="group",pt.size=0, plot_title="Mitochondrial Proportion (mtDNA%)")
p4<- QC_Plots_Complexity(seurat_object = cortex.merged, high_cutoff = 0.8, group.by="group",pt.size=0)
wrap_plots(p1,p2,p3,p4, ncol =4)


#Grouped by: Samples (8).
p5<- QC_Plots_Genes(seurat_object = cortex.merged,high_cutoff = 8000, group.by="SampleID", pt.size=0)
p6<- QC_Plots_UMIs(seurat_object = cortex.merged, high_cutoff = 30000,group.by="SampleID",pt.size=0)
p7<- QC_Plots_Mito(seurat_object = cortex.merged, high_cutoff = 2.5, group.by="SampleID", pt.size=0,plot_title="Mitochondrial Proportion (mtDNA%)")
p8<- QC_Plots_Complexity(seurat_object = cortex.merged, high_cutoff = 0.8, group.by="SampleID",pt.size=0)
wrap_plots(p5,p6,p7,p8, ncol =4)

cortex.merged <- subset(cortex.merged, nCount_RNA >= 200 & nCount_RNA <= 55000 & mitoPercent <=3.5)
median_stats_post<- Median_Stats(seurat_object = cortex.merged, group_by_var="SampleID",
                                 median_var = "SampleID")
df<-data.frame(median_stats_post)
write.csv(df, '/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/median_stats.csv')
###---Additional quality control---###
# Sample sizes.


############################################################
###---Cell Cycle Score identification---###
library(Seurat)
# Read in expression matrix. The first row is a header row, the first column is rownames.
# $ unzip /project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/input_dir/cell_cycle_files.zip

#exp.mat <- read.table(file = "./nestorawa_forcellcycle_expressionMatrix.txt", header=TRUE,
#as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
#s.genes<- str_to_title((tolower(cc.genes$s.genes)))
#g2m.genes<- str_to_title((tolower(cc.genes$g2m.genes)))

#cortex.merged<-CellCycleScoring(cortex.merged, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
#cortex.merged$CC.Difference<- (cortex.merged$S.Score)-(cortex.merged$G2M.Score)
#head(cortex.merged[[]]) # View cell cycle scores and phase assignments
#cortex.merged_final_filt<-ScaleData(cortex.merged,vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(cortex.merged))

###---Save filtered Seurat object---###
saveRDS(cortex.merged, '/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/input_dir/merged_snRNA_filtered.rds')
