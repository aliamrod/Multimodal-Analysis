# snRNA-seq UMAP
rm(list=ls())
library(Seurat)
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
library(clustree)

setwd("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir/")
cortex.merged<-readRDS('/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir/merged_snRNA_filtered.rds')

# Run the standard workflow for visualization and clustering
resolutions<-c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2)

cortex.merged<- ScaleData(cortex.merged, verbose = FALSE)
cortex.merged <- RunPCA(cortex.merged, npcs = 30, verbose = FALSE)
cortex.merged <- RunUMAP(cortex.merged, reduction = "pca", dims = 1:30)
cortex.merged <- FindNeighbors(cortex.merged, reduction = "pca", dims = 1:30)
cortex.merged <- FindClusters(cortex.merged, resolution = resolutions)
#------------------
###---Cell Cycle Score identification---###
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

# Determine percent of variation associated with each PC
pct<- cortex.merged[["pca"]]@stdev/sum(cortex.merged[["pca"]]@stdev)*100

# Calculate cumulative percents for each PC.
cumu<-cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and %variation associated with the PC as less than 5.
co1<-which(cumu>90&pct<5)[1]

# Determine the difference between variation of PC and subsequent PC.
# co2 = last point where change of % of variation is more than 0.1%.
co2<- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# Minimum of co1 versus co2.
pcs<-min(co1,co2)

# Dataframe with values.
pcs_df<-data.frame(pct=pct,
                   cumu=cumu,
                   rank=1:length(pct))

#Elbow plot to visualize pcs_dataframe.
pcs_elbow<- ggplot(pcs_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

# Change based on PCElbowPlot.
library(Matrix)
library(reshape2)
library(igraph)
library(RANN)
numPCs=pcs

print("---Make tSNE.---")
cortex.merged<-RunTSNE(object=cortex.merged,dims.use=1:numPCs, check_duplicates=F)
tsne<- TSNEPlot(object=cortex.merged)
write.table(cortex.merged[["pca"]]@feature.loadings,'./VST_PCA.txt',quote=F,row.names=T,col.names=T,sep='\t')

# Cluster selection (ie., Clustree)
head(cortex.merged[[]]) # Clustering information contained in 'meta.data' slot accessed using the [[]] operator.

# Clustering tree colored by clustering resolution.
clustree(cortex.merged)

# Clustering tree colored by the SC3 stability metric.
clustree(cortex.merged, 
         node_colour="sc3_stability")+ggtitle("SC3 Stability Index")
# Clustering tree colored by the expression of marker genes.

# Grouped by: cKO vs. FLOX groups.
qc_1<- QC_Plots_Genes(seurat_object = cortex.merged, low_cutoff = 0, high_cutoff = 6000, group.by = "group", pt.size = 0.1)
qc_2<- QC_Plots_UMIs(seurat_object = cortex.merged, low_cutoff = 0, high_cutoff = 30000, group.by = "group", pt.size = 0.1)
ggarrange(qc_1,qc_2,
          labels = c("A","B"))

#Grouped by: CKO vs. FLX groups.
# Also, implement mitochondrial and ribosomal gene percentages.
cortex.merged <- Add_Mito_Ribo_Seurat(seurat_object = cortex.merged, species = "Mouse")

p1<- QC_Plots_Genes(seurat_object = cortex.merged, low_cutoff = 0,high_cutoff = 6000, group.by="group", pt.size=0.1)
p2<- QC_Plots_UMIs(seurat_object = cortex.merged, high_cutoff = 30000,group.by="group",pt.size=0.1)
p3<- QC_Plots_Mito(seurat_object = cortex.merged, pt.size = 0.1, group.by="group")
p4<- QC_Plots_Complexity(seurat_object = cortex.merged, high_cutoff = 0.82, group.by="group",pt.size=0)+geom_violin(trim=FALSE, fill='black', color="darkred")
qc_by_group<-ggarrange(p1,p2,p3,p4,
                       ncol=4, nrow=1)

#Grouped by: Samples (8).
p5<- QC_Plots_Genes(seurat_object = cortex.merged,high_cutoff = 6000, group.by="SampleID", pt.size=0.1)
p6<- QC_Plots_UMIs(seurat_object = cortex.merged, high_cutoff = 28000,group.by="SampleID",pt.size=0.1)
p7<- QC_Plots_Mito(seurat_object = cortex.merged, pt.size = 0.1, group.by="SampleID")
p8<- QC_Plots_Complexity(seurat_object = cortex.merged, high_cutoff = 0.8, group.by="SampleID",pt.size=0)+geom_violin(trim=FALSE, fill='black', color="darkred")
qc_by_sample<-ggarrange(p5,p6,p7,p8,
                        ncol=4, nrow=1)

p9 <- QC_Plot_UMIvsGene(seurat_object=cortex.merged, low_cutoff_gene=0, high_cutoff_gene =6000, low_cutoff_UMI = 500,
                        high_cutoff_UMI=30000, group.by="SampleID")
p10 <- QC_Plot_GenevsFeature(seurat_object=cortex.merged, feature1="percent_mito", low_cutoff_gene=0,
                             high_cutoff_gene = 6000, high_cutoff_feature = 20, group.by="SampleID")

ggarrange(p9, p10,
          ncol=2,nrow=1)

############################################################
# Re-run UMAP visualization with selected clustering resolution.
# Selected resolution = 1.0.
selected_resolution<-1.0
cortex.merged<-FindClusters(cortex.merged,resolution=selected_resolution,n.start=18, n.iter=20,algorithm=1)
cortex.merged<-RunUMAP(cortex.merged,dims=1:20,min.dist=0.5,spread=1.5,components=2L)

dimplot_01<-Seurat::DimPlot(cortex.merged, reduction="umap", label=T, label.box=T, repel=TRUE)+ggtitle("Clusters in UMAP Space")
dimplot_02<-Seurat::DimPlot(cortex.merged, reduction="umap", label=T, label.box=T, repel=TRUE,group.by="Genotype", cols=c("red", "grey"))

dimplot_01+dimplot_02


# Visualization
png('/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/output_dir/umap_0403.png', width=7,height=7, res=300, units='in')
print(dimplot_01+dimplot_02)
dev.off()

# Cluster validation.
# Total number of cells in each cluster.
CKO <- sum(cortex.merged@meta.data$Genotype == "CKO")
FLOX <- sum(cortex.merged@meta.data$Genotype == "FLOX")

plot.data <- cortex.merged@meta.data %>%
  select(Genotype, cluster = paste0("RNA_snn_res.", 0.8)) %>%
  mutate(cluster = factor(as.numeric(cluster))) %>%
  group_by(cluster, Genotype) %>%
  summarise(count = n()) %>%
  mutate(clust_total = sum(count)) %>%
  mutate(clust_prop = count / clust_total) %>%
  ggplot(plot.data, aes(x = cluster, y = count, fill = Genotype)) +
  geom_col()+ggtitle("Genotypic Proportion per Cluster")+DarkTheme()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

genotypic_prop<-ggplot(plot.data, aes(x = cluster, y = clust_prop, fill = Genotype)) +
  geom_col()+scale_fill_manual(values=c(
    "red",
    "lightblue"))+DarkTheme()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+labs(x="Cluster", y="Proportional Composition")+ggtitle("Genotypic Proportion per Cluster ")



# Total counts per cluster
cluster_qc<-QC_Plots_Combined_Vln(seurat_object = cortex.merged, group.by="seurat_clusters", pt.size=0)
png(cluster_qc, '/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/Plots/cluster_qc.png')
dev.off()

# Heatmap
library(DESeq2)
library(presto)
library(dplyr)
library(presto)

DefaultAssay(cortex.merged)<-'RNA'

cl_markers <- FindAllMarkers(cortex.merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
##cl_markers %>% group_by(cluster %>% top_n(n = 2, wt = avg_log2FC))

cl_markers %>%
  filter(avg_log2FC > log(1.2) & pct.1>20 & p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group=T) %>%
  top_n(n=2, wt=avg_log2FC) %>%
  print(n=40, width=Inf)

top10_cl_markers<-cl_markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(cortex.merged, features=top10_cl_markers$gene) + NoLegend()



# SingleR
library(SingleR)
ref <- celldex::MouseRNAseqData()
results<- SingleR(test = as.SingleCellExperiment(cortex.merged), ref = ref, labels = ref$label.fine)
cortex.merged$singler_labels <- results$labels

singler<-DimPlot_scCustom(cortex.merged, reduction='umap', group.by = 'singler_labels',
                          label = TRUE, label.box=T,
                          label.size=1,
                          label.color="black")+ggtitle(
                            "SingleR Labelling")

# Tabula Muris Dataset.
library(scRNAseq)
library(TabulaMurisData)
library(ExperimentHub)
library(scuttle)
eh<- ExperimentHub()
query(eh, "TabulaMurisData")
cortex_ref<- eh[['EH1617']]
cortex_ref<- cortex_ref[, !is.na(cortex_ref$cell_ontology_class)]

cortex_ref<- logNormCounts(cortex_ref)
results_scuttle<- SingleR(test = as.SingleCellExperiment(cortex.merged), ref = cortex_ref, labels = 
                            cortex_ref$cell_ontology_id)

corrtex.merged$scuttle_labels<- results_scuttle$labels

DimPlot(cortex.merged, reduction = "umap", group.by = 'scuttle_labels', label = TRUE)

# Feature Plot

features<- c("Mki67", "Rorb", "Gad2", "Foxp1",
             "Mpeg1", "Sox10", "Reln", "Bcl11b", "Cux1", "Emx1")
DefaultAssay(cortex.merged)<-"RNA"


plot1 <- UMAPPlot(cortex.merged, group.by="Genotype")
plot2 <- UMAPPlot(cortex.merged, label = T, label.box=T)
a<-plot1
b<-plot2
c<-FeaturePlot(cortex.merged, features, ncol=5, reduction="umap", cols=c("white", "red"))

(a|b) 

# Save object.
saveRDS(cortex.merged, "./cortex_merged_0403.rds")
