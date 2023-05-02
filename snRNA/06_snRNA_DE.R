library(dplyr)
library(Seurat)
library(ggplot2)
library(MAST)


DE_Method='MAST'
top_DE <- 40
resolution <- 0.8

##### Initialize objects #####
setwd("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir")
cortex <- readRDS("./cortex_merged_0403.rds")

###### Store Initial cluster Plot #####################
pdf('/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/output_dir',height=8,width=8)
DimPlot(cortex, reduction = "umap",label = T,repel=T,pt.size = 0.2)
print(p)
dev.off()

DefaultAssay(cortex) <- 'RNA'
cortex.markers <- FindAllMarkers(cortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = DE_Method)
saveRDS(cortex.markers, "./snRNA_DE_markers.rds")

top40 <- cortex.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)
top20 <- cortex.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.table(top40,file='./snRNA_top40_DEmarkers.tsv',quote=F,col.names=T,row.names=F,sep='\t')
DefaultAssay(cortex) <- 'RNA'
pdf('./DE_top20.pdf',height=40,width=20)

DoHeatmap(cortex, features = top20$gene,raster = T) + NoLegend()
dev.off()

pdf('./DE_top40.pdf',height=40,width=20)
DoHeatmap(cortex, features = top40$gene,raster = T) + NoLegend()
dev.off()

pdf('./CellCycle_VlnPLot.pdf',height=10,width=20)
print(VlnPlot(cortex, c("S.Score","G2M.Score"),pt.size = 0.001))
dev.off()        

DimPlot(cortex, reduction = "umap",label = T,repel=T,pt.size = 0.2) + NoLegend()
dev.off()
saveRDS(cortex,'data/merged_scRNA_unfiltered_IDs.RDS') 

DefaultAssay(cortex) <- 'RNA'

cortex.markers <- FindAllMarkers(cortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = DE_Method)
saveRDS(cortex.markers,'data/scRNA_DE_markers_IDs.RDS')

top40 <- cortex.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)
top20 <- cortex.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.table(top40,file='results/scRNA_top40_DEmarkers_IDs.tsv',quote=F,col.names=T,row.names=F,sep='\t')
DefaultAssay(cortex) <- 'RNA'
pdf('plots/scRNA/DE_top20.pdf',height=40,width=20)
DoHeatmap(cortex, features = top20$gene,raster = T) + NoLegend()
dev.off()
pdf('plots/scRNA/DE_top40.pdf',height=40,width=20)
DoHeatmap(cortex, features = top40$gene,raster = T) + NoLegend()
dev.off()


saveRDS(cortex_filtered,'data/merged_scRNA_filtered_IDs.RDS')  
