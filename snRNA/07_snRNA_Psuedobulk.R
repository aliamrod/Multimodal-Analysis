# snRNA_pseudoBulk.R
library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyverse)
setwd("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir")
seur_rna= readRDS("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir/SEUR_RNA_LABELLED_0414.RDS")

a<-DimPlot(seur_rna, group.by="Genotype")
b<-DimPlot(seur_rna, group.by="cellClass_02", reduction="umap", label=T, label.size=2.5,label.box=T)+NoLegend()

# Pseudo-bulk workflow.
# Acquiring necessary metrics for aggregation across cells in a sample.
# 1. Counts matrix - sample level.
# Counts aggregate to sample level.


# Aggreate counts across cells at sample level.
View(seur_rna@meta.data)
seur_rna$sample_info<-paste0(seur_rna$Genotype, seur_rna$SampleID)
#seur_rna$SAMPLES<-rownames(seur_rna@meta.data)

split<-SplitObject(seur_rna, split.by="SampleID")

# Aggregate counts across the cells.
# AO1
cts_ao1<-AggregateExpression(split$AO1, 
                             group.by = c("cellClass", "SAMPLES"),
                             assays="RNA",
                             slot="counts",
                             return.seurat=FALSE,
                             verbose=FALSE
)

# AO2
cts_ao2<-AggregateExpression(split$AO2, 
                             group.by = c("cellClass", "SAMPLES"),
                             assays="RNA",
                             slot="counts",
                             return.seurat=FALSE,
                             verbose=FALSE
)

# AO3
cts_ao3<-AggregateExpression(split$AO3, 
                             group.by = c("cellClass", "SAMPLES"),
                             assays="RNA",
                             slot="counts",
                             return.seurat=FALSE,
                             verbose=FALSE
)
# AO4
cts_ao4<-AggregateExpression(split$AO4, 
                             group.by = c("cellClass", "SAMPLES"),
                             assays="RNA",
                             slot="counts",
                             return.seurat=FALSE,
                             verbose=FALSE
)
# AO5
cts_ao5<-AggregateExpression(split$AO5, 
                             group.by = c("cellClass", "SAMPLES"),
                             assays="RNA",
                             slot="counts",
                             return.seurat=FALSE,
                             verbose=FALSE
)

# AO6
cts_ao6<-AggregateExpression(split$AO6, 
                             group.by = c("cellClass", "SAMPLES"),
                             assays="RNA",
                             slot="counts",
                             return.seurat=FALSE,
                             verbose=FALSE
)

# AO7
cts_ao7<-AggregateExpression(split$AO7, 
                             group.by = c("cellClass", "SAMPLES"),
                             assays="RNA",
                             slot="counts",
                             return.seurat=FALSE,
                             verbose=FALSE
)

# AO8
cts_ao8<-AggregateExpression(split$AO8, 
                             group.by = c("cellClass", "SAMPLES"),
                             assays="RNA",
                             slot="counts",
                             return.seurat=FALSE,
                             verbose=FALSE
)

# Concatenate set of lists.
cts_aggregate<-c(cts_ao1,
                 cts_ao2,
                 cts_ao3,
                 cts_ao4,
                 cts_ao5,
                 cts_ao6,
                 cts_ao7)

cts<-cts_aggregate$RNA

# Transpose + convert to dataframe.
cts.t<-t(cts)

# Inspect dataframe.
cts.t[1:5, 1:5]


# Acquire values where to split.
splitRows<-gsub("_.*", "", rownames(cts.t))

meta = seur_rna@meta.data
meta2 = strsplit(meta$ID, '_') %>% do.call(rbind, .) %>% as.data.frame %>% head
meta_ids = meta2$V1

new.cts.t<-t(cts)
rownames(new.cts.t)<-gsub('.*_(.*_.*)_.*', '\\1', rownames(new.cts.t))

# Split dataframe.
cts.split<-split.data.frame(new.cts.t, 
                            f = factor(splitRows))
#gsub("(.*_.*)_.*", '\\1', string) ==> "FLOX_AO1
#gsub(".*_(.*)_.*", '\\1', string) ==> "AO1"

cts.split.modified<- lapply(cts.split, function(x){
  rownames(x)<-gsub('.*_(.*_.*)_.*', '\\1', rownames(x))
  t(x)
})

# Re-inspect list 'cts.split.modified'
#cts.split.modified$Astrocytes[1:3, 1:3]

# Get counts matrix.
counts_exc<- cts.split.modified$`Excitatory Neurons`

colData<-data.frame(samples=colnames(counts_exc))


##### ----- MAST ----- #####
# Load libraries
rm(list = ls())
library(tidyverse)
library(lme4)
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(MAST)
library(scRNAseq)
library(data.table)

#Subset clusters, save normalized counts and metadata for DEG analysis
seur_rna=readRDS("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir/SEUR_RNA_LABELLED_0414.RDS")

# Subset excitatory population.
exc<- subset(x=seur_rna, subset=cellClass_02==as.character("Excitatory Neurons"))
DefaultAssay(exc)<-'RNA'
exc<-NormalizeData(object = exc)
df<-as.data.frame(exc[["RNA"]]@data)
metadata<-as.data.frame(exc@meta.data)
meta=metadata

#meta<-meta[,c(1,2,3,4,5,6,7,8,9,10,11,21,25)];
df <- df[,match(rownames(metadata),colnames(df))];
data1<-as.matrix(df);
#save(data1,meta,file=paste("DGE/S_",test.cluster,".RData",sep=""));

#Make SingleCellAssay object to input in MAST
group1="CKO"
group2="FLOX"

test.tsne=meta
#unique_substrings<-gsub("(.*)_.*_.*", "\\1"  , rownames(test))
#rownames(test)<-unique_substrings
groups = as.character(test.tsne[,9])
axis<-factor(groups)
axis<-relevel(axis, group1)

test.data=data1[,rownames(test.tsne)] # select rows that have the same names as the rows in dataframe 'test.tsne' 
rows=rowSums(as.matrix(test.data)> 0) > ncol(test.data)*0.10 
test.data = test.data[rows,]
print(paste("Testing",nrow(test.data),"genes"))
sca=FromMatrix(as.matrix(test.data),as.data.frame(colnames(test.data)),as.data.frame(rownames(test.data)))

colData(sca)$cngeneson <- scale(cdr)
ind = as.character(test.tsne$SampleID)
cdr = colSums(assay(sca)>0)
batch=as.character(test.tsne$library.batch)
mito_perc=test.tsne$mitoPercent
colData(sca)$batch<-batch
colData(sca)$mito_perc<-scale(mito_perc)
colData(sca)$axis<-axis


#LMM using MAST to detect DEGs across the axis in a given cluster (Exc. cluster).
form=as.formula("~ axis +  cngeneson + mito_perc + batch  + (1|ind)") #(1|ind)
zlmCond = zlm(form, sca, method = "glmer", ebayes = F, silent=T)



#######
summaryCond <- summary(zlmCond, doLRT='axisCKO')
summaryDt <- summaryCond$datatable
save(summaryDt,file="DGE/LMM_C1.R")
fcHurdle <- merge(summaryDt[contrast=='axisCKO' & component=='H',.(primerid,`Pr(>Chisq)`)],summaryDt[contrast=='axisCKO' & component=='logFC', .(primerid, coef,ci.hi, ci.lo)], by='primerid');
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')];
fcHurdleSig <- merge(fcHurdle[fdr<0.05& abs(coef)>0.3], as.data.table(mcols(sca)),by='primerid');
setorder(fcHurdleSig, fdr);
gene.IDs<-as.character(fcHurdleSig$primerid);
FCs<-as.character(fcHurdleSig$coef);
FDRs<-as.character(fcHurdleSig$fdr);
out=cbind(gene.IDs,FCs,FDRs);
write.table(out,file="DGE/LMM_C1.txt",sep="\t")



##### ----- edgeR ----- #####
# snRNA_pseudoBulk.R
library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyverse)
setwd("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir")
seur_rna= readRDS("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir/SEUR_RNA_LABELLED_0414.RDS")

##### ----- edgeR ----- #####
# Load required libraries
library(Seurat)
library(edgeR)
library(scRNAseq)
library(scater)
library(scran)
library(Glimma)

# Subset excitatory population.
exc<- subset(x=seur_rna, subset=cellClass_02==as.character("Excitatory Neurons"))
DefaultAssay(exc)<-'RNA'
exc<-NormalizeData(object = exc)
df<-as.data.frame(exc[["RNA"]]@data)
metadata<-as.data.frame(exc@meta.data)
meta=metadata

# Set up the design matrix
design <- model.matrix(~Genotype, data = exc@meta.data)

# Extract the counts matrix
sce = exc
sce <- as.SingleCellExperiment(sce)
sce <- logNormCounts(sce)
var_mod <- modelGeneVar(sce)
hvg_genes <- getTopHVGs(var_mod, n=500)
hvg_sce <- sce[hvg_genes, ]
hvg_sce <- logNormCounts(hvg_sce)
#counts 
sce_counts<-counts(sce)

#counts <- assay(exc, "RNA")@data

# Convert counts to DGEList object
dge <- DGEList(counts = sce_counts, group = exc@meta.data$Genotype)

# Filter lowly expressed genes
keep <- rowSums(cpm(dge) >= 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Estimate library sizes and normalization factors
dge <- calcNormFactors(dge)

# Estimate dispersions and fit GLMs
dge <- estimateGLMCommonDisp(dge)
dge <- estimateGLMTagwiseDisp(dge)

# Perform differential expression analysis
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2) # coef = 2 specifies the comparison of interest (e.g. Mutant vs Wildtype)

# Filter differentially expressed genes by FDR and logFC
de_genes <- decideTestsDGE(lrt, p.value = 0.05, adjust.method = "BH")
de_genes <- rownames(dge)[de_genes$FDR < 0.05 & abs(lrt$coefficients[, 2]) > log2(1.5)]

# Add differential expression results to Seurat object
my_seurat@assays$RNA@data$edgeR_FDR <- NA
my_seurat@assays$RNA@data$edgeR_FDR[de_genes] <- de_genes$FDR[match(de_genes, rownames(de_genes$FDR))]
my_seurat@assays$RNA@data$edgeR_logFC <- NA
my_seurat@assays$RNA@data$edgeR_logFC[de_genes] <- lrt$coefficients[match(de_genes, rownames(lrt$coefficients)), 2]

# Save the Seurat object with differential expression results
saveRDS(my_seurat, "my_seurat_object_with_edgeR_results.rds")










####gsub("(.*)_.*_.*", "\\1"  ,'CKO_AO2_AAAGATGCACCTCGGA-1')
