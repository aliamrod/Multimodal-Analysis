# snATAC.R
rm(list=ls())
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(devtools)
library(EnsDb.Mmusculus.v79)
library(Ens)
library(ggplot2)
library(writexl)
library(patchwork)
library(BiocManager)
library(BiocParallel)
library(reshape2)
library(ComplexHeatmap)
library(ArchR)

# Adjust working directory.
setwd("/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/ANA/input_dir")
options(future.globals.maxSize = 8000*1024^2)
plan("multicore", workers = 10)
# plan()

# Convert MACS2 output table into GRanges object.
# A1
peaks.a1 = read.table('/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/MACS/A1_MACS/SAMPLE_A1_peaks.narrowPeak')
colnames(peaks.a1) = c('chr', 'start', 'end', 'name', 'score', 'strand', 'signal', 'pval',
                       'qval', 'peak')
peaks.a1 = makeGRangesFromDataFrame(peaks.a1, keep.extra.columns = TRUE)
peaks.a1 = keepStandardChromosomes(peaks.a1, pruning.mode = "coarse")

# A2
peaks.a2 = read.table('/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/MACS/A2_MACS/SAMPLE_A2_peaks.narrowPeak')
colnames(peaks.a2) = c('chr', 'start', 'end', 'name', 'score', 'strand', 'signal', 'pval',
                       'qval', 'peak')
peaks.a2 = makeGRangesFromDataFrame(peaks.a2, keep.extra.columns = TRUE)
peaks.a2 = keepStandardChromosomes(peaks.a2, pruning.mode = "coarse")

# A3
peaks.a3 = read.table('/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/MACS/A3_MACS/SAMPLE_A3_peaks.narrowPeak')
colnames(peaks.a3) = c('chr', 'start', 'end', 'name', 'score', 'strand', 'signal', 'pval',
                       'qval', 'peak')
peaks.a3 = makeGRangesFromDataFrame(peaks.a3, keep.extra.columns = TRUE)
peaks.a3 = keepStandardChromosomes(peaks.a3, pruning.mode = "coarse")


print("---Pre-processing data.---")

# CellRanger-A1
cr_counts_a1<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A1/outs/filtered_peak_bc_matrix.h5")
cr_chrassay_a1<- CreateChromatinAssay(counts = cr_counts_a1,
                                      sep = c(":", "-"),
                                      genome = "mm10", 
                                      fragments = "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A1/outs/fragments.tsv.gz")
cr_metadata_a1<- read.csv(
  file = "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A1/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

cr_seur_a1<- CreateSeuratObject(counts = cr_chrassay_a1,
                                assay = 'peaks', 
                                project = 'ATAC',
                                meta.data = cr_metadata_a1)


a1_fragpath<- "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A1/outs/fragments.tsv.gz"

#MACS2-A1
macs_counts_a1<- FeatureMatrix(
  fragments = Fragments(cr_seur_a1),
  features = peaks.a1,
  cells = colnames(cr_seur_a1)
)

# Extract gene annotations from EnsDb
annotations<- GetGRangesFromEnsDb(ensdb=EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations)<-'UCSC'
genome(annotations)<-"mm10"

# Create chromatin assay.
macs_chrassay_a1<- CreateChromatinAssay(counts=macs_counts_a1, genome="mm10",
                                        fragments=Fragments(cr_seur_a1),
                                        annotation=annotations)

# Seurat object (MACS2-based chromatin assay).
macs_seur_a1<- CreateSeuratObject(macs_chrassay_a1, assay="peaks",
                                  meta.data=cr_metadata_a1)

# Add the gene information to seur object.
Annotation(macs_seur_a1)<-annotations

##### ----- SAMPLE A2 ----- #####
# Counts
cr_counts_a2<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A2/outs/filtered_peak_bc_matrix.h5")
cr_chrassay_a2<- CreateChromatinAssay(counts = cr_counts_a2,
                                      sep = c(":", "-"),
                                      genome = "mm10", 
                                      fragments = "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A2/outs/fragments.tsv.gz")
cr_metadata_a2<- read.csv(
  file = "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A2/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

cr_seur_a2<- CreateSeuratObject(counts = cr_chrassay_a2,
                                assay = 'peaks', 
                                project = 'ATAC',
                                meta.data = cr_metadata_a2)


a2_fragpath<- "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A2/outs/fragments.tsv.gz"

macs_counts_a2<- FeatureMatrix(
  fragments = Fragments(cr_seur_a2),
  features = peaks.a2,
  cells = colnames(cr_seur_a2)
)

macs_chrassay_a2<- CreateChromatinAssay(counts=macs_counts_a2, genome="mm10",
                                        fragments=Fragments(cr_seur_a2),
                                        annotation=annotations)

macs_seur_a2<-CreateSeuratObject(macs_chrassay_a2, assay="peaks",
                                 meta.data=cr_metadata_a2)

Annotation(macs_seur_a2)<-annotations

##### ----- SAMPLE A3 ----- #####
# Counts
cr_counts_a3<- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A3/outs/filtered_peak_bc_matrix.h5")
cr_chrassay_a3<- CreateChromatinAssay(counts = cr_counts_a3,
                                      sep = c(":", "-"),
                                      genome = "mm10", 
                                      fragments = "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A3/outs/fragments.tsv.gz")
cr_metadata_a3<- read.csv(
  file = "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A3/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

cr_seur_a3<- CreateSeuratObject(counts = cr_chrassay_a3,
                                assay = 'peaks', 
                                project = 'ATAC',
                                meta.data = cr_metadata_a3)


a3_fragpath<- "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A3/outs/fragments.tsv.gz"

macs_counts_a3<- FeatureMatrix(
  fragments = Fragments(cr_seur_a3),
  features = peaks.a3,
  cells = colnames(cr_seur_a3)
)

macs_chrassay_a3<- CreateChromatinAssay(counts=macs_counts_a3, genome="mm10",
                                        fragments=Fragments(cr_seur_a3),
                                        annotation=annotations)

macs_seur_a3<-CreateSeuratObject(macs_chrassay_a3, assay="peaks",
                                 meta.data=cr_metadata_a3)

Annotation(macs_seur_a3)<-annotations

saveRDS(macs_seur_a1, "./macs_seur_a1_0402.rds")
saveRDS(macs_seur_a2, "./macs_seur_a2_0402.rds")
saveRDS(macs_seur_a3, "./macs_seur_a3_0402.rds")

#Reload objects.
macs_seur_a1<-readRDS("./macs_seur_a1_0402.rds")
macs_seur_a2<-readRDS("./macs_seur_a2_0402.rds")
macs_seur_a3<-readRDS("./macs_seur_a3_0402.rds")

##### ----- Doublet Inference ----- #####
# ArchR, v1.0.1//Dependencies.
library(ArchR)
set.seed(1)

addArchRThreads(threads=21)
addArchRGenome("mm10")

# Input data from 10X CellRanger-ATAC output.
a1.fragpath<- "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A1/outs/fragments.tsv.gz"
a2.fragpath<- "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A2/outs/fragments.tsv.gz"
a3.fragpath<- "/project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC/A3/outs/fragments.tsv.gz"

sampleName_A1=c("A1")
sampleName_A2=c("A2")
sampleName_A3=c("A3")

# DoubScores.
# A1
ArrowFiles<- createArrowFiles(
  inputFiles=a1.fragpath,
  sampleNames="A1",
  minTSS=4,
  minFrags=1000,
  addTileMat=TRUE,
  subThreading=FALSE,
  force=FALSE,
  addGeneScoreMat=TRUE)

addArchRThreads(threads=1) #adjust num. of threads to mitigate parallelity issue.
doubScores.A1<- addDoubletScores(
  input=ArrowFiles,
  k=10,
  knnMethod="UMAP",
  LSIMethod=1,
  force=TRUE
)

proj.A1<-ArchRProject(
  ArrowFiles=ArrowFiles,
  outputDirectory="/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC",
  copyArrows=TRUE
)

proj.A1<-filterDoublets(proj.A1)
saveArchRProject(ArchRProj = proj.A1, outputDirectory="/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/")

#Save and write dataframe into file.
proj.a1.df<- as.data.frame(proj.A1@cellColData)
write.table(proj.a1.df, '/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/input_dir/proj.a1.df.txt')

# A2
ArrowFiles.A2<- createArrowFiles(
  inputFiles=a2.fragpath,
  sampleNames="A2",
  minTSS=4,
  minFrags=1000,
  addTileMat=TRUE,
  subThreading=FALSE,
  force=FALSE,
  addGeneScoreMat=TRUE)

addArchRThreads(threads=1)
doubScores.A2<- addDoubletScores(
  input=ArrowFiles.A2,
  k=10,
  knnMethod="UMAP",
  LSIMethod=1,
  force=TRUE
)

proj.A2<- ArchRProject(
  ArrowFiles=ArrowFiles.A2,
  outputDirectory="/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC",
  copyArrows=TRUE
)

proj.A2<-filterDoublets(proj.A2)
saveArchRProject(ArchRProj=proj.A2, outputDirectory="/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/")

# Save and write dataframe into file (A2).
proj.a2.df<- as.data.frame(proj.A2@cellColData)
write.table(proj.a2.df, '/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/input_dir/proj.a2.df.txt')


# A3
ArrowFiles.A3<- createArrowFiles(
  inputFiles=a3.fragpath,
  sampleNames="A3",
  minTSS=4,
  minFrags=1000,
  addTileMat=TRUE,
  subThreading=FALSE,
  force=FALSE,
  addGeneScoreMat=TRUE,
)

addArchRThreads(threads=1)
doubScores.A3<- addDoubletScores(
  input=ArrowFiles.A3,
  k=10,
  knnMethod="UMAP",
  LSIMethod=1,
  force=TRUE
)

proj.A3<-ArchRProject(
  ArrowFiles=ArrowFiles.A3,
  outputDirectory="/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC",
  copyArrows=TRUE
)

proj.A3<-filterDoublets(proj.A3)
saveArchRProject(ArchRProj=proj.A3, outputDirectory="/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC")

# Save and write dataframe into file (A3).
proj.a3.df<- as.data.frame(proj.A3@cellColData)
write.table(proj.a3.df, '/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/input_dir/proj.a3.df.txt')


##### ----- Doublet Inference: Remove doubScored cellIDs ----- #####
library(dplyr)

# Read unfiltered ArchR project paths.
a1.filtered.path<-"/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/ANA/input_dir/proj.a1.df.txt" 
a2.filtered.path<-"/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/ANA/input_dir/proj.a2.df.txt"
a3.filtered.path<-"/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/ANA/input_dir/proj.a3.df.txt"

##A1
a1.filt.tab<- read.delim(file=a1.filtered.path, header=TRUE, sep=" ")

a1.filt.df1<- as.data.frame(a1.filt.tab)
a1.filt.df2<- a1.filt.df1[c(13)]

# Remove string from rows.
rownames(a1.filt.df2)<- gsub("A1#", "", as.character(rownames(a1.filt.df2)))

# Verify length of list of singlets.
length(rownames(a1.filt.df2)) #[1] 10432

##A2
a2.filt.tab<- read.delim(file=a2.filtered.path, header=TRUE, sep=" ")

a2.filt.df1<- as.data.frame(a2.filt.tab)
a2.filt.df2<- a2.filt.df1[c(13)]

# Remove string from rows.
rownames(a2.filt.df2)<- gsub("A2#", "", as.character(rownames(a2.filt.df2)))
length(rownames(a2.filt.df2)) #[1] 10187

##A3
a3.filt.tab<- read.delim(file=a3.filtered.path, header=TRUE, sep=" ")

a3.filt.df1<- as.data.frame(a3.filt.tab)
a3.filt.df2<- a3.filt.df1[c(13)]

# Remove string from rows.
rownames(a3.filt.df2)<- gsub("A3#", "", as.character(rownames(a3.filt.df2)))
length(rownames(a3.filt.df2)) #[1] 4013

# Plot unfiltered vs. filtered doublets. 
library(ggplot2)
doublet_plot <- read.table(
  header=TRUE, text='Category        SampleID NucleiCount
1   a_CellRanger      A1      11831
2   a_CellRanger      A2      11512
3   a_CellRanger      A3      4188
4   ArchR_Doublet_Inference    A1      10432
5   ArchR_Doublet_Inference                  A2       10187
6   ArchR_Doublet_Inference                A3       4013')

ggplot(doublet_plot, aes(SampleID, NucleiCount, fill = Category)) + 
  geom_bar(stat="identity", position = "dodge")+scale_fill_brewer(palette="Set1")+ggtitle("snATAC Filtration")+DarkTheme()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  geom_text(aes(label=NucleiCount), colour="white", vjust=1.5, position=position_dodge(.9))


# Create CellIDs as metadata column.
macs_seur_a1[["CellID"]]<- colnames(macs_seur_a1)
macs_seur_a2[["CellID"]]<- colnames(macs_seur_a2)
macs_seur_a3[["CellID"]]<- colnames(macs_seur_a3)


## List cells to retain.
CellsToRetain.A1<-rownames(a1.filt.df2)
CellsToRetain.A2<-rownames(a2.filt.df2)
CellsToRetain.A3<-rownames(a3.filt.df2)

## Retain list of cells.
#macs_seur_a1_filt_list<- colnames(macs_seur_a1)[colnames(macs_seur_a1) %in% toRetain]

##### ----- QC ----- #####
# Compute nucleosome signal score per cell.
macs_seur_a1<-NucleosomeSignal(object=macs_seur_a1)
macs_seur_a2<-NucleosomeSignal(object=macs_seur_a2)
macs_seur_a3<-NucleosomeSignal(object=macs_seur_a3)

# Compute TSS enrichment score per cell.
macs_seur_a1<- TSSEnrichment(object=macs_seur_a1, fast=FALSE)
macs_seur_a2<- TSSEnrichment(object=macs_seur_a2, fast=FALSE)
macs_seur_a3<- TSSEnrichment(object=macs_seur_a3, fast=FALSE)

# Add blacklist ratio and fraction of reads in peaks.
macs_seur_a1$pct_reads_in_peaks<- macs_seur_a1$peak_region_fragments / macs_seur_a1$passed_filters * 100
macs_seur_a2$pct_reads_in_peaks<- macs_seur_a2$peak_region_fragments / macs_seur_a2$passed_filters * 100
macs_seur_a3$pct_reads_in_peaks<- macs_seur_a3$peak_region_fragments / macs_seur_a3$passed_filters * 100

macs_seur_a1$blacklist_ratio<- macs_seur_a1$blacklist_region_fragments / macs_seur_a1$peak_region_fragments
macs_seur_a2$blacklist_ratio<- macs_seur_a2$blacklist_region_fragments / macs_seur_a2$peak_region_fragments
macs_seur_a3$blacklist_ratio<- macs_seur_a3$blacklist_region_fragments / macs_seur_a3$peak_region_fragments

macs_seur_a1$high.tss<- ifelse(macs_seur_a1$TSS.enrichment > 2, 'High', 'Low')
macs_seur_a2$high.tss<- ifelse(macs_seur_a2$TSS.enrichment > 2, 'High', 'Low')
macs_seur_a3$high.tss<- ifelse(macs_seur_a3$TSS.enrichment > 2, 'High', 'Low')

a<-TSSPlot(macs_seur_a1, group.by = 'high.tss')+NoLegend()+ggtitle("Sample_A1")
b<-TSSPlot(macs_seur_a2, group.by = 'high.tss')+NoLegend()+ggtitle("Sample_A2")
c<-TSSPlot(macs_seur_a3, group.by = 'high.tss')+NoLegend()+ggtitle("Sample_A3")

macs_seur_a1$nucleosome_group<- ifelse(macs_seur_a1$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
macs_seur_a2$nucleosome_group<- ifelse(macs_seur_a2$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
macs_seur_a3$nucleosome_group<- ifelse(macs_seur_a3$nucleosome_signal > 2, 'NS > 2', 'NS < 2')


#FragmentHistogram(object = macs_seur_a1, group.by = 'nucleosome_group')
VlnPlot(macs_seur_a2,features = c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
        pt.size = 0.1,ncol = 5)

# Set QC thresholds based on above and CellRanger summaries.
macs_seur_a1_subset <- subset(macs_seur_a1, subset = peak_region_fragments>1500 &
                                peak_region_fragments<50000 &
                                pct_reads_in_peaks<40 &
                                blacklist_ratio<0.05 &
                                nucleosome_signal<4 &
                                TSS.enrichment>2.5, 
                              cells=CellsToRetain.A1)

macs_seur_a2_subset<-subset(macs_seur_a2, #subset=peak_region_fragments>1500 &
                            peak_region_fragments<50000 &
                              pct_reads_in_peaks<40 &
                              blacklist_ratio<0.05 &
                              nucleosome_signal<4 &
                              TSS.enrichment>2, 
                            cells=CellsToRetain.A2)

macs_seur_a3_subset<-subset(macs_seur_a3, subset=peak_region_fragments>1500 &
                              peak_region_fragments<50000 &
                              pct_reads_in_peaks<40 &
                              blacklist_ratio<0.05 &
                              nucleosome_signal<4 &
                              TSS.enrichment>2.5, 
                            cells=CellsToRetain.A3)

# Normalization and Dim Reduction.
macs_seur_a1_subset<-RunTFIDF(macs_seur_a1_subset)
macs_seur_a1_subset<-FindTopFeatures(object=macs_seur_a1_subset)
macs_seur_a1_subset<-RunSVD(object=macs_seur_a1_subset)

macs_seur_a2_subset<-RunTFIDF(macs_seur_a2_subset)
macs_seur_a2_subset<-FindTopFeatures(object=macs_seur_a2_subset)
macs_seur_a2_subset<-RunSVD(object=macs_seur_a2_subset)

macs_seur_a3_subset<-RunTFIDF(macs_seur_a3_subset)
macs_seur_a3_subset<-FindTopFeatures(object=macs_seur_a3_subset)
macs_seur_a3_subset<-RunSVD(object=macs_seur_a3_subset)

# Non-linear dimension reduction and clustering.
macs_seur_a1_subset<-RunUMAP(
  object=macs_seur_a1_subset,
  reduction='lsi', 
  dims=2:20
)

macs_seur_a2_subset<-RunUMAP(
  object=macs_seur_a2_subset,
  reduction='lsi', 
  dims=2:20
)

macs_seur_a3_subset<-RunUMAP(
  object=macs_seur_a3_subset,
  reduction='lsi', 
  dims=2:20
)

# Add identifiers to integrated datasets.
macs_seur_a1_subset$Dataset<-"A1_CKO"
macs_seur_a2_subset$Dataset<-"A2_CKO"
macs_seur_a3_subset$Dataset<-"A3_FLOX"

macs_seur_a1_subset$Genotype<-"CKO"
macs_seur_a2_subset$Genotype<-"CKO"
macs_seur_a3_subset$Genotype<-"FLX"

saveRDS(macs_seur_a1_subset, "./macs_seur_a1_subset_0412.rds")
saveRDS(macs_seur_a2_subset, "./macs_seur_a2_subset_0412.rds")
saveRDS(macs_seur_a3_subset, "./macs_seur_a3_subset_0412.rds")

# Reload objects.
#macs_seur_a1_subset<-readRDS("./macs_seur_a1_subset_0412.rds")
#macs_seur_a2_subset<-readRDS("./macs_seur_a2_subset_0412.rds")
#macs_seur_a3_subset<-readRDS("./macs_seur_a3_subset_0412.rds")

# Integration
DefaultAssay(macs_seur_a1_subset)<-"peaks"
DefaultAssay(macs_seur_a2_subset)<-"peaks"
DefaultAssay(macs_seur_a3_subset)<-"peaks"

macs_seur_a1_subset <- RunTFIDF(macs_seur_a1_subset)
macs_seur_a2_subset <- RunTFIDF(macs_seur_a2_subset)
macs_seur_a3_subset <- RunTFIDF(macs_seur_a3_subset)

# Save objects.
saveRDS(macs_seur_a1_subset, "./macs_seur_a1_subset_0413.rds")
saveRDS(macs_seur_a2_subset, "./macs_seur_a2_subset_0413.rds")
saveRDS(macs_seur_a3_subset, "./macs_seur_a3_subset_0413.rds")

# Merge datasets pertaining to same biological condition/genotype, adding  cell ID to make sure cell names are unique
cko<- merge(
  x = macs_seur_a1_subset, 
  y = macs_seur_a2_subset,
  add.cell.ids = c("A1_CKO", "A2_CKO")
)

cko[["peaks"]]
cko<-RunTFIDF(cko)
cko<-FindTopFeatures(cko, min.cutoff = 20)
cko<-RunSVD(cko)
cko<-RunUMAP(cko, dims=2:50, reduction='lsi')

cko<-readRDS("./mer")
# Harmonize.
library(harmony)

cko.integrated <- RunHarmony(
  object = cko,
  group.by.vars = 'Dataset',
  reduction = 'lsi',
  assay.use = "peaks",
  project.dim = FALSE
)

# Recompute the UMAP using corrected LSI embeddings. 
cko.integrated <- RunUMAP(cko.integrated, dims = 2:50, reduction = 'harmony')

cko.umap <- DimPlot(cko, group.by = 'Dataset', pt.size = 0.1) + ggplot2::ggtitle("Non-Harmonized Set")
cko.integrated.umap <- DimPlot(cko.integrated, group.by = 'Dataset', pt.size = 0.1) + ggplot2::ggtitle("Harmony Integration")
cko.umap+cko.integrated.umap

saveRDS(cko.integrated, "./cko.integrated_0412.rds")

# Rename flox dataset.
flx <- macs_seur_a3_subset

# Reload objects.
setwd('/project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/ANA/input_dir/')
cko.integrated<-readRDS("./cko.integrated_0412.rds")
flx<-readRDS("./flx.rds")

#####----- Merge across sets.-----#####
merged<-merge(
  x = cko.integrated, 
  y = flx, 
  add.cell.ids = c("CKO", "FLX")
)

merged[['peaks']]
merged<-RunTFIDF(merged)
merged<-FindTopFeatures(merged, min.cutoff = 20)
merged<-RunSVD(merged)
merged<-RunUMAP(merged, dims = 2:50, reduction = "lsi")

merged.dimplot<-DimPlot(merged, group.by = "Dataset", pt.size=1)+ggplot2::ggtitle("Colored By Sample Before Corrective Integration")
saveRDS(merged, "./merged_0412.rds")
merged<-readRDS("./merged_0412.rds")

##### ----- Split Object. ----- #####
split<- SplitObject(object = merged, 
                    split.by = "Dataset")
saveRDS(split, "./split.rds")

##### ----- Integrate. ----- #####
integration.anchors<-FindIntegrationAnchors(
  object.list = split, 
  anchor.features = row.names(merged), 
  assay = c("peaks", "peaks", "peaks"),
  k.filter = NA, 
  reduction = "rlsi",
  dims = 2:30
)

saveRDS(integration.anchors, "./integration.anchors.merged_0413.rds")

# Integrate LSI embeddings.
integrated<-IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merged[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 2:30
)

saveRDS(integrated, './integrated_0413.rds')

# Create new UMAP using the integrated embeddings.
integrated<-RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated.dimplot<- DimPlot(integrated, group.by="Dataset", pt.size=1)+ggplot2::ggtitle("Colored By Sample After Corrective Integration")

a<-merged.dimplot+integrated.dimplot

##### ----- Gene activity ----- #####
seur_atac = integrated
gene.activities <- GeneActivity(seur_atac)
seur_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
seur_atac <- NormalizeData(seur_atac, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = median(seur_atac$nCount_RNA))

DefaultAssay(seur_atac) <- 'RNA'
g.list = c("Foxp1","Snap25","Bcl11b","Satb2","Gad1","Dlx2","Vim","Pax6","Hopx","Eomes","Top2a", "Rorb")
png('geneActivity.feature.plot.umap_0410.png', width = 1500, height=400*ceiling(length(g.list)/3), res=180, cols = c("white", "red")))
FeaturePlot(seur_atac,features = g.list, pt.size = .1, max.cutoff = 'q95', ncol = 3, cols=c("grey", "red"))
dev.off()

saveRDS(seur_atac, "seur_atac_0413.rds")
# Read in seur_atac object.
#seur_atac<-readRDS("./seur_atac_0413.rds")

##### ----- Integrate with RNA ----- #####
seur_rna_broad = readRDS("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir/labelled_query_02_0410.rds")
seur_rna_fine = readRDS("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir/labelled_query_01_0410.rds")
seur_atac_fine = readRDS("./seur_atac_0413.rds")
seur_atac_broad = seur_atac # Drive 2 separate objects--fine and broad labelled object.
# Load the pre-processed scRNA-seq data
seur_rna_fine <- FindVariableFeatures(
  object = seur_rna_fine,
  nfeatures = 5000
)

seur_rna_broad<-FindVariableFeatures(
  object=seur_rna_broad,
  nfeatures=5000
)

transfer.anchors.fine <- FindTransferAnchors(
  reference = seur_rna_fine,
  query = seur_atac_fine,
  reduction = 'cca',
  dims = 1:40
)

transfer.anchors.broad <- FindTransferAnchors(
  reference = seur_rna_broad,
  query = seur_atac_broad,
  reduction = 'cca',
  dims = 1:40
)

saveRDS(transfer.anchors.fine, "./transfer.anchors.fine_0413.rds")
saveRDS(transfer.anchors.broad, "./transfer.anchors.broad_0413.rds")

predicted.labels.fine <- TransferData(
  anchorset = transfer.anchors.fine,
  refdata = seur_rna_fine$manual_CellType,
  weight.reduction = seur_atac_fine[['integrated_lsi']],
  dims = 2:30
)

# Collapse subtypes into primary celltypes.
predicted.labels.broad<-TransferData(
  anchorset = transfer.anchors.broad, 
  refdata = seur_rna_broad$cellClass,
  weight.reduction = seur_atac_broad[['integrated_lsi']],
  dims = 2:30
)

seur_atac_fine<-AddMetaData(object = seur_atac_fine, metadata = predicted.labels.fine)
seur_atac_broad<-AddMetaData(object=seur_atac_broad, metadata=predicted.labels.broad)

p_00<-DimPlot(seur_rna_fine, group.by="manual_CellType",label.size=2, label.box=T, label=T, repel=T)+ggtitle('snRNA-seq')+NoLegend()
p_01<-DimPlot(seur_rna_broad, group.by="cellClass",reduction="umap", label=T,label.size=2,label.box=T, repel=T)+ggtitle('snRNA-seq')+NoLegend()
p1<-DimPlot(seur_atac_fine, pt.size=1,group.by = 'predicted.id',label = TRUE,repel = TRUE) + ggtitle('Predicted Cell Types (Fine)')+NoLegend()
p2<-DimPlot(seur_atac_broad,pt.size=1,group.by = 'predicted.id', label=TRUE, repel=TRUE)+ggtitle('Predicted Cell Types (Broad)')+NoLegend()

p3<-DimPlot(seur_atac, split.by = "Dataset", label.box=T, repel=T)
p4<-DimPlot(seur_atac, group.by = "Dataset", label.box=T, repel=T)
p3+p4

library(ggpubr)
ggarrange(p_00, p_01, p3, p4, p1, p2)



saveRDS(seur_atac_fine, "./seur_atac_labelled_fine_0413.rds")
saveRDS(seur_atac_broad, "./seur_atac_collapsed_labelled_broad_0413.rds")

# You can see that the RNA-based classifications are consistent with the UMAP visualization, 
# computed only on the ATAC-seq data. We can now easily annotate our snATAC-seq derived clusters (alternatively, we used the 
# RNA classifications themselves). Here, we transfer the cluster label (which shows finer distinctions) to annotate them.
# replace each label with its most likely prediction

# This step assigns each cluster to the most common predicted cell type. [1] indexes the vector returned by names(sort(.., 
# ...taking the first element (most common cell type)))
for(i in levels(seur_atac_broad)) {
  cells_to_reid <- WhichCells(seur_atac_broad, idents = i)
  newid <- names(sort(table(seur_atac_broad$predicted.id[cells_to_reid]),decreasing=TRUE))[1]
  Idents(seur_atac_broad, cells = cells_to_reid) <- newid
}

##### ----- Differentially accessible peaks ----- #####
seur_atac<-readRDS("./seur_atac_collapsed_labelled_broad_0413.rds")
DefaultAssay(seur_atac)<-'peaks'

# Subset excitatory population.
excitatory_subset<- subset(x=seur_atac, subset=predicted.id==as.character("Excitatory Neurons"))
DefaultAssay(excitatory_subset)<-'peaks'

sum(is.na(excitatory_subset$Genotype)) #[1] 2283
length(excitatory_subset$Genotype) # [1] 6665
excitatory_subset$Genotype[2283:6665]<-"FLOX"

da_peaks<-FindMarkers(
  object = excitatory_subset, 
  group.by = "Genotype",
  ident.1 = c("CKO"),
  ident.2 = c("FLOX"),
  min.pct=0.05,
  test.use='LR',
  latent.vars='peak_region_fragments'
)

head(da_peaks)

p11<-VlnPlot(
  object = exc_subset,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CKO", "FLX")
)

p12<- FeaturePlot(
  object = exc_subset, 
  features = rownames(da_peaks)[1],
  pt.size = 0.1, 
  max.cutoff = 'q95'
)

p11+p12
