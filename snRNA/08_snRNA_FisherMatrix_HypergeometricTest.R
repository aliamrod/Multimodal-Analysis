# Fisher Exact Test
library(Seurat)
library(scCustomize)
library(ggplot2)
library(dplyr)
library(MAST)

setwd('/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir/')

reference.path<-"/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/input_dir/GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5"
reference<-readRDS("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/input_dir/arlotta_ref_0404.rds")

query<-readRDS("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir/query_0404.rds")

# Calculate DEGs per cluster. Use this list to assign a cell type to each cluster.
query_markers<-FindAllMarkers(object=query, only.pos=T, test.use="MAST", logfc.threshold=0.25)
query_top5_markers<-Extract_Top_Markers(marker_dataframe=query_markers, num_genes=3, named_vector=FALSE,
                                        make_unique=TRUE)

reference_markers<-FindAllMarkers(object=reference, only.pos=T, test.use="MAST", logfc.threshold = 0.25)

write.table(query_markers, file="./query_markers_0410.txt")
write.table(reference_markers, file="./reference_markers_0406.txt")

# Query Dot Plot
Clustered_DotPlot(seurat_object=query, features=query_top5_markers, k=9)

#---------------------
# Fisher Exact Test.
setwd('/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/ANA/input_dir/')
reference.h5.path<- ('/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/input_dir/GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5')
reference<- readRDS("./ARLOTTA_REF.rds")
query<- readRDS("./label_transfer_ana_0404.rds")

# Calculate DEGs per cluster. Use this list to assign a cell type to each cluster.
query_markers<-FindAllMarkers(object=query, only.pos=T, test.use="MAST", logfc.threshold = 0.25)
reference_markers<-FindAllMarkers(object=reference, only.pos=T, test.use="MAST", logfc.threshold = 0.25)

write.table(query_markers, file="./query_markers.txt")
write.table(reference_markers, file="./reference_markers.txt")
#---------------------
# Fisher Exact Test.
# Load tables.
ref.path<-"./reference_markers_0404.txt"
query.path<- "./query_markers_0406.txt"

# Read cluster markers.
degtable.query<-query.path
tab<-read.table(degtable.query, sep= "", header=TRUE)
tab<- tab[tab$p_val_adj<=0.05 & tab$pct.1>=0.5, ]
#pct.1=> the percentage of cells where the feature is detected in the first group.
#pct.2=> the percentage of cells where the feature is detected in the second group
tab<- tab[c(7,6)]
tab$cluster<-as.factor(paste("Cluster_", sprintf("%02d", tab$cluster), sep=""))
colnames(tab)<-c("Gene", "DEFINITION")
tab$DEFINITION<-as.factor(tab$DEFINITION)
Genes<-as.data.frame(table(tab$DEFINITION))

# Read reference cell-type markers.
degtable.reference<-ref.path
ref_markers<-read.table(degtable.reference, sep=" ", header=TRUE)
ref_markers1<-ref_markers[ref_markers$p_val_adj<=0.05, ]
ref_markers2<- ref_markers1[ref_markers1$pct.1>=0.5,]
ref_markers3<-ref_markers2[c(7,6)]


ref_meta<-as.data.frame(reference@meta.data)
ref_meta_sel<- ref_meta[c("seurat_clusters", "newCellType")]
row.names(ref_meta_sel) <- NULL
ref_meta_sel2 <- ref_meta_sel[!duplicated(ref_meta_sel$newCellType),]

colnames(ref_meta_sel2)<- c("Cluster", "CellType")

new_ref<- merge(x = ref_markers, y = ref_meta_sel2,
                by.x = "cluster",
                by.y = "Cluster")
new_ref_mod <- new_ref[c(7, 8)]

new_ref_mod2 <- new_ref_mod[!is.na(new_ref_mod),]

new_ref_mod2<-na.omit(new_ref_mod2)
vcGenes <- list(Cajal_Retzius_Cells = new_ref_mod2[new_ref_mod2$CellType == "Cajal Retzius cells",],
                Neurons = new_ref_mod2[new_ref_mod2$CellType == "Neurons",],
                Astrocytes = new_ref_mod2[new_ref_mod2$CellType == "Astrocytes",],
                Pericytes = new_ref_mod2[new_ref_mod2$CellType == "Pericytes",],
                Microglia = new_ref_mod2[new_ref_mod2$CellType == "Microglia",],
                Corticofugal_Projection_Neurons = new_ref_mod2[new_ref_mod2$CellType == "Corticofugal Projection Neurons (CFuPNs)",],
                Inhibitory_Interneurons = new_ref_mod2[new_ref_mod2$CellType == "Interneurons",],
                Intermediate_Progenitors = new_ref_mod2[new_ref_mod2$CellType == "Intermediate progenitors",],
                Cycling_Glial_Cells = new_ref_mod2[new_ref_mod2$CellType == "Cycling glial cells",],
                Migrating_Neurons = new_ref_mod2[new_ref_mod2$CellType == "Migrating neurons", ],
                Layer_6b = new_ref_mod2[new_ref_mod2$CellType == "Layer 6b", ],
                CThPN = new_ref_mod2[new_ref_mod2$CellType == "CThPN", ],
                ULCPN=new_ref_mod2[new_ref_mod2$CellType == "UL CPN", ],
                Layer_4 = new_ref_mod2[new_ref_mod2$CellType == "Layer 4", ],
                Oligodendrocytes=new_ref_mod2[new_ref_mod2$CellType == "Oligodendrocytes", ],
                SCPN = new_ref_mod2[new_ref_mod2$CellType == "SCPN", ],
                Ependymocytes = new_ref_mod2[new_ref_mod2$CellType == "Ependymocytes", ],
                Neural_Progenitors = new_ref_mod2[new_ref_mod2$CellType == "NP", ],
                DL_CPN=new_ref_mod2[new_ref_mod2$CellType == "DL CPN", ],
                Endothelial_Cells = new_ref_mod2[new_ref_mod2$CellType == "Endothelial cells",],
                Doublet = new_ref_mod2[new_ref_mod2$CellType =="Doublet", ])
GeneSets <- vcGenes
######################### Re-arrange data
for(i in 1:length(GeneSets))
{
  colnames(GeneSets[[i]])[1] <- "Gene"
}
ln <- length(GeneSets)
cl <- length(Genes$Var1)
TEMP <- list()
INT <- list()
for (i in 1:ln)
{
  TEMP[[i]] <- tab[tab$Gene %in% GeneSets[[i]]$Gene,]
  INT[[i]] <- as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT) <- names(GeneSets)
names(TEMP) <- names(GeneSets)
######################### Create the matrix for each GeneSet
NROWS <- sapply(GeneSets, nrow)
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+
for (i in 1:length(INT))
{
  INT[[i]]$b <- NROWS[[i]] - INT[[i]]$Freq
  INT[[i]]$c <- Genes$Freq - INT[[i]]$Freq
  INT[[i]]$d <- 13500 - Genes$Freq - nrow(GeneSets[[i]])
}
######################### Function for Fisher's exact test
RunFisher <- function(row, alt = 'greater', cnf = 0.85)
{
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value),
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}
######################### Run Fisher's exact test
FisherMat=list()
for (i in 1:length(INT))
{
  FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
  rownames(FisherMat[[i]]) <- INT[[i]]$Var1
  FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat) <- names(INT)
save(FisherMat, TEMP, file = "/work/Neuroinformatics_Core/s204365/snRNA_scripts/fishers_mat.RData")
######################### Arrange a matrix for P-val
tmp <- list()
FisherP <- matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
  tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
  FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames
######################### Arrange a matrix for OR
tmp <- list()
FisherOR <- matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
  tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
  FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames
######################### Compute adjusted P-val
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames
FisherAdj[FisherAdj > 0.05] <- 1
FisherOR[FisherOR < 1] <- 0
df <- -log10(FisherAdj)
######################### Plot
library(WGCNA)
pdf("/project/Neuroinformatics_Core/Konopka_lab/s204365/RNA/output_dir/heatmap_fisher.pdf",width=16, height=16, pointsize=15)
par(mar = c(11, 7, 2, 2))
LabelMat <- paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df,
               xLabels = colnames(df),
               yLabels = rownames(df),
               colorLabels =FALSE,
               colors=colorRampPalette(c("white", "red"))(50),
               textMatrix=LabelMat,
               setStdMargins = FALSE,
               cex.text = 0.5,
               xLabelsAngle = 90)
dev.off()
######################### Plot with different color scale
low0 <- rgb(246, 248, 251, maxColorValue = 255)
low1 <- rgb(201, 230, 234, maxColorValue = 255)
low2 <- rgb(89, 158, 193, maxColorValue = 255)
mid <- rgb(67, 121, 180, maxColorValue = 255)
high1 <- rgb(65, 90, 158, maxColorValue = 255)
high2 <- rgb(52, 52, 106, maxColorValue = 255)
high3 <- rgb(15, 15, 25, maxColorValue = 255)
pdf("/project/Neuroinformatics_Core/Konopka_lab/s204365/2heatmap.pdf",width=12, height=16, pointsize=15)
par(mar = c(11, 7, 2, 2))
labeledHeatmap(Matrix = df,
               xLabels = colnames(df),
               yLabels = rownames(df),
               colorLabels =FALSE,
               colors=colorRampPalette(c(low0, low1, low2, mid, high1, high2, high3))(100),
               setStdMargins = FALSE,
               cex.text = 0.5,
               xLabelsAngle = 90)
dev.off()
