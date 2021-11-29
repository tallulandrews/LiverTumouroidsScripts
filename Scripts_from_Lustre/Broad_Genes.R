require("M3Drop")
require("scater")
require("matrixStats")
require("RColorBrewer")

set.seed(142)

CCA1 <- readRDS("CCA1_SC3_Prolif.rds")
HCC6 <- readRDS("HCC6_SC3_Prolif.rds")
HCC10 <- readRDS("HCC10_SC3_Prolif.rds")

exprs(CCA1) <- get_log_norm(CCA1)
exprs(HCC6) <- get_log_norm(HCC6)
exprs(HCC10) <- get_log_norm(HCC10)

anno_cols <- c("sc3_6_clusters", "CC_state", "Type")
CCA1_pData <- pData(CCA1)[,anno_cols]
HCC6_pData <- pData(HCC6)[,anno_cols]
HCC10_pData <- pData(HCC10)[,anno_cols]

Combined_counts <- cbind(counts(CCA1), counts(HCC6), counts(HCC10))
Combined_pData <- rbind(CCA1_pData, HCC6_pData, HCC10_pData)
Combined_fData <- fData(CCA1)[,c("Length", "feature_symbol")]

gd <- new("AnnotatedDataFrame", data=Combined_fData)
pd <- new("AnnotatedDataFrame", data=Combined_pData)
CombinedSCE <- newSCESet(countData = Combined_counts, phenoData = pd, featureData = gd)

CombinedSCE <- calculateQCMetrics(CombinedSCE)

markers <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Laura_SC3_k6_Markers.rds")

pData(CCA1)$cell_type1 <- paste("CCA1-",pData(CCA1)$sc3_6_clusters, sep="")
pData(HCC6)$cell_type1 <- paste("HCC6-",pData(HCC6)$sc3_6_clusters, sep="")
pData(HCC10)$cell_type1 <- paste("HCC10-",pData(HCC10)$sc3_6_clusters, sep="")

map = read.table("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/Hsap_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
         new = as.character(map[match(x, map[,1]),2])
         new[is.na(new)] = as.character(x[is.na(new)])
         new[duplicated(new)] = x[duplicated(new)]
         return(new)
}
fData(CCA1)$feature_symbol <- ensg2symbol(rownames(CCA1))
fData(HCC6)$feature_symbol <- ensg2symbol(rownames(HCC6))
fData(HCC10)$feature_symbol <- ensg2symbol(rownames(HCC10))

cca1_marks <- markers$CCA1$Gene[markers$CCA1$auroc > 0.8]
hcc6_marks <- markers$HCC6$Gene[markers$HCC6$auroc > 0.8]
hcc10_marks <- markers$HCC10$Gene[markers$HCC10$auroc > 0.8]

features <- unique(c( cca1_marks, hcc6_marks, hcc10_marks ))

require("scran")
source("/lustre/scratch117/cellgen/team218/TA/R-packages/new_scran/scran/R/mnnCorrect.R")
source("/lustre/scratch117/cellgen/team218/TA/R-packages/new_scran/scran/R/utils.R")
require("Matrix")
require("FNN")

laura_mnn_out <- mnnCorrect(list(CCA1=exprs(CCA1), HCC6=exprs(HCC6), HCC10=exprs(HCC10)), hvg.genes=which(rownames(exprs(CCA1)) %in% features))
Combined_corrected <- cbind(laura_mnn_out$corrected$CCA1, laura_mnn_out$corrected$HCC6, laura_mnn_out$corrected$HCC10)


Broad_genes <- c("UCK2", "SFPQ", "CLEC3B", "GTPBP4", "CAD", "RBP4", "G6PD", "CBX2", "TRMT6", "KPNA2", "MOSC2", "CEP55", "PPM1G", "RBM28", "CDC20", "SRL", "PSD4", "CDCA8", "YARS", "PPARGC1A", "APOH", "MAN2C1", "C7ORF68", "SOCS2", "SLC7A11", "GTF3C2", "BTNL9", "ADAMTS5", "YBX1", "HN1")
Broad_genes_ensg <- rownames(CCA1)[fData(CCA1)$feature_symbol %in% Broad_genes]

dat <- exprs(CombinedSCE)
rownames(dat) <- fData(CombinedSCE)$feature_symbol
png("Broad_Genes1.png", width=7, height=7, units="in", res=300)
M3DropExpressionHeatmap( Broad_genes, dat, CombinedSCE$Type )
dev.off()
png("Broad_Genes2.png", width=7, height=7, units="in", res=300)
M3DropExpressionHeatmap( Broad_genes, dat, CombinedSCE$CC_state )
dev.off()

png("Broad_Genes3.png", width=7, height=7, units="in", res=300)
M3DropExpressionHeatmap( Broad_genes_ensg, Combined_corrected, CombinedSCE$CC_state )
dev.off()

