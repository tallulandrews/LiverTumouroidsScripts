require("scater")

map = read.table("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/Hsap_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
        new = as.character(map[match(x, map[,1]),2])
        new[is.na(new)] = as.character(x[is.na(new)])
        new[duplicated(new)] = x[duplicated(new)]
        return(new)
}

set.seed(123)

sceset <- readRDS("QCed_SC_ScaterObj.rds")

sceset <- calculateQCMetrics(sceset)
CCA1 <- sceset[,pData(sceset)$Type=="CCA1"]
HCC6 <- sceset[,pData(sceset)$Type=="HCC6"]
HCC10 <- sceset[,pData(sceset)$Type=="HCC10"]

# SC3 - farm job 4 cores
require("SC3")

# CCA1 Clustering
CCA1 <- sc3_prepare(CCA1)
CCA1 <- sc3_estimate_k(CCA1)
CCA1 <- sc3(CCA1, ks=2:10, n_cores=4)
png("CCA1_k6_consensus.png")
sc3_plot_consensus(CCA1, k=6, show_pdata=c("Plate"))
dev.off()

#CCA1 <- sc3_summarise_results(CCA1,k=6)

png("CCA1_k6_PCA.png")
plotPCA(CCA1, colour_by = "sc3_6_clusters")
dev.off()

CCA1_k6_markers <- get_marker_genes(exprs(CCA1), pData(CCA1)$sc3_6_clusters)
CCA1_k6_markers$Gene = rownames(fData(CCA1))[as.numeric(rownames(CCA1_k6_markers))]
CCA1_k6_markers[is.na(CCA1_k6_markers$pvalue),]$pvalue <- 1
sig <- CCA1_k6_markers[CCA1_k6_markers$pvalue < 0.05/25507/6,]
sig[sig$clusts == 1,]$Gene
sig = sig[order(-sig$auroc),]

png("CCA1_k9_consensus.png")
sc3_plot_consensus(CCA1, k=9, show_pdata=c("Plate"))
dev.off()

#CCA1 <- sc3_summarise_results(CCA1,k=9)

png("CCA1_k9_PCA.png")
plotPCA(CCA1, colour_by = "sc3_9_clusters")
dev.off()

#CCA1_k9_markers <- get_marker_genes(exprs(CCA1), pData(CCA1)$sc3_9_clusters)

# HCC6 Clustering
HCC6 <- sc3_prepare(HCC6)
HCC6 <- sc3_estimate_k(HCC6)
HCC6@sc3$k_estimation
HCC6 <- sc3(HCC6, ks=2:(HCC6@sc3$k_estimation+1), n_cores=4)
png("HCC6_k6_consensus.png")
sc3_plot_consensus(HCC6, k=6, show_pdata=c("Plate"))
dev.off()

#HCC6 <- sc3_summarise_results(HCC6, k=11)

png("HCC6_k6_PCA.png")
plotPCA(HCC6, colour_by = "sc3_6_clusters")
dev.off()

HCC6_k6_markers <- get_marker_genes(exprs(HCC6), pData(HCC6)$sc3_6_clusters)
HCC6_k6_markers$Gene = rownames(fData(HCC6))[as.numeric(rownames(HCC6_k6_markers))]
HCC6_k6_markers[is.na(HCC6_k6_markers$pvalue),]$pvalue <- 1
HCC6_k6_markers$Symbol = ensg2symbol(HCC6_k6_markers$Gene)
sig <- HCC6_k6_markers[HCC6_k6_markers$pvalue < 0.05/25507/6,]
sig[sig$clusts == 1,]$Gene
sig = sig[order(-sig$auroc),]

# HCC10 Clustering
HCC10 <- sc3_prepare(HCC10)
HCC10 <- sc3_estimate_k(HCC10)
HCC10@sc3$k_estimation
HCC10 <- sc3(HCC10, ks=2:10, n_cores=4)

png("HCC10_k6_consensus.png")
sc3_plot_consensus(HCC10, k=6, show_pdata=c("Plate"))
dev.off()

png("HCC10_k6_PCA.png")
plotPCA(HCC10, colour_by = "sc3_6_clusters")
dev.off()

HCC10_k6_markers <- get_marker_genes(exprs(HCC10), pData(HCC10)$sc3_6_clusters)
HCC10_k6_markers$Gene = rownames(fData(HCC10))[as.numeric(rownames(HCC10_k6_markers))]
HCC10_k6_markers[is.na(HCC10_k6_markers$pvalue),]$pvalue <- 1
HCC10_k6_markers$Symbol = ensg2symbol(HCC10_k6_markers$Gene)
sig <- HCC10_k6_markers[HCC10_k6_markers$pvalue < 0.05/25507/6,]
sig[sig$clusts == 1,]$Gene
sig = sig[order(-sig$auroc),]

saveRDS(CCA1, "CCA1_SC3.rds")
saveRDS(HCC6, "HCC6_SC3.rds")
saveRDS(HCC10, "HCC10_SC3.rds")

Markers <- list(CCA1 = CCA1_k6_markers, HCC6 = HCC6_k6_markers, HCC10 = HCC10_k6_markers)
saveRDS(Markers, "Laura_SC3_k6_Markers.rds")
