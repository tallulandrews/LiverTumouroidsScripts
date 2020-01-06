require("destiny")
require("scater")
require("RColorBrewer")
source("~/Collaborations/LiverOrganoids/Laura_Scripts/Pascal_Colour_Scheme.R")

set.seed(32819)

CCA1 <- readRDS("CCA1_SC3_Prolif.rds")
HCC10 <- readRDS("HCC10_SC3_Prolif.rds")

CCA1 <- toSingleCellExperiment(CCA1)
HCC10 <- toSingleCellExperiment(HCC10)

get_log_norm <- function(SCE) {
	a <- counts(SCE);
	t <- colSums(a);
	n <- t(t(a)/t*median(t))
	return(log(n+1)/log(2))
}

exprs(CCA1) <- get_log_norm(CCA1)
exprs(HCC10) <- get_log_norm(HCC10)

do_fs <- function(SCE) {
	KW_res <- apply(exprs(SCE), 1, function(x) {kruskal.test(x, colData(SCE)$sc3_6_clusters)$p.value})
	KW_res[is.na(KW_res)] = 1
	sig <- p.adjust(KW_res, method="bon") < 0.05
	return(names(sig[sig]))
}

CCA1_fs <- do_fs(CCA1)
HCC10_fs <- do_fs(HCC10)

plot_DM <- function(SCE, fs=NULL) {
	require("matrixStats")
	keep = rowSums(exprs(SCE)) > 0 & rowVars(exprs(SCE)) > 0
	if (!is.null(fs)) {
		keep = keep & (rownames(exprs(SCE)) %in% fs);
	}
	norm <- exprs(SCE)[keep,]
	dm <- DiffusionMap(t(norm))
	dm_dims <- eigenvectors(dm)
	plot(dm_dims[,1], dm_dims[,2], col=c("black","red")[pData(SCE)$Proliferating+1], xlab="Dimension 1", ylab="Dimension 2", pch=16)
	plot(dm_dims[,1], dm_dims[,2], col=group_cols[pData(SCE)$sc3_6_clusters], xlab="Dimension 1", ylab="Dimension 2", pch=16)

}

require("matrixStats")
# CCA1 Diffusion Map

png("ForPascal_CCA1_Pseudotime.png", width=5, height=5, units="in", res=300)
plot_DM(CCA1, CCA1_fs)
keep <- rowSums(exprs(CCA1)) > 0 & rowVars(exprs(CCA1)) > 0 & rownames(exprs(CCA1)) %in% CCA1_fs
norm <- exprs(CCA1)[keep, ]
dm <- DiffusionMap(t(norm))
dm_dims <- eigenvectors(dm)
colData(CCA1)$DM1 <- dm_dims[,1]
colData(CCA1)$DM2 <- dm_dims[,2]
cluster_lab <- as.character(colData(CCA1)$sc3_6_clusters)
cluster_lab[!cluster_lab %in% names(CCA1_colours)] <- "Outliers"
cluster_lab <- factor(cluster_lab, levels=names(CCA1_colours));
plot(dm_dims[,1], dm_dims[,2], bg=CCA1_colours[cluster_lab], xlab="Dimension 1", ylab="Dimension 2", pch=21);
rect(min(dm_dims[,1])-2, min(dm_dims[,2])-2, max(dm_dims[,1])+2, max(dm_dims[,2])+2, col=grey_bg)
par(new=TRUE)
plot(dm_dims[,1], dm_dims[,2], bg=CCA1_colours[cluster_lab], xlab="Dimension 1", ylab="Dimension 2", pch=21);
par(new=FALSE)
dev.off()


png("ForPascal_HCC10_Pseudotime.png", width=5, height=5, units="in", res=300)
plot_DM(HCC10, HCC10_fs)
keep <- rowSums(exprs(HCC10)) > 0 & rowVars(exprs(HCC10)) > 0 & rownames(exprs(HCC10)) %in% HCC10_fs
norm <- exprs(HCC10)[keep, ]
dm <- DiffusionMap(t(norm))
dm_dims <- eigenvectors(dm)
colData(HCC10)$DM1 <- dm_dims[,1]
colData(HCC10)$DM2 <- dm_dims[,2]
cluster_lab <- as.character(colData(HCC10)$sc3_6_clusters)
cluster_lab[!cluster_lab %in% names(HCC10_colours)] <- "Outliers"
cluster_lab <- factor(cluster_lab, levels=names(HCC10_colours));
plot(dm_dims[,1], dm_dims[,2], bg=HCC10_colours[cluster_lab], xlab="Dimension 1", ylab="Dimension 2", pch=21);
rect(min(dm_dims[,1])-2, min(dm_dims[,2])-2, max(dm_dims[,1])+2, max(dm_dims[,2])+2, col=grey_bg)
par(new=TRUE)
plot(dm_dims[,1], dm_dims[,2], bg=HCC10_colours[cluster_lab], xlab="Dimension 1", ylab="Dimension 2", pch=21);
par(new=FALSE)
dev.off()

