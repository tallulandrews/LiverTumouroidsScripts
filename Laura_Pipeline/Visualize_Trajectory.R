source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")

source("/nfs/users/nfs_t/ta6/R-Scripts/Blank_plot.R")
source("~/Collaborations/LiverOrganoids/Laura_Pipeline/99_Functions.R")


# N features, plotting obj
args <- commandArgs(trailingOnly=TRUE)
n_marks <- as.numeric(args[1]);

#file <- "CCA1_PlottingObj_Alt.rds"
file <- args[2];
prefix <- unlist(strsplit(file, "_"))[1]

# Combine Stuff
require("scater")
SCE <- readRDS(file)
if (class(SCE) == "SCESet") {
	SCE <- toSingleCellExperiment(SCE)
}
CC_normed_mat <- assays(SCE)[["norm_exprs"]]
palette <- cluster_col(max(SCE$Manual_Clusters))
cell_colours <- SCE@metadata$palette[SCE$Manual_Clusters]
SCE@metadata$C_names[SCE@metadata$C_names == "Hypoxic"] <- "Stress"
Clusters <- SCE@metadata$C_names[SCE$Manual_Clusters]
SCE$cell_col <- cell_colours
SCE$cell_name <- Clusters;

exclude <- is.na(Clusters);
SCE <- SCE[,!exclude]
Clusters <- Clusters[!exclude]


# Markers
require("CellTypeProfiles")
source("/nfs/users/nfs_t/ta6/My_R_packages/CellTypeProfiles/R/Markers.R")
calc_markers <- function(expr_mat, clusters){
		clusters <- as.character(clusters)
		exclude <- clusters=="Outliers" | clusters=="Stress"
		expr_mat <- expr_mat[,!exclude]
		clusters <- clusters[!exclude]
		expr_mat <- expr_mat[rowSums(expr_mat > 0) > 5,]
		out <- complex_markers(expr_mat, factor(clusters))
		return(out);
}

markers_normCC <- calc_markers(assays(SCE)[["norm_exprs"]], Clusters)
markers_normCC<- markers_normCC[markers_normCC$q.value < 0.05 & markers_normCC$q.value > 0,]
markers_normCC$Symbol <- rowData(SCE)$Symbol[rownames(SCE) %in% rownames(markers_normCC)]
markers_normCC <- markers_normCC[order(-markers_normCC$AUC),]

# Visualization
tsne_plot <- function(pca_rot, SCE, perplexity=30) {
	require("Rtsne")
	set.seed(134)
	TSNE <- Rtsne(pca_rot, dims=2, perplexity=perplexity, pca=FALSE)
	cycle_pch=rev(c(1, 12, 15, 17))

	plot(TSNE$Y[,1], TSNE$Y[,2], col=SCE$cell_col, pch=18, 
		xlab=paste("tSNE 1 (", ncol(pca_rot), ",", perplexity, ")"), 
		ylab="tSNE 2")
	return(list(x=TSNE$Y[,1], y=TSNE$Y[,2]))
}
pca_plot <- function(pca, SCE) {
	sdev <- round(pca_n$sdev/sum(pca_n$sdev)*100, digits=1)
	cycle_pch=rev(c(1, 12, 15, 17))
        plot(pca$rotation[,1], pca$rotation[,2], 
		col=SCE$cell_col,
                pch=16, 
		xlab=paste("PC1 (",sdev[1],"%)",sep=""), 
		ylab=paste("PC2 (",sdev[2],"%)",sep=""))
	return(list(x=pca$rotation[,1], y=pca$rotation[,2]))

}
dim_plot <- function(dim_mat, SCE, name="Dim") {
	cycle_pch=rev(c(1, 12, 15, 17))
        plot(dim_mat[,1], dim_mat[,2], 
		col=SCE$cell_col,
                pch=16, 
		xlab=paste(name, "1", sep=" "), 
		ylab=paste(name, "2", sep=""))

	return(list(x=dim_mat[,1], y=dim_mat[,2]))
}
plot_legend <- function() {
        par(mar=c(0,0,0,0))
        blank_plot()
        legend("left", 
		c("Cluster", SCE@metadata$C_names),
                col=c("white", SCE@metadata$palette),
                pch=16,
                bty="n")
}

png(paste(prefix, n_marks, "Visualizations.png", sep="_"), 
width=9, height=6, units="in", res=300)
par(mfrow=c(2,3))
par(mar=c(4,4,1,1));


outdims <- list();
# PCA
set.seed(134)
pca_n <- prcomp(assays(SCE)[["norm_exprs"]][rownames(SCE) %in% rownames(markers_normCC)[1:n_marks],])
outdims[["pca"]]<-pca_plot(pca_n, SCE) 

SCE$PC1_normCC <- pca_n$rotation[,1]
SCE$PC2_normCC <- pca_n$rotation[,2]

# tSNE
outdims[["tsne"]]<- tsne_plot(pca_n$rotation[,1:6], SCE) #cca1: dim=5, per=30 | hcc10 dim=8, per=30 | hcc6 dim=6, per=30 | cca5 dim=6 per=30 | hcc23 dim=4 per=30 | hcc24 dim=5 

# ICA
require("fastICA")
icas_normCC <- fastICA(t(assays(SCE)[["norm_exprs"]][rownames(SCE) %in% rownames(markers_normCC)[1:n_marks],]), n.comp = 2, method="C")
ica_dims <- icas_normCC$S;

outdims[["ica"]] <- dim_plot(ica_dims, SCE, "ICA")

# DM
require("destiny")
dm <- DiffusionMap(t(assays(SCE)[["norm_exprs"]][rownames(SCE) %in% rownames(markers_normCC)[1:n_marks],]));
outdims[["dm"]] <- dim_plot(dm@eigenvectors, SCE, "DM")

# UMAP
require("umap");
umap.out <- umap(t(assays(SCE)[["norm_exprs"]][rownames(SCE) %in% rownames(markers_normCC)[1:n_marks],]));
outdims[["umap"]] <- dim_plot(umap.out$layout, SCE, "UMAP")

plot_legend()
dev.off()

saveRDS(outdims, paste(prefix, n_marks, "Visualizations_dims.rds", sep="_")) 
