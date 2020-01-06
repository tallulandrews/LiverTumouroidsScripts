source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("/nfs/users/nfs_t/ta6/R-Scripts/Blank_plot.R")
source("~/Collaborations/LiverOrganoids/Laura_Pipeline/99_Functions.R")

args <- commandArgs(trailingOnly=TRUE)
# testing args <- c("D3EM_noCC_SC3.rds", "D3EM_n_SC3.rds", "D3EM")
#noCC_SCE_file <- args[1]
#nCC_SCE_file <- args[2]
#prefix <- args[3]


# Combine Stuff
require("scater")
SCE <- readRDS("D3EM_noCC_SC3.rds")
SCE <- toSingleCellExperiment(SCE)
SCE_2 <- readRDS("D3EM_n_SC3.rds")

SCE$Clusters <- SCE$clusters_clean
SCE$Clusters_noCC <- SCE$clusters_noCC_clean
SCE$Clusters_normCC <- SCE_2$clusters_clean
SCE$Cycle <- SCE_2$CC_state
assays(SCE)[["norm_exprs"]] <- assays(SCE_2)[["norm_exprs"]]

SCE$Clusters_col <- group_cols_vs[match(SCE$Clusters, names(group_cols_vs))]
SCE$Clusters_noCC_col <- group_cols_vs[match(SCE$Clusters_noCC, names(group_cols_vs))]
SCE$Clusters_normCC_col <- group_cols_vs[match(SCE$Clusters_normCC, names(group_cols_vs))]
SCE$Cycle_col <- CC_col[SCE$CC_state]

# Make CC-removed expression matrix

CC_genes <- load_CC(set="cycling")
CC_genes <- c(as.character(CC_genes$Whitfield[,2]), as.character(CC_genes$Tirosh[,1]))

GO <- read.delim("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/GO/hsapiens_80_GO_Annoations_Emsembl.out", sep="\t", header=F)
GO_cc <- GO[GO[,3] == "cell cycle",1]


expr_mat <- assays(SCE)[["lognorm"]]
expr_mat[rownames(expr_mat) %in% GO_cc,] <- 0
expr_mat[rowData(SCE)$Symbol %in% CC_genes,] <- 0
assays(SCE)[["noCC"]] <- expr_mat

require("Polychrome")

# Olaps between clusters
table(SCE$Clusters, SCE$Clusters_noCC)
table(SCE$Clusters, SCE$Clusters_normCC)
table(SCE$Clusters_noCC, SCE$Clusters_normCC)

saveRDS(SCE, "D3EM_merged_SC3.rds");

# Markers
require("CellTypeProfiles")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/My_R_packages/CellTypeProfiles/R/Markers.R")
calc_markers <- function(expr_mat, clusters){
		clusters <- as.character(clusters)
		exclude <- clusters=="Outliers"
		expr_mat <- expr_mat[,!exclude]
		clusters <- clusters[!exclude]
		expr_mat <- expr_mat[rowSums(expr_mat > 0) > 5,]
		out <- complex_markers(expr_mat, factor(clusters))
		return(out);
}

markers_normCC <- calc_markers(assays(SCE)[["norm_exprs"]], SCE$Clusters_normCC)
markers_noCC <- calc_markers(assays(SCE)[["noCC"]], SCE$Clusters_noCC)
markers <- calc_markers(assays(SCE)[["lognorm"]], SCE$Clusters)

markers<- markers[markers$q.value < 0.05 & markers$q.value > 0,]
markers_noCC<- markers_noCC[markers_noCC$q.value < 0.05 & markers_noCC$q.value > 0,]
markers_normCC<- markers_normCC[markers_normCC$q.value < 0.05 & markers_normCC$q.value > 0,]

markers_normCC$Symbol <- rowData(SCE)$Symbol[rownames(SCE) %in% rownames(markers_normCC)]
markers_noCC$Symbol <- rowData(SCE)$Symbol[rownames(SCE) %in% rownames(markers_noCC)]
markers$Symbol <- rowData(SCE)$Symbol[rownames(SCE) %in% rownames(markers)]

markers <- markers[order(-markers$AUC),]
markers_normCC <- markers_normCC[order(-markers_normCC$AUC),]
markers_noCC <- markers_noCC[order(-markers_noCC$AUC),]

# Marker Tables
# Some Relevant Annotations:
schwalie_StemCell <- c("MFAP2", "HSPB6", "GBP3", "MYC", "CCND2", "RAC2", "ZFP36L2", "VGLL4", "ADGRG1", "IDH2", "CDK6", "SH3BGRL", "IFITM3", "LCP1", "ETV6", "RPGRIP1", "ERI3", "SDCBP", "DAP", "JTB", "GSTO1", "DAPP1", "CTSZ")
GO_telo <- unique(GO[GO[,3] == "telomere maintenance",1])
LigRec <- read.delim("~/Data/LigandReceptorPairs.csv", sep=",", header=T)
GO_tf <- unique(GO[grep("transcription factor activity", GO[,3]),1])


expr_table <- my_row_mean_aggregate(assays(SCE)[["norm_exprs"]], SCE$Clusters_normCC)
markers_normCC <- cbind(markers_normCC, expr_table[match( rownames(markers_normCC), rownames(expr_table)),])

expr_table <- my_row_mean_aggregate(assays(SCE)[["noCC"]], SCE$Clusters_noCC)
markers_noCC <- cbind(markers_noCC, expr_table[match( rownames(markers_noCC), rownames(expr_table)),])

expr_table <- my_row_mean_aggregate(assays(SCE)[["lognorm"]], SCE$Clusters)
markers <- cbind(markers, expr_table[match( rownames(markers), rownames(expr_table)),])

add_anno <- function(marker_tab) {
	marker_tab$is.cycle <- rownames(marker_tab) %in% GO_cc
	marker_tab$is.telo <- rownames(marker_tab) %in% GO_telo
	marker_tab$is.tf <- rownames(marker_tab) %in% GO_tf
	marker_tab$ReceptorOf <- LigRec[match(marker_tab$Symbol,LigRec$Receptor.ApprovedSymbol),"Ligand.ApprovedSymbol"]
	marker_tab$LigandTo <- LigRec[match(marker_tab$Symbol,LigRec$Ligand.ApprovedSymbol),"Receptor.ApprovedSymbol"]
	return(marker_tab)
}

markers_normCC <- add_anno(markers_normCC)
markers_noCC <- add_anno(markers_noCC)
markers <- add_anno(markers)

write.table(markers_normCC, "D3EM_n_Markers.csv", sep=",")
write.table(markers_noCC, "D3EM_noCC_Markers.csv", sep=",")
write.table(markers, "D3EM_Markers.csv", sep=",")


background <- rownames(SCE)[rowData(SCE)$pct_dropout < 95]
get_richments <- function(group, this_markers, rm.cc = FALSE) {
	require("gProfileR")
	require("proxy")

	group_col <- which(colnames(this_markers) == group)[1]
	assign_cols <- 2:(which(colnames(this_markers) == "p.values") -1)
	specificity <- rowSums(this_markers[,assign_cols])

	good <- specificity <=2 & this_markers$AUC > 0.75 & this_markers[,group_col]==1
	if (rm.cc) {
		good <- good & this_markers[,"is.cycle"] == FALSE
	}
	gene_list <- rownames(this_markers)[good]
	

        enrichments <- gprofiler(gene_list, organism="hsapiens",
          ordered_query=T, significant=T, custom_bg=background,
          hier_filtering="moderate", max_set_size=10000,
          src_filter=c("GO:BP", "KEGG", "REAC", "HPA"),
          correction_method="fdr", min_isect_size=3, min_set_size=10)
        enrichments<-enrichments[order(enrichments$p.value),]
        # enrichment filtering
        # remove hpa low
        exclude <- enrichments$domain == "hpa" & grepl("Low", enrichments$term.name)
        exclude <- exclude | enrichments$domain == "hpa" & grepl("Not detected", enrichments$term.name)
        exclude <- exclude | enrichments$domain == "hpa" & grepl("Uncertain", enrichments$term.name)
        enrichments <- enrichments[!exclude,]
#        enrichments <- enrichments[1:min(nrow(enrichments), n_top_rich),]
        enrichments$GroupID = rep(group, times=nrow(enrichments))
        return(enrichments)
}



# GO enrichments
#require("gProfileR")
#thing <- rowSums(markers_normCC[,2:6])

# Visualization
tsne_plot <- function(pca_rot, clusters, perplexity=30) {
	require("Rtsne")
	set.seed(134)
	TSNE <- Rtsne(pca_rot, dims=ncol(pca_rot), perplexity=perplexity, pca=FALSE)
	cycle_pch=rev(c(1, 12, 15, 17))

	layout(matrix(c(1,2), nrow=1), widths=c(8,1))
	par(mar=c(4,4,1,1))
	plot(TSNE$Y[,1], TSNE$Y[,2], col=group_cols_vs[match(clusters, names(group_cols_vs))],
		pch=cycle_pch[SCE$Cycle], xlab=paste("tSNE 1 (", ncol(pca_rot), ",", perplexity, ")"), 
		ylab="tSNE 2")
	par(mar=c(0,0,0,0))
	blank_plot()
	legend("left", c("Cluster", levels(clusters), "", "Cycle", levels(SCE$Cycle)), 
		col=c("white", group_cols_vs[match(levels(clusters), names(group_cols_vs))], 
			"white", "white", rep("black", times=length(levels(SCE$Cycle)))),
		pch=c(16, rep(16, times=length(levels(clusters))),16, 16, cycle_pch),
		bty="n")
	return(TSNE$Y)
}
pca_plot <- function(pca, clusters) {
	sdev <- round(pca_n$sdev/sum(pca_n$sdev)*100, digits=1)
	layout(matrix(c(1,2), nrow=1), widths=c(8,1))
	cycle_pch=rev(c(1, 12, 15, 17))
        par(mar=c(4,4,1,1))
        plot(pca$rotation[,1], pca$rotation[,2], 
		col=group_cols_vs[match(clusters, names(group_cols_vs))],
                pch=cycle_pch[SCE$Cycle], xlab=paste("PC1 (",sdev[1],"%)",sep=""), 
		ylab=paste("PC2 (",sdev[2],"%)",sep=""))
        par(mar=c(0,0,0,0))
        blank_plot()
        legend("left", c("Cluster", levels(clusters), "", "Cycle", levels(SCE$Cycle)),
                col=c("white", group_cols_vs[match(levels(clusters), names(group_cols_vs))],
                        "white", "white", rep("black", times=length(levels(SCE$Cycle)))),
                pch=c(16, rep(16, times=length(levels(clusters))),16, 16, cycle_pch),
                bty="n")

}




# PCA
set.seed(134)
pca_n <- prcomp(assays(SCE)[["norm_exprs"]][rownames(SCE) %in% rownames(markers_normCC)[1:1000],])
pca_no <- prcomp(assays(SCE)[["noCC"]][rownames(SCE) %in% rownames(markers_noCC)[1:1000],])
pca_base <- prcomp(assays(SCE)[["lognorm"]][rownames(SCE) %in% rownames(markers)[1:1000],])

#plot(pca$rotation[,1], pca$rotation[,2], col=SCE$Clusters_normCC_col, pch=16)
#plot(pca$rotation[,1], pca$rotation[,2], col=SCE$Cycle_col, pch=16)

png("~/Tmp_MovingFiles/D3EM_n_Cluster_PCA.png", width=9*0.5, height=8*0.5, units="in", res=300)
pca_plot(pca_n, SCE$Clusters_normCC) 
dev.off()

png("~/Tmp_MovingFiles/D3EM_no_Cluster_PCA.png", width=9*0.5, height=8*0.5, units="in", res=300)
pca_plot(pca_no, SCE$Clusters_noCC) 
dev.off()

png("~/Tmp_MovingFiles/D3EM_base_Cluster_PCA.png", width=9*0.5, height=8*0.5, units="in", res=300)
pca_plot(pca_base, SCE$Clusters) 
dev.off()

SCE$PC1_base <- pca_base$rotation[,1]
SCE$PC2_base <- pca_base$rotation[,2]
SCE$PC1_noCC <- pca_no$rotation[,1]
SCE$PC2_noCC <- pca_no$rotation[,2]
SCE$PC1_normCC <- pca_n$rotation[,1]
SCE$PC2_normCC <- pca_n$rotation[,2]

# tSNE
png("~/Tmp_MovingFiles/D3EM_n_Cluster_tSNE.png", width=9*0.5, height=8*0.5, units="in", res=300)
coords <- tsne_plot(pca_n$rotation[,1:6], SCE$Clusters_normCC) #cca1: dim=5, per=30 | hcc10 dim=8, per=30 | hcc6 dim=6, per=30 | cca5 dim=6 per=30 | hcc23 dim=4 per=30 | hcc24 dim=5 
dev.off()
SCE$tsne_x_normCC <- coords[,1]
SCE$tsne_y_normCC <- coords[,2]

png("~/Tmp_MovingFiles/D3EM_no_Cluster_tSNE.png", width=9*0.5, height=8*0.5, units="in", res=300)
coords<-tsne_plot(pca_no$rotation[,1:2], SCE$Clusters_noCC) #cca1: dim=8, per=30 | hcc10: dim=5, per=30 | hcc6: dim=6, per=30 | cca5 dim=2 per=30 | hcc23 dim=5 | hcc24 dim=3
dev.off()
SCE$tsne_x_noCC <- coords[,1]
SCE$tsne_y_noCC <- coords[,2]

png("~/Tmp_MovingFiles/D3EM_base_Cluster_tSNE.png", width=9*0.5, height=8*0.5, units="in", res=300)
coords <- tsne_plot(pca_base$rotation[,1:6], SCE$Clusters) #cca1 : dim=8, per=30 | hcc10:dim=7, per=30 | hcc6 dim=6, per=30 | cca5 dim=6 per=30 | hcc23 dim=5, per=30 | hcc24 dim=6
dev.off()
SCE$tsne_x_base <- coords[,1]
SCE$tsne_y_base <- coords[,2]

saveRDS(SCE, "D3EM_merged_SC3.rds");

