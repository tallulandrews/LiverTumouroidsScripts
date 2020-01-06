
### Marker Gene Heatmaps ###
STUFF<- readRDS(file="All_Plotting_stuff_Alt.rds")
nn_graph <- STUFF$nn_graphs
cell_colors <- STUFF$cell_colors

require("gplots")
require("scater")
require("RColorBrewer")
require("CellTypeProfiles")
clustered_rds <- c("CCA1_manual_SC3.rds", "CCA5_manual_SC3.rds", "HCC6_manual_SC3.rds", "HCC23_manual_SC3.rds", "HCC10_manual_SC3.rds", "HCC24_manual_SC3.rds", "D3DM_manual_SC3.rds", "D3EM_manual_SC3.rds", "D9DM_manual_SC3.rds", "D9EM_manual_SC3.rds");
lineIDs <- c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24", "D3DM", "D3EM", "D9DM", "D9EM");
for (i in 1:length(lineIDs)) {

	palette <- V(nn_graph[[i]]$graph)$col
	keep_cells <- cell_colors[[i]] %in% palette
	dat <- readRDS(clustered_rds[i])
	dat <- dat[,keep_cells]


	marker_file <- paste(lineIDs[i], "ManualClustering_MarkerTable.csv", sep="_");
	markers <- read.table(marker_file, sep=",")
	
	markers <- markers[!markers$Symbol == "" & markers$GeneType=="protein_coding" & markers$is.GoodMarker,]
	markers <- markers[order(markers$AUC, decreasing=T),]

	pval_col <- which(colnames(markers)=="p.values")
	cluster_cols <- 2:(pval_col-1)
	genes <- vector()
	for (j in cluster_cols) {
		top <- head(markers[markers[,j]==1,], 5)
		genes <- c(genes, as.character(top$Symbol));
	}
	genes<-unique(genes)
	heat_dat <- assays(dat)[["norm_exprs"]][rowData(dat)$feature_symbol %in% genes,]
	#heat_dat <- my_row_mean_aggregate(heat_dat, dat$Manual_Clusters)
	clabs <- dat$Manual_Clusters
	heat_dat <- heat_dat[order(genes),order(clabs)]
	genes <- genes[order(genes)]
	clabs <- clabs[order(clabs)]
	png(paste(lineIDs[i], "MarkerHeatmap_Alt.png", sep="_"), width=6, height=6, units="in", res=300)
	heatmap.2(heat_dat, col=rev(brewer.pal(8, "RdBu")), Colv=FALSE, trace="n", dendrogram="row", ColSideColors=palette[clabs])
	dev.off()
}
