require("RColorBrewer")

plate_pch = c(16, 7, 25, 3, 4, 23, 1)
type_col = brewer.pal(8, "Dark2")

tmp_cols <- brewer.pal(8, "Set1");
tmp_cols <- t(col2rgb(tmp_cols))
tmp_cols <- tmp_cols[order(tmp_cols[,1]/rowSums(tmp_cols)),]
tmp_cols <- apply(tmp_cols, 1, function(a) {rgb(a[1], a[2], a[3], maxColorValue=255)})
cluster_col = colorRampPalette(tmp_cols)
rm(tmp_cols)
prolif_col = c("black", "red")
CC_col = c("grey50", "black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("None", "G0", "G1/S", "G2/M");



require("Polychrome")
group_cols_vs <- Polychrome::alphabet.colors(26)
set.seed(1)
group_cols_vs <- sample(group_cols_vs, 20)
group_cols_vs <- c(group_cols_vs[-8], "grey50")
names(group_cols_vs) <- c(as.character(1:19), "Outliers")



get_group_cols <- function(SCE, column=NULL) {
	if (!is.null(column)) {
		return(group_cols_vs);
	}
	if (class(SCE) == "SingleCellExperiment") {
		require("SingleCellExperiment")
		group_cols <- c(cluster_col(max(as.numeric(colData(SCE)$clusters_fine))), "grey50")
		names(group_cols) <- c(as.character(1:max(as.numeric(colData(SCE)$clusters_fine))),"Outliers")
		return(group_cols);
	} else {
		require("scater")
		group_cols <- c(cluster_col(max(as.numeric(pData(SCE)$clusters_fine))), "grey50")
		names(group_cols) <- c(as.character(1:max(as.numeric(pData(SCE)$clusters_fine))),"Outliers")
		return(group_cols);
	}
}

