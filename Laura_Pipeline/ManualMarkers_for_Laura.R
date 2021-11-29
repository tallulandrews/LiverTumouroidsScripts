
#args <- commandArgs(trailingOnly=TRUE) # SCE RDSs

source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("/nfs/users/nfs_t/ta6/R-Scripts/Blank_plot.R")

#nSCEs <- length(args)
expr_type <- "lognorm"

type_col <- type_col[1:nSCEs]

#SCE_list <- list();
#keep_genes <- c()
#consistent_genes <- c();
#for (f in args) {
#	require("scater")
#	obj <- readRDS(f);
#	keep_genes <- c(keep_genes, as.character(rownames(obj)[ rowData(obj)$pct_dropout < 90 ]));
#	if (length(consistent_genes) == 0) {
#		consistent_genes <- rownames(obj)
#	} else {
#		consistent_genes <- consistent_genes[consistent_genes %in% rownames(obj)]
#	}
#	tmp <- unlist(strsplit(f, "\\."))
#	SCE_list[[tmp[1]]] <- obj;
#}
keep_genes <- sort(unique(keep_genes));
keep_genes <- keep_genes[keep_genes %in% consistent_genes]

require("CellTypeProfiles")
source("/nfs/users/nfs_t/ta6/My_R_packages/CellTypeProfiles/R/Markers.R")
require("proxy")

for (i in 1:nSCEs) {
	obj <- SCE_list[[i]]
	obj_name <- sub("_SC3", "", names(SCE_list[i]))
	obj <- obj[match(keep_genes, rownames(obj)),]

	profiles_manual <- my_row_mean_aggregate(assays(obj)[["lognorm"]], obj$Manual_Clusters)
	colnames(profiles_manual) <- paste(obj_name, "Expr", colnames(profiles_manual), sep="_")
	marker_table <- rowData(obj)[,grep("manual_marker", colnames(rowData(obj)))]
	colnames(marker_table) <- sub("manual_marker", paste(obj_name, "marker", sep="_"), colnames(marker_table));
	colnames(marker_table) <- sub("Feature", "Good", colnames(marker_table));
	profiles_manual <- cbind(profiles_manual, marker_table)
	profiles_manual$Symbol <- rowData(obj)$Symbol
	write.table(profiles_manual, file=paste("ForLaura",names(SCE_list)[i], "MegaMarkers_manual.csv", sep="_"), row.names=TRUE, col.names=TRUE, sep=",")


	profiles_coarse <- my_row_mean_aggregate(assays(obj)[["lognorm"]], obj$clusters_coarse)
	colnames(profiles_coarse) <- paste(obj_name, "Expr", colnames(profiles_coarse), sep="_")
	marker_table <- rowData(obj)[,grep("coarse_marker", colnames(rowData(obj)))]
	colnames(marker_table) <- sub("coarse_marker", paste(obj_name, "marker", sep="_"), colnames(marker_table));
	colnames(marker_table) <- sub("Feature", "Good", colnames(marker_table));
	profiles_coarse <- cbind(profiles_coarse, marker_table)
	profiles_coarse$Symbol <- rowData(obj)$Symbol
	write.table(profiles_coarse, file=paste("ForLaura",names(SCE_list)[i], "MegaMarkers_coarse.csv", sep="_"), row.names=TRUE, col.names=TRUE, sep=",")

	write.table( table(obj$clusters_coarse, obj$Manual_Clusters), file=paste(obj_name, "Coarse_vs_Manual.txt", sep="_") )

	# Clean Markers (not already stored?)
	profiles_clean <- my_row_mean_aggregate(assays(obj)[["lognorm"]], obj$clusters_clean)
	colnames(profiles_clean) <- paste(obj_name, "Expr", colnames(profiles_clean), sep="_")

	markers_clean <- complex_markers(assays(obj)[["lognorm"]], colData(obj)$clusters_clean)
	colnames(markers_clean) <- paste("clean_marker", colnames(markers_clean), sep="_")
	identical(rownames(markers_clean), rownames(rowData(obj)))
	rowData(obj) <- cbind(rowData(obj), markers_clean)
	saveRDS(obj, file=paste(obj_name, "Markers.rds", sep="_"))
	
	marker_table <- rowData(obj)[,grep("clean_marker", colnames(rowData(obj)))]
        colnames(marker_table) <- sub("clean_marker", paste(obj_name, "marker", sep="_"), colnames(marker_table));
        colnames(marker_table) <- sub("Feature", "Good", colnames(marker_table));
        profiles_clean <- cbind(profiles_clean, marker_table)
        profiles_clean$Symbol <- rowData(obj)$Symbol
        write.table(profiles_clean, file=paste("ForLaura",names(SCE_list)[i], "MegaMarkers_clean.csv", sep="_"), row.names=TRUE, col.names=TRUE, sep=",")
}


