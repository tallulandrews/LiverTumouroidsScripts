args <- commandArgs(trailingOnly=TRUE) # SCE RDSs

source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/99_Functions.R")
source("/nfs/users/nfs_t/ta6/R-Scripts/Blank_plot.R")

CC_genes <- load_CC("all")
G2M_genes <- CC_genes$Tirosh[CC_genes$Tirosh[,2] =="G2M",1]
G1S_genes <- CC_genes$Tirosh[CC_genes$Tirosh[,2] =="G1S",1]

nSCEs <- length(args)
for (f in args) {
        require("scater")
	require("CellTypeProfiles")
        obj <- readRDS(f);

	cluster_expr <- my_row_mean_aggregate(obj, obj$clusters_clean)
	G2M_expr <- colMeans( cluster_expr[fData(obj)$Symbol %in% G2M_genes,] )
	G1S_expr <- colMeans( cluster_expr[fData(obj)$Symbol %in% G1S_genes,] )
	total_expr <- colMeans(cluster_expr)

	out <- factor_counts(obj$clusters_clean)
	out <- cbind(out, table(obj$clusters_clean, obj$CC_state_new))

	out <- as.data.frame(out)
	out$G2M_expr_fc <- G2M_expr/total_expr
	out$G1S_expr_fc <- G1S_expr/total_expr
	out$global_mean_expr <- total_expr

	out <- cbind(as.matrix(out), table(obj$clusters_clean, obj$clusters_noCC_clean))
	cnames <- colnames(out)
	cnames[1] <- "total cells"
	cnames[2:5] <- paste("NewCC", cnames[2:5], sep="_")
	cnames[9:length(cnames)] <- paste("noCC_cluster", cnames[9:length(cnames)], sep="_")
	colnames(out) <- cnames
	tmp <- unlist(strsplit(f, "_"))

	write.table(out, file=paste(tmp[1], "cluster_table.csv", sep="_"), sep=",")
}

