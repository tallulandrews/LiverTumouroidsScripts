
args <- commandArgs(trailingOnly=TRUE) # SCE RDSs

source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("/nfs/users/nfs_t/ta6/R-Scripts/Blank_plot.R")

nSCEs <- length(args)
expr_type <- "lognorm"

type_col <- type_col[1:nSCEs]

SCE_list <- list();
keep_genes <- c()
consistent_genes <- c();
for (f in args) {
	require("scater")
	obj <- readRDS(f);
	keep_genes <- c(keep_genes, as.character(rownames(obj)[ fData(obj)$pct_dropout < 90 ]));
	if (length(consistent_genes) == 0) {
		consistent_genes <- rownames(obj)
	} else {
		consistent_genes <- consistent_genes[consistent_genes %in% rownames(obj)]
	}
	tmp <- unlist(strsplit(f, "\\."))
	SCE_list[[tmp[1]]] <- obj;
}
keep_genes <- sort(unique(keep_genes));
keep_genes <- keep_genes[keep_genes %in% consistent_genes]

require("CellTypeProfiles")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/My_R_packages/CellTypeProfiles/R/Markers.R")
require("proxy")

for (i in 1:nSCEs) {
	obj <- SCE_list[[i]]
	obj_name <- sub("_SC3", "", names(SCE_list[i]))
	obj <- obj[match(keep_genes, rownames(obj)),]

	profiles_fine <- my_row_mean_aggregate(get_exprs(obj, "lognorm"), obj$clusters_fine)
	colnames(profiles_fine) <- paste(obj_name, "Expr", colnames(profiles_fine), sep="_")
	marker_table <- fData(obj)[,grep("fine_marker", colnames(fData(obj)))]
	colnames(marker_table) <- sub("fine_marker", paste(obj_name, "marker", sep="_"), colnames(marker_table));
	colnames(marker_table) <- sub("Feature", "Good", colnames(marker_table));
	profiles_fine <- cbind(profiles_fine, marker_table)
	profiles_fine$Symbol <- fData(obj)$Symbol
	write.table(profiles_fine, file=paste("ForLaura",names(SCE_list)[i], "MegaMarkers_fine.csv", sep="_"), row.names=TRUE, col.names=TRUE, sep=",")


	profiles_coarse <- my_row_mean_aggregate(get_exprs(obj, "lognorm"), obj$clusters_coarse)
	colnames(profiles_coarse) <- paste(obj_name, "Expr", colnames(profiles_coarse), sep="_")
	marker_table <- fData(obj)[,grep("coarse_marker", colnames(fData(obj)))]
	colnames(marker_table) <- sub("coarse_marker", paste(obj_name, "marker", sep="_"), colnames(marker_table));
	colnames(marker_table) <- sub("Feature", "Good", colnames(marker_table));
	profiles_coarse <- cbind(profiles_coarse, marker_table)
	profiles_coarse$Symbol <- fData(obj)$Symbol
	write.table(profiles_coarse, file=paste("ForLaura",names(SCE_list)[i], "MegaMarkers_coarse.csv", sep="_"), row.names=TRUE, col.names=TRUE, sep=",")

	png(paste(obj_name, "Cluster_Legend.png", sep="_"), width=3, height=3, units="in", res=300)
	blank_plot()
	legend("left", levels(obj$clusters_fine), fill=cluster_col(max(as.numeric(obj$clusters_fine))), bty="n")
	legend("right", levels(obj$clusters_coarse), fill=cluster_col(max(as.numeric(obj$clusters_coarse))), bty="n")
	dev.off()
	
	write.table( table(obj$clusters_coarse, obj$clusters_fine), file=paste(obj_name, "Coarse_vs_Fine.txt", sep="_") )

	# Clean Markers (not already stored?)
	profiles_clean <- my_row_mean_aggregate(get_exprs(obj, "lognorm"), obj$clusters_clean)
	colnames(profiles_clean) <- paste(obj_name, "Expr", colnames(profiles_clean), sep="_")

	markers_clean <- complex_markers(get_exprs(obj, "lognorm"), pData(obj)$clusters_clean)
	colnames(markers_clean) <- paste("clean_marker", colnames(markers_clean), sep="_")
	identical(rownames(markers_clean), rownames(fData(obj)))
	fData(obj) <- cbind(fData(obj), markers_clean)
	saveRDS(obj, file=paste(obj_name, "Markers.rds", sep="_"))
	
	marker_table <- fData(obj)[,grep("clean_marker", colnames(fData(obj)))]
        colnames(marker_table) <- sub("clean_marker", paste(obj_name, "marker", sep="_"), colnames(marker_table));
        colnames(marker_table) <- sub("Feature", "Good", colnames(marker_table));
        profiles_clean <- cbind(profiles_clean, marker_table)
        profiles_clean$Symbol <- fData(obj)$Symbol
        write.table(profiles_clean, file=paste("ForLaura",names(SCE_list)[i], "MegaMarkers_clean.csv", sep="_"), row.names=TRUE, col.names=TRUE, sep=",")
}


# Cross References
cca1 <- read.delim("ForLaura_CCA1_SC3_MegaMarkers_clean.csv", sep=",")
hcc6 <- read.delim("ForLaura_HCC6_SC3_MegaMarkers_clean.csv", sep=",")
hcc10 <- read.delim("ForLaura_HCC10_SC3_MegaMarkers_clean.csv", sep=",")
d3dm <- read.delim("ForLaura_D3DM_SC3_MegaMarkers_clean.csv", sep=",")
d3em <- read.delim("ForLaura_D3EM_SC3_MegaMarkers_clean.csv", sep=",")
cca1 <- cca1[cca1$CCA1_marker_q.value < 0.05 & cca1$CCA1_marker_AUC > 0.7,]
hcc6 <- hcc6[hcc6$HCC6_marker_q.value < 0.05 & hcc6$HCC6_marker_AUC > 0.7,]
hcc10 <- hcc10[hcc10$HCC10_marker_q.value < 0.05 & hcc10$HCC10_marker_AUC > 0.7,]
d3dm <- d3dm[d3dm$D3DM_marker_q.value < 0.05 & d3dm$D3DM_marker_AUC > 0.7,]
d3em <- d3em[d3em$D3EM_marker_q.value < 0.05 & d3em$D3EM_marker_AUC > 0.7,]

cca1_a <- cca1[,c(18, 9:15)]
hcc6_a <- hcc6[,c(18, 9:15)]
hcc10_a <- hcc10[,c(18, 9:15)]
d3dm_a <- d3dm[,c(14, 7:11)]
d3em_a <- d3em[,c(16, 8:13)]

make_intersect_table <- function(assign1, assign2) {

	M <- matrix(0, ncol=ncol(assign2)-1, nrow=ncol(assign1)-1)
	for (i in 2:ncol(assign1)) {
		g1 <- assign1[assign1[,i]==1,1]
		for (j in 2:ncol(assign2)) {
			g2 <- assign2[assign2[,j]==1,1]
			M[i-1,j-1] <- length(intersect(g1,g2))/length(union(g1, g2));
		}
	}
	rownames(M) <- colnames(assign1)[2:ncol(assign1)]
	colnames(M) <- colnames(assign2)[2:ncol(assign2)]
	return(M)
}

a_list <- list(cca1_a, hcc6_a, hcc10_a, d3dm_a, d3em_a)

out_tab <- matrix(0, nrow=1, ncol=length(a_list))
colnames(out_tab) <- c("CCA1", "HCC6", "HCC10", "D3DM", "D3EM")
rownames_out <- c("Non")

for (a in 1:length(a_list)) {
out <- matrix(0, nrow=ncol(a_list[[a]])-1, ncol=length(a_list))
rownames_out <- c(rownames_out, colnames(a_list[[a]])[2:length(a_list[[a]])])
for (b in 1:length(a_list)) {

M <- make_intersect_table(a_list[[a]], a_list[[b]])
c2h <- apply(M, 1, function(x){colnames(M)[which(x == max(x))]})
out[,b] <- c2h

}
out_tab <- rbind(out_tab, out)
}
rownames(out_tab) <- rownames_out
out_tab <- out_tab[-1,]

check_recip <- function(c) {
	tmp <- unlist(strsplit(rownames(out_tab)[c], "_"))
	out_tab[out_tab[c,],tmp[1]] == rownames(out_tab)[c]
}

recip <- sapply(1:nrow(out_tab), check_recip)
colnames(recip) <- rownames(out_tab)
recip <- t(recip)


Reciprocal <- out_tab; Reciprocal[!recip] <- "-";
write.table(Reciprocal, file="Reciprocal_match.txt")
