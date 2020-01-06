# Edited 19 Feb PM 2018 to CC_state_new rather than CC_state
args <- commandArgs(trailingOnly=TRUE) # rds file for data, prefix for output

source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
outprefix=args[2];
expr_type <- "lognorm"
cluster_column <- "clusters_clean"
if (length(args) > 2) {
	expr_type <- args[3]
}
if (length(args) > 3) {
	cluster_column <- args[4]
}


set.seed(74468)

require("scater")
require("SingleCellExperiment")
require("matrixStats")

SCE <- readRDS(args[1])
cluster_col_set <- get_group_cols(SCE)
if (length(args) > 3) {
	cluster_col_set <- get_group_cols(SCE, column=cluster_column)
}
color_for_legend <- cluster_col_set[names(cluster_col_set) %in% colData(SCE)[,cluster_column]]

cluster_names_for_legend <- names(cluster_col_set)[names(cluster_col_set) %in% colData(SCE)[,cluster_column]];
cluster_names_for_legend[cluster_names_for_legend == "Outliers"] <- "O";

# Convert to SingleCellExperiment as necessary
if (class(SCE)[1] == "SCESet") {
	SCE <- toSingleCellExperiment(SCE)
}
# Use updated cell-cycle if available
if (sum(colnames(colData(SCE)) == "Cycle") == 1) {
	colData(SCE)$CC_state_new <- colData(SCE)$Cycle
}

print("loaded")

if ( !( "fine_marker_is.Feature" %in% colnames(rowData(SCE)) ) ) {
	rowData(SCE)$fine_marker_is.Feature <- abs(rowData(SCE)$fine_marker_q.value) < 0.05 & rowData(SCE)$fine_marker_AUC > 0.7
}

keep <- colData(SCE)[,cluster_column] != "Outliers"
#SCE_clean <- SCE[,keep]
SCE_clean <- SCE[,keep]
keep <- rowSums(assays(SCE_clean)[[expr_type]]) > 0 & rowVars(assays(SCE_clean)[[expr_type]]) > 0
SCE_clean <- SCE_clean[keep,]

print("filtered")

# KW genes - Feature Selection
KW_res <- apply(assays(SCE_clean)[[ expr_type ]], 1, function(x) {
		kruskal.test(x, colData(SCE_clean)[,cluster_column])$p.value})
KW_res[is.na(KW_res)] = 1
sig <- p.adjust(KW_res, method="bon") < 0.05
rowData(SCE_clean)$KW_p.value <- KW_res
rowData(SCE_clean)$KW_is.Feature <- sig

rowData(SCE)$KW_p.value <- rep(1, times=nrow(SCE))
rowData(SCE)[ match(names(KW_res), rownames(SCE)), "KW_p.value"] <- KW_res;
rowData(SCE)$KW_is.Feature <- p.adjust(rowData(SCE)$KW_p.value, method="bon") < 0.05

print("KW_DE Done")


## Diffusion Map ##
#SCE_FS <- SCE_clean[rowData(SCE_clean)$KW_is.Feature | rowData(SCE_clean)$fine_marker_is.Feature,]
SCE_FS <- SCE[rowData(SCE)$KW_is.Feature | rowData(SCE)$fine_marker_is.Feature,]

require("destiny")
dm <- DiffusionMap(t(assays(SCE_FS)[[expr_type]]));
dm_dims <- eigenvectors(dm)
colData(SCE)$KW_DM1 <- dm_dims[,1]
colData(SCE)$KW_DM2 <- dm_dims[,2]
colData(SCE)$KW_DM3 <- dm_dims[,3]

dpt <- DPT(dm)

pseudo_tips <- destiny::tips(dpt)
tip_data <- colData(SCE)[pseudo_tips,] 
pos_score <- abs(tip_data$KW_DM1)

start <- which(pos_score == min(pos_score[tip_data$CC_state_new %in% c("G1S", "G2M")]))
ends <- which(tip_data$CC_state_new %in% c("G0", "None"))
pseudo_time <- dpt$DPT1
branch <- dpt$Branch
if (length(start) > 0 & length(ends) > 0) {
if(mean(pseudo_time[pseudo_tips[start]]) > mean(pseudo_time[pseudo_tips[ends]])){
	pseudo_time <- max(pseudo_time)-pseudo_time
}
}
colData(SCE)$DPT_time <- pseudo_time
colData(SCE)$DPT_branch <- branch

print("diffusion map")

# Plots
require("scatterplot3d")
#scatterplot3d::scatterplot3d(colData(SCE)$KW_DM1, colData(SCE)$KW_DM2 , colData(SCE)$KW_DM3, 
#	color=cluster_col_set[ factor(colData(SCE)$clusters_fine, levels=names(cluster_col_set)) ],
#	pch=16, xlab="DM1", ylab="DM2", zlab="DM3")
#legend("topleft", pch=16, col=cluster_col_set, cluster_names_for_legend, bty="n", ncol=2)

png(paste(outprefix, "_DM_Pseudotime1.png", sep=""), width=6, height=6, units="in", res=300)
scatterplot3d::scatterplot3d(colData(SCE)$KW_DM1, colData(SCE)$KW_DM2 , colData(SCE)$KW_DM3, 
	color=cluster_col_set[ factor(colData(SCE)[,cluster_column], levels=names(cluster_col_set)) ], 
	pch=16, xlab="DM1", ylab="DM2", zlab="DM3")
legend("topleft", pch=16, col=color_for_legend, cluster_names_for_legend, bty="n", ncol=2)
dev.off()

#scatterplot3d::scatterplot3d(colData(SCE)$KW_DM1, colData(SCE)$KW_DM2 , colData(SCE)$KW_DM3, 
#	color=cluster_col_set[ factor(colData(SCE)$clusters_coarse, levels=names(cluster_col_set)) ], 
#	pch=16, xlab="DM1", ylab="DM2", zlab="DM3")
#legend("topleft", pch=16, col=cluster_col_set, cluster_names_for_legend, bty="n", ncol=2)

png(paste(outprefix, "_DM_Pseudotime2.png", sep=""), width=6, height=6, units="in", res=300)
scatterplot3d::scatterplot3d(colData(SCE)$KW_DM1, colData(SCE)$KW_DM2 , colData(SCE)$KW_DM3, 
	color=CC_col[colData(SCE)$CC_state_new], 
	pch=16, xlab="DM1", ylab="DM2", zlab="DM3")
legend("topleft", pch=16, col=CC_col, names(CC_col), bty="n")
dev.off()


## TSCAN ##
#require("TSCAN")
#
#dat <- counts(SCE_FS)
#procdat <- TSCAN::preprocess(t(t(dat)/colSums(dat)*1000000), minexpr_percent = 0.05, cvcutoff = 0.01)
#colnames(procdat) <- 1:ncol(dat)
#datclust <- TSCAN::exprmclust(procdat, clusternum = max(as.numeric(SCE_FS$clusters_coarse)):max(as.numeric(SCE_FS$clusters_fine)))
#png(paste(outprefix, "TSCANPseudo.png", sep="_"), width=6, height=6, units="in", res=300)
#TSCAN::plotmclust(datclust)
#dev.off()

## Monocle ##
require("monocle")
SCE_FS2 <- SCE_FS[,colData(SCE_FS)[,cluster_column] != "Outliers"]
pd <- as.data.frame(colData(SCE_FS2))
pd <- new("AnnotatedDataFrame", data=pd)
fd <- as.data.frame(rowData(SCE_FS2))
fd <- new("AnnotatedDataFrame", data=fd)
fd$gene_short_name <- fd$Symbol
dat <- counts(SCE_FS2)
rownames(fd) <- rownames(dat)
dCellData <- monocle::newCellDataSet(dat, phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())
dCellData <- monocle::setOrderingFilter(dCellData, rownames(dat))
#dCellData <- estimateSizeFactors(dCellData)
sizeFactors(dCellData) <- colSums(dat);
dCellData <- estimateDispersions(dCellData)
dCellData <- monocle::reduceDimension(dCellData, pseudo_expr = 1, reduction_method="DDRTree", norm_method="log", scaling=TRUE)
dCellData <- monocle::orderCells(dCellData, reverse = TRUE)
png(paste(outprefix, "_Monocle_Pseudotime1.png", sep=""), width=6, height=6, units="in", res=300)
monocle::plot_cell_trajectory(dCellData, color_by="CC_state_new",  show_branch_points=F)
dev.off()
png(paste(outprefix, "_Monocle_Pseudotime2.png", sep=""), width=6, height=6, units="in", res=300)
monocle::plot_cell_trajectory(dCellData, color_by=cluster_column, show_branch_points=F)
dev.off()

if (mean(dCellData$Pseudotime[dCellData$CC_state_new %in% c("G1S", "G2M")]) > mean(dCellData$Pseudotime[dCellData$CC_state_new %in% c("None", "G0")]) ) {
	dCellData$Pseudotime <- max(dCellData$Pseudotime) - dCellData$Pseudotime + min(dCellData$Pseudotime)
}

column_match <- match(colnames(dCellData), colnames(SCE))
row_match <- match(rownames(dCellData), rownames(SCE))

colData(SCE)$Monocle_time <- NA
colData(SCE)$Monocle_time[column_match] <- dCellData$Pseudotime
colData(SCE)$Monocle_branch <- NA
colData(SCE)$Monocle_branch[column_match] <- dCellData$State

monocle_Dims_cells <- dCellData@reducedDimS
monocle_Dims_genes <- dCellData@reducedDimW

colData(SCE)$Monocle_D1 <- NA
colData(SCE)$Monocle_D2 <- NA
rowData(SCE)$Monocle_D1 <- NA
rowData(SCE)$Monocle_D2 <- NA
colData(SCE)$Monocle_D1[column_match] <- monocle_Dims_cells[1,]
colData(SCE)$Monocle_D2[column_match] <- monocle_Dims_cells[2,]
rowData(SCE)$Monocle_D1[row_match] <- monocle_Dims_genes[,1]
rowData(SCE)$Monocle_D2[row_match] <- monocle_Dims_genes[,2]

# Big Monocle
pd <- as.data.frame(colData(SCE))
pd <- new("AnnotatedDataFrame", data=pd)
fd <- as.data.frame(rowData(SCE))
fd <- new("AnnotatedDataFrame", data=fd)
fd$gene_short_name <- fd$Symbol
dat <- counts(SCE)
rownames(fd) <- rownames(dat)
dCellData <- monocle::newCellDataSet(dat, phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())
dCellData <- monocle::setOrderingFilter(dCellData, rownames(SCE_FS))
sizeFactors(dCellData) <- colSums(dat);
dCellData <- estimateDispersions(dCellData)
dCellData <- monocle::reduceDimension(dCellData, pseudo_expr = 1, reduction_method="DDRTree", norm_method="log", scaling=TRUE)
dCellData <- monocle::orderCells(dCellData, reverse = TRUE)

column_match <- match(colnames(dCellData), colnames(SCE))
row_match <- match(rownames(dCellData), rownames(SCE))

colData(SCE)$Monocle_time_all <- NA
colData(SCE)$Monocle_time_all[column_match] <- dCellData$Pseudotime
colData(SCE)$Monocle_branch_all <- NA
colData(SCE)$Monocle_branch_all[column_match] <- dCellData$State

saveRDS(SCE, file=paste(outprefix, "_Pseudo.rds", sep=""))

# BEAM?

print("Monocle")


## mpath ##

#### Mpath crazyness ####
require("Mpath")

fixed_landmark_designation <- function(lognorm, groupID, distMethod = "euclidean", 
    method = "kmeans", numcluster = NULL, diversity_cut = 0.6, size_cut = 0.05, suppress.plot=TRUE) {
	require("igraph")
	require("vegan")
	if (distMethod == "euclidean") {
		hc <- hclust(dist(t(lognorm)), method = "ward.D")
	} else {
		data.cor <- cor(lognorm, method = distMethod)
		data.cor <- as.dist(1 - data.cor)
		hc <- hclust(data.cor, method = "ward.D")
	}
	nmi_res <- data.frame(cluster_num = c(1:ncol(lognorm)), 
        nmi = vector(length = ncol(lognorm)))
	for (k in 1:ncol(lognorm)) {
		ct <- cutree(hc, k = k)
		if (identical(names(ct), names(groupID))) {
			nmi_res[k, "nmi"] <- compare(as.factor(ct), as.factor(groupID), method = "nmi")
		} else {
			print("Error: names of groupID does not match column names of lognorm.\n")
		}
	}
	nmi_res_sort <- nmi_res[order(nmi_res$nmi, decreasing = TRUE),]
	optimal_cluster_num <- nmi_res_sort[1, "cluster_num"]
	optimal_nmi <- nmi_res_sort[1, "nmi"]
	if (is.null(numcluster)) {
		ct <- cutree(hc, k = optimal_cluster_num)
	} else {
		ct <- cutree(hc, k = numcluster)
	}
	clusters <- data.frame(table(ct))
	clusters$diversity <- diversity(t(table(groupID, ct)), index = "shannon")

	if (method == "kmeans") {
		km <- kmeans(clusters[, 2:3], 2)
		km_c <- km$cluster
		km_center <- km$center
		if (km_center[1, 1] > km_center[2, 1] & km_center[1, 2] > km_center[2, 2]) {
			km_c[km_c == 1] <- "Landmark clusters"
			km_c[km_c == 2] <- "Nonlandmark clusters"
		} else {
			km_c[km_c == 2] <- "Landmark clusters"
			km_c[km_c == 1] <- "Nonlandmark clusters"
		}
	} else if (method == "diversity") {
		km_c <- vector(length = nrow(clusters))
		km_c[clusters$diversity <= diversity_cut] <- "Landmark clusters"
		km_c[clusters$diversity > diversity_cut] <- "Nonlandmark clusters"
	} else if (method == "size") {
		km_c <- vector(length = nrow(clusters))
		km_c[clusters$Freq >= ceiling(size_cut * ncol(lognorm)) + 1] <- "Landmark clusters"
		km_c[clusters$Freq < ceiling(size_cut * ncol(lognorm)) + 1] <- "Nonlandmark clusters"
	} else if (method == "diversity_size") {
		km_c <- vector(length = nrow(clusters))
		km_c[clusters$diversity <= diversity_cut & clusters$Freq >= 
			ceiling(size_cut * ncol(lognorm)) + 1] <- "Landmark clusters"
		km_c[clusters$diversity > diversity_cut | clusters$Freq < 
			ceiling(size_cut * ncol(lognorm)) + 1] <- "Nonlandmark clusters"
	}

	cc_num <- as.character(clusters[km_c == "Landmark clusters", "ct"])
	ct_cc <- ct[ct %in% cc_num]
	cc <- data.frame(cell = names(ct_cc), cluster = ct_cc)
	cc$true_group <- groupID[as.character(cc$cell)]
	count <- table(cc$cluster, cc$true_group)
	num_name <- data.frame(num = vector(length = length(cc_num)), 
        name = vector(length = length(cc_num)))
	for (i in 1:nrow(count)) {
		num_name$num[i] <- rownames(count)[i]
		num_name$name[i] <- colnames(count)[order(count[i, ], decreasing = TRUE)[1]]
	}
	row.names(num_name) <- num_name[, "num"]
	cc$name <- num_name[as.character(cc$cluster), "name"]
	cc$landmark_cluster <- paste(cc$name, cc$cluster, sep = "_")
	cc_all <- data.frame(cell = names(ct), cluster = ct)
	cc_all$true_group <- groupID[as.character(cc_all$cell)]
	count_all <- table(cc_all$cluster, cc_all$true_group)
	num_name <- data.frame(num = vector(length = nrow(count_all)), 
        name = vector(length = nrow(count_all)))
	for (i in 1:nrow(count_all)) {
		num_name$num[i] <- rownames(count_all)[i]
		num_name$name[i] <- colnames(count_all)[order(count_all[i, ], decreasing = TRUE)[1]]
	}
	num_name$name <- paste(num_name$name, "_", num_name$num, sep = "")
	row.names(num_name) <- num_name$num
	clusters$name <- num_name[as.character(clusters$ct), "name"]
	col_palatte <- c("red", "black", height = 2, width = 2)

	if (!suppress.plot) {
		plot(clusters$Freq, clusters$diversity, xlab = "Cluster size", 
      	      ylab = "Cluster diversity", col = col_palatte[as.factor(km_c)], 
      	      pch = 16, xlim = c(min(clusters$Freq) - 1, max(clusters$Freq) + 
      	          1), ylim = c(min(clusters$diversity) - 0.1, max(clusters$diversity) + 
      	          0.1))
      	text(clusters$Freq, clusters$diversity - 0.02, labels = clusters$name, 
      	      cex = 1, adj = c(0.5, 0.9))
      	if (method == "diversity") {
      	      abline(h = diversity_cut, lty = 2)
      	} else if (method == "size") {
      	      abline(v = ceiling(size_cut * ncol(lognorm)) + 1, lty = 2)
      	} else if (method == "diversity_size") {
      	      abline(h = diversity_cut, lty = 2)
      	      abline(v = ceiling(size_cut * ncol(lognorm)) + 1, lty = 2)
      	}
		legend("topright", legend = unique(km_c), pch = 16, col = col_palatte[as.factor(unique(km_c))], 
      	      bty = "n")
	}
	return(cc[, c("cell", "landmark_cluster")])
}

make_mpath_network <- function(SCE_mpath, fs=rownames(counts(SCE_mpath)), distMethod="pearson") {
	require("igraph")
	norm <- 2^assays(SCE_mpath)[[ "lognorm"]]-1
	norm <- norm[rownames(norm) %in% fs,]
	cellID <- rownames(colData(SCE_mpath))
	groupID <- colData(SCE_mpath)$clusters_fine
	Ng <- length(levels(groupID))
	names(groupID) <- cellID

	#mpath_landmarks <- fixed_landmark_designation(norm, groupID, distMethod="pearson", method="diversity_size")
	faux_landmarks <- cbind(cellID, groupID); colnames(faux_landmarks) <- c("cell", "landmark_cluster")
	mpath_network <- build_network(norm, faux_landmarks, distMethod=distMethod)
#	mpath_network <- round(mpath_network/summary(factor(groupID)), digits=2)
	undir_net = mpath_network+t(mpath_network)
#	graph <- graph_from_adjacency_matrix(undir_net, mode="undirected", weighted=TRUE)
	graph <- graph_from_adjacency_matrix(mpath_network, mode="directed", weighted=NULL)
	set.seed(123)

	# Rescale edge weights to have reasonable line-thicknesses
	adj_mat <- cbind(rep(1:Ng, times=Ng), rep(1:Ng, each=Ng), as.vector(undir_net))
	rescale <- function(x) {
		ceiling(x/max(x)*5)
	}
	adj_mat <- cbind(adj_mat,rescale(adj_mat[,3]))


	my_layout <- layout_with_fr(graph);
	plot_edge <- function(adj_mat_row) {
		if(adj_mat_row[3] > 0) {
			xcoords = c(my_layout[adj_mat_row[1],1],my_layout[adj_mat_row[2],1]) 
			ycoords = c(my_layout[adj_mat_row[1],2],my_layout[adj_mat_row[2],2])
			lines(xcoords, ycoords, col="black", lwd=adj_mat_row[4])
			text(mean(xcoords), mean(ycoords), as.character(adj_mat_row[3]))
		}
	}

	par(xpd=TRUE)
	cluster_labs <- 1:max(as.numeric(colData(SCE_mpath)$clusters_fine))
	plot(my_layout, pch=16, col=cluster_col_set[cluster_labs], xaxt="n", yaxt="n", bty="n", cex=5, xlab="",ylab="")
	apply(adj_mat, 1, plot_edge)
	points(my_layout, pch=16, col=cluster_col_set[cluster_labs], cex=5)
	text(my_layout[,1], my_layout[,2], labels=as.character(cluster_labs), adj=c(0.5,0.5))
#	points(my_layout, pch=c("1","2","3","4","5","6"))
	invisible(list(network = mpath_network, plot_layout = my_layout))
}

#png(paste(outprefix, "_mpath_Pseudotime.png", sep=""), width=6, height=6, units="in", res=300)
#make_mpath_network(SCE_FS, distMethod="spearman")
#dev.off()

#print("Mpath")
