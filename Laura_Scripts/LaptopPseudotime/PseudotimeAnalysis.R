require("destiny")
require("scater")
require("RColorBrewer")

group_cols <- brewer.pal(n=6, "Set2")

CCA1 <- readRDS("CCA1_SC3_Prolif.rds")
HCC6 <- readRDS("HCC6_SC3_Prolif.rds")
HCC10 <- readRDS("HCC10_SC3_Prolif.rds")

get_log_norm <- function(SCE) {
	a <- counts(SCE);
	t <- colSums(a);
	n <- t(t(a)/t*median(t))
	return(log(n+1)/log(2))
}

exprs(CCA1) <- get_log_norm(CCA1)
exprs(HCC6) <- get_log_norm(HCC6)
exprs(HCC10) <- get_log_norm(HCC10)

do_fs <- function(SCE) {
	KW_res <- apply(exprs(SCE), 1, function(x) {kruskal.test(x, pData(SCE)$sc3_6_clusters)$p.value})
	KW_res[is.na(KW_res)] = 1
	sig <- p.adjust(KW_res, method="bon") < 0.05
	return(names(sig[sig]))
}

CCA1_fs <- do_fs(CCA1)
HCC6_fs <- do_fs(HCC6)
HCC10_fs <- do_fs(HCC10)

plot_DM <- function(SCE, fs=NULL) {
	require("matrixStats")
	keep = rowSums(exprs(SCE)) > 0 & rowVars(exprs(SCE)) > 0
	if (!is.null(fs)) {
		keep = keep & (rownames(exprs(SCE)) %in% fs);
	}
	norm <- exprs(SCE)[keep,]
#	fs <- M3Drop_Differential_Expression(2^norm -1)
#	dm <- DiffusionMap(t(norm[rownames(norm) %in% fs$Gene,]))
	dm <- DiffusionMap(t(norm))
	dm_dims <- eigenvectors(dm)
	plot(dm_dims[,1], dm_dims[,2], col=c("black","red")[pData(SCE)$Proliferating+1], xlab="Dimension 1", ylab="Dimension 2", pch=16)
	plot(dm_dims[,1], dm_dims[,2], col=group_cols[pData(SCE)$sc3_6_clusters], xlab="Dimension 1", ylab="Dimension 2", pch=16)

}

png("CCA1_DM_Pseudotime.png", width=6*2, height=6, units="in", res=300)
par(mfrow=c(1,2))
plot_DM(CCA1, CCA1_fs)
dev.off()
png("HCC6_DM_Pseudotime.png", width=6*2, height=6, units="in", res=300)
par(mfrow=c(1,2))
plot_DM(HCC6, HCC6_fs)
dev.off()
png("HCC10_DM_Pseudotime.png", width=6*2, height=6, units="in", res=300)
par(mfrow=c(1,2))
plot_DM(HCC10, HCC10_fs)
dev.off()

do_TSCAN <- function(SCE, fs=NULL) {
	require("TSCAN")
	dat <- counts(SCE)
	if (!is.null(fs)) {
		dat <- dat[rownames(dat) %in% fs,]
	}
	procdat <- TSCAN::preprocess(dat)
	colnames(procdat) <- 1:ncol(dat)
	datclust <- TSCAN::exprmclust(procdat, clusternum = 6)
	TSCAN::plotmclust(datclust)
}

do_monocle <- function(SCE, fs, name="Data") {
	require("monocle")
	pd <- pData(SCE)
	pd$sc3_6_clusters <- factor(pd$sc3_6_clusters)
	pd <- new("AnnotatedDataFrame", data=pd)
	fd <- fData(SCE)
	fd <- new("AnnotatedDataFrame", data=fd)

	dCellData <- monocle::newCellDataSet(counts(SCE), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())
	dCellData <- monocle::setOrderingFilter(dCellData, fs)
	dCellData <- estimateSizeFactors(dCellData)
	dCellData <- estimateDispersions(dCellData)
	dCellData <- monocle::reduceDimension(dCellData, pseudo_expr = 1)
	dCellData <- monocle::orderCells(dCellData, reverse = TRUE)
	png(paste(name, "_Monocle_Pseudotime1.png", sep=""), width=6, height=6, units="in", res=300)
	monocle::plot_cell_trajectory(dCellData, color_by="Proliferating")
	dev.off()
	png(paste(name, "_Monocle_Pseudotime2.png", sep=""), width=6, height=6, units="in", res=300)
	monocle::plot_cell_trajectory(dCellData, color_by="sc3_6_clusters")
	dev.off()
}

do_monocle(CCA1, CCA1_fs, name="CCA1")
do_monocle(HCC6, HCC6_fs, name="HCC6")
do_monocle(HCC10, HCC10_fs, name="HCC10")

#### Mpath crazyness ####

fixed_landmark_designation <- function(norm, groupID, distMethod = "euclidean", 
    method = "kmeans", numcluster = NULL, diversity_cut = 0.6, size_cut = 0.05, suppress.plot=TRUE) {
	require("igraph")
	require("vegan")
	lognorm <- log(norm+1)/log(2);
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
		if (sum(names(ct) != names(groupID)) == 0) {
			nmi_res[k, "nmi"] <- compare(as.factor(ct), as.factor(groupID), method = "nmi")
		} else {
			print("Error: names of groupID does not match column names of norm.\n")
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
		km_c[clusters$Freq >= ceiling(size_cut * ncol(norm)) + 1] <- "Landmark clusters"
		km_c[clusters$Freq < ceiling(size_cut * ncol(norm)) + 1] <- "Nonlandmark clusters"
	} else if (method == "diversity_size") {
		km_c <- vector(length = nrow(clusters))
		km_c[clusters$diversity <= diversity_cut & clusters$Freq >= 
			ceiling(size_cut * ncol(norm)) + 1] <- "Landmark clusters"
		km_c[clusters$diversity > diversity_cut | clusters$Freq < 
			ceiling(size_cut * ncol(norm)) + 1] <- "Nonlandmark clusters"
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
      	      abline(v = ceiling(size_cut * ncol(norm)) + 1, lty = 2)
      	} else if (method == "diversity_size") {
      	      abline(h = diversity_cut, lty = 2)
      	      abline(v = ceiling(size_cut * ncol(norm)) + 1, lty = 2)
      	}
		legend("topright", legend = unique(km_c), pch = 16, col = col_palatte[as.factor(unique(km_c))], 
      	      bty = "n")
	}
	return(cc[, c("cell", "landmark_cluster")])
}

require("RColorBrewer")

group_cols <- brewer.pal(n=6, "Set2")


make_mpath_network <- function(SCE, fs=rownames(counts(SCE)), distMethod="pearson") {
	cnts <- counts(SCE)
	sf <- colSums(cnts)
	norm <- t(t(cnts)/sf*1000000)
	norm <- norm[rownames(norm) %in% fs,]
	cellID <- rownames(pData(SCE))
	groupID <- pData(SCE)$sc3_6_clusters
	names(groupID) <- cellID

	#mpath_landmarks <- fixed_landmark_designation(norm, groupID, distMethod="pearson", method="diversity_size")
	faux_landmarks <- cbind(cellID, groupID); colnames(faux_landmarks) <- c("cell", "landmark_cluster")
	mpath_network <- build_network(norm, faux_landmarks, distMethod="pearson")
#	mpath_network <- round(mpath_network/summary(factor(groupID)), digits=2)
	undir_net = mpath_network+t(mpath_network)
#	graph <- graph_from_adjacency_matrix(undir_net, mode="undirected", weighted=TRUE)
	graph <- graph_from_adjacency_matrix(mpath_network, mode="directed", weighted=NULL)
	set.seed(123)
	adj_mat <- cbind(rep(1:6, times=6), rep(1:6, each=6), as.vector(undir_net))
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
	plot(my_layout, pch=16, col=group_cols, xaxt="n", yaxt="n", bty="n", cex=5, xlab="",ylab="")
	apply(adj_mat, 1, plot_edge)
	points(my_layout, pch=16, col=group_cols, cex=5)
	points(my_layout, pch=c("1","2","3","4","5","6"))
	invisible(list(network = mpath_network, plot_layout = my_layout))
}

png("CCA1_mpath_network.png", width=6, height=6, units="in", res=300)
make_mpath_network(CCA1, CCA1_fs)
dev.off()
png("HCC6_mpath_network.png", width=6, height=6, units="in", res=300)
make_mpath_network(HCC6, HCC6_fs)
dev.off()
png("HCC10_mpath_network.png", width=6, height=6, units="in", res=300)
make_mpath_network(HCC10, HCC10_fs)
dev.off()