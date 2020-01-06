require("M3Drop")
require("scater")
require("matrixStats")
require("RColorBrewer")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/DiffExpr/DE_functions.R")

get_log_norm <- function(SCE) {
	dat <- counts(SCE);
	fac <- colSums(dat);
	norm <- t(t(dat)/fac*median(fac))
	return(log(norm +1)/log(2));
}

get_cluster_profile <- function(expr_mat, cluster_bool) {
	# This is garbage
	profile <- rowMeans(expr_mat[,cluster_bool]) - rowMeans(expr_mat)
	profileD <- rowMeans(expr_mat[,cluster_bool] > 0) - rowMeans(expr_mat > 0)
	stddev <- sqrt( rowVars(expr_mat[,cluster_bool])/sum(cluster_bool) )
	stddev_all <- sqrt( rowVars(expr_mat)/length(expr_mat[1,]) )
	stddev[stddev==0] = Inf
	return(list(diff = profile, diffD = profileD, stddev = stddev, sdev_all = stddev_all))
}

cluster_relative_heatmap <- function(expr_mat, clusters, distfun=dist, hclustfun=function(x){hclust(x,method="complete")}, gene_list = rownames(expr_mat), npermute=0) {
	require("gplots")
	require("permute")
	cluster_labs = unique(clusters)
#	my_profiles <- matrix(0, ncol=length(cluster_labs), nrow=length(expr_mat[,1]))
#	for(i in 1:length(cluster_labs)) {
#		my_profiles[,i] = get_cluster_profile(expr_mat, clusters==cluster_labs[i])$diff
#	}
	cluster_means <- my_row_mean_aggregate(expr_mat, clusters)
#	my_profiles <- cluster_means-rowMeans(expr_mat)
	my_profiles <- cluster_means-rowMeans(cluster_means)
	my_profiles <- my_profiles[rownames(my_profiles) %in% gene_list,]

	D <- vector();
	M <- my_profiles;
	ncol <- length(my_profiles[1,])
	nrow <- length(my_profiles[,1])
	if (npermute > 0) {
		set.seed(101); # reproduciblity
		for(rep in 1:npermute) {
			P <- shuffleSet(ncol, nset=nrow, quiet=TRUE);
			perm <- t(sapply(seq_len(nrow(P)), function(i, P, M) M[i,P[i,]],M=M,P=P))
			D <- c(D,as.vector(distfun(t(perm))))
		}
		threshold <- quantile(D,probs=0.05)
	
		my_dists <- distfun(t(my_profiles))
		my_hclust <- hclustfun(my_dists)
		my_sig <- cutree(my_hclust, h=threshold)
		my_col_vec <- rainbow(n=max(my_sig))[my_sig]
	}

#	colnames(my_profiles) <- cluster_labs
	heatcols <- colorRampPalette(c("blue","white","red"))(255)
	if (npermute == 0) {
		heatout <- heatmap.2(my_profiles, trace="n", col=heatcols, symbreaks=TRUE, key.title="", key.xlab="Relative Expression", hclustfun=hclustfun, distfun=distfun)
	} else {
		heatout <- heatmap.2(my_profiles, trace="n", col=heatcols, symbreaks=TRUE, key.title="", key.xlab="Relative Expression", hclustfun=hclustfun, distfun=distfun, ColSideColors=my_col_vec)
	}
	invisible(heatout)
}

my_row_mean_aggregate <- function(mat, groups) {
#	out <-t( aggregate(t(mat), by=list(groups), mean) )
#	colnames(out) <- out[1,]
#	out <- out[-1,]
#	storage.mode(out) <- "numeric"

	# Much faster version of the above!
	MAT <- as.matrix(mat)
	x <- split(seq(ncol(MAT)), groups)
	result <- sapply(x, function(a) rowMeans(MAT[,a]))

	return(result);
}

cluster_heatmap <- function(expr_mat, clusters, DEonly=TRUE) {
	require("gplots")
	cluster_labs = unique(clusters)
	my_profiles <- my_row_mean_aggregate(expr_mat, cluster)
#	my_profiles <-t( aggregate(t(expr_mat), by=list(cluster), mean) )
#	colnames(my_profiles) <- my_profiles[1,]
#	my_profiles <- my_profiles[-1,]
#	storage.mode(my_profiles) <- "numeric"

	rownames(my_profiles) <- rownames(expr_mat);
	to_plot = rownames(expr_mat);
	if (DEonly) {
		kw <- apply(expr_mat, 1, function(x) {kruskal.test(x ~ factor(clusters))$p.value})
		to_plot = rownames(expr_mat[p.adjust(kw) < 0.05,])
	}
	heatcols <- colorRampPalette(c("black","blue","white","red","yellow"))(255)
	heatmap.2(my_profiles[rownames(my_profiles) %in% to_plot,], symbreaks=TRUE, trace="n", col=heatcols, key.title="", key.xlab="Relative Expression", hclustfun=function(x){hclust(x, method="average")})
}

cosine_dist <- function(x,y) {x %*% y / sqrt(x%*%x * y%*%y)}

match_clusters <- function(SCE1, SCE2) {
	profiles1 <- list();
	for (i in 1:6) {
		profiles1[[i]] = get_cluster_profile(exprs(SCE1), pData(SCE1)$sc3_6_clusters == i)
	}
	
	profiles2 <- list();
	for (i in 1:6) {
		profiles2[[i]] = get_cluster_profile(exprs(SCE2), pData(SCE2)$sc3_6_clusters == i)
	}
	matches1 =rep(0, times=length(profiles1))
	matches1val =rep(0, times=length(profiles1))
	matches2 = rep(0, times=length(profiles2))
	matches2val = rep(0, times=length(profiles2))
	for (i in 1:6) {
		for (j in 1:6) {
			cor_val <- cor(profiles1[[i]]$diff, profiles2[[j]]$diff)
			if (cor_val > matches1val[i]) {
				matches1val[i] = cor_val;
				matches1[i] = j;
			}
			if (cor_val > matches2val[j]) {
				matches2val[j] = cor_val;
				matches2[j] = i;
			}

		}
	}
	recip = which(1:length(matches1) == matches2[matches1])
	recip_table = cbind(recip,matches1[recip])
	colnames(recip_table) = c("SCE1", "SCE2")
	match_table = cbind(matches1, matches2)
	colnames(match_table) = c("SCE1->SCE2", "SCE2->SCE1")
	return(list(matches=match_table, recip=recip_table))
}

fancy_heatmap <- function(SCE, KW_res, sce_markers) {
	require("gplots");
#	KW_res <- apply(exprs(SCE), 1, function(x) {kruskal.test(x ~ pData(SCE)$sc3_6_clusters)$p.value})
	group_cols <- brewer.pal(n=6, "Set2")
	

	toplot = p.adjust(KW_res, method="bon") < 0.05
	toplot[is.na(toplot)] = FALSE
	toplot = unique(c( names(toplot[toplot]), sce_markers$Gene[sce_markers$pvalue < 0.05/length(sce_markers[,1])/6]))
	row_col = rep(0, times=length(toplot))
	names(row_col) <- toplot
	marker_dat <- sce_markers[(match(names(row_col), sce_markers$Gene)),]
	row_col = group_cols[marker_dat[,2]]
	row_col[marker_dat$pvalue > 0.05/length(sce_markers[,1])/6] = "white"
	col_col <- group_cols[pData(SCE)$sc3_6_clusters]

	heatcols <- colorRampPalette(c("black","blue","cyan","green","yellow", "white"))(255)

	heatmap.2(exprs(SCE)[rownames(exprs(SCE)) %in% toplot,], trace="n", RowSideColors=row_col, ColSideColors=col_col, col=heatcols, key.title="", key.xlab="Log2 Expression", hclustfun=function(x){hclust(x, method="ward.D2")})


}

do_Combat<- function(norm, batch) {
	require("sva")
	mod0 = model.matrix(~1, data=norm)
	combat_data = ComBat(dat=Combined_data, batch=factor(batch), mod=mod0, par.prior=TRUE, prior.plots=FALSE)
	return(combat_data)
}

do_MeanRemoval<- function(norm, batch, clusters=NULL) {
	batch <- as.factor(batch);
	batch_labs <- levels(batch);
	corrected <- norm;
	for( b in 1:length(batch_labs)) {

		if (!is.null(clusters)) {
			cluster_means <- my_row_mean_aggregate(norm[,batch==batch_labs[b]], clusters[batch==batch_labs[b]])

		        corrected[,batch==batch_labs[b]] <- corrected[,batch==batch_labs[b]]-rowMeans(cluster_means)
		} else {
			corrected[,batch==batch_labs[b]] <- corrected[,batch==batch_labs[b]]-rowMeans(corrected[,batch==batch_labs[b]])
		}
	}
	return(corrected);
}
