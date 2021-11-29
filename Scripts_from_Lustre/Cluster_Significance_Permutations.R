require("M3Drop")
require("scater")
require("matrixStats")
require("RColorBrewer")

get_log_norm <- function(SCE) {
	dat <- counts(SCE);
	fac <- colSums(dat);
	norm <- t(t(dat)/fac*median(fac))
	return(log(norm +1)/log(2));
}

CCA1 <- readRDS("CCA1_SC3.rds")
HCC6 <- readRDS("HCC6_SC3.rds")
HCC10 <- readRDS("HCC10_SC3.rds")

exprs(CCA1) <- get_log_norm(CCA1)
exprs(HCC6) <- get_log_norm(HCC6)
exprs(HCC10) <- get_log_norm(HCC10)

my_row_mean_aggregate <- function(mat, groups) {
        # Much faster version
        MAT <- as.matrix(mat)
        x <- split(seq(ncol(MAT)), groups)
        result <- sapply(x, function(a) rowMeans(MAT[,a]))

        return(result);
}

#CCA1_dist_mat <- 1-cor(my_row_mean_aggregate(exprs(CCA1), pData(CCA1)$sc3_6_clusters))


permute_within_groups <- function(expr_mat, clusters) {
	clusters <- as.numeric(factor(clusters))
#	reorder<- order(clusters)
	new_mat <- expr_mat#[,reorder]
	clust <- clusters#[reorder]
	require("permute")
	
	nrow=length(new_mat[,1])
	p_mat = new_mat;

	
	for (i in 1:max(clust)) {
		p <- permute::shuffleSet(sum(clust==i), nset=nrow, quiet=TRUE, check=FALSE)
		perm <- t(sapply(seq_len(nrow(p)), function(j,P) p_mat[j,clust==i][P[j,]], P=p))
		p_mat[,clust==i]<-perm
	}
	return(p_mat)
}

#thing <- permute_within_groups(exprs(CCA1), pData(CCA1)$sc3_6_clusters) 
#dist_fun <- function(x) {as.dist(1-cor(x, method="spearman"))}


#p_dists <- as.matrix(dist_fun(thing))
#diag(p_dists) <- 2

#closest <- apply(p_dists, 1, min)
#c_closest <- sapply(seq_len(length(closest)), function(x){pData(CCA1)$sc3_6_clusters[which(p_dists[x,] == closest[x])]} )
#table(c_closest, pData(CCA1)$sc3_6_clusters)

#### FORGET ALL THE PERMUTATION STUFF JUST USE THE BELOW ####
dist_fun <- function(x) {as.dist(1-cor(x, method="spearman"))}

NN_check <- function(SCE, dist_fun=function(x){dist(t(x))}) {
	o_dists <- as.matrix(dist_fun(exprs(SCE)))
	diag(o_dists) <- 2

	closest <- apply(o_dists, 1, min)
	c_closest <- sapply(seq_len(length(closest)), function(x){pData(SCE)$sc3_6_clusters[which(o_dists[x,] == closest[x])]} )
	factor(c_closest, levels=1:6)
	t(t(table(c_closest, pData(SCE)$sc3_6_clusters))/summary(factor(pData(SCE)$sc3_6_clusters)))
}

#CCA1
Scor_NN_check <- NN_check(CCA1, function(x) {as.dist(1-cor(x, method="spearman"))})
Pcor_NN_check <- NN_check(CCA1, function(x) {as.dist(1-cor(x, method="pearson"))})
Eucl_NN_check <- NN_check(CCA1, function(x){dist(t(x))})
Pcor_NN_check<- rbind(Pcor_NN_check, rep(0, times=6))
Pcor_NN_check<- Pcor_NN_check[c(1,6,2,3,4,5),]
CCA1_overall <- (Scor_NN_check+Pcor_NN_check+Eucl_NN_check)/3

#HCC6
Scor_NN_check <- NN_check(HCC6, function(x) {as.dist(1-cor(x, method="spearman"))})
Pcor_NN_check <- NN_check(HCC6, function(x) {as.dist(1-cor(x, method="pearson"))})
Eucl_NN_check <- NN_check(HCC6, function(x){dist(t(x))})
Scor_NN_check<- rbind(Scor_NN_check, rep(0, times=6))
Pcor_NN_check<- rbind(Pcor_NN_check, rep(0, times=6))
HCC6_overall <- (Scor_NN_check+Pcor_NN_check+Eucl_NN_check)/3

#HCC10
Scor_NN_check <- NN_check(HCC10, function(x) {as.dist(1-cor(x, method="spearman"))})
Pcor_NN_check <- NN_check(HCC10, function(x) {as.dist(1-cor(x, method="pearson"))})
Eucl_NN_check <- NN_check(HCC10, function(x){dist(t(x))})
HCC10_overall <- (Scor_NN_check+Pcor_NN_check+Eucl_NN_check)/3


o_dists <- as.matrix(dist_fun(exprs(HCC6)))
diag(o_dists) <- 2

closest <- apply(o_dists, 1, min)
c_closest <- sapply(seq_len(length(closest)), function(x){pData(HCC6)$sc3_6_clusters[which(o_dists[x,] == closest[x])]} )
table(c_closest, pData(HCC6)$sc3_6_clusters)



o_dists <- as.matrix(dist_fun(exprs(HCC10)))
diag(o_dists) <- 2

closest <- apply(o_dists, 1, min)
c_closest <- sapply(seq_len(length(closest)), function(x){pData(HCC10)$sc3_6_clusters[which(o_dists[x,] == closest[x])]} )
table(c_closest, pData(HCC10)$sc3_6_clusters)



