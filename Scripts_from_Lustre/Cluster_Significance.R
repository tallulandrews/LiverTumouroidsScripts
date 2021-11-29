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

## DE genes
CCA1_kw <- apply(exprs(CCA1), 1, function(x) {kruskal.test(x ~ pData(CCA1)$sc3_6_clusters)$p.value})
HCC6_kw <- apply(exprs(HCC6), 1, function(x) {kruskal.test(x ~ pData(HCC6)$sc3_6_clusters)$p.value})
HCC10_kw <- apply(exprs(HCC10), 1, function(x) {kruskal.test(x ~ pData(HCC10)$sc3_6_clusters)$p.value})


get_params <- function(SCE, kw.out) {
	require("M3Drop")
	expr_mat = 2^exprs(SCE)-1;
	expr_mat <- expr_mat[!is.na(kw.out) | kw.out > 0.05,]
	vals <- bg__calc_variables(expr_mat)
	fit <- bg__fit_MM(vals$p, vals$s)
	mean2disp <- bg__get_mean2disp(expr_mat);
	return(list(K=fit$K, mean2disp=mean2disp));
}

add_dropouts <- function(x, mu, K) {
#	p_drop <- 1- mu/(mu+K);
	expect_pos <- mu*length(x)/sum(x>0);
	p_drop <- K/expect_pos
	toss <- runif(n=length(x));
	x[toss < p_drop] <- 0;
	return(x);
}
my_row_mean_aggregate <- function(mat, groups) {
	MAT <- as.matrix(mat)
	x <- split(seq(ncol(MAT)), groups)
	results <- sapply(x, function(a) rowMeans(MAT[,a]))
	return(result);
}

do_sim <- function(expr_mat, clusters, params, n=100, n.cores=1) {
	cluster_means <- my_row_mean_aggregate(expr_mat, clusters);
	tmp <- cluster_means[,clusters]	
	expr_mat[expr_mat==0] = tmp[expr_mat==0]; # correct zeros if possible

	for (s in 1:n) {
		sim = apply(expr_mat, 2, function(x) {
			sapply(x, function(mu){add_dropouts(rnbinom(1, mu=mu, size=1/params$mean2disp(mu)), my, params$K)})
			}) 
		out <- do_SC3(sim, n.cores=n.cores)
	}
	return(sim);
}

do_SC3<-function(expr_mat, n.cores=1) {
	require("scater")
	require("SC3")
	pd <- new("AnnotatedDataFrame", data=data.frame(num=1:length(expr_mat[,1])))
	sceset <- newSCESet(countData=sim, phenoData=pd)
	sceset <- sc3_prepare(sceset)
	sceset <- sc3(sceset, ks=6, n_cores=n.cores)
	return(CCA1@sc3$consensus$"6")
#	return(pData(sceset)$sc3_6_clusters);	
}
