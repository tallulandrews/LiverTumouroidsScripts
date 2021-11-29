type_col = c("purple", "blue", "forestgreen", "black")
plate_pch = c(16, 7, 25)
library("RColorBrewer")
group_cols <- brewer.pal(n=6, "Set2")

map = read.table("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/Hsap_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
        new = as.character(map[match(x, map[,1]),2])
        new[is.na(new)] = as.character(x[is.na(new)])
        new[duplicated(new)] = x[duplicated(new)]
        return(new)
}


pca_plot <- function(SCE) {
	set.seed(1)
	size_factors <- colSums(counts(SCE))
	norm <- t(t(counts(SCE))/size_factors*median(size_factors))
	norm <- log(norm+1)/log(2)
#	norm <- exprs(SCE)

	plate <- factor(pData(SCE)[,2], levels=c("868", "869", "870"))
#	type  <- factor(pData(SCE)[,4], levels=c("CCA1","HCC6","HCC10","no cells"))
	type  <- pData(SCE)$sc3_6_clusters

	PCA <- prcomp(norm, scale=TRUE)
	plot(PCA$rotation[,1], PCA$rotation[,2], col=group_cols[type], pch=plate_pch[plate], xlab=paste("PC1 (",round(PCA$sdev[1]/sum(PCA$sdev)*100, digits=2)," %)",sep=""), ylab=paste("PC2 (",round(PCA$sdev[2]/sum(PCA$sdev)*100, digits=2)," %)",sep=""))
	legend("topright",as.character(1:6),fill=group_cols, bty="n")
	return(PCA)
}


tsne_plot <- function(SCE, perplexity=10) {
	require("Rtsne")
	set.seed(123)
	size_factors <- colSums(counts(SCE))
	norm <- t(t(counts(SCE))/size_factors*median(size_factors))
	norm <- log(norm+1)/log(2)
	colnames(norm) = 1:length(norm[1,])
	rownames(norm) = 1:length(norm[,1])

	plate <- factor(pData(SCE)[,2], levels=c("868", "869", "870"))
#	type  <- factor(pData(SCE)[,4], levels=c("CCA1","HCC6","HCC10","no cells"))
	type  <- pData(SCE)$sc3_6_clusters

	tsne <- Rtsne(t(norm), perplexity=perplexity)
	plot(tsne$Y[,1], tsne$Y[,2], col=groups_col[type], pch=plate_pch[plate], xlab="Component 1", ylab="Component 2")
	return(tsne$Y)
}


destiny_plot <- function(SCE) {
	require("destiny")
	set.seed(1)
	size_factors <- colSums(counts(SCE))
	norm <- t(t(counts(SCE))/size_factors*median(size_factors))
	norm <- log(norm+1)/log(2)
#	norm <- exprs(SCE)

	plate <- factor(pData(SCE)[,2], levels=c("868", "869", "870"))
#	type  <- factor(pData(SCE)[,4], levels=c("CCA1","HCC6","HCC10","no cells"))
	type  <- pData(SCE)$sc3_6_clusters

	dm <- DiffusionMap(norm)
	dm_dims <- eigenvectors(dm)
	plot(dm_dims[,1], dm_dims[,2], col=group_cols[type], pch=plate_pch[plate], xlab="Dimension 1", ylab="Dimension 2", main="")
	legend("topright",as.character(1:6),fill=group_cols, bty="n")
	return(dm)
}

monocle_pseudotime <- function(SCE) {
	require("M3Drop")
	require("monocle")
	
	pd <- new("AnnotatedDataFrame", data=pData(SCE))
	fd <- new("AnnotatedDataFrame", data=fData(SCE))
	CDS <- newCellDataSet(counts(SCE), phenoData=pd, featureData = fd)

	size_factors <- colSums(counts(SCE))
        norm <- t(t(counts(SCE))/size_factors*median(size_factors))
        norm <- log(norm+1)/log(2)
	m3dgenes <- M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=0.05)

	CDS <- setOrderingFilter(CDS, which(rownames(norm) %in% rownames(m3dgenes)))
	CDS <- reduceDimension(CDS, max_components = 2, reduction_method="DDRTree", norm_method="log", pseudo_expr = 1)
	plot(CDS@reducedDimS[1,], CDS@reducedDimS[2,], col=group_cols[type], pch=plate_pch[plate], xlab="Dimension 1", ylab="Dimension 2", main="")
	legend("topright",as.character(1:6),fill=group_cols, bty="n")

	return(t(CDS@reducedDimS));
}

require("scater")

CCA1 <- readRDS("CCA1_SC3.rds")
fData(CCA1)$Symbol = ensg2symbol(rownames(fData(CCA1)))
HCC6 <- readRDS("HCC6_SC3.rds")
fData(HCC6)$Symbol = ensg2symbol(rownames(fData(HCC6)))
HCC10 <- readRDS("HCC10_SC3.rds")
fData(HCC10)$Symbol = ensg2symbol(rownames(fData(HCC10)))

cellcycle <- read.table("~/Data/Whitfield_CC.txt")
cellcycle_simple <- as.matrix(cellcycle[cellcycle[,1] != "CC",])
cellcycle_simple[cellcycle_simple[,1] == "G2",] = "G2M";
cellcycle_simple[cellcycle_simple[,1] == "S",] = "G1S";
cellcycle_simple = cellcycle_simple[cellcycle_simple[,1] != "MG1",];
G0_genes <- read.table("~/Data/Reactome_G0.txt", header=F)

get_prolif <- function(CC, SCE, suppress.plot=TRUE) {
	require(mixtools)
	cc_g <- fData(SCE)$Symbol %in% CC
	score <- colSums(exprs(SCE)[cc_g,])
	mix <- normalmixEM(score)
	if (!suppress.plot) {
		plot(mix, which=2)
	}
	p1 <- dnorm(score, mean = mix$mu[1], sd=mix$sigma[1]) 
	p2 <- dnorm(score, mean = mix$mu[2], sd=mix$sigma[2])
	if (mix$mu[1] < mix$mu[2]) {
		assign <- p2 > p1
	} else {
		assign <- p1 > p2
	}
	return(list(score=score, assign=assign))
}

get_prolif2 <- function(CC, SCE, name="cycle genes", suppress.plot=TRUE) {
	require(mixtools)
        require(matrixStats)
	norm <- counts(SCE)
	nfactor <- colSums(norm);
	norm <- t(t(norm)/nfactor*median(nfactor))
        loggednorm <- log(norm+1)/log(2)
        cc_g <- fData(SCE)$Symbol %in% CC
	loggednorm <- loggednorm[cc_g,]
	loggednorm <- loggednorm[rowVars(loggednorm) > 0,]
        loggednorm <- (loggednorm-rowMeans(loggednorm))/rowVars(loggednorm)
        score <- colSums(loggednorm)
        mix <- normalmixEM(score)
        p1 <- dnorm(score, mean = mix$mu[1], sd=mix$sigma[1])
        p2 <- dnorm(score, mean = mix$mu[2], sd=mix$sigma[2])
	
	#means_test <- test.equality(score, arbmean=TRUE)$p.value # This is not reliable running in multiple times on same data gets very different answers!
	
	m1 = mix$mu[1]
	m2 = mix$mu[2]
#	means_test <- (m1-m2)/(sqrt(mix$sigma[1]^2+mix$sigma[2]^2)/sqrt(length(norm[1,])))
#	means_test <- pnorm(abs(means_test), lower.tail=F) < 0.05
	means_test <- abs(m1-m2) > min(mix$sigma)
	

        if (m1 < m2 & means_test) {
                assign <- p2 > p1
        } else if (m2 < m1 & means_test) {
                assign <- p1 > p2
        } else {
		assign <- rep(FALSE, times=length(p1))
		name=paste(name, "- No Difference")
	}
        if (!suppress.plot) {
                plot(mix, which=2, main2=name, xlab2="Score")
        }
        return(list(score=score, assign=assign))

}

#prolif_analysis2 <- function(SCE, perplexity=10) {
#	require("M3Drop")
#	norm <- counts(SCE)
#	nfactor <- colSums(norm);
#	norm <- t(t(norm)/nfactor*median(nfactor))
#	Features <- M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=0.05)
#	tSNE <- tsne_plot(SCE[rownames(norm) %in% rownames(Features),], perplexity=perplexity)
#
#	sets <- levels(factor(cellcycle[,1]))
#	out <- matrix(0, nrow=length(pData(SCE)[,1]), ncol=length(sets))
#	for( i in 1:length(sets) ) {
#		s = sets[i]
#		assign<-get_prolif2(cellcycle[cellcycle[,1]==s,2], SCE)$assign
#		out[,i] = assign
#	}
#	CC_score <- factor(rowSums(out))
#	CC_col = c("grey80","grey75","grey50","grey25","grey15","black")
#	plot(tSNE[,1], tSNE[,2], col=CC_col[CC_score], pch=16)
#}

prolif_analysis<- function(SCE, name="Data") {
	
	# Get Prolif
	png(paste(name, "Proliferation_MixtureModels.png", sep="_"), width=8, height=8*2/3, units="in", res=300)
	par(mfrow=c(2,3))
	sets <- levels(factor(cellcycle[,1]))
	out <- matrix(0, nrow=length(pData(SCE)[,1]), ncol=length(sets))
	for( i in 1:length(sets) ) {
		s = sets[i]
		assign<-get_prolif2(cellcycle[cellcycle[,1]==s,2], SCE, name=s, suppress.plot=F)$assign
		out[,i] = assign
	}
	dev.off()

	summary(factor(rowSums(out)))
	prolif = rowSums(out) >= 3

#	plot(t[,1], t[,2], pch=plate_pch[plate], col=c("black", "red")[prolif+1], xlab="Component 1", ylab="Component 2")
	# Prolif vs Clusters
	png(paste(name, "Proliferation_Figure.png", sep="_"), width=8, height=4, units="in", res=300)
	par(mfrow=c(1,2))
	par(mar=c(4,4,1,1))

#	require("M3Drop")
#	FS <- M3DropFeatureSelection(counts(SCE), mt_method="fdr", mt_threshold=0.05, suppress.plot=TRUE)
	pca <- pca_plot(SCE)
#	t <- tsne_plot(SCE, perplexity=perplexity)

	plate <- factor(pData(SCE)[,2], levels=c("868", "869", "870"))
	type  <- factor(pData(SCE)[,4], levels=c("CCA1","HCC6","HCC10","no cells"))

	plot(pca$rotation[,1], pca$rotation[,2], col=c("black", "red")[prolif+1], pch=plate_pch[plate], xlab=paste("PC1 (",round(pca$sdev[1]/sum(pca$sdev)*100, digits=2)," %)",sep=""), ylab=paste("PC2 (",round(pca$sdev[2]/sum(pca$sdev)*100, digits=2)," %)",sep=""))
	dev.off()

	# Output
	pData(SCE)$Proliferating <- prolif

	output <- cbind(aggregate(prolif, by=list(pData(SCE)$sc3_6_clusters), sum),summary(factor(pData(SCE)$sc3_6_clusters)))
	colnames(output) = c("sc3","prolif", "ncells")
#	return(output)
	return(prolif)
}

library("M3Drop")

cca1_prolif_out <- prolif_analysis(CCA1, "CCA1")
pData(CCA1)$Proliferating <- cca1_prolif_out

hcc6_prolif_out <- prolif_analysis(HCC6, "HCC6")
pData(HCC6)$Proliferating <- hcc6_prolif_out

hcc10_prolif_out <- prolif_analysis(HCC10, "HCC10")
pData(HCC10)$Proliferating <- hcc10_prolif_out

prolif_analysis2 <- function(SCE, name="Data") {
	
#	g0 <- get_prolif2(G0_genes[,1], SCE, name="G0", suppress.plot=F)$assign
	png(paste(name, "MixtureModels2_Figure.png", sep="_"), width=8, height=4, units="in", res=300)
	par(mfrow=c(1,2))
	par(mar=c(4,4,1,1))
	g1 <- get_prolif2(cellcycle_simple[cellcycle_simple[,1] == "G1S",2], SCE, name="G1S", suppress.plot=F)$assign
	g2 <- get_prolif2(cellcycle_simple[cellcycle_simple[,1] == "G2M",2], SCE, name="G2M", suppress.plot=F)$assign
	dev.off()
	state <- rep("Non", times=length(g1))
#	state[g0] = "G0";
	state[g1] = "G1";
	state[g2] = "G2";
#	state[g0 & g1 | g0 & g2] = "G0+G1/G2";
	state[g1 & g2] = "Both";
	
	png(paste(name, "Proliferation2_Figure.png", sep="_"), width=8, height=4, units="in", res=300)
	par(mfrow=c(1,2))
	par(mar=c(4,4,1,1))
	pca <- pca_plot(SCE)

	plate <- factor(pData(SCE)[,2], levels=c("868", "869", "870"))
	type  <- factor(pData(SCE)[,4], levels=c("CCA1","HCC6","HCC10","no cells"))
	state <- factor(state, levels=c("Non", "G1","G2","Both"))

	plot(pca$rotation[,1], pca$rotation[,2], col=c("black","goldenrod1", "darkorange", "red")[state], pch=plate_pch[plate], xlab=paste("PC1 (",round(pca$sdev[1]/sum(pca$sdev)*100, digits=2)," %)",sep=""), ylab=paste("PC2 (",round(pca$sdev[2]/sum(pca$sdev)*100, digits=2)," %)",sep=""))
	legend("topright", c("Non", "G1", "G2", "Both"), col=c("black","goldenrod1", "darkorange", "red"), pch=16)
	dev.off()

	pData(SCE)$CC_state <- state
	# Output
	output <- cbind(aggregate(state, by=list(pData(SCE)$sc3_6_clusters), summary),summary(factor(pData(SCE)$sc3_6_clusters)))
	colnames(output) = c("sc3","state", "ncells")
#	return(output)
	return(state)
}

cca1_prolif2_out <- prolif_analysis2(CCA1, "CCA1")
pData(CCA1)$CC_state <- cca1_prolif2_out

hcc6_prolif2_out <- prolif_analysis2(HCC6, "HCC6")
pData(HCC6)$CC_state <- hcc6_prolif2_out

hcc10_prolif2_out <- prolif_analysis2(HCC10, "HCC10")
pData(HCC10)$CC_state <- hcc10_prolif2_out

saveRDS(CCA1, file="CCA1_SC3_Prolif.rds")
saveRDS(HCC6, file="HCC6_SC3_Prolif.rds")
saveRDS(HCC10, file="HCC10_SC3_Prolif.rds")
