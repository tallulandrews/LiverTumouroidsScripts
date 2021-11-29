require("monocle")
require("M3Drop")
require("scater")
require("matrixStats")
require("RColorBrewer")

CCA1 <- readRDS("CCA1_SC3.rds")
HCC6 <- readRDS("HCC6_SC3.rds")
HCC10 <- readRDS("HCC10_SC3.rds")

markers <- readRDS("Laura_SC3_k6_Markers.rds")

fancy_heatmap <- function(SCE, markers) {
	

}

monocle_setup <- function(SCE) {
	pd <- new("AnnotatedDataFrame", data=pData(SCE))
        fd <- new("AnnotatedDataFrame", data=fData(SCE[keep,]))
	
        CDS <- newCellDataSet(counts(SCE)[keep,], phenoData=pd, featureData = fd, expressionFamily=negbinomial.size())
	CDS <- estimateSizeFactors(CDS)
	return(CDS)
}

















SCE = HCC10

monocle_pseudotime <- function(SCE) {
	thing <- counts(SCE);
#	keep <- rowSums(thing > 5) > 2;
#	keep <- keep & rowVars(thing) > 1;
	pd <- new("AnnotatedDataFrame", data=pData(SCE))
        fd <- new("AnnotatedDataFrame", data=fData(SCE[keep,]))
	
        CDS <- newCellDataSet(counts(SCE)[keep,], phenoData=pd, featureData = fd, expressionFamily=negbinomial.size())
	CDS <- estimateSizeFactors(CDS)

        size_factors <- colSums(counts(SCE))
        norm <- t(t(counts(SCE))/size_factors*median(size_factors))
        m3dgenes <- M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=0.05)
        norm <- log(norm+1)/log(2)
	sig_marks <- markers$HCC10$Gene[markers$HCC10$pvalue < 0.05/length(markers$HCC10[,1])/6]

        CDS <- setOrderingFilter(CDS, which(rownames(norm) %in% rownames(sig_marks)))
        CDS <- reduceDimension(CDS, max_components = 2, reduction_method="DDRTree", norm_method="log", pseudo_expr = 1)
	CDS <- orderCells(CDS, reverse=TRUE)

	type  <- pData(SCE)$sc3_6_clusters
	plate <- factor(pData(SCE)[,2], levels=c("868", "869", "870"))

	plate_pch = c(16, 7, 25)
	type_col = c("purple", "blue", "forestgreen", "black")
	group_cols <- brewer.pal(n=6, "Set2")

	plot_cell_trajectory(CDS)

        plot(CDS@reducedDimS[1,], CDS@reducedDimS[2,], col=group_cols[type], pch=plate_pch[plate], xlab="Dimension 1", ylab="Dimension 2", main="")
        legend("topright",as.character(1:6),fill=group_cols, bty="n")

        return(t(CDS@reducedDimS));
}

plot_genes_in_pseudotime(HCC10_Monocle[Laura_markers,])

Laura_markers = c("ENSG00000088325","ENSG00000122952", "ENSG00000130649", "ENSG00000119888", "ENSG00000139292")
