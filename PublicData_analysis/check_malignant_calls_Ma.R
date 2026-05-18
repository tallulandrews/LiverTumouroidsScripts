if (FALSE) {
	source("/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Merixtell/MerixtellLoadData.R")
	require(Seurat)
	require(infercnv)
	require(mclust)
	seur_obj <- load_MaTumors(min.cells=20, min.features=500)
	out_tag <- "MaTumors"

	infercnv_out <- readRDS("MaTumors_infercnv_obj.rds")

	aneuploidy <- Matrix(round(infercnv_out@expr.data, digits=2)-1)
	#aneuploidy[abs(aneuploidy) <= 0.1] <- 0
	total_cnv_signal <- colSums(abs(aneuploidy))

	n_tumours <- length(unique(seur_obj@meta.data$Sample))
	classification <- Mclust(log2(total_cnv_signal+1), G=2)
	compare <- data.frame(seur_obj@meta.data[match(names(classification$classification), rownames(seur_obj@meta.data)),], 
				classification$classification, 
				total_cnv_signal)
	tmp<- table(compare$Type, compare$classification.classification)
	tmp["Malignant cell",]/colSums(tmp)

	TP <- tmp["Malignant cell",2]
	FP <- sum(tmp[rownames(tmp) != "Malignant cell",2])
	TN <- sum(tmp[rownames(tmp) != "Malignant cell",1])
	FN <- tmp["Malignant cell",1]

	png(paste(out_tag, "inferCNV_mclust_classification.png", sep="_"), width=6, height=6, units="in", res=150)
	par(mar=c(7,4,1,1))
	barplot(t(tmp), las=2, col=c("grey75", "black"), ylab="N Cells")
	legend("topleft", fill=c("grey75", "black"), c("normal", "tumour"), title="inferCNV classification", bty="n")
	legend("topright", paste("Accuracy:", round((TP+TN)/(TP+FP+TN+FN), digits=2)), bty="n")
	dev.off()
}



call_malignant <- function(infercnv_obj) {
	require(mclust)
	require(infercnv)
	aneuploidy <- Matrix(round(infercnv_obj@expr.data, digits=2)-1)
	total_cnv_signal <- Matrix::colSums(abs(aneuploidy))
	classification <- Mclust(log2(total_cnv_signal+1), G=2)
	categories <- as.character(classification$classification)
	categories[categories == "1"] <- "Non-Malignant"
        categories[categories == "2"] <- "Malignant"
	names(categories) <- names(total_cnv_signal)
	return(categories)
}
