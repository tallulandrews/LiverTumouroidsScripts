set.seed(8932)

source("/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Merixtell/MerixtellLoadData.R")
source("/home/tandrew6/projects/def-tandrew6/tandrew6/Meritxell/check_malignant_calls_Ma.R")
final_TIC_Diff_markers <- read.table("Final_Markers.txt", header=FALSE)
metabolic_markers <- read.table("B27_components_genesets_Harmonizome_aggregate_gene_list.csv", sep=",", header=TRUE)
colnames(final_TIC_Diff_markers) <- c("Gene", "Cluster")
metabolic_markers <- data.frame(Gene=metabolic_markers[,2], Pathway=metabolic_markers[,1])
require(Seurat)
require(CycleMix)

seurat_analysis_pipeline <- function(seur_obj) {
	set.seed(1891)
#	seur_obj <- Seurat::SCTransform(seur_obj)
	seur_obj <- Seurat::NormalizeData(seur_obj)
	seur_obj <- Seurat::ScaleData(seur_obj)
	seur_obj <- Seurat::FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)
	hvgs <- VariableFeatures(seur_obj); hvgs <- hvgs[!grepl("^MT-", hvgs)];
	VariableFeatures(seur_obj) <- hvgs;
	seur_obj <- Seurat::RunPCA(seur_obj, features = VariableFeatures(object = seur_obj))
	seur_obj <- Seurat::FindNeighbors(seur_obj, dims = 1:20)
	seur_obj <- Seurat::FindClusters(seur_obj)
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes
	seur_obj <- Seurat::CellCycleScoring(seur_obj, s.features = s.genes, g2m.features=g2m.genes, set.ident=TRUE)
	seur_obj <- Seurat::RunUMAP(seur_obj, dims = 1:20, parallel=FALSE)
	return(seur_obj)
}

visualize_CC <- function(seur_obj, outname="project") {
	require(RColorBrewer)
	CC_col = c("grey50", "black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("None", "G0", "G1/S", "G2/M");


	# My CC Phase Assignment
	require(CycleMix)
	require(SingleCellExperiment)
	require(mclust)
	require(ggplot2)
	SCE <- as.SingleCellExperiment(seur_obj)
	logcounts(SCE) <- as.matrix(logcounts(SCE))
	rowData(SCE)$feature_symbol <- rownames(seur_obj)
	CC_genes <- HGeneSets$Tirosh
	CC_genes <- CC_genes[CC_genes[,1] %in% rowData(SCE)$feature_symbol,]
	cc_assignments <- CycleMix::classifyCells(SCE, CC_genes)
	seur_obj@meta.data$CycleMix_phase <- cc_assignments$phase
	
	tmp_col <- CC_col[2:4]; names(tmp_col) <- c("None", "G1S", "G2M")
	png(paste(outname, "CycleMix_phase.png", sep="_"), width=8, height=8, units="in", res=150)
	print(DimPlot(seur_obj, group.by="CycleMix_phase") + 
		scale_color_manual(values=tmp_col))
	dev.off()
	tmp_col <- CC_col[2:4]; names(tmp_col) <- c("G1", "S", "G2M")
	png(paste(outname, "Seurat_Phase.png", sep="_"), width=8, height=8, units="in", res=150)
	print(DimPlot(seur_obj, group.by="Phase")+ scale_color_manual(values=tmp_col))
	dev.off()
	png(paste(outname, "SeuratPhase_bycluster.png", sep="_"), width=8, height=8, units="in", res=150)
	barplot(table(factor(seur_obj$Phase, levels=c("G1", "S", "G2M")), seur_obj@meta.data$seurat_clusters), col=CC_col[2:4] )
	dev.off()
	
	png(paste(outname, "CycleMix_bycluster.png", sep="_"), width=8, height=8, units="in", res=150)
	barplot(table(factor(seur_obj$CycleMix_phase, levels=c("None", "G1S", "G2M")), seur_obj@meta.data$seurat_clusters), col=CC_col[2:4] )
	dev.off()
	
}

visualize_key_genes <- function(seur_obj, genes=final_TIC_Diff_markers) {
	labelled_genes <- unlist(genes[!duplicated(genes[,1]),1]);
	names(labelled_genes) <- genes[!duplicated(genes[,1]),2];
	this_plot <- DotPlot(seur_obj, features=labelled_genes, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	return(this_plot)
}

gene_figures <- function(seur_obj, outname, type=c("CC", "HCC")) {
	if (type == "CC") {
		TIC_Diff <- final_TIC_Diff_markers[grep("CCA", final_TIC_Diff_markers[,2]),]
	} else if (type == "HCC") {
		TIC_Diff <- final_TIC_Diff_markers[grep("HCC", final_TIC_Diff_markers[,2]),]
	}
	png(paste(outname, type, "Markers.png", sep="_"), width=10, height=6, units="in", res=100)
	this_plot <- visualize_key_genes(seur_obj, 
			genes=TIC_Diff[grep("CCA", final_TIC_Diff_markers[,2]),])
	print(this_plot)
	dev.off()

	png(paste(outname, "Metabolism.png", sep="_"), width=16, height=6, units="in", res=100)
	this_plot <- visualize_key_genes(seur_obj, genes=metabolic_markers[metabolic_markers[,2] %in% c("Fattyacid_lipids", "Linoleic Acid", "Linolenic Acid"),])
	print(this_plot)
	dev.off()
	
	png(paste(outname, "Metabolism2.png", sep="_"), width=16, height=6, units="in", res=100)
	this_plot <- visualize_key_genes(seur_obj, genes=metabolic_markers[metabolic_markers[,2] %in% c("Mitochondria", "Beta-oxidation"),])
	print(this_plot)
	dev.off()
}

general_figures <- function(seur_obj, outname, sample_col="Sample") {
	png(paste(outname, "UMAP.png", sep="_"), width=6, height=6, units="in", res=100)
	print(DimPlot(seur_obj, group.by="seurat_clusters", label=TRUE))
	dev.off()
	
	png(paste(outname, "SampleUMAP.png", sep="_"), width=6, height=6, units="in", res=100)
	print(DimPlot(seur_obj, group.by=sample_col, label=TRUE))
	dev.off()
}

seur_obj <- load_MaTumors(min.cells=20, min.features=500)
out_tag <- "MaTumors"

# Load Tumor / non-tumor annotations

# Subset the tumor cells
Ma_tumour_obj <- seur_obj[,seur_obj@meta.data$Type == "Malignant cell"]
HCC_obj <- Ma_tumour_obj[,Ma_tumour_obj@meta.data$tumour_type == "Hepatocellularcarcinoma"]
CC_obj <- Ma_tumour_obj[,Ma_tumour_obj@meta.data$tumour_type == "Intrahepaticcholangiocarcinoma"]
#cluster the tumour cells
HCC_obj <- seurat_analysis_pipeline(HCC_obj)
CC_obj <- seurat_analysis_pipeline(CC_obj)

visualize_CC(HCC_obj, "MaHCC")
visualize_CC(CC_obj, "MaCC")

saveRDS(HCC_obj, "MaHCC_clustered.rds")
saveRDS(CC_obj, "MaCC_clustered.rds")

### MA CC ###
png(paste(out_tag, "CC", "markerUMAP.png", sep="_"), width=6*2, height=6, units="in", res=100)
print(FeaturePlot(CC_obj, features=c("CA9", "ATP1B3")))
dev.off()

cor(CC_obj@assays$RNA@data["CA9",], CC_obj@assays$RNA@data["ATP1B3",])

gene_figures(CC_obj, paste(out_tag, "CC", sep="_"), type="CC")
general_figures(CC_obj, paste(out_tag, "CC", sep="_"), sample_col="Sample")
### MA HCC ###
gene_figures(HCC_obj, paste(out_tag, "HCC", sep="_"), type="HCC")
general_figures(HCC_obj, paste(out_tag, "HCC", sep="_"), sample_col="Sample")


###### Zhang #####
set.seed(5927)
ZhangMats <- load_ZhangCC_aslist()
out_tag <- "ZhangCC"

seur_obj <- convert_MatrixList2Seurat(ZhangMats, min.cells=20, min.features=500)
infercnv_out <- readRDS("ZhangCC_infercnv_obj.rds")

is.tumourcell <- call_malignant(infercnv_out)
seur_obj <- seur_obj[,colnames(seur_obj) %in% names(is.tumourcell)]
seur_obj@meta.data$CNV_calls <- is.tumourcell[match(colnames(seur_obj), names(is.tumourcell))]

Zhang_tumour_obj <- seur_obj[,seur_obj@meta.data$CNV_calls == "Malignant"]

Zhang_tumour_obj <- seurat_analysis_pipeline(Zhang_tumour_obj)
visualize_CC(Zhang_tumour_obj, out_tag)
gene_figures(Zhang_tumour_obj, out_tag, type="CC")
general_figures(Zhang_tumour_obj, out_tag, sample_col="origin")
saveRDS(Zhang_tumour_obj, "Zhang_clustered.rds")

###### Su #####
set.seed(1920)
seur_obj <- load_SuHCC(min.cells=20, min.features=500)
out_tag <- "SuHCC"

infercnv_out <- readRDS("SuHCC_infercnv_obj.rds")

is.tumourcell <- call_malignant(infercnv_out)
seur_obj <- seur_obj[,colnames(seur_obj) %in% names(is.tumourcell)]
seur_obj@meta.data$CNV_calls <- is.tumourcell[match(colnames(seur_obj), names(is.tumourcell))]

Su_tumour_obj <- seur_obj[,seur_obj@meta.data$CNV_calls == "Malignant"]

Su_tumour_obj <- seurat_analysis_pipeline(Su_tumour_obj)
visualize_CC(Su_tumour_obj, out_tag)
gene_figures(Su_tumour_obj, out_tag, type="HCC")
general_figures(Su_tumour_obj, out_tag, sample_col="patient")

saveRDS(Su_tumour_obj, "Su_clustered.rds")
