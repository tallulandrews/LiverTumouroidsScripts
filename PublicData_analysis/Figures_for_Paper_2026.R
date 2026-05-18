set.seed(8932)
dir = "/home/tandrew6/projects/def-tandrew6/tandrew6/Meritxell"

source("/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Merixtell/MerixtellLoadData.R")
source("/home/tandrew6/projects/def-tandrew6/tandrew6/Meritxell/check_malignant_calls_Ma.R")
final_TIC_Diff_markers <- read.table("Final_Markers.txt", header=FALSE)
metabolic_markers <- read.table("B27_components_genesets_Harmonizome_aggregate_gene_list.csv", sep=",", header=TRUE)
colnames(final_TIC_Diff_markers) <- c("Gene", "Cluster")
metabolic_markers <- data.frame(Gene=metabolic_markers[,2], Pathway=metabolic_markers[,1])
Atlas_markers <- read.delim("~/projects/def-tandrew6/tandrew6/UHN/LiverMap2.0/Livemap2.0_Markers.csv", sep=",")
require(Seurat)
require(CycleMix)
require(mclust)
require(ggplot2)

calc_CC <- function(seur_obj, outname="project") {
	set.seed(2901)
	require(RColorBrewer)
	CC_col = c("grey50", "black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("None", "G0", "G1/S", "G2/M");


	# My CC Phase Assignment
	require(CycleMix)
	require(SingleCellExperiment)
	SCE <- as.SingleCellExperiment(seur_obj)
	logcounts(SCE) <- as.matrix(logcounts(SCE))
	rowData(SCE)$feature_symbol <- rownames(seur_obj)
	CC_genes <- HGeneSets$Tirosh
	CC_genes <- CC_genes[CC_genes[,1] %in% rowData(SCE)$feature_symbol,]
	cc_assignments <- CycleMix::classifyCells(SCE, CC_genes)
	seur_obj@meta.data$CycleMix_phase <- cc_assignments$phase
	return(seur_obj)
}

plot_CC <- function(seur_obj, outname="project") {
	CC_col = c("grey50", "black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("None", "G0", "G1/S", "G2/M");
	
	tmp_col <- CC_col[2:4]; names(tmp_col) <- c("None", "G1S", "G2M")
	png(paste(outname, "CycleMix_phase_2026.png", sep="_"), width=8, height=8, units="in", res=300)
	print(DimPlot(seur_obj, group.by="CycleMix_phase") + 
		scale_color_manual(values=tmp_col))
	dev.off()
	tmp_col <- CC_col[2:4]; names(tmp_col) <- c("G1", "S", "G2M")
	png(paste(outname, "Seurat_Phase_2026.png", sep="_"), width=8, height=8, units="in", res=300)
	print(DimPlot(seur_obj, group.by="Phase")+ scale_color_manual(values=tmp_col))
	dev.off()

	png(paste(outname, "SeuratPhase_bycluster_2026.png", sep="_"), width=8, height=8, units="in", res=150)
	barplot(table(factor(seur_obj$Phase, levels=c("G1", "S", "G2M")), seur_obj@meta.data$seurat_clusters), col=CC_col[2:4] )
	dev.off()
	
	png(paste(outname, "CycleMix_bycluster_2026.png", sep="_"), width=8, height=8, units="in", res=150)
	barplot(table(factor(seur_obj$CycleMix_phase, levels=c("None", "G1S", "G2M")), seur_obj@meta.data$seurat_clusters), col=CC_col[2:4] )
	dev.off()
	
}

visualize_key_genes <- function(seur_obj, outname, type=c("CC", "HCC")) {
	if (type == "CC") {
		genes <- final_TIC_Diff_markers[grep("CCA", final_TIC_Diff_markers[,2]),]
		genes2 <- Atlas_markers[Atlas_markers$General.Type=="Cholangiocyte", c("Gene", "Specific.Type")]
		genes <- data.frame(Gene=c(genes[,1], genes2[,1]), Type=c(genes[,2], genes2[,2]))
	} else if (type == "HCC") {
		genes <- final_TIC_Diff_markers[grep("HCC", final_TIC_Diff_markers[,2]),]
		genes <- genes[genes[,2] %in% c("HCC10_Diffs", "HCC10_TICs", "HCC24_TICs", "HCC24_Diffs"),]
		genes2 <- Atlas_markers[Atlas_markers$General.Type=="Hepatocyte", c("Gene", "Specific.Type")]
		genes <- data.frame(Gene=c(genes[,1], genes2[,1]), Type=c(genes[,2], genes2[,2]))
	}
	labelled_genes <- unlist(genes[!duplicated(genes[,1]),1]);
	names(labelled_genes) <- genes[!duplicated(genes[,1]),2];
	png(paste(outname, type, "Markers_2026.png", sep="_"), width=14, height=6, units="in", res=100)
	this_plot <- DotPlot(seur_obj, features=labelled_genes, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	print(this_plot)
	dev.off()
	return(this_plot)
}

test_DE_key_genes <- function(seur_obj, group1, group2, type=c("CC", "HCC"), features=rownames(seur_obj)) {
	DE <- FindMarkers(seur_obj, ident.1=group1, ident.2=group2, min.pct=0, logfc.threshold=-Inf,
				features=intersect(features, rownames(seur_obj)), group.by="seurat_clusters")

	return(DE)
}

general_figures <- function(seur_obj, outname, sample_col="Sample") {
	png(paste(outname, "UMAP_2026.png", sep="_"), width=6, height=6, units="in", res=100)
	print(DimPlot(seur_obj, group.by="seurat_clusters", label=TRUE))
	dev.off()
	
	png(paste(outname, "SampleUMAP_2026.png", sep="_"), width=6, height=6, units="in", res=100)
	print(DimPlot(seur_obj, group.by=sample_col, label=TRUE))
	dev.off()
}

###### Ma ########

#Ma_HCC <- readRDS("MaHCC_clustered.rds")
Ma_CC <- readRDS("MaCC_clustered.rds")

Ma_CC <- calc_CC(Ma_CC)

Ma_CC_samples <- c("S358_P16_LCP46")
Ma_CC <- Ma_CC[,Ma_CC@meta.data$Sample %in% Ma_CC_samples & Ma_CC@meta.data$seurat_clusters %in% c("0", "4", "5")]
Ma_CC@meta.data$seurat_clusters <- paste("MaCC", Ma_CC@meta.data$seurat_clusters, sep="_")

Ma_CC@meta.data$Sample <- paste("Ma", as.character(Ma_CC@meta.data$Sample))
visualize_key_genes(Ma_CC, type="CC", outname="Figure_MaCC")

Ma_CC@meta.data$UMAP1 <- Ma_CC@reductions$umap@cell.embeddings[,1]
Ma_CC@meta.data$UMAP2 <- Ma_CC@reductions$umap@cell.embeddings[,2]

#Ma_DE <- test_DE_key_genes(Ma_CC, group1="MaCC_5", group2="MaCC_4")


###### Zhang ########

Zhang_CC <- readRDS("Zhang_clustered.rds")
Zhang_CC <- calc_CC(Zhang_CC)
Zhang_CC@meta.data$Sample <- sub("_UMI.*gz", "", sub("^GSM.*ICC_", "", Zhang_CC@meta.data$origin))
Zhang_CC_samples <- c("18_Tumor","23_Tumor","24_Tumor1", "24_Tumor2")

Zhang_CC <- Zhang_CC[, Zhang_CC@meta.data$Sample %in% Zhang_CC_samples & 
			Zhang_CC@meta.data$seurat_clusters %in% c("1", "7", "6", "19", "13", "3")]
Zhang_CC@meta.data$Sample <- sub("Tumor.*", "Tumor", Zhang_CC@meta.data$Sample)
Zhang_CC@meta.data$Sample <- paste("Zh", as.character(Zhang_CC@meta.data$Sample))
Zhang_CC@meta.data$Sample <- factor(Zhang_CC@meta.data$Sample)
Zhang_CC@meta.data$seurat_clusters <- factor(Zhang_CC@meta.data$seurat_clusters)
Zhang_CC@meta.data$seurat_clusters <- paste("Zhang", as.character(Zhang_CC@meta.data$seurat_clusters), sep="_")

visualize_key_genes(Zhang_CC, type="CC", outname="Figure_ZhangCC")
visualize_key_genes(Zhang_CC, type="HCC", outname="Figure_ZhangHCC")

Zhang_CC@meta.data$UMAP1 <- Zhang_CC@reductions$umap@cell.embeddings[,1]
Zhang_CC@meta.data$UMAP2 <- Zhang_CC@reductions$umap@cell.embeddings[,2]

#Zhang_DE18 <- test_DE_key_genes(Zhang_CC, group1="Zhang_7", group2="Zhang_1") # Tumour 18
#Zhang_DE23 <- test_DE_key_genes(Zhang_CC, group1="Zhang_19", group2="Zhang_6") # Tumour 23
#Zhang_DE24 <- test_DE_key_genes(Zhang_CC, group1="Zhang_13", group2="Zhang_3") # Tumour 24

###### Su #####
Su_HCC <- readRDS("Su_clustered.rds")
Su_HCC <- calc_CC(Su_HCC)
Su_HCC_samples <- c("Patient 2", "Patient 5", "Patient 9")

Su_HCC <- Su_HCC[,Su_HCC@meta.data$patient %in% Su_HCC_samples & 
		  Su_HCC@meta.data$seurat_clusters %in% c("1", "2", "3", "4", "5", "6", "8", "9")]
Su_HCC@meta.data$seurat_clusters <- as.character(Su_HCC@meta.data$seurat_clusters)
Su_HCC@meta.data$seurat_clusters[Su_HCC@meta.data$seurat_clusters=="6" & Su_HCC@meta.data$patient == "Patient 2"] <- "6a"
Su_HCC@meta.data$seurat_clusters[Su_HCC@meta.data$seurat_clusters=="6" & Su_HCC@meta.data$patient == "Patient 9"] <- "6b"
Su_HCC@meta.data$seurat_clusters <- paste("Su", Su_HCC@meta.data$seurat_clusters, sep="_")
Su_HCC@meta.data$Sample <- paste("Su", Su_HCC@meta.data$patient, sep="_")

visualize_key_genes(Su_HCC, type="HCC", outname="Figure_SuHCC")

Su_HCC@meta.data$UMAP1 <- Su_HCC@reductions$umap@cell.embeddings[,1]
Su_HCC@meta.data$UMAP2 <- Su_HCC@reductions$umap@cell.embeddings[,2]
#Su_DE2 <- test_DE_key_genes(Su_HCC, group1="Su_6a", group2="Su_1") # Patient 2
#Su_DE5 <- test_DE_key_genes(Su_HCC, group1="Su_3", group2="Su_8") # Patient 5
#Su_DE9 <- test_DE_key_genes(Su_HCC, group1="Su_6b", group2="Su_4") # Patient 9

##### Combine Samples for One Map ######
Colour_Scheme <- list(
	"Su Patient 2"=c("#b2e2e2", "#66c2a4", "#238b45"),
	"Su Patient 5"=c("#cbc9e2", "#9e9ac8", "#6a51a3"),
	"Su Patient 9"=c("#fcae91", "#fb6a4a", "#cb181d"),
	"Ma P16_LCP46"=c("#bdd7e7", "#6baed6", "#2171b5"),
	"Zh 18_Tumor"=c("#fbb4b9", "#f768a1", "#ae017e"),
	"Zh 23_Tumor"=c("#fed98e", "#fe9929", "#cc4c02"),
	"Zh 24_Tumor"=c("#cccccc", "#969696", "#525252")
	)

## UMAP Coords
Su_x <- Su_HCC@meta.data$UMAP1 - max(Su_HCC@meta.data$UMAP1)
Ma_x <- (Ma_CC@meta.data$UMAP1-max(Ma_CC@meta.data$UMAP1)) - min(Ma_CC@meta.data$UMAP1)
Zhang_x <- abs(Zhang_CC@meta.data$UMAP1-max(Zhang_CC@meta.data$UMAP1))
rescale_f <- max(abs(Zhang_x), abs(Su_x), abs(Ma_x))
Su_x <- Su_x/rescale_f
Ma_x <- Ma_x/rescale_f
Zhang_x <- Zhang_x/rescale_f

Su_y <- Su_HCC@meta.data$UMAP2
Ma_y <- Ma_CC@meta.data$UMAP2 + max(Zhang_CC@meta.data$UMAP2[Zhang_CC@meta.data$Sample=="Zh 18_Tumor"])+2
Zhang_y <- Zhang_CC@meta.data$UMAP2; Zhang_y[Zhang_CC@meta.data$Sample=="Zh 23_Tumor"] <- -1*Zhang_y[Zhang_CC@meta.data$Sample=="Zh 23_Tumor"]

all <- merge(Ma_CC, merge(Zhang_CC, Su_HCC))
all <- ScaleData(all)
all <- Seurat::FindVariableFeatures(all)
all <- RunPCA(all, features=VariableFeatures(all))
all <- Seurat::RunUMAP(all, dims = 1:20, parallel=FALSE)
all@reductions$umap@cell.embeddings[,1] <- c(Ma_x, Zhang_x, Su_x)
all@reductions$umap@cell.embeddings[,2] <- c(Ma_y, Zhang_y, Su_y)


plot_CC(all, outname="All_UMAP")
general_figures(all, outname="All_Plots")


png( "Figure_All_SampleUMAP_2026.png", width=7, height=6, units="in", res=600)
print(DimPlot(all, group.by="Sample", label=FALSE, cols=c(
		Colour_Scheme[["Ma P16_LCP46"]][3], 
		Colour_Scheme[["Su Patient 2"]][3], 
		Colour_Scheme[["Su Patient 5"]][3], 
		Colour_Scheme[["Su Patient 9"]][3], 
		Colour_Scheme[["Zh 18_Tumor"]][3], 
		Colour_Scheme[["Zh 23_Tumor"]][3], 
		Colour_Scheme[["Zh 24_Tumor"]][3] 
		)
	))
dev.off()

require("RColorBrewer")
CC_col = c("grey50", "black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("None", "G0", "G1/S", "G2/M");
CC_col = c("black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("None", "G1S", "G2M");

png( "Figure_All_CycleMixUMAP_2026.png", width=7, height=6, units="in", res=600)
print(DimPlot(all, group.by="CycleMix_phase", label=FALSE, cols=CC_col))
dev.off()

CC_col = c("black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("G1", "S", "G2M");
png( "Figure_All_CellCycleUMAP_2026.png", width=7, height=6, units="in", res=600)
print(DimPlot(all, group.by="Phase", label=FALSE, cols=CC_col))
dev.off()


png( "Figure_All_ClusterUMAP_2026.png", width=7, height=6, units="in", res=600)
# 3 = TIC, 2 = Prog, 1=Diff

all@meta.data$seurat_clusters <- factor(all@meta.data$seurat_clusters, levels=c( # Updated these levels for new Diff/Mid for Su on Mar 3 2026
		"MaCC_4", "MaCC_0","MaCC_5",
		"Su_1", "Su_5", "Su_6a",
		"Su_8", "Su_9", "Su_3",
		"Su_4", "Su_2", "Su_6b",
		"Zhang_1", "Zhang_7", 
		"Zhang_3", "Zhang_13",
		"Zhang_6", "Zhang_19") )
cluster2fancy <- c("MaCC_4"="Ma1_Diff", "MaCC_0"="Ma1_Mid", "MaCC_5"="Ma1_CPCs",
			"Su_1" = "Su2_Diff", "Su2_5"="Su2_Mid", "Su_6a"="Su2_CPCs",
			"Su_8" = "Su3_Diff", "Su_9"="Su3_Mid", "Su_3"="Su3_CPCs",
			"Su_4"="Su4_Diff","Su_2"="Su4_Mid","Su_6b"="Su4_CPCs",
			"Zhang_1"="Zhang5_Diff",  "Zhang_7"="Zhang5_CPCs",
			"Zhang_3"="Zhang7_Diff", "Zhang_13"="Zhang7_CPCs",
			"Zhang_6"="Zhang6_Diff","Zhang_19"="Zhang6_CPCs")
all@meta.data$fancy_names <- cluster2fancy[all@meta.data$seurat_clusters]
all@meta.data$fancy_names <- factor(all@meta.data$fancy_names, levels=c("Ma1_Diff", "Ma1_Mid", "Ma1_CPCs","Su2_Diff", "Su2_Mid", "Su2_CPCs", "Su3_Diff", "Su3_Mid", "Su3_CPCs", "Su4_Diff", "Su4_Mid", "Su4_CPCs", "Zhang5_Diff", "Zhang5_CPCs", "Zhang7_Diff", "Zhang7_CPCs", "Zhang6_Diff", "Zhang6_CPCs"))

print(DimPlot(all, group.by="seurat_clusters", label=FALSE, cols=c(
                Colour_Scheme[["Ma P16_LCP46"]][1:3],
                Colour_Scheme[["Su Patient 2"]][1:3],
		Colour_Scheme[["Su Patient 5"]][1:3],
		Colour_Scheme[["Su Patient 9"]][1:3],
		Colour_Scheme[["Zh 18_Tumor"]][c(1,3)], 
		Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
		Colour_Scheme[["Zh 23_Tumor"]][c(1,3)]
                )
        ))
dev.off()

tmp <- final_TIC_Diff_markers[grepl("TICs", final_TIC_Diff_markers[,2]),1]
tmp <- tmp[duplicated(tmp)]
detect<- rowMeans(all@assays$RNA@counts[rownames(all) %in% tmp,] > 0)
TIC_genes <- names(detect)[detect > 0.1]

tmp <- final_TIC_Diff_markers[grepl("^CC.*Diffs", final_TIC_Diff_markers[,2]),1]
detect<- rowMeans(all@assays$RNA@counts[rownames(all) %in% tmp,] > 0)
CC_genes <- c("CA9",names(detect)[detect > 0.1])
# Add Liver2.0 Markers
CC_genes <- c(CC_genes, "KRT7", "KRT19", "EPCAM", "CLDN1","MUC1", "MUC5B")

tmp <- final_TIC_Diff_markers[grepl("^HCC.*Diffs", final_TIC_Diff_markers[,2]),1]
detect<- rowMeans(all@assays$RNA@counts[rownames(all) %in% tmp,] > 0)
HCC_genes <- names(detect)[detect > 0.1]
# Add Liver2.0 Markers
HCC_genes <- c(HCC_genes, "ALB", "SDS", "OAT", "CYP2A1", "CYP3A4", "CYP2E1", "ADH1B", "HNF4A")

Diff_genes <- CC_genes[CC_genes %in% HCC_genes]

CC_genes <- CC_genes[!CC_genes %in% Diff_genes]
HCC_genes <- HCC_genes[!HCC_genes %in% Diff_genes]
HCC_genes <- c(HCC_genes, "ADH4", "ABAT", "SLC29A4")
HCC_genes <- HCC_genes[! HCC_genes %in% c("CD44", "GLI4")]
names(CC_genes) <- rep("Chol", length(CC_genes))
names(HCC_genes) <- rep("Hep", length(HCC_genes))
names(Diff_genes) <- rep("Diff", length(Diff_genes))
names(TIC_genes) <- rep("TICs",length(TIC_genes))

# DE
full_set <- c(TIC_genes, Diff_genes, CC_genes, HCC_genes)

Ma_DE <- test_DE_key_genes(Ma_CC, group1="MaCC_5", group2="MaCC_4", features=full_set)
Zhang_DE18 <- test_DE_key_genes(Zhang_CC, group1="Zhang_7", group2="Zhang_1", features=full_set) # Tumour 18
Zhang_DE23 <- test_DE_key_genes(Zhang_CC, group1="Zhang_19", group2="Zhang_6", features=full_set) # Tumour 23
Zhang_DE24 <- test_DE_key_genes(Zhang_CC, group1="Zhang_13", group2="Zhang_3", features=full_set) # Tumour 24
Su_DE2 <- test_DE_key_genes(Su_HCC, group1="Su_6a", group2="Su_1", features=full_set) # Patient 2
Su_DE5 <- test_DE_key_genes(Su_HCC, group1="Su_3", group2="Su_8", features=full_set) # Patient 5
Su_DE9 <- test_DE_key_genes(Su_HCC, group1="Su_6b", group2="Su_4", features=full_set) # Patient 9

Ma_DE <- Ma_DE[match(full_set, rownames(Ma_DE)),]
Su_DE2 <- Su_DE2[match(full_set, rownames(Su_DE2)),]
Su_DE5 <-Su_DE5[match(full_set, rownames(Su_DE5)),]
Su_DE9 <- Su_DE9[match(full_set, rownames(Su_DE9)),]
Zhang_DE18 <-Zhang_DE18[match(full_set, rownames(Zhang_DE18)),]
Zhang_DE23 <-Zhang_DE23[match(full_set, rownames(Zhang_DE23)),]
Zhang_DE24 <-Zhang_DE24[match(full_set, rownames(Zhang_DE24)),]

require(pheatmap)
log2fcs <- cbind(Ma_DE$avg_log2FC, Su_DE2$avg_log2FC, Su_DE5$avg_log2FC, Su_DE9$avg_log2FC, Zhang_DE18$avg_log2FC, Zhang_DE23$avg_log2FC, Zhang_DE24$avg_log2FC)
log2fcs[is.na(log2fcs)] <- 0
sig <- cbind(Ma_DE$p_val_adj, Su_DE2$p_val_adj, Su_DE5$p_val_adj, Su_DE9$p_val_adj, Zhang_DE18$p_val_adj, Zhang_DE23$p_val_adj, Zhang_DE24$p_val_adj)
sig[is.na(sig)]<-1
colnames(log2fcs) <- c("Ma", "Su2", "Su5", "Su9", "Zhang18", "Zhang23", "Zhang24")
rownames(log2fcs) <- full_set
colnames(sig) <- c("Ma", "Su2", "Su5", "Su9", "Zhang18", "Zhang23", "Zhang24")
rownames(sig) <- full_set
gene_lab <- data.frame(type=names(full_set))
rownames(gene_lab) <- rownames(log2fcs)

sig_plot <- matrix("", nrow=nrow(sig), ncol=ncol(sig))
sig_plot[sig < 0.05] <- "*"
sig_plot[sig < 10^-5] <- "**"
sig_plot[sig < 10^-10] <- "***"
png("Patient_data_DE_heatmap_with_sig_stars_2026.png", width=6, height=12, units="in", res=100)
pheatmap(log2fcs, annotation_row=gene_lab, display_numbers=sig_plot,cluster_rows=FALSE)
dev.off()


#png("Figure_DotPlot_Combined_newgenes_2026.png", width=11, height=6, units="in", res=300)
pdf("Figure_DotPlot_Combined_newgenes_2026.pdf", width=11, height=6)
this_plot <- DotPlot(all, features=c(TIC_genes, Diff_genes, CC_genes, HCC_genes), group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(this_plot)
dev.off()
pdf("Figure_DotPlot_newnames_newgenes_2026.pdf", width=11, height=6)
this_plot <- DotPlot(all, features=c(TIC_genes, Diff_genes, CC_genes, HCC_genes), group.by="fancy_names") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(this_plot)
dev.off()


saveRDS(all, "Combined_tumour_data_2026.rds")
all <- readRDS("Combined_tumour_data_2026.rds")
### Custom Tailored Dotplot ###
gene_list <- read.table("Gene_lists_DotPlot_from_Patricia_TO_send_TA_toimport.csv", header=T, sep=",")
genes <- as.character(gene_list[,1]); names(genes) <- gene_list[,2]
names(genes)[genes %in% genes[duplicated(genes)]] <- "Diffs"
genes <- genes[!duplicated(genes)]
names(genes)[names(genes) %in% c("CCs", "HCCs")] <- "Diffs"

png("Figure_DotPlot_Combined2_2026.png", width=11, height=6, units="in", res=300)
this_plot <- DotPlot(all, features=genes, group.by="fancy_names") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(this_plot)
dev.off()

sets <- list(MA1 = c("MaCC_0", "MaCC_4", "MaCC_5"),
		Su1 = c("Su_1", "Su_5", "Su_6a"),
		Su2 = c("Su_4", "Su_2", "Su_6b"),
		Su3 = c("Su_8", "Su_9", "Su_3"),
		Zh1 = c("Zhang_1", "Zhang_7"),
		Zh2 = c("Zhang_3", "Zhang_13"),
		Zh3 = c("Zhang_6", "Zhang_19")
	)

for (a in names(sets)) {
	tmp <- all[, all@meta.data$seurat_clusters %in% sets[[a]]]
	tmp@meta.data$seurat_clusters <- factor(tmp@meta.data$seurat_clusters,  sets[[a]])

	png(paste("Figure", a, "full_DotPlot_2026.png", sep="_"), width=11, height=2.25, units="in", res=300)
	this_plot <- DotPlot(tmp, features=genes, group.by="seurat_clusters", scale=FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	print(this_plot)
	dev.off()
}

genes_trimmed <- genes[genes %in% rownames(all)[rowSums(all@assays$RNA@counts > 0) > 0.05]]
for (a in names(sets)) {
	tmp <- all[, all@meta.data$seurat_clusters %in% sets[[a]]]
	tmp@meta.data$seurat_clusters <- factor(tmp@meta.data$seurat_clusters,  sets[[a]])

	png(paste("Figure", a, "DotPlot_narrow_2026.png", sep="_"), width=11, height=2.25, units="in", res=300)
	this_plot <- DotPlot(tmp, features=genes_trimmed, group.by="seurat_clusters", scale=FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	print(this_plot)
	dev.off()
}


# Figures:
# Main : CycleMix, Clusters


# Project "all" data onto PCs used for lineage scoring then calculate lineage scores for them:
pca <- readRDS("/project/6069895/tandrew6/Meritxell/Lineage_Scoring/Figure1_alt_NewGeneList_Lineage_scoring_rawPCA.rds")
require(Matrix)
norm_expr <- all@assays$RNA@data
my_rot_mat <- pca$x
common_genes <- intersect(rownames(my_rot_mat), rownames(norm_expr))
my_rot_mat <- my_rot_mat[match(common_genes, rownames(my_rot_mat)),]
norm_expr <- norm_expr[match(common_genes, rownames(norm_expr)),]
projection <- t(norm_expr) %*% my_rot_mat

thing <- apply(projection, 2, scale)
chol_score <- (thing[,1]+thing[,2])/2
hep_score <- (-thing[,2]+thing[,3])/2
stem_score <- (thing[,5]-thing[,3])/2

all@meta.data$chol_score = chol_score
all@meta.data$hep_score = hep_score
all@meta.data$stem_score = stem_score

aggregate(all@meta.data$chol_score, by=list(all@meta.data$seurat_clusters), mean)
aggregate(all@meta.data$hep_score, by=list(all@meta.data$seurat_clusters), mean)
aggregate(all@meta.data$stem_score, by=list(all@meta.data$seurat_clusters), mean)

table(all@meta.data$CycleMix_phase, all@meta.data$seurat_clusters)

png("Projected_2026.png", width=8, height=8, units="in", res=300)
plot(all@meta.data$hep_score-all@meta.data$chol_score, all@meta.data$stem_score, pch=16, 
		col=c(
                Colour_Scheme[["Ma P16_LCP46"]][1:3],
                Colour_Scheme[["Su Patient 2"]][1:3],
                Colour_Scheme[["Su Patient 5"]][1:3],
                Colour_Scheme[["Su Patient 9"]][1:3],
                Colour_Scheme[["Zh 18_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 23_Tumor"]][c(1,3)]
                )[all@meta.data$seurat_clusters], xlab="Cholangiocyte ---- Hepatocyte", ylab="Stemness")
dev.off()

x_vals <- all@meta.data$hep_score-all@meta.data$chol_score
y_vals <-  all@meta.data$stem_score
colour_scheme <- c(
                Colour_Scheme[["Ma P16_LCP46"]][1:3],
                Colour_Scheme[["Su Patient 2"]][1:3],
                Colour_Scheme[["Su Patient 5"]][1:3],
                Colour_Scheme[["Su Patient 9"]][1:3],
                Colour_Scheme[["Zh 18_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 23_Tumor"]][c(1,3)]
                )
colour_scheme_greyscale <- c("grey85","grey65","grey35",
				"grey85","grey65","grey35",
				"grey85","grey65","grey35",
				"grey85","grey65","grey35",
				"grey85","grey35",
				"grey85","grey35",
				"grey85","grey35")
is.cc <- all@meta.data$orig.ident %in% c("Ma_tumors","X18T","X23T","X24T1","X24T2")
png("Projected_2026_greyscale.png", width=8, height=8, units="in", res=300)
plot(all@meta.data$hep_score-all@meta.data$chol_score, all@meta.data$stem_score, pch=c(16,1)[is.cc+1], 
		col=colour_scheme_greyscale[all@meta.data$seurat_clusters], 
		xlab="Cholangiocyte ---- Hepatocyte", ylab="Stemness")
legend("topright", c("Diff","Mid","CPC","CC","HCC"), 
		col=c("grey85","grey65","grey35", "grey35","grey35"), 
		pch=c(16,16,16,1,16))
dev.off()
png("Projected_2026_CConly.png", width=8, height=8, units="in", res=300)
cc_cells <- all@meta.data$orig.ident %in% c("Ma_tumors","X18T","X23T","X24T1","X24T2")
plot(all@meta.data$hep_score[cc_cells]-all@meta.data$chol_score[cc_cells], all@meta.data$stem_score[cc_cells], pch=1, cex=0.5,
                col=colour_scheme[all@meta.data$seurat_clusters[cc_cells]], 
		xlab="Cholangiocyte ---- Hepatocyte", ylab="Stemness")
dev.off()

tumour2cluster <- list("Ma1"=c("Ma1_Diff", "Ma1_Mid", "Ma1_CPCs"),
			"Su2"=c("Su2_Diff", "Su2_Mid", "Su2_CPCs"),
			"Su3"=c("Su3_Diff","Su3_Mid","Su3_CPCs"),
			"Su4"=c("Su4_Diff","Su4_Mid","Su4_CPCs"),
			"Zhang5"=c("Zhang5_Diff","Zhang5_CPCs"),
			"Zhang7"=c("Zhang7_Diff","Zhang7_CPCs"),
			"Zhang6"=c("Zhang6_Diff","Zhang6_CPCs"))
for(patient in names(tumour2cluster)) {
	patient_for_files <- gsub(" ", "_", patient)
	png(paste0("Projected_2026_",patient_for_files,"grey.png"), width=5.5, height=5, units="in", res=300)
	par(mar=c(4,4,1,1))
	these_cells <- all@meta.data$fancy_names %in% tumour2cluster[[patient]]
	plot(x_vals[these_cells],y_vals[these_cells],
		pch=16, cex=1,
                #col=colour_scheme[all@meta.data$seurat_clusters[these_cells]],
                col=colour_scheme_greyscale[all@meta.data$seurat_clusters[these_cells]],
                xlab="Cholangiocyte ---- Hepatocyte", ylab="Stemness")
	legend("topright",  levels(factor(all@meta.data$fancy_names[these_cells])),
		pch=16,cex=0.75,
		#col=colour_scheme[sort(unique(as.numeric(all@meta.data$seurat_clusters[these_cells])))]
		col=colour_scheme_greyscale[sort(unique(as.numeric(all@meta.data$seurat_clusters[these_cells])))]
		)
	dev.off()
}


##### New Format: Diff vs Stem for HCC and CC separately #####

Su2_hep_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su2_Diff"],
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su2_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su2_hep_test2 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su2_CPCs"],
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su2_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su2_hep_test3 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su2_Diff"],
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su2_CPCs"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su3_hep_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su3_Diff"],
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su3_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su3_hep_test2 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su3_CPCs"],
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su3_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su3_hep_test3 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su3_CPCs"],
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su3_Diff"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su4_hep_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su4_Diff"],
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su4_Mid"]) # Su 3 vs Su 8 sig a p < 10^-20, Su 3 vs Su 9 sig at p < 10-10
Su4_hep_test2 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su4_CPCs"],
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su4_Mid"]) # Su 3 vs Su 8 sig a p < 10^-20, Su 3 vs Su 9 sig at p < 10-10
Su4_hep_test3 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su4_CPCs"],
        all@meta.data$hep_score[all@meta.data$fancy_names == "Su4_Diff"]) # Su 3 vs Su 8 sig a p < 10^-20, Su 3 vs Su 9 sig at p < 10-10


Su2_stem_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su2_Diff"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su2_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su2_stem_test2 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su2_CPCs"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su2_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su2_stem_test3 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su2_Diff"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su2_CPCs"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su3_stem_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su3_Diff"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su3_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su3_stem_test2 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su3_CPCs"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su3_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su3_stem_test3 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su3_CPCs"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su3_Diff"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su4_stem_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su4_Diff"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su4_Mid"]) # Su 3 vs Su 8 sig a p < 10^-20, Su 3 vs Su 9 sig at p < 10-10
Su4_stem_test2 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su4_CPCs"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su4_Mid"]) # Su 3 vs Su 8 sig a p < 10^-20, Su 3 vs Su 9 sig at p < 10-10
Su4_stem_test3 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su4_CPCs"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Su4_Diff"]) # Su 3 vs Su 8 sig a p < 10^-20, Su 3 vs Su 9 sig at p < 10-10


error_bar <- function(x1, x2, pval, loc, vstep=0.2) {
	pval = signif(pval, digits=1)
	if (pval < 0.05) {
		y = max(loc$stats[4,x1:x2])
		lines(c(x1,x2), c(y+vstep, y+vstep))
#		text(mean(c(x1,x2)), y+vstep, pos=3, paste("p =", pval), cex=0.5)
		if (pval < 0.000005) {
			text(mean(c(x1,x2)), y+vstep, pos=3, "***", cex=0.5)
		} else if (pval < 0.0005) {
			text(mean(c(x1,x2)), y+vstep, pos=3, "**", cex=0.5)
		} else if (pval < 0.05) {
			text(mean(c(x1,x2)), y+vstep, pos=3, "*", cex=0.5)
		}
	}
}

my_wilcox <- function(x,y){
	res <- wilcox.test(x,y)
	out <- list(stat=res$statistic, p.value=round(res$p.value, digits=1))
	return(out)
}

require(ggplot2); require(ggpubr); require(ggsignif)
pdf("Lineage_score_violinplots_SuHCC_2026_pval_equalaxes.pdf", width=8, height=10)
hep_samples_only <- all[,all@meta.data$orig.ident %in% c("HCC2", "HCC5", "HCC9")]
tmp <- levels(hep_samples_only@meta.data$fancy_names)
tmp <- tmp[tmp %in% hep_samples_only@meta.data$fancy_names]
hep_samples_only@meta.data$fancy_names <- factor(hep_samples_only@meta.data$fancy_names, levels=tmp)

colour_scheme <- c(Colour_Scheme[["Su Patient 2"]][1:3],
                        Colour_Scheme[["Su Patient 5"]][1:3],
                        Colour_Scheme[["Su Patient 9"]][1:3])
names(colour_scheme) <- levels(hep_samples_only@meta.data$fancy_names)

hep_plot <- ggplot(hep_samples_only@meta.data, aes(x=fancy_names, y=hep_score, fill=fancy_names)) + 
	geom_violin(width=1.4) +geom_boxplot(width=0.1, color="black", fill="white") +
	geom_signif(comparisons=list(c("Su2_Diff", "Su2_CPCs"), 
					c("Su3_Diff", "Su3_CPCs"),
					c("Su4_Diff", "Su4_CPCs")),
		     annotations = paste("p =",signif(c(Su2_hep_test3$p.value,
							Su3_hep_test3$p.value,
							Su4_hep_test3$p.value), digits=1)),
		    tip_length=0
		    ) +
	scale_fill_manual(values=colour_scheme) + 
	theme_classic()+theme(plot.title=element_text(color="black", face="bold", hjust=0.5)) + 
	labs(title="Hepatocyte Score",x ="Cluster", y = "Hepatocyte Score")

stem_plot <- ggplot(hep_samples_only@meta.data, aes(x=fancy_names, y=stem_score, fill=fancy_names)) + 
	geom_violin(width=1.4) +geom_boxplot(width=0.1, color="black", fill="white")+
	geom_signif(comparisons=list(c("Su2_Diff", "Su2_CPCs"), 
					c("Su3_Diff", "Su3_CPCs"),
					c("Su4_Diff", "Su4_CPCs")),
		     annotations = paste("p =",signif(c(Su2_stem_test3$p.value,
                                                        Su3_stem_test3$p.value,
                                                        Su4_stem_test3$p.value), digits=1)),
		    tip_length=0,y_position=c(0.1,2.7,1.0)
		    ) +
	scale_fill_manual(values=colour_scheme) + 
	theme_classic()+theme(plot.title=element_text(color="black", face="bold", hjust=0.5)) + 
        labs(title="TIC/CPC Score",x ="Cluster", y = "TIC/CPC Score")

ggarrange(stem_plot+ylim(-4.5,3.5), hep_plot+ylim(-2,5.5), ncol=1, nrow=2)
#ggarrange(stem_plot, hep_plot, ncol=1, nrow=2)

dev.off()

pdf("Projected_boxes_hep_samples_2026_stars.pdf", width=8, height=10)
hep_samples_only <- all[,all@meta.data$orig.ident %in% c("HCC2", "HCC5", "HCC9")]
#tmp <- levels(hep_samples_only@meta.data$seurat_clusters)
#tmp <- tmp[tmp %in% hep_samples_only@meta.data$seurat_clusters]
tmp <- levels(hep_samples_only@meta.data$fancy_names)
tmp <- tmp[tmp %in% hep_samples_only@meta.data$fancy_names]
#tmp <- tmp[c(2,1,3,5,4,6,7,8,9)]
hep_samples_only@meta.data$fancy_names <- factor(hep_samples_only@meta.data$fancy_names, levels=tmp)
box_names <- paste(c(rep("HCC2",3), rep("HCC5", 3), rep("HCC9",3)),tmp, sep="_")
par(mfrow=c(2,1))
par(mar=c(6,4,4,1))
loc <- boxplot(hep_samples_only@meta.data$hep_score~hep_samples_only@meta.data$fancy_names,
		col=c(Colour_Scheme[["Su Patient 2"]][1:3],
	                Colour_Scheme[["Su Patient 5"]][1:3],
	                Colour_Scheme[["Su Patient 9"]][1:3]),
			names=box_names, outpch=16, outcex=0.5, notch=TRUE,
		 border=c(Colour_Scheme[["Su Patient 2"]][1:3],
                        Colour_Scheme[["Su Patient 5"]][1:3],
                        Colour_Scheme[["Su Patient 9"]][1:3]),
		 main="Hepatocyte", ylab="Hepatocyte Lineage Score", las=2, xlab="")
points(1:ncol(loc$stats),loc$stats[3,], pch=18, col="black")

### Significance Bars ###
error_bar(1,2,Su2_hep_test1$p.value, loc)
error_bar(2,3,Su2_hep_test2$p.value, loc)
error_bar(1,3,Su2_hep_test3$p.value, loc, vstep=0.6)

error_bar(4,5,Su3_hep_test1$p.value, loc)
error_bar(5,6,Su3_hep_test2$p.value, loc)
error_bar(4,6,Su3_hep_test3$p.value, loc, vstep=0.6)

error_bar(7,8,Su4_hep_test1$p.value, loc)
error_bar(8,9,Su4_hep_test2$p.value, loc)
error_bar(7,9,Su4_hep_test3$p.value, loc, vstep=0.6)

loc <- boxplot(hep_samples_only@meta.data$stem_score~hep_samples_only@meta.data$fancy_names,
		col=c(Colour_Scheme[["Su Patient 2"]][1:3],
	                Colour_Scheme[["Su Patient 5"]][1:3],
	                Colour_Scheme[["Su Patient 9"]][1:3]),
			names=box_names, outpch=16, outcex=0.5, notch=TRUE,
 		border=c(Colour_Scheme[["Su Patient 2"]][1:3],
                        Colour_Scheme[["Su Patient 5"]][1:3],
                        Colour_Scheme[["Su Patient 9"]][1:3]),
		 main="Stemness", ylab="Stem Lineage Score", las=2, xlab="")
points(1:ncol(loc$stats),loc$stats[3,], pch=18, col="black")

error_bar(1,2,Su2_stem_test1$p.value, loc)
error_bar(2,3,Su2_stem_test2$p.value, loc)
error_bar(1,3,Su2_stem_test3$p.value, loc, vstep=0.7)

error_bar(4,5,Su3_stem_test1$p.value, loc)
error_bar(5,6,Su3_stem_test2$p.value, loc)
error_bar(4,6,Su3_stem_test3$p.value, loc, vstep=0.7)

error_bar(7,8,Su4_stem_test1$p.value, loc)
error_bar(8,9,Su4_stem_test2$p.value, loc)
error_bar(7,9,Su4_stem_test3$p.value, loc, vstep=0.7)

dev.off()

# Chol boxplots

Ma1_chol_test1 <- wilcox.test(
        all@meta.data$chol_score[all@meta.data$fancy_names == "Ma1_Diff"],
        all@meta.data$chol_score[all@meta.data$fancy_names == "Ma1_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Ma1_chol_test2 <- wilcox.test(
        all@meta.data$chol_score[all@meta.data$fancy_names == "Ma1_CPCs"],
        all@meta.data$chol_score[all@meta.data$fancy_names == "Ma1_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Ma1_chol_test3 <- wilcox.test(
        all@meta.data$chol_score[all@meta.data$fancy_names == "Ma1_Diff"],
        all@meta.data$chol_score[all@meta.data$fancy_names == "Ma1_CPCs"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Zhang5_chol_test1 <- wilcox.test(
        all@meta.data$chol_score[all@meta.data$fancy_names == "Zhang5_Diff"],
        all@meta.data$chol_score[all@meta.data$fancy_names == "Zhang5_CPCs"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Zhang7_chol_test1 <- wilcox.test(
        all@meta.data$chol_score[all@meta.data$fancy_names == "Zhang7_CPCs"],
        all@meta.data$chol_score[all@meta.data$fancy_names == "Zhang7_Diff"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Zhang6_chol_test1 <- wilcox.test(
        all@meta.data$chol_score[all@meta.data$fancy_names == "Zhang6_CPCs"],
        all@meta.data$chol_score[all@meta.data$fancy_names == "Zhang6_Diff"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05


Ma1_stem_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Ma1_Diff"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Ma1_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Ma1_stem_test2 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Ma1_CPCs"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Ma1_Mid"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Ma1_stem_test3 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Ma1_Diff"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Ma1_CPCs"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Zhang5_stem_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Zhang5_Diff"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Zhang5_CPCs"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Zhang7_stem_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Zhang7_CPCs"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Zhang7_Diff"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Zhang6_stem_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$fancy_names == "Zhang6_CPCs"],
        all@meta.data$stem_score[all@meta.data$fancy_names == "Zhang6_Diff"])

chol_samples_only <- all[,all@meta.data$orig.ident %in% c("Ma_tumors", "X18T", "X23T", "X24T1", "X24T2")]
pdf("Lineage_score_violinplots_MaZhangCC_2026_pval_equalaxes.pdf", width=10, height=10)
tmp <- levels(chol_samples_only@meta.data$fancy_names)
tmp <- tmp[tmp %in% chol_samples_only@meta.data$fancy_names]
chol_samples_only@meta.data$fancy_names <- factor(chol_samples_only@meta.data$fancy_names, levels=tmp)

colour_scheme <- c(Colour_Scheme[["Ma P16_LCP46"]][1:3],
                        Colour_Scheme[["Zh 18_Tumor"]][c(1,3)],
                        Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
                        Colour_Scheme[["Zh 23_Tumor"]][c(1,3)])
names(colour_scheme) <- levels(chol_samples_only@meta.data$fancy_names)

chol_plot <- ggplot(chol_samples_only@meta.data, aes(x=fancy_names, y=chol_score, fill=fancy_names)) +
        geom_violin(width=1.4) +geom_boxplot(width=0.1, color="black", fill="white") +
        geom_signif(comparisons=list(c("Ma1_Diff", "Ma1_CPCs"),
                                        c("Zhang5_Diff", "Zhang5_CPCs"),
                                        c("Zhang7_Diff", "Zhang7_CPCs"),
					c("Zhang6_Diff", "Zhang6_CPCs")),
                     annotations = paste("p =",signif(c(Ma1_chol_test3$p.value,
                                                        Zhang5_chol_test1$p.value,
                                                        Zhang7_chol_test1$p.value,
							Zhang6_chol_test1$p.value), digits=1)),
                    tip_length=0, y_position=c(3.2,2.1,1.9, 0.5)
                    ) +
        scale_fill_manual(values=colour_scheme) +
        theme_classic()+theme(plot.title=element_text(color="black", face="bold", hjust=0.5)) +
        labs(title="Cholangiocyte Score",x ="Cluster", y = "Cholangiocyte Score")

stem_plot <- ggplot(chol_samples_only@meta.data, aes(x=fancy_names, y=stem_score, fill=fancy_names)) +
        geom_violin(width=1.4) +geom_boxplot(width=0.1, color="black", fill="white")+
        geom_signif(comparisons=list(c("Ma1_Diff", "Ma1_CPCs"),
                                        c("Zhang5_Diff", "Zhang5_CPCs"),
                                        c("Zhang7_Diff", "Zhang7_CPCs"),
                                        c("Zhang6_Diff", "Zhang6_CPCs")),
                     annotations = paste("p =",signif(c(Ma1_stem_test3$p.value,
                                                        Zhang5_stem_test1$p.value,
                                                        Zhang7_stem_test1$p.value,
							Zhang6_stem_test1$p.value), digits=1)),
                    tip_length=0, y_position=c(1.9,2.1,3.1, 2.4)
                    ) +
        scale_fill_manual(values=colour_scheme) +
        theme_classic()+theme(plot.title=element_text(color="black", face="bold", hjust=0.5)) +
        labs(title="TIC/CPC Score",x ="Cluster", y = "TIC/CPC Score")

ggarrange(stem_plot+ylim(-4.5,3.5), chol_plot+ylim(-2,5.5), ncol=1, nrow=2)

dev.off()


pdf("Projected_boxes_chol_samples_2026_stars.pdf", width=8, height=10)
chol_samples_only <- all[,!(all@meta.data$orig.ident %in% c("HCC2", "HCC5", "HCC9"))]
tmp <- levels(chol_samples_only@meta.data$fancy_names)
tmp <- tmp[tmp %in% chol_samples_only@meta.data$fancy_names]
chol_samples_only@meta.data$fancy_names <- factor(chol_samples_only@meta.data$fancy_names, levels=tmp)
box_names <- paste(c(rep("Ma1",3), rep("Zh_18", 2), rep("Zh_24",2), rep("Zh_23",2)),tmp, sep="_")
par(mfrow=c(2,1))
par(mar=c(6,4,4,1))
loc <- boxplot(chol_samples_only@meta.data$chol_score~chol_samples_only@meta.data$fancy_names,
                col=c(
		Colour_Scheme[["Ma P16_LCP46"]][1:3],
                Colour_Scheme[["Zh 18_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 23_Tumor"]][c(1,3)]),
                        names=box_names, outpch=16, outcex=0.5, notch=TRUE,
		border=c(
		Colour_Scheme[["Ma P16_LCP46"]][1:3],
                Colour_Scheme[["Zh 18_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 23_Tumor"]][c(1,3)]),
                 main="Cholangiocyte", ylab="Cholangiocyte Lineage Score", las=2, xlab="")

points(1:ncol(loc$stats),loc$stats[3,], pch=18, col="black")

error_bar(1,2,Ma1_chol_test1$p.value, loc)
error_bar(2,3,Ma1_chol_test2$p.value, loc, vstep=0.1)
error_bar(1,3,Ma1_chol_test3$p.value, loc, vstep=0.6)

error_bar(4,5,Zhang5_chol_test1$p.value, loc)
error_bar(6,7,Zhang7_chol_test1$p.value, loc)
error_bar(8,9,Zhang6_chol_test1$p.value, loc)

loc <- boxplot(chol_samples_only@meta.data$stem_score~chol_samples_only@meta.data$fancy_names,
                col=c(
		Colour_Scheme[["Ma P16_LCP46"]][1:3],
                Colour_Scheme[["Zh 18_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 23_Tumor"]][c(1,3)]),
                        names=box_names, outpch=16, outcex=0.5, notch=TRUE,
		border=c(
		Colour_Scheme[["Ma P16_LCP46"]][1:3],
                Colour_Scheme[["Zh 18_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 23_Tumor"]][c(1,3)]),
                 main="Stemness", ylab="Stem Lineage Score", las=2, xlab="")
points(1:ncol(loc$stats),loc$stats[3,], pch=18, col="black")

error_bar(1,2,Ma1_stem_test1$p.value, loc)
error_bar(2,3,Ma1_stem_test2$p.value, loc)
error_bar(1,3,Ma1_stem_test3$p.value, loc, vstep=0.6)

error_bar(4,5,Zhang5_stem_test1$p.value, loc)
error_bar(6,7,Zhang7_stem_test1$p.value, loc)
error_bar(8,9,Zhang6_stem_test1$p.value, loc)

dev.off()




# Calculate Significance. 
Ma1_test1 <- wilcox.test(
	all@meta.data$hep_score[all@meta.data$seurat_clusters == "MaCC_5"],  
	all@meta.data$hep_score[all@meta.data$seurat_clusters == "MaCC_4"]) # or Ma_0 # Not-sig any combo


Su1_test1 <- wilcox.test(
	all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_6a"],  
	all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_5"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05


Su2_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_3"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_9"]) # Su 3 vs Su 8 sig a p < 10^-20, Su 3 vs Su 9 sig at p < 10-10


Su3_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_6b"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_4"]) # Su6b vs Su2 = ns, Su6b vs Su4 p < 10^-7


Zhang1_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Zhang_7"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Zhang_1"]) #  p < 10^-7

Zhang2_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Zhang_13"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Zhang_3"]) # p < 10^-35

Zhang3_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Zhang_19"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Zhang_6"]) # ns


png("Projected_box_hep_2026.png", width=8, height=8, units="in", res=100)
loc <- boxplot(all@meta.data$hep_score~all@meta.data$seurat_clusters, col=c(
                Colour_Scheme[["Ma P16_LCP46"]][1:3],
                Colour_Scheme[["Su Patient 2"]][1:3],
                Colour_Scheme[["Su Patient 5"]][1:3],
                Colour_Scheme[["Su Patient 9"]][1:3],
                Colour_Scheme[["Zh 18_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 23_Tumor"]][c(1,3)]
                ), main="", ylab="Hepatocyte Lineage Score", las=2, xlab="")
#lines(c(1,3), c(0.75, 0.75)) # Ma
#text(2,0.75, "ns", pos=3)
lines(c(5,6), c(4.5, 4.5)) # Su 1
text(5.5, 4.5, "p < 0.05", pos=3) 
lines(x=c(7,7,9,9), y=c(0.8, -0.5, -0.5,-0.3)) # Su2
lines(x=c(8,8,9), y=c(1.2, -0.5, -0.5))
text(8, -0.5, "p < 1e-10", pos=1)
lines(x=c(10,10,12,12), y=c(4.5, 4.8, 4.8, 3.8)) #Su3
text(11, 4.8, "p < 1e-7", pos=3)
#lines(c(13,14), c(0.5, 0.5)) # Zhang1
#text(13.5, 0.5, "p < 1e-7", pos=3)
#lines(c(15,16), c(0.75, 0.75)) # Zhang2
#text(15.5, 0.75, "p < 1e-35", pos=3)
#lines(c(17,18), c(1, 1)) # Zhang2
#text(17.5, 1, "ns", pos=3)

mtext("Ma 1", 3, 1, at=2)
mtext("Su 2", 3, 1, at=5)
mtext("Su 3", 3, 1, at=8)
mtext("Su 4", 3, 1, at=11)
mtext("Zhang 5", 3, 2, at=13.5)
mtext("Zhang 6", 3, 1, at=15.5)
mtext("Zhang 7", 3, 2, at=17.5)
abline(v=c(3.5, 6.5, 9.5, 12.5, 14.5, 16.5), lty=2, col="grey85")

dev.off()

# Calculate Significance. 
Ma1_test1 <- wilcox.test(
	all@meta.data$chol_score[all@meta.data$seurat_clusters == "MaCC_5"],  
	all@meta.data$chol_score[all@meta.data$seurat_clusters == "MaCC_0"]) # Ma5 vs Ma4 p < 10^-5, Ma5 vs Ma0 p < 0.05


Zhang1_test1 <- wilcox.test(
        all@meta.data$chol_score[all@meta.data$seurat_clusters == "Zhang_7"],
        all@meta.data$chol_score[all@meta.data$seurat_clusters == "Zhang_1"]) #  p < 10^-7

Zhang2_test1 <- wilcox.test(
        all@meta.data$chol_score[all@meta.data$seurat_clusters == "Zhang_13"],
        all@meta.data$chol_score[all@meta.data$seurat_clusters == "Zhang_3"]) # p < 10^-35

Zhang3_test1 <- wilcox.test(
        all@meta.data$chol_score[all@meta.data$seurat_clusters == "Zhang_19"],
        all@meta.data$chol_score[all@meta.data$seurat_clusters == "Zhang_6"]) # ns

png("Projected_box_chol_2026.png", width=8, height=8, units="in", res=100)
boxplot(all@meta.data$chol_score~all@meta.data$seurat_clusters, col=c(
                Colour_Scheme[["Ma P16_LCP46"]][1:3],
                Colour_Scheme[["Su Patient 2"]][1:3],
                Colour_Scheme[["Su Patient 5"]][1:3],
                Colour_Scheme[["Su Patient 9"]][1:3],
                Colour_Scheme[["Zh 18_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 23_Tumor"]][c(1,3)]
                ), main="", ylab="Cholangiocyte Lineage Score", las=2, xlab="")
lines(c(1,1,3,3), c(1.5,2.5, 2.5,1.5)) # Ma1
text(2, 2.5, "p < 1e-5", pos=3)
lines(c(2,2,3,3), c(1.5,1.85, 1.85,1.5)) # Ma1
text(2.5, 1.85, "p < 0.05", pos=3)
#lines(c(4,4,6,6), c(-2, -1.5, -1.5, -2)) # Su1
#text(5, -1.5, "p < 1e-4", pos=3)
lines(c(13,13,14,14), c(1.5, 2, 2, 1.25)) # Zhang 1
text(13.5, 2, "p < 10^-35", pos=3)
lines(c(15,15,16,16), c(0.25, 1, 1, 0.5)) # Zhang 2
text(15.5, 1, "p < 10^-60", pos=3)
lines(c(17,17,18,18), c(-0, 0.35, 0.35, -0)) # Zhang 3
text(17.5, 0.35, "p < 10^-35", pos=3)

mtext("Ma 1", 3, 1, at=2)
mtext("Su 2", 3, 1, at=5)
mtext("Su 3", 3, 1, at=8)
mtext("Su 4", 3, 1, at=11)
mtext("Zhang 5", 3, 2, at=13.5)
mtext("Zhang 6", 3, 1, at=15.5)
mtext("Zhang 7", 3, 2, at=17.5)
abline(v=c(3.5, 6.5, 9.5, 12.5, 14.5, 16.5), lty=2, col="grey85")

dev.off()




# Calculate Significance. 
Ma1_test1 <- wilcox.test( 
	all@meta.data$stem_score[all@meta.data$seurat_clusters == "MaCC_5"],  
	all@meta.data$stem_score[all@meta.data$seurat_clusters == "MaCC_0"]) # 10-10


Su1_test1 <- wilcox.test(
	all@meta.data$stem_score[all@meta.data$seurat_clusters == "Su_6a"],  
	all@meta.data$stem_score[all@meta.data$seurat_clusters == "Su_1"]) # 10-12


Su2_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$seurat_clusters == "Su_3"],
        all@meta.data$stem_score[all@meta.data$seurat_clusters == "Su_8"]) # 10-10


Su3_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$seurat_clusters == "Su_6b"],
        all@meta.data$stem_score[all@meta.data$seurat_clusters == "Su_2"]) # 10-17

Zhang1_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$seurat_clusters == "Zhang_7"],
        all@meta.data$stem_score[all@meta.data$seurat_clusters == "Zhang_1"]) #  p < 10^-100

Zhang2_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$seurat_clusters == "Zhang_13"],
        all@meta.data$stem_score[all@meta.data$seurat_clusters == "Zhang_3"]) # p < 10^-63

Zhang3_test1 <- wilcox.test(
        all@meta.data$stem_score[all@meta.data$seurat_clusters == "Zhang_19"],
        all@meta.data$stem_score[all@meta.data$seurat_clusters == "Zhang_6"]) # p < 10-100


png("Projected_box_stem_2026.png", width=8, height=8, units="in", res=100)
boxplot(all@meta.data$stem_score~all@meta.data$seurat_clusters, col=c(
                Colour_Scheme[["Ma P16_LCP46"]][1:3],
                Colour_Scheme[["Su Patient 2"]][1:3],
                Colour_Scheme[["Su Patient 5"]][1:3],
                Colour_Scheme[["Su Patient 9"]][1:3],
                Colour_Scheme[["Zh 18_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 24_Tumor"]][c(1,3)],
                Colour_Scheme[["Zh 23_Tumor"]][c(1,3)]
                ), main="", ylab="Stemness Lineage Score", las=2, xlab="")


lines(c(1,3), c(0.5,0.5)) 
text(2, 0.5, "p < 1e-10", pos=3)
lines(c(4,6), c(-3.5,-3.5)) 
text(5, -3.5, "p < 1e-12", pos=1)
lines(c(7,9), c(2,2)) 
text(8, 2, "p < 1e-10", pos=3)
lines(c(10,12), c(-0.5,-0.5)) 
text(11, -0.5, "p < 1e-10", pos=3)
lines(c(13,14), c(1.5,1.5)) 
text(13.5, 1.5, "p < 1e-100,", pos=3)
lines(c(15,16), c(2,2)) 
text(15.5, 2, "p < 1e-63,", pos=3)
lines(c(17,18), c(1.5,1.5)) 
text(17.5, 1.5, "p < 1e-100,", pos=3)

mtext("Ma 1", 3, 1, at=2)
mtext("Su 2", 3, 1, at=5)
mtext("Su 3", 3, 1, at=8)
mtext("Su 4", 3, 1, at=11)
mtext("Zhang 5", 3, 2, at=13.5)
mtext("Zhang 6", 3, 1, at=15.5)
mtext("Zhang 7", 3, 2, at=17.5)
abline(v=c(3.5, 6.5, 9.5, 12.5, 14.5, 16.5), lty=2, col="grey85")
dev.off()





