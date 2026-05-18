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
	png(paste(outname, "CycleMix_phase.png", sep="_"), width=8, height=8, units="in", res=300)
	print(DimPlot(seur_obj, group.by="CycleMix_phase") + 
		scale_color_manual(values=tmp_col))
	dev.off()
	tmp_col <- CC_col[2:4]; names(tmp_col) <- c("G1", "S", "G2M")
	png(paste(outname, "Seurat_Phase.png", sep="_"), width=8, height=8, units="in", res=300)
	print(DimPlot(seur_obj, group.by="Phase")+ scale_color_manual(values=tmp_col))
	dev.off()

	png(paste(outname, "SeuratPhase_bycluster.png", sep="_"), width=8, height=8, units="in", res=150)
	barplot(table(factor(seur_obj$Phase, levels=c("G1", "S", "G2M")), seur_obj@meta.data$seurat_clusters), col=CC_col[2:4] )
	dev.off()
	
	png(paste(outname, "CycleMix_bycluster.png", sep="_"), width=8, height=8, units="in", res=150)
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
	png(paste(outname, type, "Markers.png", sep="_"), width=14, height=6, units="in", res=100)
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
	png(paste(outname, "UMAP.png", sep="_"), width=6, height=6, units="in", res=100)
	print(DimPlot(seur_obj, group.by="seurat_clusters", label=TRUE))
	dev.off()
	
	png(paste(outname, "SampleUMAP.png", sep="_"), width=6, height=6, units="in", res=100)
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


png( "Figure_All_SampleUMAP.png", width=7, height=6, units="in", res=600)
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

png( "Figure_All_CycleMixUMAP.png", width=7, height=6, units="in", res=600)
print(DimPlot(all, group.by="CycleMix_phase", label=FALSE, cols=CC_col))
dev.off()

CC_col = c("black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("G1", "S", "G2M");
png( "Figure_All_CellCycleUMAP.png", width=7, height=6, units="in", res=600)
print(DimPlot(all, group.by="Phase", label=FALSE, cols=CC_col))
dev.off()


png( "Figure_All_ClusterUMAP.png", width=7, height=6, units="in", res=600)
# 3 = TIC, 2 = Prog, 1=Diff

all@meta.data$seurat_clusters <- factor(all@meta.data$seurat_clusters, levels=c(
		"MaCC_4", "MaCC_0","MaCC_5",
		"Su_5", "Su_1", "Su_6a",
		"Su_9", "Su_8", "Su_3",
		"Su_4", "Su_2", "Su_6b",
		"Zhang_1", "Zhang_7", 
		"Zhang_3", "Zhang_13",
		"Zhang_6", "Zhang_19") )

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
png("Patient_data_DE_heatmap_with_sig_stars.png", width=6, height=12, units="in", res=100)
pheatmap(log2fcs, annotation_row=gene_lab, display_numbers=sig_plot,cluster_rows=FALSE)
dev.off()


#png("Figure_DotPlot_Combined_newgenes.png", width=11, height=6, units="in", res=300)
pdf("Figure_DotPlot_Combined_newgenes.pdf", width=11, height=6)
this_plot <- DotPlot(all, features=c(TIC_genes, Diff_genes, CC_genes, HCC_genes), group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(this_plot)
dev.off()


saveRDS(all, "Combined_tumour_data.rds")
all <- readRDS("Combined_tumour_data.rds")
### Custom Tailored Dotplot ###
gene_list <- read.table("Gene_lists_DotPlot_from_Patricia_TO_send_TA_toimport.csv", header=T, sep=",")
genes <- as.character(gene_list[,1]); names(genes) <- gene_list[,2]
names(genes)[genes %in% genes[duplicated(genes)]] <- "Diffs"
genes <- genes[!duplicated(genes)]
names(genes)[names(genes) %in% c("CCs", "HCCs")] <- "Diffs"

png("Figure_DotPlot_Combined2.png", width=11, height=6, units="in", res=300)
this_plot <- DotPlot(all, features=genes, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(this_plot)
dev.off()

sets <- list(MA1 = c("MaCC_0", "MaCC_4", "MaCC_5"),
		Su1 = c("Su 1", "Su 5", "Su 6a"),
		Su2 = c("Su 4", "Su 2", "Su 6b"),
		Su3 = c("Su 8", "Su 9", "Su 3"),
		Zh1 = c("Zhang_1", "Zhang_7"),
		Zh2 = c("Zhang_3", "Zhang_13"),
		Zh3 = c("Zhang_6", "Zhang_19")
	)

for (a in names(sets)) {
	tmp <- all[, all@meta.data$seurat_clusters %in% sets[[a]]]
	tmp@meta.data$seurat_clusters <- factor(tmp@meta.data$seurat_clusters,  sets[[a]])

	png(paste("Figure", a, "full_DotPlot.png", sep="_"), width=11, height=2.25, units="in", res=300)
	this_plot <- DotPlot(tmp, features=genes, group.by="seurat_clusters", scale=FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	print(this_plot)
	dev.off()
}

genes_trimmed <- genes[genes %in% rownames(all)[rowSums(all@assays$RNA@counts > 0) > 0.05]]
for (a in names(sets)) {
	tmp <- all[, all@meta.data$seurat_clusters %in% sets[[a]]]
	tmp@meta.data$seurat_clusters <- factor(tmp@meta.data$seurat_clusters,  sets[[a]])

	png(paste("Figure", a, "DotPlot_narrow.png", sep="_"), width=11, height=2.25, units="in", res=300)
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

png("Projected.png", width=8, height=8, units="in", res=100)
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


##### New Format: Diff vs Stem for HCC and CC separately #####

Su1_hep_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_6a"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_5"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su1_hep_test2 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_6a"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_1"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su1_hep_test3 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_5"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_1"]) # or Su_5 # Su_6a vs Su1 sig at p < 0.05

Su2_hep_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_3"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_9"]) # Su 3 vs Su 8 sig a p < 10^-20, Su 3 vs Su 9 sig at p < 10-10
Su2_hep_test2 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_3"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_8"]) # Su 3 vs Su 8 sig a p < 10^-20, Su 3 vs Su 9 sig at p < 10-10
Su2_hep_test3 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_8"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_9"]) # Su 3 vs Su 8 sig a p < 10^-20, Su 3 vs Su 9 sig at p < 10-10


Su3_hep_test1 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_6b"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_4"]) # Su6b vs Su2 = ns, Su6b vs Su4 p < 10^-7

Su3_hep_test2 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_2"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_4"]) # Su6b vs Su2 = ns, Su6b vs Su4 p < 10^-7

Su3_hep_test3 <- wilcox.test(
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_6b"],
        all@meta.data$hep_score[all@meta.data$seurat_clusters == "Su_2"]) # Su6b vs Su2 = ns, Su6b vs Su4 p < 10^-7

pdf("Projected_boxes_hep_samples.pdf", width=8, height=10)
hep_samples_only <- all[,all@meta.data$orig.ident %in% c("HCC2", "HCC5", "HCC9")]
tmp <- levels(hep_samples_only@meta.data$seurat_clusters)
tmp <- tmp[tmp %in% hep_samples_only@meta.data$seurat_clusters]
tmp <- tmp[c(2,1,3,5,4,6,7,8,9)]
hep_samples_only@meta.data$seurat_clusters <- factor(hep_samples_only@meta.data$seurat_clusters, levels=tmp)
box_names <- paste(c(rep("HCC2",3), rep("HCC5", 3), rep("HCC9",3)),tmp, sep="_")
par(mfrow=c(2,1))
par(mar=c(6,4,4,1))
loc <- boxplot(hep_samples_only@meta.data$hep_score~hep_samples_only@meta.data$seurat_clusters,
		col=c(Colour_Scheme[["Su Patient 2"]][1:3],
	                Colour_Scheme[["Su Patient 5"]][1:3],
	                Colour_Scheme[["Su Patient 9"]][1:3]),
			names=box_names, outpch=16, outcex=0.5, notch=TRUE,
		 border=c(Colour_Scheme[["Su Patient 2"]][1:3],
                        Colour_Scheme[["Su Patient 5"]][1:3],
                        Colour_Scheme[["Su Patient 9"]][1:3]),
		 main="Hepatocyte", ylab="Hepatocyte Lineage Score", las=2, xlab="")
points(1:ncol(loc$stats),loc$stats[3,], pch=18, col="black")


loc <- boxplot(hep_samples_only@meta.data$stem_score~hep_samples_only@meta.data$seurat_clusters,
		col=c(Colour_Scheme[["Su Patient 2"]][1:3],
	                Colour_Scheme[["Su Patient 5"]][1:3],
	                Colour_Scheme[["Su Patient 9"]][1:3]),
			names=box_names, outpch=16, outcex=0.5, notch=TRUE,
 		border=c(Colour_Scheme[["Su Patient 2"]][1:3],
                        Colour_Scheme[["Su Patient 5"]][1:3],
                        Colour_Scheme[["Su Patient 9"]][1:3]),
		 main="Stemness", ylab="Stem Lineage Score", las=2, xlab="")
points(1:ncol(loc$stats),loc$stats[3,], pch=18, col="black")
dev.off()


pdf("Projected_boxes_chol_samples.pdf", width=8, height=10)
chol_samples_only <- all[,!(all@meta.data$orig.ident %in% c("HCC2", "HCC5", "HCC9"))]
tmp <- levels(chol_samples_only@meta.data$seurat_clusters)
tmp <- tmp[tmp %in% chol_samples_only@meta.data$seurat_clusters]
chol_samples_only@meta.data$seurat_clusters <- factor(chol_samples_only@meta.data$seurat_clusters, levels=tmp)
box_names <- paste(c(rep("Ma1",3), rep("Zh_18", 2), rep("Zh_24",2), rep("Zh_23",2)),tmp, sep="_")
par(mfrow=c(2,1))
par(mar=c(6,4,4,1))
loc <- boxplot(chol_samples_only@meta.data$chol_score~chol_samples_only@meta.data$seurat_clusters,
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

loc <- boxplot(chol_samples_only@meta.data$stem_score~chol_samples_only@meta.data$seurat_clusters,
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


png("Projected_box_hep.png", width=8, height=8, units="in", res=100)
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

png("Projected_box_chol.png", width=8, height=8, units="in", res=100)
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


png("Projected_box_stem.png", width=8, height=8, units="in", res=100)
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





