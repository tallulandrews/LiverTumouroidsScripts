set.seed(8932)
dir = "/home/tandrew6/projects/def-tandrew6/tandrew6/Meritxell"

source("/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Merixtell/MerixtellLoadData.R")
source("/home/tandrew6/projects/def-tandrew6/tandrew6/Meritxell/check_malignant_calls_Ma.R")
require(Seurat)
require(CycleMix)
require(mclust)
require(ggplot2)

#### -------------------------- ######
#### Colour Schemes ####
Colour_Scheme <- list(
	"Su Patient 2"=c("#b2e2e2", "#66c2a4", "#238b45"),
	"Su Patient 5"=c("#cbc9e2", "#9e9ac8", "#6a51a3"),
	"Su Patient 9"=c("#fcae91", "#fb6a4a", "#cb181d"),
	"Ma P16_LCP46"=c("#bdd7e7", "#6baed6", "#2171b5"),
	"Zh 18_Tumor"=c("#fbb4b9", "#f768a1", "#ae017e"),
	"Zh 23_Tumor"=c("#fed98e", "#fe9929", "#cc4c02"),
	"Zh 24_Tumor"=c("#cccccc", "#969696", "#525252")
	)

require("RColorBrewer")
CC_col = c("grey50", "black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("None", "G0", "G1/S", "G2/M");
CC_col = c("black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("None", "G1S", "G2M");
CC_col = c("black", brewer.pal(3, "PuOr")[-2]); names(CC_col) = c("G1", "S", "G2M");

#### -------------------------- ######
########### READ IN DATA ##############

###### Ma ########
Ma_CC <- readRDS("MaCC_clustered.rds")

Ma_CC <- calc_CC(Ma_CC)

Ma_CC_samples <- c("S358_P16_LCP46")
Ma_CC <- Ma_CC[,Ma_CC@meta.data$Sample %in% Ma_CC_samples & Ma_CC@meta.data$seurat_clusters %in% c("0", "4", "5")]
Ma_CC@meta.data$seurat_clusters <- paste("MaCC", Ma_CC@meta.data$seurat_clusters, sep="_")

Ma_CC@meta.data$Sample <- paste("Ma", as.character(Ma_CC@meta.data$Sample))

Ma_CC@meta.data$UMAP1 <- Ma_CC@reductions$umap@cell.embeddings[,1]
Ma_CC@meta.data$UMAP2 <- Ma_CC@reductions$umap@cell.embeddings[,2]

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

Zhang_CC@meta.data$UMAP1 <- Zhang_CC@reductions$umap@cell.embeddings[,1]
Zhang_CC@meta.data$UMAP2 <- Zhang_CC@reductions$umap@cell.embeddings[,2]

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

Su_HCC@meta.data$UMAP1 <- Su_HCC@reductions$umap@cell.embeddings[,1]
Su_HCC@meta.data$UMAP2 <- Su_HCC@reductions$umap@cell.embeddings[,2]

all <- merge(Ma_CC, merge(Zhang_CC, Su_HCC))
all <- ScaleData(all)
all <- Seurat::FindVariableFeatures(all)
all <- RunPCA(all, features=VariableFeatures(all))
all <- Seurat::RunUMAP(all, dims = 1:20, parallel=FALSE)

#### -------------------------- ######

##### Rename Clusters #####

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

#### -------------------------- ######
##### Make Heatmap #####

genes_4 <- read.table("/home/tandrew6/scripts/Meritxell/Genes_from_Heatmap_04_2026.txt",sep="\t")
require(pheatmap)

sub_sample_cells <- c()
for (i in unique(all@meta.data$fancy_names)) {
	n_per_Cell <- all@meta.data$nFeature_RNA[all@meta.data$fancy_names==i]
	names(n_per_Cell) <- colnames(all)[all@meta.data$fancy_names==i]
	sub_sample_cells <- c(sub_sample_cells,names(head(sort(n_per_Cell),10)))
}


genes_4 <- genes_4[genes_4[,1] %in% rownames(all),]
genes_for_heatmap <- genes_4[,1]
for_heatmap <- all@assays$RNA@data[genes_for_heatmap,sub_sample_cells]
colnames(for_heatmap) <- sub_sample_cells

orig_dataset <- as.character(all@meta.data[sub_sample_cells,"fancy_names"])
orig_dataset[grepl("^Ma",orig_dataset)] <- "Ma"
orig_dataset[grepl("^Su",orig_dataset)] <- "Su"
orig_dataset[grepl("^Zhang",orig_dataset)] <- "Zhang"
cell_anno <- data.frame(fancy_names=all@meta.data[sub_sample_cells,"fancy_names"], orig_dataset=orig_dataset)
rownames(cell_anno) <- sub_sample_cells

gene_anno <- data.frame(type=genes_4[,2])
rownames(gene_anno) <- genes_4[,1]

keep_genes <- rowSums(for_heatmap) > 0
for_heatmap <- for_heatmap[keep_genes,]
gene_anno <- data.frame(type=gene_anno[keep_genes,1])
rownames(gene_anno) <- rownames(for_heatmap)

col <- c("TIC"="orange", "HCC"="firebrick", "CC"="forestgreen")
png("test.png", width=8, height=6, units="in", res=150)
pheatmap(for_heatmap, annotation_row = gene_anno,annotation_col = cell_anno, 
		show_rownames=F, show_colnames=F,clustering_method="ward.D2",cluster_rows=FALSE,
		annotation_colors=list(type=c("CC"=as.character(col["CC"]), 
						"HCC"=as.character(col["HCC"]),
						"TIC"=as.character(col["TIC"])),
			fancy_names=c(
				"Ma1_Diff"=as.character(col["CC"]),
				"Ma1_Mid"=as.character(col["CC"]),
				"Ma1_CPCs"=as.character(col["TIC"]),
				"Su2_Diff"=as.character(col["HCC"]),
				"Su2_Mid"=as.character(col["HCC"]),
				"Su2_CPCs"=as.character(col["TIC"]),
				"Su3_Diff"=as.character(col["HCC"]),
				"Su3_Mid"=as.character(col["HCC"]),
				"Su3_CPCs"=as.character(col["TIC"]),
				"Su4_Diff"=as.character(col["HCC"]),
				"Su4_Mid"=as.character(col["HCC"]),
				"Su4_CPCs"=as.character(col["TIC"]),
				"Zhang5_Diff"=as.character(col["CC"]),
				"Zhang5_CPCs"=as.character(col["TIC"]),
				"Zhang7_Diff"=as.character(col["CC"]),
				"Zhang7_CPCs"=as.character(col["TIC"]),
				"Zhang6_Diff"=as.character(col["CC"]),
				"Zhang6_CPCs"=as.character(col["TIC"])),
			orig_dataset=c("Ma"="dodgerblue", "Su"="purple", "Zhang"="grey35")
		)
	)
dev.off()

require(MultiPath)
aggregated <- group_rowmeans(all@assays$RNA@data[genes_for_heatmap,], all@meta.data$fancy_names)

gene_anno <- data.frame(type=genes_4[,2])
rownames(gene_anno) <- genes_4[,1]

cell_anno <- data.frame(type2=c("CC", "CC", "TIC", 
				"HCC", "HCC", "TIC",
				"HCC", "HCC", "TIC",
				"HCC", "HCC", "TIC",
				"CC", "TIC", "CC", "TIC", "CC", "TIC"),
			orig_dataset=c("Ma","Ma","Ma",
					"Su", "Su", "Su",
					"Su", "Su", "Su",
					"Su", "Su", "Su",
					"Zhang","Zhang","Zhang","Zhang","Zhang","Zhang"
					))
rownames(cell_anno) <- colnames(aggregated)

png("Heatmap_lineage_genes.png", width=8, height=6, units="in", res=150)
pheatmap(aggregated, annotation_row = gene_anno,annotation_col = cell_anno,
                show_rownames=F, show_colnames=F,clustering_method="ward.D2",cluster_rows=FALSE,
                annotation_colors=list(
					type=c("CC"=as.character(col["CC"]),
                                                "HCC"=as.character(col["HCC"]),
                                                "TIC"=as.character(col["TIC"])),
					type2=c("CC"=as.character(col["CC"]),
                                                "HCC"=as.character(col["HCC"]),
                                                "TIC"=as.character(col["TIC"])),
					orig_dataset=c("Ma"="dodgerblue", "Su"="purple", "Zhang"="grey35")
                ),
		border_color=NA
        )
dev.off()




#### -------------------------- ######
################# DotPlot ##################

genes_1 <- read.table("/home/tandrew6/scripts/Meritxell/Genes_from_Figure1_04_2026.txt", sep="\t")
genes_2 <- read.table("/home/tandrew6/scripts/Meritxell/Genes_from_ExtFigure3_04_2026.txt", sep="\t")
genes_3 <- read.table("/home/tandrew6/scripts/Meritxell/Genes_from_Figure2d_04_2026.txt",sep="\t")

pdf("MoreDotPlots.pdf", width=15, height=5)

## genes1 ##
these_genes <- genes_1[,1]
tmp_names <- genes_1[,2]
tmp_names[grepl("TIC",tmp_names)] <- "TIC"
names(these_genes) <- tmp_names
these_genes <- these_genes[!duplicated(these_genes)]

require(ggplot2)
DotPlot(all, group.by="fancy_names", features=these_genes)+ggtitle("Figure1 Genes")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## genes2 ##
these_genes <- genes_2[,1]
tmp_names <- genes_2[,2]
tmp_names[grepl("TIC",tmp_names)] <- "TIC"
names(these_genes) <- tmp_names
these_genes <- these_genes[!duplicated(these_genes)]

require(ggplot2)
DotPlot(all, group.by="fancy_names", features=these_genes)+ggtitle("ExtFigure3 Genes")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

pdf("MoreDotPlots2.pdf", width=10, height=5)
## genes3 ##
these_genes <- genes_3[,1]
tmp_names <- genes_3[,2]
tmp_names[grepl("TIC",tmp_names)] <- "TIC"
names(these_genes) <- tmp_names
these_genes <- these_genes[!duplicated(these_genes)]

require(ggplot2)
DotPlot(all, group.by="fancy_names", features=these_genes)+ggtitle("Figure2d Genes")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()










