require("scater")
Camp_QC = readRDS("Camp_QC.rds")
# Shorten names so combination names aren't ridiculously long
my_shortened_celltype <- as.character(pData(Camp_QC)$cell_type1)
my_shortened_celltype[my_shortened_celltype == "definitive endoderm"] <- "de"
my_shortened_celltype[my_shortened_celltype == "hepatic endoderm"] <- "he"
my_shortened_celltype[my_shortened_celltype == "mesenchymal stem cell"] <- "msc"
my_shortened_celltype[my_shortened_celltype == "immature hepatoblast"] <- "ih"
my_shortened_celltype[my_shortened_celltype == "mature hepatocyte"] <- "m-hepato"
my_shortened_celltype[my_shortened_celltype == "erythroblasts"] <- "erythro"
my_shortened_celltype[my_shortened_celltype == "lymphoblasts"] <- "lympho"
my_shortened_celltype[my_shortened_celltype == "fetal hepatocytes"] <- "f-hepato"
my_shortened_celltype[my_shortened_celltype == "stellate"] <- "stellate"
pData(Camp_QC)$cell_type2 <- factor(my_shortened_celltype);
my_shortened_source <- as.character(pData(Camp_QC)$Source)
my_shortened_source[my_shortened_source=="iPSC line TkDA3-4"] <- "iPSC line"
my_shortened_source[my_shortened_source=="Mesenchymal stem cell"] <- "MSC"
my_shortened_source[my_shortened_source=="liver bud"] <- "Liver bud"
pData(Camp_QC)$Source2 <- factor(my_shortened_source);


# Get M3Drop features
M3D_features <- rownames(fData(Camp_QC))[fData(Camp_QC)$M3D_qval < 0.05]

# Heatmap setup
source("~/R-Scripts/heatmap.3.R")
require("RColorBrewer")
heatcolours <- rev(brewer.pal(11, "RdBu"))
col_breaks <- c(-100, seq(-2, 2, length = 10), 100)

# cell annotation colours set up
type_col_palette = colorRampPalette(c("grey50","red","orange","goldenrod","forestgreen","cornflowerblue","navy","purple"))
#type_col = brewer.pal(length(levels(pData(Camp)$cell_type1)), "Set2")
type_col = type_col_palette(length(levels(pData(Camp_QC)$cell_type1)))
#source_col = c("#fb8072","#80b1d3","#fdb462","#b3de69")
source_col = c("#377eb8","#fdc086","#4daf4a","#a65628","black")

CellType <- type_col[pData(Camp_QC)$cell_type1]
Source <- source_col[pData(Camp_QC)$Source]

# Simple heatmap
#heatout <- heatmap.3(exprs(Camp_QC)[featureNames(Camp_QC) %in% M3D_features,], scale="row", ColSideColors=cbind(CellType,Source), ColSideColorsSize=2, breaks=col_breaks, col=heatcolours, hclustfun=function(x){hclust(x, method="ward.D2")})

#### Marker Genes take 1 ######
# Assign genes to group they are a marker of source or cell-type which ever is better
require("M3Drop")

trunc_data <- exprs(Camp_QC)[featureNames(Camp_QC) %in% M3D_features,]
markers_type <- M3DropGetMarkers(trunc_data, pData(Camp_QC)$cell_type1)
markers_source <- M3DropGetMarkers(trunc_data, pData(Camp_QC)$Source)

trunc_data <- trunc_data[sort(rownames(trunc_data)),]
markers_type <- markers_type[sort(rownames(markers_type)),]
markers_source <- markers_source[sort(rownames(markers_source)),]

Gene_assign <- rep("None", times=length(markers_type[,1]))
min_AUC = 0.7
for(i in 1:length(markers_type[,1])) {
	if (markers_type[i,1] > markers_source[i,1] & markers_type[i,1] > min_AUC) {
		Gene_assign[i] <- as.character(markers_type[i,2])
	}
	if (markers_type[i,1] < markers_source[i,1] & markers_source[i,1] > min_AUC) {
#		Gene_assign[i] <- as.character(markers_source[i,2])
		Gene_assign[i] <- "None"
	}
}

assign_col <- c("white", type_col, source_col)
assign_levels <- c("None", levels(pData(Camp_QC)$cell_type1), levels(pData(Camp_QC)$Source))

# Remake heatmap coloring genes
keep <- !(Gene_assign =="None")

toplot <- trunc_data[keep,]
gene_col <- matrix(assign_col[factor(Gene_assign[keep], levels=assign_levels)], nrow=1)
sort_order <- order(Gene_assign[keep], rowMeans(toplot))
toplot<-toplot[sort_order,]
gene_col<- matrix(gene_col[,sort_order], nrow=1)

png("Camp_orig_markers_heatmap.png", width=6, height=6, units="in", res=300)
heatout <- heatmap.3(toplot, scale="row", breaks=col_breaks, col=heatcolours, 
			dendrogram="column", Rowv=FALSE,
			ColSideColors=cbind(CellType,Source), ColSideColorsSize=2, 
			RowSideColors=gene_col, RowSideColorsSize=1,
			hclustfun=function(x){hclust(x, method="ward.D2")}) 
dev.off();

##### Marker Genes take 2 ######
# Allow gene to be marker of >1 group:
# using pROC
# Whole dataset not just M3Drop features

trunc_data <- exprs(Camp_QC); # remove for testing purposes
source("~/R-Scripts/New_Marker_Genes.R")
require("pROC")

combo_labels = paste(pData(Camp_QC)$cell_type2,pData(Camp_QC)$Source2)
newMarkers_both <- New_Complex_Markers(trunc_data, combo_labels, strict_only=FALSE, n_max=length(unique(combo_labels))-1)
Gene_assign <- get_combo_names(newMarkers_both)
# AUC threshold
min_AUC = 0.7
# Eliminate sets of groups with few marker genes.
other_cats <- names(summary(factor(Gene_assign)))[summary(factor(Gene_assign))<=10]
Gene_assign[Gene_assign %in% other_cats] <- "Other"

# Get rid of uninteresting markers & non-markers
keep <- (rowSums(newMarkers_both[,2:(length(newMarkers_both[1,])-2)]) > 0) & newMarkers_both[,1] > min_AUC & Gene_assign != "Other"
toplot <- trunc_data[keep,]
gene_colour_bar = factor(Gene_assign[keep])
gene_colour_palette = rainbow(length(levels(gene_colour_bar)));
sort_order <- order(Gene_assign[keep], rowMeans(toplot))
toplot<- toplot[sort_order,]
gene_colour_bar <- gene_colour_bar[sort_order]

png("Camp_Marker_heatmap.png", width=6, height=6, units="in", res=300)
heatout <- heatmap.3(toplot, scale="row", breaks=col_breaks, col=heatcolours, Rowv=FALSE,
                        ColSideColors=cbind(CellType,Source), ColSideColorsSize=2,
                        RowSideColors=matrix(gene_colour_palette[gene_colour_bar], nrow=1), RowSideColorsSize=1,
                        hclustfun=function(x){hclust(x, method="ward.D2")})
dev.off()

png("Camp_Marker_heatmap_legend.png", width=4, height=6, units="in", res=300)
plot(1:10, col="white", axes=FALSE, xlab="", ylab="")
legend("topleft", levels(gene_colour_bar), fill=gene_colour_palette, bty="n")
legend("bottomright", levels(pData(Camp_QC)$cell_type1), fill=type_col, bty="n")
dev.off()
fData(Camp_QC)$MarkerWeak <- Gene_assign

#### Strict Markers only
trunc_data <- exprs(Camp_QC); # remove for testing
# This should cover cell-type & source specific markers through diff combinations of the intersection tags
newMarkers_both <- New_Complex_Markers(trunc_data, combo_labels, strict_only=TRUE, n_max=length(unique(combo_labels))-1)
Gene_assign <- get_combo_names(newMarkers_both)
# AUC threshold
min_AUC = 0.7
# Eliminate sets of groups with few marker genes.
other_cats <- names(summary(factor(Gene_assign)))[summary(factor(Gene_assign))<=10]
Gene_assign[Gene_assign %in% other_cats] <- "Other"


# Get rid of uninteresting markers & non-markers
keep <- (rowSums(newMarkers_both[,2:(length(newMarkers_both[1,])-2)]) > 0) & newMarkers_both[,1] > min_AUC & Gene_assign != "Other"
toplot <- trunc_data[keep,]
gene_colour_bar = factor(Gene_assign[keep])
gene_colour_palette = rainbow(length(levels(gene_colour_bar)));
sort_order <- order(Gene_assign[keep], rowMeans(toplot))
toplot<- toplot[sort_order,]
gene_colour_bar <- gene_colour_bar[sort_order]

png("Camp_MarkerStrict_heatmap.png", width=6, height=6, units="in", res=300)
heatout <- heatmap.3(toplot, scale="row", breaks=col_breaks, col=heatcolours, Rowv=FALSE,
                        ColSideColors=cbind(CellType,Source), ColSideColorsSize=2,
                        RowSideColors=matrix(gene_colour_palette[gene_colour_bar], nrow=1), RowSideColorsSize=1,
                        hclustfun=function(x){hclust(x, method="ward.D2")})
dev.off()
png("Camp_MarkerStrict_heatmap_legend.png", width=4, height=6, units="in", res=300)
plot(1:10, col="white", axes=FALSE, xlab="", ylab="")
legend("topleft", levels(gene_colour_bar), fill=gene_colour_palette, bty="n")
legend("bottomright", levels(pData(Camp_QC)$cell_type1), fill=type_col, bty="n")
dev.off()

fData(Camp_QC)$MarkerStrict <- Gene_assign

saveRDS(Camp_QC, "Camp_human_w_Markers.rds")
