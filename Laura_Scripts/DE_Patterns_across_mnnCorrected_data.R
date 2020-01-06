require("scater")
require("RColorBrewer")


### Add Laura Liver
require("M3Drop")
require("scater")
require("matrixStats")
require("RColorBrewer")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/DiffExpr/DE_functions.R")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/CellProfiles/0_Cluster_Profiles_Functions.R")

CCA1 <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/CCA1_SC3.rds")
HCC6 <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/HCC6_SC3.rds")
HCC10 <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/HCC10_SC3.rds")

exprs(CCA1) <- get_log_norm(CCA1)
exprs(HCC6) <- get_log_norm(HCC6)
exprs(HCC10) <- get_log_norm(HCC10)

pData(CCA1)$cell_type1 <- paste("CCA1-",pData(CCA1)$sc3_6_clusters, sep="")
pData(HCC6)$cell_type1 <- paste("HCC6-",pData(HCC6)$sc3_6_clusters, sep="")
pData(HCC10)$cell_type1 <- paste("HCC10-",pData(HCC10)$sc3_6_clusters, sep="")

map = read.table("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/Hsap_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
         new = as.character(map[match(x, map[,1]),2])
         new[is.na(new)] = as.character(x[is.na(new)])
         new[duplicated(new)] = x[duplicated(new)]
         return(new)
}
fData(CCA1)$feature_symbol <- ensg2symbol(rownames(CCA1))
fData(HCC6)$feature_symbol <- ensg2symbol(rownames(HCC6))
fData(HCC10)$feature_symbol <- ensg2symbol(rownames(HCC10))

### Laura ###
# my package
require("CellTypeProfiles")
pData(CCA1)$cell_type1 <- factor(pData(CCA1)$cell_type1)
pData(HCC6)$cell_type1 <- factor(pData(HCC6)$cell_type1)
pData(HCC10)$cell_type1 <- factor(pData(HCC10)$cell_type1)

truth_cca1 <- colnames(CCA1_profile$profiles)
truth_hcc6 <- colnames(HCC6_profile$profiles)
truth_hcc10 <- colnames(HCC10_profile$profiles)

###### Laura Only at single-cell level correction ######
### mnn - Laura only ###
CCA1_profile = get_cluster_profile(exprs(CCA1), pData(CCA1)$cell_type1, is.log=2, out.log=2, feature_selection=marker.features)
HCC6_profile = get_cluster_profile(exprs(HCC6), pData(HCC6)$cell_type1, is.log=2, out.log=2, feature_selection=marker.features)
HCC10_profile = get_cluster_profile(exprs(HCC10), pData(HCC10)$cell_type1, is.log=2, out.log=2, feature_selection=marker.features)

markers = list();
markers$CCA1 <- CCA1_profile$is.feature
markers$HCC6 <- HCC6_profile$is.feature
markers$HCC10 <- HCC10_profile$is.feature

cca1_marks <- rownames(markers$CCA1)[rowSums(markers$CCA1) > 0]
hcc6_marks <- rownames(markers$HCC6)[rowSums(markers$HCC6) > 0]
hcc10_marks <- rownames(markers$HCC10)[rowSums(markers$HCC10) > 0]

features <- unique(c( cca1_marks, hcc6_marks, hcc10_marks ))

# Colours
Laura_obj <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ForReviewPaper_stuff.rds")
dataset_cols <- Laura_obj$colours
dataset_pch <- c(1, 18, 7)
require("RColorBrewer")
clust_cols <- brewer.pal(6, "Set2")
prolif_cols = c("black","red")

## MNN Correct ##
source("../4_Multi_Datasets_MNN.R")
require("scater")
source("/lustre/scratch117/cellgen/team218/TA/R-packages/new_scran/scran/R/mnnCorrect.R")
source("/lustre/scratch117/cellgen/team218/TA/R-packages/new_scran/scran/R/utils.R")
require("Matrix")
require("FNN")

laura_mnn_out <- mnnCorrect(list(CCA1=exprs(CCA1), HCC6=exprs(HCC6), HCC10=exprs(HCC10)), hvg.genes=which(rownames(exprs(CCA1)) %in% features))
Combined_corrected <- cbind(laura_mnn_out$corrected$CCA1, laura_mnn_out$corrected$HCC6, laura_mnn_out$corrected$HCC10)
saveRDS(Combined_corrected, file="/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/Laura_MNN_stuff.rds")

## CTP ##

profile_List = list(CCA1=CCA1_profile, HCC6=HCC6_profile, HCC10=HCC10_profile)

#matches = combine_and_match_clusters(profile_List, suppress.plot=FALSE)
matches = combine_and_match_clusters(profile_List, multihit=TRUE, suppress.plot=FALSE, sig.threshold=0.00001, CI.level=0.50)
corrected_profiles = glm_of_matches(matches)

source("~/R-Scripts/Ensembl_Stuff.R")
map_symbol_ensg(names(which(p.adjust(corrected_profiles$effect_p_values[3,], method="fdr") < 0.05)), is.org="Hsap", is.name="ensg")

heatout <- cluster_profile_heatmap(corrected_profiles$corrected_profiles, matches)

cca1_ctp_norm <- correct_sng_cells(CCA1_profile$norm_mat, "CCA1", corrected_profiles)
hcc6_ctp_norm <- correct_sng_cells(HCC6_profile$norm_mat, "HCC6", corrected_profiles)
hcc10_ctp_norm <- correct_sng_cells(HCC10_profile$norm_mat, "HCC10", corrected_profiles)

cca1_ctp_norm2 <- correct_sng_cells(CCA1_profile$norm_mat, "CCA1", corrected_profiles, allow.negatives=TRUE)
hcc6_ctp_norm2 <- correct_sng_cells(HCC6_profile$norm_mat, "HCC6", corrected_profiles, allow.negatives=TRUE)
hcc10_ctp_norm2 <- correct_sng_cells(HCC10_profile$norm_mat, "HCC10", corrected_profiles, allow.negatives=TRUE)

Combined_ctp_norm <- cbind(cca1_ctp_norm, hcc6_ctp_norm, hcc10_ctp_norm);

Laura_CTP_stuff <- list(profiles=profile_List, matches=matches, glm_out=corrected_profiles, sng_cell_norm=Combined_ctp_norm);
saveRDS(Laura_CTP_stuff, file="/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/Laura_CTP_stuff.rds")

# Merge
Combined_raw <- cbind(exprs(CCA1), exprs(HCC6), exprs(HCC10))
dataset <- rep(c(1,2,3), times=c(ncol(laura_mnn_out$corrected$CCA1),ncol(laura_mnn_out$corrected$HCC6), ncol(laura_mnn_out$corrected$HCC10)))
clust_ID <- c(CCA1$sc3_6_clusters, HCC6$sc3_6_clusters, HCC10$sc3_6_clusters)
prolif <- c(CCA1$sc3_6_clusters %in% c("4","1"), HCC6$sc3_6_clusters ==4,  HCC10$sc3_6_clusters==1 )

Laura_Merge_Stuff = list(combined_raw=Combined_raw, dataset_lab=dataset, clust_ID= clust_ID, prolif=prolif)
saveRDS(Laura_Merge_Stuff, file="/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/Laura_Merge_stuff.rds")

# tSNE
library("Rtsne")
set.seed(123)
tsne_corrected <- Rtsne(t(Combined_corrected[rownames(Combined_corrected) %in% features,]))
png("Laura_lines_tsne_mnn.png", width=6, height=6, units="in", res=300)
plot(tsne_corrected$Y[,1], tsne_corrected$Y[,2], col=dataset_cols[dataset], pch=16, xlab="Component 1", ylab="Component 2")
legend("topleft",c("L1","L2","L3"), pch=16, col=dataset_cols, bty="n")
dev.off()

set.seed(123)
tsne_ctp_norm <- Rtsne(t(Combined_ctp_norm[rownames(Combined_ctp_norm) %in% features,]))
png("Laura_lines_tsne_ctp.png", width=6, height=6, units="in", res=300)
plot(tsne_ctp_norm$Y[,1], tsne_ctp_norm$Y[,2], col=dataset_cols[dataset], pch=16, xlab="Component 1", ylab="Component 2")
legend("topleft",c("L1","L2","L3"), pch=16, col=dataset_cols, bty="n")
dev.off()

tsne_raw <- Rtsne(t(Combined_raw[rownames(Combined_raw) %in% features,]))
png("Laura_lines_tsne_raw.png", width=6, height=6, units="in", res=300)
plot(tsne_raw$Y[,1], tsne_raw$Y[,2], col=dataset_cols[dataset], pch=16, xlab="Component 1", ylab="Component 2")
legend("topleft",c("L1","L2","L3"), pch=16, col=dataset_cols, bty="n")
dev.off()
plot(tsne_raw$Y[,1], tsne_raw$Y[,2], col=clust_cols[clust_ID], pch=dataset_pch[dataset], xlab="Component 1", ylab="Component 2")


# PCA
pca_corrected <- prcomp(Combined_corrected[rownames(Combined_corrected) %in% features,])
png("Laura_lines_pca_mnn.png", width=6, height=6, units="in", res=300)
plot(pca_corrected$rotation[,1], pca_corrected$rotation[,2], col=dataset_cols[dataset], pch=16, xlab="Component 1", ylab="Component 2")
legend("topleft",c("L1","L2","L3"), pch=16, col=dataset_cols, bty="n")
dev.off()

pca_ctp_norm <- prcomp(Combined_ctp_norm[rownames(Combined_ctp_norm) %in% features,])
png("Laura_lines_pca_ctp.png", width=6, height=6, units="in", res=300)
plot(pca_ctp_norm$rotation[,1], pca_ctp_norm$rotation[,2], col=dataset_cols[dataset], pch=16, xlab="Component 1", ylab="Component 2")
legend("topleft",c("L1","L2","L3"), pch=16, col=dataset_cols, bty="n")
dev.off()

pca_raw <- prcomp(Combined_raw[rownames(Combined_raw) %in% features,])
png("Laura_lines_pca_raw.png", width=6, height=6, units="in", res=300)
plot(pca_raw$rotation[,1], pca_raw$rotation[,2], col=dataset_cols[dataset], pch=16, xlab="Component 1", ylab="Component 2")
legend("topleft",c("L1","L2","L3"), pch=16, col=dataset_cols, bty="n")
dev.off()

# Diffusion Maps
require("destiny")
set.seed(101)
dm_corrrected <- DiffusionMap(t(Combined_corrected[rownames(Combined_corrected) %in% features,]))
dm_ctp_norm <- DiffusionMap(t(Combined_ctp_norm[rownames(Combined_ctp_norm) %in% features,]))
dm_raw <- DiffusionMap(t(Combined_raw[rownames(Combined_raw) %in% features,]))

png("Laura_mnn_DM.png", width=6, height=6, units="in", res=300)
plot(eigenvectors(dm_corrrected)[,1], eigenvectors(dm_corrrected)[,3], col=prolif_cols[prolif+1], xlab="Dimension 1", ylab="Dimension 3", pch=dataset_pch[dataset])
legend("topleft", c("Prolif", "Diff", "L1", "L2", "L3"), col=c(rev(prolif_cols), "black", "black", "black"), pch=c(16,16, dataset_pch), bty="n")
dev.off()
png("Laura_ctp_norm_DM.png", width=6, height=6, units="in", res=300)
plot(eigenvectors(dm_ctp_norm)[,1], eigenvectors(dm_ctp_norm)[,3], col=prolif_cols[prolif+1], xlab="Dimension 1", ylab="Dimension 3", pch=dataset_pch[dataset])
legend("topleft", c("Prolif", "Diff", "L1", "L2", "L3"), col=c(rev(prolif_cols), "black", "black", "black"), pch=c(16,16, dataset_pch), bty="n")
dev.off()
png("Laura_raw_DM.png", width=6, height=6, units="in", res=300)
plot(eigenvectors(dm_raw)[,1], eigenvectors(dm_raw)[,3], col=prolif_cols[prolif+1], xlab="Dimension 1", ylab="Dimension 3", pch=dataset_pch[dataset])
legend("topleft", c("Prolif", "Diff", "L1", "L2", "L3"), col=c(rev(prolif_cols), "black", "black", "black"), pch=c(16,16, dataset_pch), bty="n")
dev.off()
png("Laura_raw_DM_datasets.png", width=6, height=6, units="in", res=300)
plot(eigenvectors(dm_raw)[,1], eigenvectors(dm_raw)[,2], col=dataset_cols[dataset], xlab="Dimension 1", ylab="Dimension 2", pch=16)
legend("topleft", c("L1", "L2", "L3"), col=c(col=dataset_cols), pch=16, bty="n")
dev.off()

# Consistent markers
features_C <- unique(cca1_marks[cca1_marks %in% hcc6_marks | cca1_marks %in% hcc10_marks], hcc6_marks[hcc6_marks %in% hcc10_marks])
dm_raw_C <- DiffusionMap(t(Combined_raw[rownames(Combined_raw) %in% features_C,]))
png("Laura_consistent_markers_DM.png", width=6, height=6, units="in", res=300)
plot(-1*eigenvectors(dm_raw_C)[,1], -1*eigenvectors(dm_raw_C)[,3], col=prolif_cols[prolif+1], xlab="Dimension 1", ylab="Dimension 3", pch=dataset_pch[dataset])
legend("bottomright", c("Prolif", "Diff", "L1", "L2", "L3"), col=c(rev(prolif_cols), "black", "black", "black"), pch=c(16,16, dataset_pch), bty="n")
dev.off()

# Are cell-types consistent across datasets?
my_match_heatmap <- function(sim_matrix, is_dist_mat=FALSE) {
        heatcols <- colorRampPalette(c("blue","khaki1","red"))(255)
	if (is_dist_mat){
		heatcols <- rev(heatcols);
        	dendro<-hclust(as.dist(sim_matrix))
		xlab_name <- "Distance"
	} else {
        	dendro<-hclust(as.dist(1-sim_matrix))
		xlab_name <- "Similarity"
	}


        require("CellTypeProfiles")

        heatout <- heatmap.3(sim_matrix, trace="n", scale="none", col=heatcols, symbreaks=FALSE,
                        key.title="", key.xlab=xlab_name, symm=TRUE,
                        Rowv=as.dendrogram(dendro), Colv=as.dendrogram(dendro))
}

cell_type_Combined_corrected <- my_row_mean_aggregate(Combined_corrected, paste("L",dataset,"-", clust_ID, sep=""))
d_corrected <- dist(t(cell_type_Combined_corrected))
sc_corrected <- cor(cell_type_Combined_corrected, method="spearman")
png("Laura_only_mnn.png", width=6, height=6, units="in", res=300)
my_match_heatmap(d_corrected, is_dist_mat=TRUE)
dev.off()

cell_type_Cfeatures_corrected <- my_row_mean_aggregate(Combined_raw[rownames(Combined_raw) %in% features_C,], paste("L",dataset,"-", clust_ID, sep=""))
d_sC <- dist(t(cell_type_Cfeatures_corrected))
sc_sC <- cor(cell_type_Cfeatures_corrected, method="spearman")
png("Laura_only_consistentmarkers.png", width=6, height=6, units="in", res=300)
my_match_heatmap(d_sC, is_dist_mat=TRUE)
dev.off()

# meta 
dataset_list <- c(L1=CCA1, L2=HCC6, L3=HCC10)
source("../3_Multi_Dataset_MetaNeighbour.R")
meta_LO_out <- do_metaneighbour(dataset_list)
colnames(meta_LO_out) <- rownames(meta_LO_out) <- paste("L",rep(c(1,2,3), each=6),"-", rep(1:6, time=3), sep="")
png("Laura_only_meta.png", width=6, height=6, units="in", res=300)
my_match_heatmap(meta_LO_out)
dev.off()

# scmap
source("../2_Muli_Dataset_scmap.R")
scmap_LO_out <- do_scmap(dataset_list, n_total=6*3)
colnames(scmap_LO_out) <- rownames(scmap_LO_out) <- paste("L",rep(c(1,2,3), each=6),"-", rep(1:6, time=3), sep="")
png("Laura_only_scmap.png", width=6, height=6, units="in", res=300)
my_match_heatmap(scmap_LO_out+t(scmap_LO_out))
dev.off()

#ctp 
colnames(CCA1_profile$profiles) = c(1:6)
colnames(HCC6_profile$profiles) = c(1:6)
colnames(HCC10_profile$profiles) = c(1:6)
profile_List = list(L1=CCA1_profile, L2=HCC6_profile, L3=HCC10_profile)

matches = combine_and_match_clusters(profile_List, features=features)
corrected_profiles = glm_of_matches(matches)

#ctp_LO_mat <- cor(corrected_profiles, method="spearman")
ctp_LO_mat <- dist(t(corrected_profiles))
png("Laura_only_ctp.png", width=6, height=6, units="in", res=300)
my_match_heatmap(ctp_LO_mat)
dev.off()

### DE across matched cell-types ####

# Wilcox.test for matched cell-type & corrected matrices
# mnnCorrect
differentiated <- (Laura_Merge_Stuff$clust_ID %in% c(5,6) & Laura_Merge_Stuff$dataset_lab == 1) |
		  (Laura_Merge_Stuff$clust_ID == 6 & Laura_Merge_Stuff$dataset_lab == 2)
intermediate <- (Laura_Merge_Stuff$clust_ID %in% c(3) & Laura_Merge_Stuff$dataset_lab == 1) |
		(Laura_Merge_Stuff$clust_ID %in% c(1,2,3) & Laura_Merge_Stuff$dataset_lab == 2) |
		(Laura_Merge_Stuff$clust_ID %in% c(4,5) & Laura_Merge_Stuff$dataset_lab == 3) 
proliferating <- (Laura_Merge_Stuff$clust_ID %in% c(1,4) & Laura_Merge_Stuff$dataset_lab == 1) |
		(Laura_Merge_Stuff$clust_ID %in% c(4) & Laura_Merge_Stuff$dataset_lab == 2) |
		(Laura_Merge_Stuff$clust_ID %in% c(1) & Laura_Merge_Stuff$dataset_lab == 3) 

p_vs_i <- apply(Combined_corrected, 1, function(x) {wilcox.test(x[proliferating], x[intermediate])$p.value})
p_vs_d <- apply(Combined_corrected, 1, function(x) {wilcox.test(x[proliferating], x[differentiated])$p.value})
i_vs_d <- apply(Combined_corrected, 1, function(x) {wilcox.test(x[differentiated], x[intermediate])$p.value})

p_vs_i_dir <- rowMeans(Combined_corrected[,proliferating])-rowMeans(Combined_corrected[,intermediate])
p_vs_d_dir <- rowMeans(Combined_corrected[,proliferating])-rowMeans(Combined_corrected[,differentiated])
i_vs_d_dir <- rowMeans(Combined_corrected[,intermediate])-rowMeans(Combined_corrected[,differentiated])

require("gProfileR")
# p-high
p_high <-rownames(Combined_corrected)[ p.adjust(p_vs_i, method="bon") < 0.01 & p_vs_i_dir > 1]
out <- gprofiler(p_high, organism="hsapiens")
out <- out[order(out$p.value),]
out <- out[out$domain %in% c("BP", "MF"),]
out$prop_accounted <- out$overlap.size/out$query.size
head(out[,c(3,6, 9, 10, 12, 15)], 50)

png("Top_Prolif_DE_GO_rich.png", width=7.5, height=6, units="in", res=300)
n=25
par(mar=c(4,20,1,1))
barplot(rev(out[1:n,15]*100), names=rev(out[1:n,"term.name"]), horiz=TRUE, las=1, xlab="Percent of DE Genes", col="black")
dev.off()
p_high_out <- out

# d-high
d_high <-rownames(Combined_corrected)[ p.adjust(i_vs_d, method="bon") < 0.01 & i_vs_d_dir < -1]
out <- gprofiler(d_high, organism="hsapiens")
out <- out[order(out$p.value),]
out <- out[out$domain %in% c("BP", "MF"),]
out$prop_accounted <- out$overlap.size/out$query.size
head(out[,c(3,6, 9, 10, 12, 15)], 50)

png("Top_Diff_DE_GO_rich.png", width=7.5, height=6, units="in", res=300)
n=25
par(mar=c(4,20,1,1))
barplot(rev(out[1:n,15]*100), names=rev(out[1:n,"term.name"]), horiz=TRUE, las=1, xlab="Percent of DE Genes", col="forestgreen")
dev.off()
d_high_out <- out


# i-high
i_high1 <-rownames(Combined_corrected)[ p.adjust(i_vs_d, method="bon") < 0.01 & i_vs_d_dir > 0.6]
i_high2 <-rownames(Combined_corrected)[ p.adjust(p_vs_i, method="bon") < 0.01 & p_vs_i_dir < -0.6]
i_high <-i_high1[i_high1 %in% i_high2]
out <- gprofiler(i_high, organism="hsapiens")
out <- out[order(out$p.value),]
out <- out[out$domain %in% c("BP", "MF"),]
out$prop_accounted <- out$overlap.size/out$query.size
head(out[,c(3,6, 9, 10, 12, 15)], 50)

png("Top_Inter_DE_GO_rich.png", width=7.5, height=6, units="in", res=300)
n=25
par(mar=c(4,20,1,1))
barplot(rev(out[1:n,15]*100), names=rev(out[1:n,"term.name"]), horiz=TRUE, las=1, xlab="Percent of DE Genes", col="dodgerblue")
dev.off()
i_high_out <- out

# Classify
thresh = 0.6;
decreasing <- rownames(Combined_corrected)[ (p.adjust(p_vs_i, method="bon") < 0.01 & p_vs_i_dir > thresh ) &
						( p.adjust(i_vs_d, method="bon") < 0.01 & i_vs_d_dir > thresh )]
increasing <- rownames(Combined_corrected)[ (p.adjust(p_vs_i, method="bon") < 0.01 & p_vs_i_dir < -thresh) &
						( p.adjust(i_vs_d, method="bon") < 0.01 & i_vs_d_dir < -thresh)]

i_high <- rownames(Combined_corrected)[ (p.adjust(p_vs_i, method="bon") < 0.01 & p_vs_i_dir < -thresh) &
						( p.adjust(i_vs_d, method="bon") < 0.01 & i_vs_d_dir > thresh ) ]
i_low <- rownames(Combined_corrected)[ (p.adjust(p_vs_i, method="bon") < 0.01 & p_vs_i_dir > thresh) &
						( p.adjust(i_vs_d, method="bon") < 0.01 & i_vs_d_dir < -thresh ) ]

p_high <- rownames(Combined_corrected)[ (p.adjust(p_vs_i, method="bon") < 0.01 & p_vs_i_dir > thresh ) &
					( p.adjust(i_vs_d, method="bon") > 0.01 & abs(i_vs_d_dir) < thresh ) ]
p_low <- rownames(Combined_corrected)[ (p.adjust(p_vs_i, method="bon") < 0.01 & p_vs_i_dir < -thresh ) &
                                        ( p.adjust(i_vs_d, method="bon") > 0.01 & abs(i_vs_d_dir) < thresh ) ]

d_high <-rownames(Combined_corrected)[ ( p.adjust(i_vs_d, method="bon") < 0.01 & i_vs_d_dir < -thresh ) &
					( p.adjust(p_vs_i, method="bon") > 0.01 & abs(p_vs_i_dir) < thresh ) ]
d_low <- rownames(Combined_corrected)[ ( p.adjust(i_vs_d, method="bon") < 0.01 & i_vs_d_dir > thresh ) &
                                        ( p.adjust(p_vs_i, method="bon") > 0.01 & abs(p_vs_i_dir) < thresh ) ]

png("Laura_Liver_Combined_DE_Classes.png", width=10, height=9, units="in", res=300)
par(mfcol=c(3,3))
par(mar=c(1,1,1,1))
plot(1:3, c(1,0,0), xaxt="n", yaxt="n", xlab="", ylab="", type="l", lwd=3); text(1,0.1,length(p_high), cex=2, font=2, pos=4)
plot(1:3, c(1,0.5,0), xaxt="n", yaxt="n", xlab="", ylab="", type="l", lwd=3); text(1,0.1,length(decreasing), cex=2, font=2, pos=4)
plot(1:3, c(1,1,0), xaxt="n", yaxt="n", xlab="", ylab="", type="l", lwd=3); text(1,0.1,length(d_low), cex=2, font=2, pos=4)
plot(1:3, c(0,1,1), xaxt="n", yaxt="n", xlab="", ylab="", type="l", lwd=3); text(3,0.1,length(p_low), cex=2, font=2, pos=2)
plot(1:3, c(0,0.5,1), xaxt="n", yaxt="n", xlab="", ylab="", type="l", lwd=3); text(3,0.1,length(increasing), cex=2, font=2, pos=2)
plot(1:3, c(0,0,1), xaxt="n", yaxt="n", xlab="", ylab="", type="l", lwd=3); text(3,0.1,length(d_high), cex=2, font=2, pos=2)
plot(1:3, c(0,1,0), xaxt="n", yaxt="n", xlab="", ylab="", type="l", lwd=3); text(2,0.1,length(i_high), cex=2, font=2)
plot(1:3, c(1,0,1), xaxt="n", yaxt="n", xlab="", ylab="", type="l", lwd=3); text(1,0.1,length(i_low), cex=2, font=2, pos=4)
dev.off()


my_do_gprof <- function(g) {
	require("gProfileR")
	out <- gprofiler(g, organism="hsapiens")
	out <- out[order(out$p.value),]
	out <- out[out$domain %in% c("BP", "MF"),]
	out$prop_accounted <- out$overlap.size/out$query.size
	return(out)
}

dec_out <- my_do_gprof(decreasing)
inc_out <- my_do_gprof(increasing)
ih_out <- my_do_gprof(i_high)
il_out <- my_do_gprof(i_low)
ph_out <- my_do_gprof(p_high)
pl_out <- my_do_gprof(p_low)
dh_out <- my_do_gprof(d_high)
dl_out <- my_do_gprof(d_low)

head(dec_out[,"term.name"], 25)

source("~/R-Scripts/Ensembl_Stuff.R")
map_symbol_ensg(unlist(strsplit(dl_out[dl_out$term.name=="translation","intersection"], ",")), is.org="Hsap", is.name="ensg")



# GLM with significant cell-type effects for CTP


### Single strain figures ###

# CCA1 clusters
set.seed(123)
#marks <- map_symbol_ensg(markers$CCA1$Gene[markers$CCA1$pvalue < sig_threshold], is.org="Hsap", is.name="ensg")
marks <- cca1_marks

tmp <- exprs(CCA1)[rownames(exprs(CCA1)) %in% marks,]
tmp <- tmp[rowSums(tmp) > 10,]
tsne_CCA1 <- Rtsne(t(tmp))
png("L1_tsne.png", width=6, height=6, units="in", res=300)
plot(tsne_CCA1$Y, col=clust_cols[CCA1$sc3_6_clusters], pch=16, xlab="tSNE Dim 1", ylab="tSNE Dim 2", main="L1 Clusters")
legend("bottomright", col=clust_cols, c(as.character(1:6)), pch=16, bty="n")
dev.off()

# HCC6 clusters
set.seed(123)
#marks <- map_symbol_ensg(markers$HCC6$Gene[markers$HCC6$pvalue < sig_threshold], is.org="Hsap", is.name="ensg")
marks <- hcc6_marks

tmp <- exprs(HCC6)[rownames(exprs(HCC6)) %in% marks,]
tmp <- tmp[rowSums(tmp) > 10,]
tsne_HCC6 <- Rtsne(t(tmp))
png("L2_tsne.png", width=6, height=6, units="in", res=300)
plot(tsne_HCC6$Y, col=clust_cols[HCC6$sc3_6_clusters], pch=16, xlab="tSNE Dim 1", ylab="tSNE Dim 2", main="L2 Clusters")
legend("bottomright", col=clust_cols, c(as.character(1:6)), pch=16, bty="n")
dev.off()

# HCC10 clusters
set.seed(123)
#marks <- map_symbol_ensg(markers$HCC10$Gene[markers$HCC10$pvalue < sig_threshold], is.org="Hsap", is.name="ensg")
marks <- hcc10_marks

tmp <- exprs(HCC10)[rownames(exprs(HCC10)) %in% marks,]
tmp <- tmp[rowSums(tmp) > 10,]
tsne_HCC10 <- Rtsne(t(tmp))
png("L3_tsne.png", width=6, height=6, units="in", res=300)
plot(tsne_HCC10$Y, col=clust_cols[HCC10$sc3_6_clusters], pch=16, xlab="tSNE Dim 1", ylab="tSNE Dim 2", main="L3 Clusters")
legend("bottomright", col=clust_cols, c(as.character(1:6)), pch=16, bty="n")
dev.off()


### Laura Liver Seurat ###
