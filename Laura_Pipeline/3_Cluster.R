# R3.4
args <- commandArgs(trailingOnly=TRUE) # rds file for data, prefix for output, number of cores

source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("~/Collaborations/LiverOrganoids/Laura_Pipeline/99_Functions.R")

n.cores = args[3];
outprefix = args[2];
k_min = 4
if (length(args) > 3) {
	expr_type = args[4];
} else {
	expr_type = "lognorm"
}

if (length(args) > 4) {
	CC_rm = as.logical(args[5]);
} else {
	CC_rm=FALSE
}


if (CC_rm) {
	outprefix <- paste(outprefix, "noCC", sep="_")
}

require("scater")
SCE <- readRDS(args[1]);
#SCE <- readRDS("Test_QCed_BatchCor.rds")
#SCE <- SCE[,SCE$Type=="CCA1"]
if (class(SCE)[1] == "SCESet") {
	SCE <- scater::toSingleCellExperiment(SCE)
}

set.seed(123)

max_ks <- ceiling(dim(SCE)[2]/25)

# SC3 - farm job 4 cores
require("SC3")
source("~/R-Scripts/Ensembl_Stuff.R")

rowData(SCE)$feature_symbol <- General_Map(rownames(SCE), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")

SCE_orig <- SCE

# Optional : remove CC
if (CC_rm) {
	SCE <- remove_cell_cycle(SCE)
}

# SCE Clustering
assays(SCE)[["logcounts"]] <- assays(SCE)[[expr_type]]
SCE <- sc3_prepare(SCE)
SCE <- sc3_estimate_k(SCE)
estimated_k <- SCE@metadata$sc3$k_estimation # Does this work?
if( estimated_k > max_ks ) {
	max_ks <- estimated_k
}
SCE <- sc3(SCE, ks=2:max_ks, n_cores=as.numeric(n.cores), biology=FALSE)

print("SC3 - finished")
## Redo clustering ##
# use consensus agreement as distance directly for clustering
# Calculate silhouette & stability & internal-vs-external consistency
# Algorithmicly select k
# save best k clusters



Ks <- get_optimal_k(SCE)
optimal_k <- Ks$optim
fine_k <- Ks$fine

print("optimal k found")

png(paste(outprefix, "k", optimal_k, "SC3_consensus.png", sep="_"), width=7, height=7, units="in", res=300)
out <- my_clustering(SCE, as.character(optimal_k), suppress.plot=FALSE)
dev.off()

require("CellTypeProfiles")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/My_R_packages/CellTypeProfiles/R/Markers.R")
if (!CC_rm) {
	colData(SCE_orig)$clusters_coarse <- factor(out$Cs); # Save optimal K clusters
	markers_coarse <- complex_markers(assays(SCE_orig)[[expr_type]], colData(SCE_orig)$clusters_coarse)
	markers_coarse$is.Feature <- markers_coarse$q.value < 0.05 & markers_coarse$q.value >= 0 & markers_coarse$AUC > 0.7
	colnames(markers_coarse) <- paste("coarse_marker", colnames(markers_coarse), sep="_");
} else {
	colData(SCE_orig)$clusters_noCC_coarse <- factor(out$Cs)
	markers_coarse <- complex_markers(assays(SCE_orig)[[expr_type]], colData(SCE_orig)$clusters_noCC_coarse)
	markers_coarse$is.Feature <- markers_coarse$q.value < 0.05 & markers_coarse$q.value >= 0 & markers_coarse$AUC > 0.7
	colnames(markers_coarse) <- paste("noCC_coarse_marker", colnames(markers_coarse), sep="_");
}

identical(rownames(markers_coarse), rownames(rowData(SCE_orig)))
rowData(SCE_orig) <- cbind(rowData(SCE_orig), markers_coarse)


png(paste(outprefix, "k", fine_k, "SC3_consensus.png", sep="_"), width=7, height=7, units="in", res=300)
out <- my_clustering(SCE, as.character(fine_k), suppress.plot=FALSE)
dev.off()

if (!CC_rm) {
	colData(SCE_orig)$clusters_fine <- factor(out$Cs); # Save optimal K clusters
	markers_fine <- complex_markers(assays(SCE_orig)[[ expr_type ]], colData(SCE_orig)$clusters_fine)
	markers_fine$is.Feature <- markers_fine$q.value < 0.05 & markers_fine$q.value >= 0 & markers_fine$AUC > 0.7
	markers_refinement <- markers_fine
	colnames(markers_fine) <- paste("fine_marker", colnames(markers_fine), sep="_");
} else {
	colData(SCE_orig)$clusters_noCC_fine <- factor(out$Cs); # Save optimal K clusters
	markers_fine <- complex_markers(assays(SCE_orig)[[ expr_type ]], colData(SCE_orig)$clusters_noCC_fine)
	markers_fine$is.Feature <- markers_fine$q.value < 0.05 & markers_fine$q.value >= 0 & markers_fine$AUC > 0.7
	markers_refinement <- markers_fine
	colnames(markers_fine) <- paste("noCC_fine_marker", colnames(markers_fine), sep="_");
}
identical(rownames(markers_fine), rownames(rowData(SCE_orig)))
rowData(SCE_orig) <- cbind(rowData(SCE_orig), markers_fine)


## Refine Clusters ##
# calculate markers & refine clusters
# save refined clusters

# Calculate Markers
#refinement <- refine_clusters(SCE_orig, SCE, expr_type)
if (!CC_rm) {
	refinement <- refine_clusters(SCE, expr_type, colData(SCE_orig)$clusters_fine, markers_refinement)
} else {
	refinement <- refine_clusters(SCE, expr_type, colData(SCE_orig)$clusters_noCC_fine, markers_refinement)
}
newCs <- refinement$newCs
#markers <- refinement$markers

print("clusters refined")

if (!CC_rm) {
	colData(SCE_orig)$clusters_clean <- factor(newCs);
} else {
	colData(SCE_orig)$clusters_noCC_clean <- factor(newCs);
}


#identical(rownames(markers), rownames(rowData(SCE_orig)))
#colnames(markers) <- paste("fine_marker", colnames(markers), sep="_");
#rowData(SCE_orig) <- cbind(rowData(SCE_orig), markers)

saveRDS(SCE_orig, file=paste(outprefix,"SC3.rds", sep="_"))
