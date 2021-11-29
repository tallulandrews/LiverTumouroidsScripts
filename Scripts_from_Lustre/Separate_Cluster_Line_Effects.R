require("M3Drop")
require("scater")
require("matrixStats")
require("RColorBrewer")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/DiffExpr/DE_functions.R")
source("Cluster_Profiles_Functions.R")
require("CellTypeProfiles")
source("complex_markers.R")

set.seed(142)

CCA1 <- readRDS("CCA1_SC3_Prolif.rds")
HCC6 <- readRDS("HCC6_SC3_Prolif.rds")
HCC10 <- readRDS("HCC10_SC3_Prolif.rds")

exprs(CCA1) <- get_log_norm(CCA1)
exprs(HCC6) <- get_log_norm(HCC6)
exprs(HCC10) <- get_log_norm(HCC10)

### Tidy Clusters ###
Tidy_Clusters <- function(SCE, labels) {
	marks <- complex_markers(exprs(SCE), labels)
	marks_sig <- marks[marks$q.value < 0.05,] # significantly DE
	marks_sig_tab <- marks_sig[,-c(1, ncol(marks_sig), ncol(marks_sig)-1)] # just assignment table
	marks_sig_tab[rowSums(marks_sig_tab) > ncol(marks_sig_tab)/2,] <- 1-marks_sig_tab[rowSums(marks_sig_tab) > ncol(marks_sig_tab)/2,] # Invert negative markers
	marks_uniqueness <- rowSums(marks_sig_tab) # number of clusters is a marker of
	cluster_sizes <- factor_counts(labels) # number of cells in each cluser
	expected_markers <- qbinom(1-0.05, size=nrow(marks_sig_tab), prob=0.05/ncol(marks_sig_tab)) 
	unique_markers <- colSums(marks_sig_tab[marks_uniqueness ==1,])

	#shared_table <- sapply(1:length(names(cluster_sizes)), function(c) {colSums(marks_sig_tab[ marks_sig_tab[,c] == 1 & marks_uniqueness >1,])})
	#diag(shared_table) = 0

	new_clusters <- as.character(labels)
	for (c in 1:length(names(cluster_sizes))) {
		c_name <- names(cluster_sizes)[c]
		unique_marks <- sum(marks_sig_tab[,c] == 1 & marks_uniqueness ==1)
		new_name <- c_name;
		if (cluster_sizes[c] <=sum(cluster_sizes)*0.05) { # suspect clusters = < 5% of total dataset
			#5% FDR equally distributed across clusters
			if (unique_marks > expected_markers) {
				new_name <- "Outliers"
			} else { # which share most DE with
				shared <- colSums(marks_sig_tab[ marks_sig_tab[,c] == 1 & marks_uniqueness >1,])
				equiv <- which(shared == max(shared[-c]))
				new_name <- names(shared)[equiv]
			}
		} else {
			shared <- colSums(marks_sig_tab[ marks_sig_tab[,c] == 1 & marks_uniqueness >1,])
			equiv <- which(shared == max(shared[-c]))
			if( shared[equiv] > unique_marks & shared[equiv]/sum(shared[-c]) > 0.5 ) {
				new_name <- names(shared)[equiv]
			}
		}
		if ( sum(new_clusters == new_name) == 0 & new_name != "Outliers") {
			new_clusters[labels == new_name] <- c_name
		} else {
			new_clusters[new_clusters==c_name] <- new_name
		}
	}
	print(table(new_clusters, labels))
	return(new_clusters)
}

pData(CCA1)$sc3_clusters <- pData(CCA1)$sc3_6_clusters
pData(HCC6)$sc3_clusters <- pData(HCC6)$sc3_6_clusters
pData(HCC10)$sc3_clusters <- pData(HCC10)$sc3_6_clusters

pData(CCA1)$tidy_clusters <- Tidy_Clusters(CCA1, pData(CCA1)$sc3_clusters)
pData(HCC6)$tidy_clusters <- Tidy_Clusters(HCC6, pData(HCC6)$sc3_clusters)
pData(HCC10)$tidy_clusters <- Tidy_Clusters(HCC10, pData(HCC10)$sc3_clusters)

Assign_Well_coords <- function(SCE) {
	x <- as.numeric(factor(substr(SCE$Well, 1, 1)))
	y <- as.numeric(substr(SCE$Well, 2, 3))
	SCE$Well_x <- x
	SCE$Well_y <- y
	return(SCE)
}

CCA1 <- Assign_Well_coords(CCA1)
HCC6 <- Assign_Well_coords(HCC6)
HCC10 <- Assign_Well_coords(HCC10)

#pData(CCA1)$Dataset <- rep("CCA1", times=ncol(CCA1))
#pData(HCC6)$Dataset <- rep("HCC6", times=ncol(HCC6))
#pData(HCC10)$Dataset <- rep("HCC10", times=ncol(HCC10))

### Combine datasets ###
require("scater")
anno_cols <- c("sc3_clusters", "tidy_clusters", "Well", "Well_x", "Well_y", "Plate", "Type", "CC_state")
CCA1_pData <- pData(CCA1)[,anno_cols]
HCC6_pData <- pData(HCC6)[,anno_cols]
HCC10_pData <- pData(HCC10)[,anno_cols]

#CCA1_pData$Final_Clusters <- paste("CCA1", pData(CCA1)$sc3_6_clusters, sep="_")
#HCC6_pData$Final_Clusters <- paste("HCC6", pData(HCC6)$sc3_6_clusters, sep="_")
#HCC10_pData$Final_Clusters <- paste("HCC10", pData(HCC10)$sc3_6_clusters, sep="_")

Combined_counts <- cbind(counts(CCA1), counts(HCC6), counts(HCC10))
Combined_pData <- rbind(CCA1_pData, HCC6_pData, HCC10_pData)
Combined_fData <- fData(CCA1)[,1:2]

gd <- new("AnnotatedDataFrame", data=Combined_fData)
pd <- new("AnnotatedDataFrame", data=Combined_pData)
CombinedSCE <- newSCESet(countData = Combined_counts, phenoData = pd, featureData = gd)

CombinedSCE <- calculateQCMetrics(CombinedSCE)
CombinedSCE$Final_Clusters <- paste(CombinedSCE$Type, CombinedSCE$sc3_clusters, sep="_")


### Weighted Means ###

batch <- factor(pData(CombinedSCE)$Type)
batch_labs <- levels(batch)

# Correct For Line Effects
corrected_line <- exprs(CombinedSCE)
for( b in 1:length(batch_labs)) {
	cluster_means <- my_row_mean_aggregate(corrected_line[,batch==batch_labs[b]], pData(CombinedSCE)$Final_Clusters[batch==batch_labs[b]])

        corrected_line[,batch==batch_labs[b]] <- corrected_line[,batch==batch_labs[b]]-rowMeans(cluster_means)
}

require("CellTypeProfiles")
source("/nfs/users/nfs_t/ta6/R-Scripts/Ensembl_Stuff.R")
set.seed(671)
png("WeightedMean_Profiles_Heatmap.png", width=6, height=6, units="in", res=300)
cluster_relative_heatmap(corrected_line, pData(CombinedSCE)$Final_Clusters, npermute=100)
dev.off()

my_profiles_line <- my_row_mean_aggregate(corrected_line,  pData(CombinedSCE)$Final_Clusters)

write.table(my_profiles_line, file="WeightedMean_Profiles.csv", sep=",")
forLaura <- data.frame(my_profiles, Symbol = map_symbol_ensg(rownames(my_profiles), is.org="Hsap", is.name="ensg"), Ensg = rownames(my_profiles))

write.table(forLaura, file="ForLaura_WeightedMean_Profiles_full.csv", sep=",")

# Correct For Cell-Type Effects
CombinedSCE_Cleaned <- CombinedSCE[,CombinedSCE$tidy_clusters != "Outliers"]


matched_celltypes <- paste(CombinedSCE_Cleaned$Type, CombinedSCE_Cleaned$tidy_clusters, sep="_")
Prolif <- c("CCA1_1", "HCC6_4", "CCA1_4", "HCC10_1")
Diff <- c("CCA1_6", "HCC6_5", "HCC10_3")
Inter <- c("HCC6_1", "HCC6_2", "CCA1_3", "HCC10_4", "HCC10_5", "HCC10_6", "HCC6_3")

matched_celltypes[matched_celltypes %in% Prolif] <- "Prolif"
matched_celltypes[matched_celltypes %in% Diff] <- "Diff"
matched_celltypes[matched_celltypes %in% Inter] <- "Inter"

matched_celltypes <- factor(matched_celltypes)
matched_labs <- levels(matched_celltypes)
corrected_ct <- exprs(CombinedSCE_Cleaned)
for( m in 1:length(matched_labs)) {
	this_subset <- matched_celltypes==matched_labs[m]
        dataset_means <- my_row_mean_aggregate(corrected_ct[,this_subset], 
				factor(pData(CombinedSCE_Cleaned)$Type[this_subset]))

        corrected_ct[,this_subset] <- corrected_ct[,this_subset]-rowMeans(dataset_means)
}


set.seed(749)
cluster_relative_heatmap(corrected_ct, pData(CombinedSCE)$Final_Clusters, npermute=100)

my_profiles_ct <- my_row_mean_aggregate(corrected_ct,  pData(CombinedSCE_Cleaned)$Final_Clusters)

### GLM ###
matched_celltypes <- paste(CombinedSCE$Type, CombinedSCE$tidy_clusters, sep="_")
Prolif <- c("CCA1_1", "HCC6_4", "CCA1_4", "HCC10_1")
Diff <- c("CCA1_6", "HCC6_5", "HCC10_3")
Inter <- c("HCC6_1", "HCC6_2", "CCA1_3", "HCC10_4", "HCC10_5", "HCC10_6", "HCC6_3")

matched_celltypes[matched_celltypes %in% Prolif] <- "Prolif"
matched_celltypes[matched_celltypes %in% Diff] <- "Diff"
matched_celltypes[matched_celltypes %in% Inter] <- "Inter"

pData(CombinedSCE)$matched_celltypes <- matched_celltypes

CombinedSCE_Cleaned <- CombinedSCE[,CombinedSCE$tidy_clusters != "Outliers"]

CombinedSCE_Cleaned$matched_celltypes <- factor(CombinedSCE_Cleaned$matched_celltypes)
CombinedSCE_Cleaned$Type <- factor(CombinedSCE_Cleaned$Type)

glm_mod <- model.matrix(~CombinedSCE_Cleaned$matched_celltypes+CombinedSCE_Cleaned$Type+CombinedSCE_Cleaned$Well_x+CombinedSCE_Cleaned$Well_y+CombinedSCE_Cleaned$Plate)

run_glm <- function(gene_expr) {
	res <- glm(gene_expr ~ CombinedSCE_Cleaned$matched_celltypes+CombinedSCE_Cleaned$Type+CombinedSCE_Cleaned$Well_x+CombinedSCE_Cleaned$Well_y+CombinedSCE_Cleaned$Plate)
	return(summary(res)$coef)
}
out <- apply(exprs(CombinedSCE_Cleaned), 1, run_glm)
rownames(out) <- c(colnames(glm_mod), paste(colnames(glm_mod), "Err"), paste(colnames(glm_mod), "t"), paste(colnames(glm_mod), "p"))

# what do we want to look at?
line_specific <- p.adjust(out[28,], method="bon") < 0.05 |  p.adjust(out[29,], method="bon") < 0.05
line_effect <- out[c(4,5),]
line_effect_max <- apply(line_effect, 2, max)

type_specific <- p.adjust(out[26,], method="bon") < 0.05 |  p.adjust(out[27,], method="bon") < 0.05
type_effect <- out[c(2,3),]
type_effect_max <- apply(type_effect, 2, max)
# classification: each line up/down x each type up/down
# HCC10+HCC6+=4, HCC10+HCC6-=2, HCC10-HCC6- = 1, HCC10-HCC6+ = 3
class_line_dir <- as.numeric(out[4,] > 0)+1
class_line_dir[out[5,] > 0] <- class_line_dir[out[5,] > 0]+2
class_line_dir[!line_specific] <- 0

# Inter+Prolif+ = 4, Inter+Prolif- = 2, Inter-Prolif+ =3, Inter-Prolif- = 1
class_type_dir <- as.numeric(out[2,] > 0)+1
class_type_dir[out[3,] > 0] <- class_type_dir[out[3,] > 0]+2
class_type_dir[!type_specific] <- 0

### Don't like the above anymore ###
# Instead let's do the three lines against each other while controling for cell-cycle

#cellcycle <- read.table("~/Data/Whitfield_CC.txt")
#cellcycle_simple <- as.matrix(cellcycle[cellcycle[,1] != "CC",])
#cellcycle_simple[cellcycle_simple[,1] == "G2",] = "G2M";
#cellcycle_simple[cellcycle_simple[,1] == "S",] = "G1S";
#cellcycle_simple = cellcycle_simple[cellcycle_simple[,1] != "MG1",];
#G0_genes <- read.table("~/Data/Reactome_G0.txt", header=F)

run_glm <- function(gene_expr) {
        res <- glm(gene_expr ~ CombinedSCE_Cleaned$CC_state+CombinedSCE_Cleaned$Type+CombinedSCE_Cleaned$Plate)
        return(summary(res)$coef)
}

out <- apply(exprs(CombinedSCE_Cleaned), 1, run_glm)
glm_mod <- model.matrix(~CombinedSCE_Cleaned$CC_state+factor(CombinedSCE_Cleaned$Type)+CombinedSCE_Cleaned$Plate)
rownames(out) <- c(colnames(glm_mod), paste(colnames(glm_mod), "Err"), paste(colnames(glm_mod), "t"), paste(colnames(glm_mod), "p"))


line_specific <- p.adjust(out[26,], method="bon") < 0.05 |  p.adjust(out[27,], method="bon") < 0.05
line_effect <- out[c(5,6),]
# Top genes each kind
line_DE <- t(out)
class_cca_l <- line_DE[,5] > 0 & line_DE[,6] > 0 & p.adjust(line_DE[,26], method="bon") < 0.05 & p.adjust(line_DE[,27], method="bon") < 0.05
class_cca_h <- line_DE[,5] < 0 & line_DE[,6] < 0 & p.adjust(line_DE[,26], method="bon") < 0.05 & p.adjust(line_DE[,27], method="bon") < 0.05
class_hcc10_h <- line_DE[,6] > 0 & p.adjust(line_DE[,26], method="bon") < 0.05 & p.adjust(line_DE[,27], method="bon") > 0.05
class_hcc10_l <- line_DE[,6] < 0 & p.adjust(line_DE[,26], method="bon") < 0.05 & p.adjust(line_DE[,27], method="bon") > 0.05
class_hcc6_h <- line_DE[,5] > 0 & p.adjust(line_DE[,26], method="bon") > 0.05 & p.adjust(line_DE[,27], method="bon") < 0.05
class_hcc6_l <- line_DE[,5] < 0 & p.adjust(line_DE[,26], method="bon") > 0.05 & p.adjust(line_DE[,27], method="bon") < 0.05

require("gProfileR")
cca1_h_func <- gprofiler(rownames(line_DE)[class_cca_h])
cca1_l_func <- gprofiler(rownames(line_DE)[class_cca_l])
hcc6_l_func <- gprofiler(rownames(line_DE)[class_hcc6_l])
hcc6_h_func <- gprofiler(rownames(line_DE)[class_hcc6_h])
hcc10_l_func <- gprofiler(rownames(line_DE)[class_hcc10_l])
hcc10_h_func <- gprofiler(rownames(line_DE)[class_hcc10_h])





## edgeR version
require("edgeR")
edgeR_obj <- DGEList(counts=counts(CombinedSCE_Cleaned))
design <- model.matrix(~ CombinedSCE_Cleaned$matched_celltypes+factor(CombinedSCE_Cleaned$Type)+CombinedSCE_Cleaned$Well_x+CombinedSCE_Cleaned$Well_y+CombinedSCE_Cleaned$Plate);
edgeR_obj <- estimateDisp(edgeR_obj, design)
fit <- glmQLFit(edgeR_obj, design);
# celltype
ct.test <- glmQLFTest(fit, coef=2:3)
line.test <- glmQLFTest(fit, coef=4:5)
pos.test <- glmQLFTest(fit, coef=6:7)
topTags(ct.test)

## MAST version
require("MAST")
Mast_pdata <- pData(CombinedSCE_Cleaned)
Mast_pdata$line <- factor(Mast_pdata$Type)
Mast_pdata$matched_celltypes <- factor(Mast_pdata$matched_celltypes)
#Mast_pdata$wellKey <- CombinedSCE_Cleaned$Well
MAST_obj <- FromMatrix(exprs(CombinedSCE_Cleaned), Mast_pdata, fData(CombinedSCE_Cleaned))
zlmCond <- zlm(~ matched_celltypes + line + Well_x + Well_y + Plate + total_features, MAST_obj)
summaryCondCCA1 <- summary(zlmCond, doLRT = "lineCCA1")
summaryCondHCC10 <- summary(zlmCond, doLRT = "lineHCC10")
summaryCondInter <- summary(zlmCond, doLRT = "matched_celltypesInter")
summaryCondProlif <- summary(zlmCond, doLRT = "matched_celltypesProlif")

MAST_out <- list(CCA1=summaryCondCCA1, HCC10=summaryCondHCC10, Inter=summaryCondInter, Prolif=summaryCondProlif)
saveRDS(MAST_out, file="MAST_output.rds")

summaryCondHCC10 <- summaryCondHCC10$datatable[component == "H",]
summaryCondCCA1  <- summaryCondCCA1$datatable[component == "H",]
summaryCondInter  <- summaryCondInter$datatable[component == "H",]
summaryCondProlif  <- summaryCondProlif$datatable[component == "H",]
