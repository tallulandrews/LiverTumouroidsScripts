Camp <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/camp_human.rds")
Camp_mouse <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/camp_mouse.rds")
require("scater")
require("RColorBrewer")


my_QC <- function(dat, min_detect=2000, min_counts=10000) {
        par(mar=c(4,4,1,1))
        plot(colSums(dat), colSums(dat > 0), xlab="total counts", ylab="total genes")

        abline(h=min_detect)
        abline(v=min_counts)
        filter_detect <- colSums(dat > 0) > min_detect
        filter_total <- colSums(dat) > min_counts

        return((filter_detect & filter_total))
}

pData(Camp)$cell_type1 <- as.character(pData(Camp)$cell_type1)
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="definitive endoderm"] <- "d-endo"
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="hepatic endoderm"] <- "h-endo"
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="immature hepatoblast"] <- "i-hepato"
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="mature hepatocyte"] <- "m-hepato"
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="mesenchymal stem cell"] <- "msc"
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="adult hepatocytes"] <- "a-hepato"
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="erythroblasts"] <- "rbc"
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="fetal hepatocytes"] <- "f-hepato"
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="lymphoblasts"] <- "lymph"
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="lymphoblasts"] <- "lymph"
pData(Camp)$cell_type1[pData(Camp)$cell_type1=="endothelial"] <- "endoth"
pData(Camp)$cell_type1 <- factor(pData(Camp)$cell_type1)

pData(Camp_mouse)$cell_type1 <- as.character(pData(Camp_mouse)$cell_type1)
pData(Camp_mouse)$cell_type1[pData(Camp_mouse)$cell_type1=="immune 1"] <- "immune"
pData(Camp_mouse)$cell_type1[pData(Camp_mouse)$cell_type1=="immune 2"] <- "immune"
pData(Camp_mouse)$cell_type1[pData(Camp_mouse)$cell_type1=="hepatoblasts 1"] <- "hepato"
pData(Camp_mouse)$cell_type1[pData(Camp_mouse)$cell_type1=="hepatoblasts 2"] <- "hepato"
pData(Camp_mouse)$cell_type1<-factor(pData(Camp_mouse)$cell_type1)


Camp_liver_bud <- Camp[,pData(Camp)$age %in% c("liver bud","late liver bud")]
Camp_liver_bud <- Camp_liver_bud[,my_QC(exprs(Camp_liver_bud))]

Camp_transplant <- Camp[,pData(Camp)$age %in% c("10 days post-transplant","15 days post-transplant","3 days post-transplant")]
Camp_transplant <- Camp_transplant[,my_QC(exprs(Camp_transplant))]

Camp_2D <- Camp[,pData(Camp)$age %in% c("0 days", "6 days","8 days","14 days","21 days")]
Camp_2D <- Camp_2D[,my_QC(exprs(Camp_2D), 4500, 20000)]

Camp_Htissue <- Camp[,pData(Camp)$Source == "Liver"]
Camp_Htissue <- Camp_Htissue[,my_QC(exprs(Camp_Htissue), 1000,5000)]

Camp_Mtissue <- Camp_mouse
Camp_Mtissue <- Camp_Mtissue[,my_QC(exprs(Camp_Mtissue), 3000,15000)]

saveRDS(Camp_liver_bud, file="/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/camp_lb_qced.rds")
saveRDS(Camp_transplant, file="/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/camp_trans_qced.rds")
saveRDS(Camp_2D, file="/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/camp_2d_qced.rds")
saveRDS(Camp_Htissue, file="/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/camp_ht_qced.rds")
saveRDS(Camp_Mtissue, file="/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/camp_mt_qced.rds")

# Mouse2Human
source("~/R-Scripts/Ensembl_Stuff.R")
m_expr_mat <- exprs(Camp_Mtissue)
mgenes <- rownames(m_expr_mat)
problems <- grep("Rik",mgenes)
mgenes[problems] <- sub("^X","", mgenes[problems])
ensmusg <- map_symbol_ensg(mgenes, is.org="Mmus", is.name="symbol")
ensg <- map_Hsap_Mmus_one2one(ensmusg, is.org="Mmus")
keep <- grepl("ENSG",ensg);

m_expr_mat<-m_expr_mat[keep,]
ensg<-ensg[keep]
hsymbol <- map_symbol_ensg(ensg, is.org="Hsap", is.name="ensg")
non_zero <- rowMeans(m_expr_mat) > 0
m_expr_mat<- m_expr_mat[non_zero,]
hsymbol<- hsymbol[non_zero]

nomatch <- hsymbol==""
m_expr_mat <- m_expr_mat[!nomatch,]
hsymbol <- hsymbol[!nomatch]

# Deal with duplicates - keep most highly expressed version
dups <- unique(hsymbol[duplicated(hsymbol)])
for (D in dups) {
	thing <- rowMeans(m_expr_mat)
	remove = hsymbol==D & thing < max(thing[hsymbol == D])
	hsymbol <- hsymbol[!remove]
	m_expr_mat<-m_expr_mat[!remove,]
}
rownames(m_expr_mat) <- hsymbol
# Build human-version scater object
PD <- new("AnnotatedDataFrame", data=data.frame(cell_type1=pData(Camp_Mtissue)$cell_type1))
rownames(PD) <- colnames(m_expr_mat)
M_Tissue<- newSCESet(fpkmData=2^m_expr_mat-1, phenoData=PD, logExprsOffset=1)
fData(M_Tissue)$feature_symbol=rownames(m_expr_mat);

saveRDS(Camp_Mtissue, file="/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/camp_mt2h_qced.rds")

##### Do Pipeline - CAMP ONLY
#source("~/NetworkInferencePipeline/Dropouts/My_R_packages/CellTypeProfiles/R/Profiles.R")

require("CellTypeProfiles")

# Get profiles
LB_profiles <- get_cluster_profile(exprs(Camp_liver_bud), factor(as.character(pData(Camp_liver_bud)$cell_type1)), is.log=2, feature_selection=marker.features)
D2_profiles <- get_cluster_profile(exprs(Camp_2D), factor(as.character(pData(Camp_2D)$cell_type1)), is.log=2, feature_selection=marker.features)
H_profiles <- get_cluster_profile(exprs(Camp_Htissue), factor(as.character(pData(Camp_Htissue)$cell_type1)), is.log=2, feature_selection=marker.features)
M_profiles <- get_cluster_profile(m_expr_mat, factor(as.character(pData(Camp_Mtissue)$cell_type1)), is.log=2, feature_selection=marker.features)

#### Plot Consistently #####
raw_truth =  c(colnames(LB_profiles$profiles), colnames(D2_profiles$profiles), colnames(H_profiles$profiles), colnames(M_profiles$profiles))
truth = raw_truth
truth[grep("-endo", truth)] <- "endoderm"
truth[grep("hepa", truth)] <- "hepatocyte"
truth[truth=="lymph"] <- "blood"
truth[truth=="rbc"] <- "blood"
truth[truth=="immune"] <- "blood"
truth[truth=="Kupffer"] <- "blood"
names(truth) <- raw_truth
dataset = rep(c("LB","D2","H","M"), times=c(ncol(LB_profiles$profiles), ncol(D2_profiles$profiles), ncol(H_profiles$profiles), ncol(M_profiles$profiles)))
names(dataset) <- raw_truth


# Get combine & match
matches = combine_and_match_clusters(list(LB=LB_profiles, D2=D2_profiles, HT=H_profiles, MT=M_profiles))

glm_out <- glm_of_matches(matches)
# Plot
png("Camp_only_CTP.png", width=6, height=6, units="in", res=300)
heatout <- cluster_profile_heatmap(glm_out$corrected_profiles, matches, npermute=0, ann=truth, dataset=dataset,features_only=TRUE)
dev.off()
#my_score <- scmap_clade_o_score(heatout, truth)
#saveRDS(scores, file="Liver_ctp_score.rds")

### scmap ###
source("../2_Muli_Dataset_scmap.R")
pData(Camp_liver_bud)$cell_type1 <- paste("LB", pData(Camp_liver_bud)$cell_type1, sep="-")
pData(Camp_2D)$cell_type1 <- paste("2D", pData(Camp_2D)$cell_type1, sep="-")
pData(Camp_Htissue)$cell_type1 <- paste("HT", pData(Camp_Htissue)$cell_type1, sep="-")
pData(M_Tissue)$cell_type1 <- paste("MT", pData(M_Tissue)$cell_type1, sep="-")

Camp_liver_bud<- calculateQCMetrics(Camp_liver_bud)
Camp_2D<- calculateQCMetrics(Camp_2D)
Camp_Htissue<- calculateQCMetrics(Camp_Htissue)
M_Tissue<- calculateQCMetrics(M_Tissue)

SCElist <- list(Camp_liver_bud, Camp_2D, Camp_Htissue, M_Tissue)
all_groups <- c(levels(factor(pData(Camp_liver_bud)$cell_type1)), 
		levels(factor(pData(Camp_2D)$cell_type1)), 
		levels(factor(pData(Camp_Htissue)$cell_type1)), 
		levels(factor(pData(M_Tissue)$cell_type1)))
names(truth) <- all_groups

scmap_out <- do_scmap(SCElist, n_total=length(all_groups))
png("Camp_only_scmap.png", width=6, height=6, units="in", res=300)
heatout <- scmap_heatmap(scmap_out, ann=truth, dataset=dataset)
dev.off()
#scmap_scores <- scmap_clade_o_score(heatout, truth)
#saveRDS(scmap_scores, file="Liver_scmap_score.rds")

### metaneighbour ###
source("../3_Multi_Dataset_MetaNeighbour.R")

#TRUTH = truth;
meta_out <- do_metaneighbour(SCElist)
png("Camp_only_meta.png", width=6, height=6, units="in", res=300)
heatout <- metaneighbour_heatmap(meta_out, truth, dataset=dataset)
dev.off()
#meta_scores <- scmap_clade_o_score(heatout, truth)
#saveRDS(meta_scores, file="Liver_meta_score.rds")

### mnn ###
source("../4_Multi_Datasets_MNN.R")
#features <- unique(c(rownames(LB_profiles)[LB_profiles[,ncol(LB_profiles)]==1], 
#	rownames(D2_profiles)[D2_profiles[,ncol(D2_profiles)]==1],
#	rownames(H_profiles)[H_profiles[,ncol(H_profiles)]==1],
#	rownames(M_profiles)[M_profiles[,ncol(M_profiles)]==1]))
corrected_profiles <- do_mnn(SCElist, features)

png("Camp_only_MNN.png", width=6, height=6, units="in", res=300)
heatout <- cluster_profile_heatmap(corrected_profiles, matches, npermute=0, ann=truth, dataset=dataset,features_only=TRUE)
dev.off()
### CCA ###
#source("../5_Multi_Datasets_CCA.R")
#cca_out <- do_cca(SCElist, n_total=length(all_groups)) # error in align

#png("Camp_only_CCA.png", width=6, height=6, units="in", res=300)
#heatout <- cca_heatmap(cca_out, ann=truth, dataset=dataset)
#dev.off()




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

markers <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Laura_SC3_k6_Markers.rds")

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

CCA1_profile = get_cluster_profile(exprs(CCA1), pData(CCA1)$cell_type1, is.log=2, out.log=2, feature_selection=marker.features)
HCC6_profile = get_cluster_profile(exprs(HCC6), pData(HCC6)$cell_type1, is.log=2, out.log=2, feature_selection=marker.features)
HCC10_profile = get_cluster_profile(exprs(HCC10), pData(HCC10)$cell_type1, is.log=2, out.log=2, feature_selection=marker.features)

rownames(CCA1_profile$profiles) <- fData(CCA1)$feature_symbol
rownames(HCC6_profile$profiles) <- fData(HCC6)$feature_symbol
rownames(HCC10_profile$profiles) <- fData(HCC10)$feature_symbol
rownames(CCA1_profile$is.feature) <- fData(CCA1)$feature_symbol
rownames(HCC6_profile$is.feature) <- fData(HCC6)$feature_symbol
rownames(HCC10_profile$is.feature) <- fData(HCC10)$feature_symbol

#### Laura + Camp ####
#truth_cca1 <- paste("CCA1-",colnames(CCA1_profile[,-ncol(CCA1_profile)]), sep="")
#truth_hcc6 <- paste("HCC6-",colnames(HCC6_profile[,-ncol(HCC6_profile)]), sep="")
#truth_hcc10 <- paste("HCC10-",colnames(HCC10_profile[,-ncol(HCC10_profile)]), sep="")
truth_cca1 <- colnames(CCA1_profile$profiles)
truth_hcc6 <- colnames(HCC6_profile$profiles)
truth_hcc10 <- colnames(HCC10_profile$profiles)
truth = c(truth_cca1, truth_hcc6, truth_hcc10,
        colnames(LB_profiles$profiles),
        colnames(D2_profiles$profiles),
        colnames(H_profiles$profiles), 
	colnames(M_profiles$profiles))
dataset = rep(c("T", "T", "T", "CL", "CL", "HT", "MT"),
                times=c(ncol(CCA1_profile$profiles), ncol(HCC6_profile$profiles), 
			ncol(HCC10_profile$profiles), ncol(LB_profiles$profiles), 
			ncol(D2_profiles$profiles), ncol(H_profiles$profiles), 
			ncol(M_profiles$profiles)))

truth1 <- truth
truth1[truth1 %in% c("CCA1-4", "CCA1-1", "HCC6-4", "HCC10-1", "CCA1-2")] = "L-Prolif"
truth1[truth1 %in% c("CCA1-6", "HCC6-5")] = "L-Diff"
truth1[grep("CCA1", truth1)] = "L-Inter"
truth1[grep("HCC6", truth1)] = "L-Inter"
truth1[grep("HCC10", truth1)] = "L-Inter"
truth1[grep("hepa", truth1)] = "hepatocyte"
truth1[truth1 %in% c("rbc", "immune", "lymph", "Kupffer")] = "blood"
truth1[grep("-endo", truth1)] = "endoderm"
truth1 <- factor(truth1, levels=c("blood", "endoderm", "hepatocyte", "L-Prolif", "L-Diff", "L-Inter", "msc", "ipsc", "stellate", "endoth"))
dataset <- factor(dataset)


truth1_col <- c("rosybrown1", "#bae4b3", "red3", "black", "green", "cyan", "#d7b5d8","#525252","#2171b5", "#ffffb2")
#dataset_col <- c("#c2e699","#df65b0","#b3cde3","#8c96c6","#88419d","#238443","black")
dataset_col <- c("#c2e699","#df65b0","black", "forestgreen")

my_heatmap <- function(sim_matrix) {
        sim_matrix <- sim_matrix/max(abs(sim_matrix)); # [-1,1]
        sim_matrix <- (sim_matrix+t(sim_matrix))/2;

        heatcols <- colorRampPalette(c("blue","khaki1","red"))(255)

        dendro<-hclust(as.dist(1-sim_matrix))

        require("CellTypeProfiles")
        ColumnCols <- matrix(truth1_col[truth1], ncol=1); colnames(ColumnCols) = c("CellType")
        RowCols <- matrix(dataset_col[dataset], nrow=1); rownames(RowCols) = c("Dataset");

        heatout <- heatmap.3(sim_matrix, trace="n", scale="none", col=heatcols, symbreaks=FALSE,
                        key.title="", key.xlab="Similarity", symm=TRUE,
                        Rowv=as.dendrogram(dendro), Colv=as.dendrogram(dendro),
                        ColSideColorsSize=length(ColumnCols[1,]), RowSideColorsSize=length(RowCols[,1]),
                        ColSideColors=as.matrix(ColumnCols), RowSideColors=RowCols)
        return(list(heatmap_out=heatout, dist_mat=sim_matrix));
}





profile_List = list(CCA1=CCA1_profile, HCC6=HCC6_profile, HCC10=HCC10_profile, LB=LB_profiles, D2=D2_profiles, HT=H_profiles, MT=M_profiles)

matches = combine_and_match_clusters(profile_List)
corrected_profiles = glm_of_matches(matches)

ctp_sim_mat <- cor(corrected_profiles$corrected_profiles, method="spearman")
png("Laura_Camp_CTP.png", width=6, height=6, units="in", res=300)
heatout <- my_heatmap(ctp_sim_mat)
#ann = truth;
#ann[grep("CCA1", truth)] = "Line"
#ann[grep("HCC6", truth)] = "Line"
#ann[grep("HCC10", truth)] = "Line"
#cluster_profile_heatmap(corrected_profiles, matches, npermute=0, ann=ann, dataset=dataset)
dev.off()

# scmap 
source("../2_Muli_Dataset_scmap.R")
all_groups = c(levels(factor(pData(Camp_liver_bud)$cell_type1)),
                levels(factor(pData(Camp_2D)$cell_type1)),
                levels(factor(pData(Camp_Htissue)$cell_type1)),
                levels(factor(pData(M_Tissue)$cell_type1)),
		levels(factor(pData(CCA1)$cell_type1)),
                levels(factor(pData(HCC6)$cell_type1)),
                levels(factor(pData(HCC10)$cell_type1)))

SCElist <- list(CCA1, HCC6, HCC10, Camp_liver_bud, Camp_2D, Camp_Htissue, M_Tissue)
scmap_out <- do_scmap(SCElist, n_total=length(all_groups))

png("Laura_Camp_scmap.png", width=6, height=6, units="in", res=300)
heatout <- my_heatmap(scmap_out)
dev.off()

# meta
raw_truth <- c(levels(factor(pData(CCA1)$cell_type1)),
		levels(factor(pData(HCC6)$cell_type1)),
		levels(factor(pData(HCC10)$cell_type1)),
		levels(factor(pData(Camp_liver_bud)$cell_type1)),
		levels(factor(pData(Camp_2D)$cell_type1)),
		levels(factor(pData(Camp_Htissue)$cell_type1)),
		levels(factor(pData(M_Tissue)$cell_type1)))

source("../3_Multi_Dataset_MetaNeighbour.R")
rownames(CCA1) <- fData(CCA1)$feature_symbol
rownames(HCC6) <- fData(HCC6)$feature_symbol
rownames(HCC10) <- fData(HCC10)$feature_symbol
SCElist <- list(CCA1, HCC6, HCC10, Camp_liver_bud, Camp_2D, Camp_Htissue, M_Tissue)
meta_out <- do_metaneighbour(SCElist)

png("Laura_Camp_meta.png", width=6, height=6, units="in", res=300)
heatout <- my_heatmap(meta_out)
dev.off()

#mnn
source("../4_Multi_Datasets_MNN.R")
mnn_out <- do_mnn(SCElist)
mnn_sim_mat <- as.matrix(dist(t(mnn_out)))
mnn_sim_mat <- 1-mnn_sim_mat/max(mnn_sim_mat)

png("Laura_Camp_mnn.png", width=6, height=6, units="in", res=300)
heatout <- my_heatmap(mnn_sim_mat)
dev.off()

#cca
source("../5_Multi_Datasets_CCA.R")
cca_out <- do_cca(SCElist, n_total=length(all_groups));

png("Laura_Camp_cca.png", width=6, height=6, units="in", res=300)
heatout <- my_heatmap(cca_out)
dev.off()

png("Laura_Camp_legends.png", width=6, height=6, units="in", res=300)
plot(1,1, col="white", bty="n", xaxt="n", yaxt="n")
legend("topleft", title="Celltype", levels(truth1), fill=truth1_col, bty="n")
legend("topright", title="Dataset", levels(dataset), fill=dataset_col, bty="n")
dev.off()



###### Laura Only at single-cell level correction ######
### mnn - Laura only ###
#CCA1_kw <- apply(exprs(CCA1), 1, function(x) {kruskal.test(x ~ pData(CCA1)$sc3_6_clusters)$p.value})
#HCC6_kw <- apply(exprs(HCC6), 1, function(x) {kruskal.test(x ~ pData(HCC6)$sc3_6_clusters)$p.value})
#HCC10_kw <- apply(exprs(HCC10), 1, function(x) {kruskal.test(x ~ pData(HCC10)$sc3_6_clusters)$p.value})

# Marker genes = features
#sig_threshold <- 0.05/(nrow(exprs(CCA1))*6)
#features <- unique(c(markers$CCA1$Gene[markers$CCA1$pvalue < sig_threshold], markers$HCC6$Gene[markers$HCC6$pvalue < sig_threshold], markers$HCC10$Gene[markers$HCC10$pvalue < sig_threshold]))

markers = list();
markers$CCA1 <- CCA1_profile$is.feature
markers$HCC6 <- HCC6_profile$is.feature
markers$HCC10 <- HCC10_profile$is.feature

cca1_marks <- rownames(markers$CCA1)[rowSums(markers$CCA1) > 0]
hcc6_marks <- rownames(markers$HCC6)[rowSums(markers$HCC6) > 0]
hcc10_marks <- rownames(markers$HCC10)[rowSums(markers$HCC10) > 0]

features <- unique(c( cca1_marks, hcc6_marks, hcc10_marks ))

# Fix gene names
source("~/R-Scripts/Ensembl_Stuff.R")
#features_symbol <- map_symbol_ensg(features, is.org="Hsap", is.name="ensg")
#features[features_symbol != ""] <- features_symbol[features_symbol != ""]

# Colours
Laura_obj <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ForReviewPaper_stuff.rds")
dataset_cols <- Laura_obj$colours
dataset_pch <- c(1, 18, 7)
require("RColorBrewer")
clust_cols <- brewer.pal(6, "Set2")
prolif_cols = c("black","red")

## MNN Correct ##
laura_mnn_out <- mnnCorrect(list(CCA1=exprs(CCA1), HCC6=exprs(HCC6), HCC10=exprs(HCC10)), hvg.genes=which(rownames(exprs(CCA1)) %in% features))
Combined_corrected <- cbind(laura_mnn_out$corrected$CCA1, laura_mnn_out$corrected$HCC6, laura_mnn_out$corrected$HCC10)

## CTP ##
CCA1_profile = get_cluster_profile(exprs(CCA1), pData(CCA1)$cell_type1, is.log=2, out.log=2, feature_selection=marker.features)
HCC6_profile = get_cluster_profile(exprs(HCC6), pData(HCC6)$cell_type1, is.log=2, out.log=2, feature_selection=marker.features)
HCC10_profile = get_cluster_profile(exprs(HCC10), pData(HCC10)$cell_type1, is.log=2, out.log=2, feature_selection=marker.features)

profile_List = list(CCA1=CCA1_profile, HCC6=HCC6_profile, HCC10=HCC10_profile)

matches = combine_and_match_clusters(profile_List)
corrected_profiles = glm_of_matches(matches)

cca1_ctp_norm <- correct_sng_cells(CCA1_profile$norm_mat, "CCA1", corrected_profiles)
hcc6_ctp_norm <- correct_sng_cells(HCC6_profile$norm_mat, "HCC6", corrected_profiles)
hcc10_ctp_norm <- correct_sng_cells(HCC10_profile$norm_mat, "HCC10", corrected_profiles)

Combined_ctp_norm <- cbind(cca1_ctp_norm, hcc6_ctp_norm, hcc10_ctp_norm);

Laura_CTP_stuff <- list(profiles=profile_List, matches=matches, glm_out=corrected_profiles, sng_cell_norm=Combined_ctp_norm);
saveRDS(Laura_CTP_stuff, file="/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/Laura_CTP_stuff.rds")

# Merge
Combined_raw <- cbind(exprs(CCA1), exprs(HCC6), exprs(HCC10))
dataset <- rep(c(1,2,3), times=c(ncol(laura_mnn_out$corrected$CCA1),ncol(laura_mnn_out$corrected$HCC6), ncol(laura_mnn_out$corrected$HCC10)))
clust_ID <- c(CCA1$sc3_6_clusters, HCC6$sc3_6_clusters, HCC10$sc3_6_clusters)
prolif <- c(CCA1$sc3_6_clusters %in% c("4","1"), HCC6$sc3_6_clusters ==4,  HCC10$sc3_6_clusters==1 )

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
#cca1_marks_ensg <- markers$CCA1$Gene[markers$CCA1$pvalue < sig_threshold]
#hcc6_marks_ensg <- markers$HCC6$Gene[markers$HCC6$pvalue < sig_threshold]
#hcc10_marks_ensg <- markers$HCC10$Gene[markers$HCC10$pvalue < sig_threshold]

#cca1_marks <- map_symbol_ensg(cca1_marks_ensg, is.org="Hsap", is.name="ensg")
#hcc6_marks <- map_symbol_ensg(hcc6_marks_ensg, is.org="Hsap", is.name="ensg")
#hcc10_marks <- map_symbol_ensg(hcc10_marks_ensg, is.org="Hsap", is.name="ensg")
#cca1_marks[cca1_marks == ""] <- cca1_marks_ensg[cca1_marks == ""]
#hcc6_marks[hcc6_marks == ""] <- hcc6_marks_ensg[hcc6_marks == ""]
#hcc10_marks[hcc10_marks == ""] <- hcc10_marks_ensg[hcc10_marks == ""]
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
   #     ColumnCols <- matrix(truth1_col[truth1], ncol=1); colnames(ColumnCols) = c("CellType")
   #     RowCols <- matrix(dataset_col[dataset], nrow=1); rownames(RowCols) = c("Dataset");

        heatout <- heatmap.3(sim_matrix, trace="n", scale="none", col=heatcols, symbreaks=FALSE,
                        key.title="", key.xlab=xlab_name, symm=TRUE,
                        Rowv=as.dendrogram(dendro), Colv=as.dendrogram(dendro))
    #                    ColSideColorsSize=length(ColumnCols[1,]), RowSideColorsSize=length(RowCols[,1]),
    #                    ColSideColors=as.matrix(ColumnCols), RowSideColors=RowCols)
    #    return(list(heatmap_out=heatout, dist_mat=sim_matrix));
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
# Same features as mnn
sig_threshold <- 0.05/(nrow(exprs(CCA1))*6)
features <- unique(c(cca1_marks, hcc6_marks, hcc10_marks))
# Fix gene names
#features_symbol <- map_symbol_ensg(features, is.org="Hsap", is.name="ensg")
#features[features_symbol != ""] <- features_symbol[features_symbol != ""]

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
