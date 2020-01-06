# GO enrichments of markers
# External markers
# diff CC phase of same functional cell-type?
source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("/nfs/users/nfs_t/ta6/R-Scripts/Blank_plot.R")

args <- commandArgs(trailingOnly=TRUE) # SCE RDSs
nSCEs <- length(args)
expr_type <- "lognorm"

type_col <- type_col[1:nSCEs]

SCE_list <- list();
keep_genes <- c()
background_genes <- c()
max_ngenes <- 0;
for (f in args) {
	require("scater")
	obj <- readRDS(f);
	if (nrow(obj) > max_ngenes) {max_ngenes <- nrow(obj)}
	keep_genes <- c(keep_genes, as.character(fData(obj)[ fData(obj)$pct_dropout < 90, "Symbol"]));
	background_genes <- c(background_genes, rownames(obj)[ fData(obj)$pct_dropout <= 99]);
	tmp <- unlist(strsplit(f, "\\."))
	SCE_list[[tmp[1]]] <- obj;
}
keep_genes <- unique(keep_genes);
background_genes <- unique(background_genes);


## External Markers ##

# Read in Putative Markers
Lit_Markers <- read.delim("~/Collaborations/LiverOrganoids/Literature_Markers.txt", sep="\t", header=TRUE)


CC_Stem <- read.delim("~/Collaborations/LiverOrganoids/Markers_CC_Stemness.txt", sep=" ", header=TRUE)
HCC_Stem <- read.delim("~/Collaborations/LiverOrganoids/Markers_HCC_Stemness.txt", sep=" ", header=FALSE)
CHC_Stem <- read.delim("~/Collaborations/LiverOrganoids/Markers_CHC_Stemness.txt", sep=" ", header=TRUE)


Hepato_MSigdb <- read.delim("~/Collaborations/LiverOrganoids/Markers_HCC_Differentiation.txt")

HPA <- readRDS("/lustre/scratch117/cellgen/team218/TA/OtherDownloadedData/Good_HPA_Hep_Chol_Markers.rds")

Samp <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExternalData/Sampaziotis_MarkersV2.rds")

Hepato_Camp <- read.table("~/Collaborations/LiverOrganoids/FromLaura_Camp_Hepatocyte_Markers.txt")
Stem_Camp <- read.table("~/Collaborations/LiverOrganoids/FromLaura_Camp_Hepatoblast_Markers.txt")
MSC_Camp <- read.table("~/Collaborations/LiverOrganoids/FromLaura_Camp_Mesenchymal_Markers.txt")

CC_1 <- read.table("~/Data/Whitfield_CC.txt")
CC_2 <- read.table("~/Collaborations/LiverOrganoids/New_CC_171117.txt", header=FALSE)
CC_3 <- read.delim("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/GO/hsapiens_80_GO_Annoations_Emsembl.out", sep="\t", header=FALSE)
CC_3 <- CC_3[CC_3[,3] == "cell cycle",1]
source("~/R-Scripts/Ensembl_Stuff.R")
CC_3 <- unique(map_symbol_ensg(as.character(CC_3), is.org="Hsap", is.name="ensg"))

All_CC <- unique(c(as.character(CC_1[,2]), as.character(CC_2[,1]), as.character(CC_3)))

All_Hepato <-c(unique(as.character(Hepato_MSigdb[,1])),
	unique(as.character(Hepato_Camp[,1])),
	unique(as.character(names(Samp$HB))),
	#as.character(HPA[HPA[,8] < HPA[,16],2]),
	unique(as.character(Lit_Markers[Lit_Markers[,2]=="hepatocyte",1]))
		)
Good_Hepato <- All_Hepato[duplicated(All_Hepato)]
All_Chol <- c(unique(as.character(HPA[HPA[,8] > HPA[,16],2])),
	unique(as.character(Lit_Markers[Lit_Markers[,2] %in% c("cholangiocyte", "mature cholangiocyte"),1])),
	unique(as.character(names(Samp$mChol)))
	)
Good_Chol <- All_Chol[duplicated(All_Chol)]
All_Stem <- c(unique(as.character(CC_Stem[CC_Stem[,2]==1,1])),
	unique(as.character(CHC_Stem[CHC_Stem[,2]==1,1])),
	unique(as.character(HCC_Stem[HCC_Stem[,2]==1,1])),
	unique(as.character(Stem_Camp[,1])), 
	unique(names(Samp$Stem)))
Good_Stem <- All_Stem[duplicated(All_Stem)]

All_MSC <- c(as.character(MSC_Camp[,1]))

All_Hepato <- unique(All_Hepato); All_Hepato <- All_Hepato[!is.na(All_Hepato)]
All_Chol <- unique(All_Chol); All_Chol <- All_Chol[!is.na(All_Chol)]
All_Stem <- unique(All_Stem); All_Stem <- All_Stem[!is.na(All_Stem)]
All_MSC <- unique(All_MSC); All_MSC <- All_MSC[!is.na(All_MSC)]

# Remove any gene that appears on conflicting lists
All_Markers <- c(All_Hepato, All_Chol, All_Stem, All_MSC, All_CC)
exclude <- All_Markers[duplicated(All_Markers)]
All_Markers <- unique(All_Markers); All_Markers <- All_Markers[!All_Markers %in% exclude];
All_Hepato <- All_Hepato[!All_Hepato %in% exclude]
All_Chol <- All_Chol[!All_Chol %in% exclude]
All_Stem <- All_Stem[!All_Stem %in% exclude]
All_MSC <- All_MSC[!All_MSC %in% exclude]

for (i in 1:nSCEs) {
	obj <- SCE_list[[i]]
	fData(obj)$MarkCellType <- rep("", times=nrow(obj));
	fData(obj)$MarkCellType[fData(obj)$Symbol %in% All_Hepato] <- "Hepato"
	fData(obj)$MarkCellType[fData(obj)$Symbol %in% All_Chol] <- "Chol"
	fData(obj)$MarkCellType[fData(obj)$Symbol %in% All_Stem] <- "Stem"
	fData(obj)$MarkCellType[fData(obj)$Symbol %in% All_MSC] <- "MSC"
	fData(obj)$MarkCellType[fData(obj)$Symbol %in% All_CC] <- "CC"

	fData(obj)$MarkCellType2[fData(obj)$Symbol %in% Good_Hepato] <- "Hepato"
	fData(obj)$MarkCellType2[fData(obj)$Symbol %in% Good_Chol] <- "Chol"
	fData(obj)$MarkCellType2[fData(obj)$Symbol %in% Good_Stem] <- "Stem"
	SCE_list[[i]] <- obj;
}

# Test if overrepresented in kept genes:
kept_test <- function(x) {
	q <- sum(x %in% keep_genes)
	m <- length(keep_genes)
	n <- max_ngenes-m
	k <- length(x)
	return(phyper(q, m, n, k, lower.tail=FALSE))
}
if (kept_test(All_Hepato) > 10^-5) {print("Warning Hepato markers not enriched among kept genes")}
if (kept_test(All_Chol) > 10^-5) {print("Warning Chol markers not enriched among kept genes")}
if (kept_test(All_Stem) > 10^-5) {print("Warning Stem markers not enriched among kept genes")}
if (kept_test(All_MSC) > 10^-5) {print("Warning MSC markers not enriched among kept genes")}

### Move the below to a separate script I think ###
# Should I look for GO enrichments for each cell-type? then weed out genes not associated with an enriched function?


# Remove low & non-DE genes from the lists
All_Hepato <- All_Hepato[All_Hepato %in% keep_genes]
All_Chol <- All_Chol[All_Chol %in% keep_genes]
All_Stem <- All_Stem[All_Stem %in% keep_genes]
All_MSC <- All_MSC[All_MSC %in% keep_genes]
All_Markers <- All_Markers[All_Markers %in% keep_genes]
All_Markers <- All_Markers[!All_Markers %in% All_CC]

require("CellTypeProfiles");

	
# gene-gene correlations of markers -> cluster gene-gene network and ID clusters enriched for markers with a particular annotation. 
# Clustering = complete, theshold = ?

mark_cols = c("forestgreen", "navy", "purple", "salmon")

Score <- rep(0, times=length(All_Markers))
names(Score) <- All_Markers
for (i in 1:nSCEs) {
	obj <- SCE_list[[i]]
	obj <- calculateQCMetrics(obj)
	obj <- obj[fData(obj)$Symbol %in% All_Markers & fData(obj)$pct_dropout < 99,]
	expr_mat <- get_exprs(obj, expr_type)
	rownames(expr_mat) <- fData(obj)$Symbol
	c_mat <- cor(t(expr_mat), method="spearman")
#	heatout <- heatmap.3(c_mat, distfun=function(x){as.dist(1-x)}, 
#		hclustfun=function(x){hclust(x, method="ward.D")},
#		ColSideColors=matrix(mark_cols[factor(fData(obj)$MType)], nrow=nrow(c_mat)), ColSideColorsSize=1)
	score <- apply(c_mat, 1, min)
	score2 <- apply(c_mat, 1, function(x){max(x[x < 1])})
	overall <- apply(cbind(rank(-1*score),rank(score2)), 1, min)
	overall <- overall[match(names(Score), names(overall))]
	overall[is.na(overall)] <- 0
	Score <- Score+overall
}

Score <- rev(sort(Score))

all_plot <- function(gene) {
	require("RColorBrewer")
	get_splits <- function(gene) {
		vals = c();
		for (obj in SCE_list) {
			mat <- get_exprs(obj, expr_type)
			vals <- c(vals, mat[fData(obj)$Symbol == gene, ]) 
		}
		splits <- seq(from=min(vals), to=max(vals), length=5)
		return(splits);
	}

	splits <- get_splits(gene);
	a <- ceiling(sqrt(length(SCE_list)))
	par(mfrow=c(a,a))
	par(mar=c(1,4,3,1))
	
	for( i in 1:length(SCE_list) ) {
		obj <- SCE_list[[i]]
		binned <- cut( get_exprs(obj, expr_type)[fData(obj)$Symbol == gene,], breaks=splits, include.lowest=TRUE )
		my_cols <- colorRampPalette(brewer.pal(6, "Blues"))(length(splits)-1)

		my_title <- names(SCE_list)[i]
		my_title <- unlist(strsplit(my_title, "_"))
		my_title <- my_title[1]
		pt_col <- my_cols[binned]
		plot(obj$DM1, obj$DM2, pch=16, main=my_title, xaxt="n", yaxt="n", ylab=gene, bty="n", col=pt_col)
	} 	
}

# Or is this more informative:
Hep_score <- vector()
Chol_score <- vector()
Stem_score <- vector()
donor_id <- vector()
for (i in 1:length(SCE_list)) {
	obj <- SCE_list[[i]]
	H_score <- colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Hepato,])
	C_score <- colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Chol,])
	S_score <- colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Stem,])
	pData(obj)$ScoreHep <- H_score
	pData(obj)$ScoreChol <- C_score
	pData(obj)$ScoreStem <- S_score
	SCE_list[[i]] <- obj;
	Hep_score <- c(Hep_score, H_score)
	Chol_score <- c(Chol_score, C_score)
	Stem_score <- c(Stem_score, S_score)
	donor_id <- c(donor_id, rep(names(SCE_list)[i], times=ncol(obj)))
}

dataset_names <- names(SCE_list); dataset_names <- matrix(unlist(strsplit(dataset_names, "_")), ncol=2, byrow=TRUE)[,1];

#plot(Stem_score, Hep_score, col=c("forestgreen", "purple", "black")[factor(donor_id)], pch=16)
#plot(Chol_score, Hep_score, col=c("forestgreen", "purple", "black")[factor(donor_id)], pch=16)
#plot(Chol_score, Stem_score, col=c("forestgreen", "purple", "black")[factor(donor_id)], pch=16)
png("All_CellType_Scores.png", width=7, height=7, units="in", res=300)
layout(matrix(c(1,2,1,2), nrow=2, byrow=TRUE), widths=c(6,1))
par(mar=c(4,4,1,1))
plot(Chol_score-Hep_score, Stem_score, col=type_col[factor(donor_id)], pch=16, xlab="hepatocyte ----- cholangiocyte", ylab="stem-ness")
reset_mar <- blank_plot()
legend("right", levels(factor(dataset_names)), pch=16, col=type_col, bty="n")
par(mar=reset_mar);
dev.off()

### Lineage Score v2 ###
Chol_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Chol.txt", header=TRUE)
Hep_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Hep.txt", header=TRUE)

# Remove overlapping markers
exclude <- Hep_lineage[(Hep_lineage[Hep_lineage[,2] != "Prog",1] %in% Chol_lineage[Chol_lineage[,2] != "Prog",1]),1]
common_stem <- Hep_lineage[(Hep_lineage[Hep_lineage[,2] == "Prog",1] %in% Chol_lineage[Chol_lineage[,2] == "Prog",1]),1]
write.table(common_stem, "Common_Stem.txt")

Chol_lineage <- Chol_lineage[!Chol_lineage[,1] %in% exclude, ]
Hep_lineage <- Hep_lineage[! Hep_lineage[,1] %in% exclude, ]
Chol_lineage[,2] <- factor(Chol_lineage[,2], levels=c("Prog", "Chol"))
Hep_lineage[,2] <- factor(Hep_lineage[,2], levels=c("Prog", "Hep"))

lin_keep="";
for (i in 1:length(SCE_list)) {
	obj <- SCE_list[[i]]
	lin_keep <- c(lin_keep, as.character(fData(obj)$Symbol[fData(obj)$fine_marker_q.value < 0.05 & fData(obj)$fine_marker_AUC > 0.6]))
}
lin_keep <- unique(lin_keep)
control_genes <- keep_genes[!keep_genes %in% lin_keep]

Chol_lineage <- Chol_lineage[Chol_lineage[,1] %in% lin_keep,]
Hep_lineage <- Hep_lineage[Hep_lineage[,1] %in% lin_keep,]

for (i in 1:length(SCE_list)) {
	require("matrixStats")
	obj <- SCE_list[[i]]
	mat <- get_exprs(obj, "lognorm")
	Cholm <- fData(obj)$Symbol %in% Chol_lineage[Chol_lineage[,2] =="Chol",1]
	Cholp <- fData(obj)$Symbol %in% Chol_lineage[Chol_lineage[,2] =="Prog",1]
	Hepm <- fData(obj)$Symbol %in% Hep_lineage[Hep_lineage[,2] =="Hep",1]
	Hepp <- fData(obj)$Symbol %in% Hep_lineage[Hep_lineage[,2] =="Prog",1]
	common_prog <- fData(obj)$Symbol %in% common_stem
	Cholm <- colMeans(mat[Cholm,]);
	Cholp <- colMeans(mat[Cholp,]);
	Hepm <- colMeans(mat[Hepm,]); 
	Hepp <- colMeans(mat[Hepp,]); 
	common_prog <- colMeans(mat[common_prog,]); 
	control <- colMeans(mat[fData(obj)$Symbol %in% control_genes,])

	pData(obj)$lin_score_Chol <- Cholm
	pData(obj)$lin_score_Hep <- Hepm
	pData(obj)$lin_score_Chol_Prog <- Cholp
	pData(obj)$lin_score_Hep_Prog <- Hepp
	pData(obj)$Lineage_x <- Cholm/Cholp; pData(obj)$Lineage_x[Hepm > Cholm] <- (-Hepm/Hepp)[Hepm > Cholm];
	pData(obj)$Lineage_y <- (Cholp+Hepp)/2

	cluster_col_set <- get_group_cols(obj)

	png(paste(names(SCE_list)[i], "Lineage.png", sep="_"), width=8, height=4.5, units="in", res=300)
	par(mfrow=c(1,2))
	par(mar=c(4,4,1,1))
	plot(pData(obj)$lin_score_Chol, pData(obj)$lin_score_Chol_Prog, xlab="Chol Mature", ylab="Chol Prog", col=cluster_col_set[obj$clusters_clean], pch=16)
	plot(pData(obj)$lin_score_Hep, pData(obj)$lin_score_Hep_Prog, xlab="Hep Mature", ylab="Hep Prog", col=cluster_col_set[obj$clusters_clean], pch=16)
	dev.off()

	SCE_list[[i]] <- obj;
}

Linx <- vector()
Liny <- vector()
donor_id <- vector()
for (i in 1:length(SCE_list)) {
        obj <- SCE_list[[i]]
        Linx <- c(Linx, pData(obj)$Lineage_x)
        Liny <- c(Liny, pData(obj)$Lineage_y)
        donor_id <- c(donor_id, rep(names(SCE_list)[i], times=ncol(obj)))
}

dataset_names <- names(SCE_list); dataset_names <- matrix(unlist(strsplit(dataset_names, "_")), ncol=2, byrow=TRUE)[,1];

png("All_Lineage.png", width=8/1.2, height=7/1.2, units="in", res=300)
layout(matrix(c(1,2,1,2), nrow=2, byrow=TRUE), widths=c(6,1))
par(mar=c(4,4,1,1))
plot(Linx, Liny, col=type_col[factor(donor_id)], pch=16, xlab="hepatocyte ----- cholangiocyte", ylab="stem-ness")
reset_mar <- blank_plot()
legend("right", levels(factor(dataset_names)), pch=16, col=type_col, bty="n")
par(mar=reset_mar);
dev.off()

	


for (i in 1:length(SCE_list)) {
	saveRDS(SCE_list[[i]], file=paste(names(SCE_list)[i], "Function.rds", sep="_"))
}

