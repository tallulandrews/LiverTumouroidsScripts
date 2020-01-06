# GO enrichments of markers
# External markers
# diff CC phase of same functional cell-type?
source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("/nfs/users/nfs_t/ta6/R-Scripts/Blank_plot.R")

require("scater")

# Set Up CCA1 and HCC10 as "Truth"?
CCA1 <- readRDS("CCA1_merged_SC3.rds")
HCC10 <- readRDS("HCC10_merged_SC3.rds")
expr_type <- "lognorm"
cycle_pch=rev(c(1, 12, 15, 17))
                                              

nSCEs <- 2
type_col <- type_col[1:nSCEs]

get_features <- function(SCE) {
	SCE <- SCE[rowData(SCE)$pct_dropout < 90,]
	#SCE <- SCE[rowData(SCE)$biotype =="protein_coding",]
	SCE <- SCE[rowData(SCE)$Symbol != "",]
	KW.base <- apply(assays(SCE)[[ "lognorm" ]], 1, function(x) {
                kruskal.test(x, colData(SCE)[,"Clusters"])$p.value})
	KW.norm <- apply(assays(SCE)[[ "norm_exprs" ]], 1, function(x) {
                kruskal.test(x, colData(SCE)[,"Clusters_normCC"])$p.value})
	KW.no <- apply(assays(SCE)[[ "noCC" ]], 1, function(x) {
                kruskal.test(x, colData(SCE)[,"Clusters_noCC"])$p.value})

	sig <- p.adjust(KW.base, method="bon") < 0.05 |  p.adjust(KW.norm, method="bon") < 0.05 |  p.adjust(KW.no, method="bon") < 0.05 | rowData(SCE)$noCC_fine_marker_is.Feature | rowData(SCE)$fine_marker_is.Feature

	sig[is.na(sig)] <- FALSE
	return(as.character(rowData(SCE)[sig,"Symbol"]))
}


CCA1 <- CCA1[rowData(CCA1)[,"Symbol"] %in% get_features(CCA1),]
HCC10 <- HCC10[rowData(HCC10)[,"Symbol"] %in% get_features(HCC10),]

### Clean Up Lineage Genes ###
Chol_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Chol.txt", header=TRUE)
Hep_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Hep.txt", header=TRUE)

# Remove overlapping markers
Hep_both <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & grepl("Hep", Hep_lineage[,2]), 1]
Chol_both <- Chol_lineage[ grepl("Prog", Chol_lineage[,2]) & grepl("Chol", Chol_lineage[,2]), 1]
Prog_both <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Prog",1], 1]

Conflict1 <- Hep_lineage[ grepl("Hep", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Chol",1], 1]
Conflict2 <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Chol",1], 1]
Conflict3 <- Hep_lineage[ grepl("Hep", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Prog",1], 1]

Conflicts <- c(as.character(Conflict1), as.character(Conflict2), as.character(Conflict3), Chol_both, Hep_both)

#Chol_lineage <- Chol_lineage[Chol_lineage[,1] %in% marker_genes & Chol_lineage[,1] %in% keep_genes,]
#Hep_lineage <- Hep_lineage[Hep_lineage[,1] %in% marker_genes & Hep_lineage[,1] %in% keep_genes,]

Chol_lineage[,2] <- as.character(Chol_lineage[,2])
Chol_lineage[Chol_lineage[,2] == "Prog",2] <- "Chol-Prog" 
Chol_lineage[Chol_lineage[,2] == "Chol",2] <- "Chol-Mature" 
Chol_lineage <- Chol_lineage[!(Chol_lineage[,1] %in% Conflicts) & !(Chol_lineage[,1] %in% Chol_both),]

Hep_lineage[,2] <- as.character(Hep_lineage[,2])
Hep_lineage[Hep_lineage[,2] == "Prog",2] <- "Hep-Prog" 
Hep_lineage[Hep_lineage[,2] == "Hep",2] <- "Hep-Mature" 
Hep_lineage <- Hep_lineage[!(Hep_lineage[,1] %in% Conflicts),]

Lineage <- rbind(Chol_lineage,Hep_lineage)
Lineage[Lineage[,1] %in% Prog_both,2] <- "Common-Prog"
Lineage <- Lineage[!duplicated(Lineage[,1]),]


### ICA-based ####
require("fastICA")

do_ICA <- function(SCE, n, expr_type) {
	mat <- assays(SCE)[[expr_type]]
	ica <- fastICA(mat, n.comp=n, method="C")
	gene_symbols <- rowData(SCE)$Symbol
	rownames(ica$S) <- rowData(SCE)$Symbol

	tab <- rbind(
		colMeans(ica$S[gene_symbols %in% Lineage[ Lineage[,2] == "Chol-Prog" ,1],]),
		colMeans(ica$S[gene_symbols %in% Lineage[ Lineage[,2] == "Chol-Mature" ,1],]),
		colMeans(ica$S[gene_symbols %in% Lineage[ Lineage[,2] == "Hep-Prog" ,1],]),
		colMeans(ica$S[gene_symbols %in% Lineage[ Lineage[,2] == "Hep-Mature" ,1],]),
		colMeans(ica$S[gene_symbols %in% Lineage[ Lineage[,2] == "Common-Prog" ,1],]))
	rownames(tab) <- c("C-Prog", "C", "H-Prog", "H", "Prog")
	print(tab);

	# Which dimensions correspond to stemness vs lineage?
	score_prog <- (tab[1,]-tab[2,] + tab[3,]-tab[4,] + tab[5,]-tab[2,] + tab[5,]-tab[4,])/4
	score_lin <- (tab[2,]-tab[4,])
	threshold <- mean(abs(tab))
	rel_prog <- abs(score_prog)/threshold
	rel_lin <- abs(score_lin)/threshold
	prog_dim <- which(rel_prog == max(rel_prog) & rel_prog > 1)
	lin_dim <- which(rel_lin == max(rel_lin) & rel_lin > 1)
	print(paste("Prog:", prog_dim, collapse=" "))
	print(paste("Lin:", lin_dim, collapse=" "))
	if (score_prog[prog_dim] < 0) {prog_dim <- -1*prog_dim}	
	if (score_lin[lin_dim] > 0) {lin_dim <- -1*lin_dim}	

	return(list(tab=tab, ica=ica, lin=lin_dim, prog=prog_dim));
}

ica_plot <- function(ica, cols, pchs=16) {
	plot(ica$ica$A[abs(ica$lin),]*sign(ica$lin),         
		ica$ica$A[abs(ica$prog),]*sign(ica$prog),          
		col=cols, pch=pchs, xlab="Chol -------- Hep", ylab="Diff ------- Stem")
}

### CCA1 ###

# Base
set.seed(20391)
icas_cca1_base <- do_ICA(CCA1, 2, "lognorm")
ica_plot(icas_cca1_base, get_group_cols(CCA1, "Clusters")[CCA1$Clusters], cycle_pch[CCA1$Cycle])

# Norm
set.seed(20391)
icas_cca1_nCC<- do_ICA(CCA1, 2, "norm_exprs")
ica_plot(icas_cca1_nCC, get_group_cols(CCA1, "Clusters_normCC")[CCA1$Clusters_normCC], cycle_pch[CCA1$Cycle])

# noCC
set.seed(20391)
icas_cca1_noCC <- do_ICA(CCA1, 3, "noCC")
ica_plot(icas_cca1_noCC, get_group_cols(CCA1, "Clusters_noCC")[CCA1$Clusters_noCC], cycle_pch[CCA1$Cycle])

### HCC10 ###

# Base
set.seed(39203)
icas_hcc10_base <- do_ICA(HCC10, 4, "lognorm")
icas_hcc10_base$prog <- -3
ica_plot(icas_hcc10_base, get_group_cols(HCC10, "Clusters")[HCC10$Clusters], cycle_pch[HCC10$Cycle])

# Norm
set.seed(39203)
icas_hcc10_nCC<- do_ICA(HCC10, 2, "norm_exprs")
ica_plot(icas_hcc10_nCC, get_group_cols(HCC10, "Clusters_normCC")[HCC10$Clusters_normCC], cycle_pch[HCC10$Cycle])

# noCC
set.seed(39203)
icas_hcc10_noCC <- do_ICA(HCC10, 5, "noCC")
ica_plot(icas_hcc10_noCC, get_group_cols(HCC10, "Clusters_noCC")[HCC10$Clusters_noCC], cycle_pch[HCC10$Cycle])


png("CCA1_CCcorrect_Lineage_ICA.png", width=5*3, height=5, units="in", res=300)
par(mfrow=c(1,3))
par(mar=c(4,4,1,1))
ica_plot(icas_cca1_base, get_group_cols(CCA1, "Clusters")[CCA1$Clusters], cycle_pch[CCA1$Cycle])
ica_plot(icas_cca1_noCC, get_group_cols(CCA1, "Clusters_noCC")[CCA1$Clusters_noCC], cycle_pch[CCA1$Cycle])
ica_plot(icas_cca1_nCC, get_group_cols(CCA1, "Clusters_normCC")[CCA1$Clusters_normCC], cycle_pch[CCA1$Cycle])
dev.off()

png("HCC10_CCcorrect_Lineage_ICA.png", width=5*3, height=5, units="in", res=300)
par(mfrow=c(1,3))
par(mar=c(4,4,1,1))
ica_plot(icas_hcc10_base, get_group_cols(HCC10, "Clusters")[HCC10$Clusters], cycle_pch[HCC10$Cycle])
ica_plot(icas_hcc10_noCC, get_group_cols(HCC10, "Clusters_noCC")[HCC10$Clusters_noCC], cycle_pch[HCC10$Cycle])
ica_plot(icas_hcc10_nCC, get_group_cols(HCC10, "Clusters_normCC")[HCC10$Clusters_normCC], cycle_pch[HCC10$Cycle])
dev.off()

# Get Top genes each dimension:
top_genes <- function(ica, lineage_table) {
	weight_prog <- ica$ica$S[,abs(ica$prog)]*sign(ica$prog)
	weight_lin  <- ica$ica$S[,abs(ica$lin)]*sign(ica$lin)
	stemness <- weight_prog
	cholness <- -1*weight_lin
	hepness <- weight_lin
	# Add annotations
	stuff <- data.frame(Stem=stemness, Chol=cholness, Hep=hepness, Ann=lineage_table[match(names(stemness), lineage_table[,1]), 2])
	rownames(stuff) <- rownames(ica$ica$S)

	prog <- rownames(stuff)[stuff$Stem > quantile(stuff$Stem, probs=0.75) &
			#stuff$Chol < quantile(stuff$Chol, probs=0.5) &
			#stuff$Hep < quantile(stuff$Hep, probs=0.5) &
			grepl("Prog", stuff$Ann) ]
	hep <- rownames(stuff)[stuff$Hep > quantile(stuff$Hep, probs=0.75) & 
			stuff$Stem < quantile(stuff$Stem, probs=0.25) &
			stuff$Ann == "Hep-Mature" ]
	chol <- rownames(stuff)[stuff$Chol > quantile(stuff$Chol, probs=0.75) & 
			stuff$Stem < quantile(stuff$Stem, probs=0.25) &
			stuff$Ann == "Chol-Mature" ]
	return(list(weights=stuff, prog=prog, hep=hep, chol=chol))
}

#### HERE HERE HERE HERE HERE HERE HERE ####

cca1_genes_base <- top_genes(icas_cca1_base, Lineage)
cca1_genes_nCC <- top_genes(icas_cca1_nCC, Lineage)
cca1_genes_noCC <- top_genes(icas_cca1_noCC, Lineage)

hcc10_genes_base <- top_genes(icas_hcc10_base, Lineage)
hcc10_genes_nCC <- top_genes(icas_hcc10_nCC, Lineage)
hcc10_genes_noCC <- top_genes(icas_hcc10_noCC, Lineage)








## Novel Markers 
Lit_Markers <- read.delim("~/Collaborations/LiverOrganoids/Literature_Markers.txt", sep="\t", header=TRUE)
CC_Stem <- read.delim("~/Collaborations/LiverOrganoids/Markers_CC_Stemness.txt", sep=" ", header=TRUE)
HCC_Stem <- read.delim("~/Collaborations/LiverOrganoids/Markers_HCC_Stemness.txt", sep=" ", header=FALSE)
CHC_Stem <- read.delim("~/Collaborations/LiverOrganoids/Markers_CHC_Stemness.txt", sep=" ", header=TRUE)
Cancer_Stem <- c(as.character(CC_Stem[CC_Stem[,2]>0,1]), as.character(HCC_Stem[HCC_Stem[,2]>0,1]), as.character(CHC_Stem[CHC_Stem[,2]>0,1]))
Hepato_MSigdb <- read.delim("~/Collaborations/LiverOrganoids/Markers_HCC_Differentiation.txt")

HPA <- readRDS("/lustre/scratch117/cellgen/team218/TA/OtherDownloadedData/Good_HPA_Hep_Chol_Markers.rds")
HPA_clean <- cbind(as.character(HPA[,2]), as.character(HPA[,4]))
HPA_clean[as.numeric(HPA[,5]) < as.numeric(HPA[,13]),2] <- "hepatocyte"
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



Common_Stem_unsup <- data.frame(Common_Stem_unsup, 
	LineageMarker=as.character(Lineage[match(rownames(Common_Stem_unsup), Lineage[,1]),2]), 
	CC=rownames(Common_Stem_unsup) %in% All_CC, 
	SampPaper=as.character(Lit_Markers[ match(rownames(Common_Stem_unsup), Lit_Markers[,1]) ,2]),
	CancerStem=rownames(Common_Stem_unsup) %in% Cancer_Stem,
	HPA = as.character(HPA[ match(rownames(Common_Stem_unsup), HPA[,1]),2])
	)
	
Mature_HCC10_unsup <- data.frame(Mature_HCC10_unsup, 
	LineageMarker=as.character(Lineage[match(rownames(Mature_HCC10_unsup), Lineage[,1]),2]), 
	CC=rownames(Mature_HCC10_unsup) %in% All_CC, 
	SampPaper=as.character(Lit_Markers[ match(rownames(Mature_HCC10_unsup), Lit_Markers[,1]) ,2]),
	CancerStem=rownames(Mature_HCC10_unsup) %in% Cancer_Stem,
	HPA = as.character(HPA[ match(rownames(Mature_HCC10_unsup), HPA[,1]),2])
	)
	
Mature_CCA1_unsup <- data.frame(Mature_CCA1_unsup, 
	LineageMarker=as.character(Lineage[match(rownames(Mature_CCA1_unsup), Lineage[,1]),2]), 
	CC=rownames(Mature_CCA1_unsup) %in% All_CC, 
	SampPaper=as.character(Lit_Markers[ match(rownames(Mature_CCA1_unsup), Lit_Markers[,1]) ,2]),
	CancerStem=rownames(Mature_CCA1_unsup) %in% Cancer_Stem,
	HPA = as.character(HPA[ match(rownames(Mature_CCA1_unsup), HPA[,1]),2])
	)
	

write.table(Common_Stem_unsup, file="ICA_Markers_Stem.csv", sep=",", row.names=T, col.names=T)
write.table(Mature_HCC10_unsup, file="ICA_Markers_Hep.csv", sep=",", row.names=T, col.names=T)
write.table(Mature_CCA1_unsup, file="ICA_Markers_Chol.csv", sep=",", row.names=T, col.names=T)

## New Markers
#Super_Hep <- HCC10_Hep[rownames(HCC10_Hep) %in% rownames(CCA1_Hep),]
#Super_Chol <- HCC10_Chol[rownames(HCC10_Chol) %in% rownames(CCA1_Chol),]
#Super_Stem <-  HCC10_Prog[rownames(HCC10_Prog) %in% rownames(CCA1_Prog),]

intensity_plot <- function(SCE1, SCE2, genes, name1="", name2="", ylab1="", ylab2="") {
        vals1 <- exprs(SCE1)[rowData(SCE1)$Symbol %in% genes,];
	if (!is.null(dim(vals1))) {
		vals1 <- colMeans(vals1)
	}

        vals2 <- exprs(SCE2)[rowData(SCE2)$Symbol %in% genes,];
	if (!is.null(dim(vals2))) {
		vals2 <- colMeans(vals2)
	}
        splits <- quantile(c(vals1, vals2), probs=seq(from=0, to=1, length=8))
        splits <- unique(splits);
#        splits <- seq(from=0, to=max(splits), length=7)

        bins <- cut(vals1, breaks=splits, include.lowest=TRUE);
        my_cols <- colorRampPalette(brewer.pal(8, "Blues"))(length(levels(bins)));
        plot(SCE1$coords1, SCE1$coords2, col=my_cols[bins], pch=16, main=name1, xaxt="n", yaxt="n", ylab=ylab1, xlab="", bty="n")

        bins <- cut(vals2, breaks=splits, include.lowest=TRUE);
        my_cols <- colorRampPalette(brewer.pal(8, "Blues"))(length(levels(bins)));
        plot(SCE2$coords1, SCE2$coords2, col=my_cols[bins], pch=16, main=name2, xaxt="n", yaxt="n", ylab=ylab2, xlab="", bty="n")
	return(splits)
}

png("Cheating_Lineage_plot.png", width=4*3, height=4*2, units="in", res=300)
par(mfcol=c(2,3))
stem_splits <- intensity_plot(CCA1, HCC10, rownames(Common_Stem), name1="Stem", ylab1="CCA1", ylab2="HCC10")
hepa_splits <- intensity_plot(CCA1, HCC10, rownames(Mature_Hep_HCC10), name1="Hepatocyte")
chol_splits <- intensity_plot(CCA1, HCC10, rownames(Mature_Chol_CCA1), name1="Cholangiocyte")
dev.off()

my_legend <- function(splits) {
	my_cols <- colorRampPalette(brewer.pal(8, "Blues"))(length(splits)-1);
	a<-diff(splits)/2 + splits[-length(splits)]
	legend("bottomleft", as.character(c(round(a, digits=1))), fill=my_cols, bty="n", ncol=3)
}


### Markers for Laura
marker_plot_for_laura <- function(genes, name) {
#	png(name, width=4*length(genes), height=2*4, units="in", res=300)
	par(mfcol=c(2,length(genes)))
	par(mar=c(0,2,2,0))
	for (g in genes) {
		splits <- intensity_plot(CCA1, HCC10, g, name1=g, ylab1="CCA1", ylab2="HCC10")
		my_legend(splits)
	}
#	dev.off()
}
list_HCC10_CSC <-  c("EPCAM", "EZR", "PSMD14", "KRT18")
list_HCC10_HEP <- c("CYP2E1", "AFM", "ALDH6A1", "G6PC", "ADH6", "SLC51A")
list_CCA1_CSC <- c("C1QBP", "EXOSC9", "ODC1", "TUBB2B")
list_CCA1_Chol <- c("NDRG1", "MSLN", "CD164", "CA9")

marker_plot_for_laura(list_HCC10_CSC, "ForLaura_ExprPlots_HCC10-CSC.png")
marker_plot_for_laura(list_HCC10_HEP, "ForLaura_ExprPlots_HCC10-HEP.png")
marker_plot_for_laura(list_CCA1_CSC, "ForLaura_ExprPlots_CCA1-CSC.png")
marker_plot_for_laura(list_CCA1_Chol, "ForLaura_ExprPlots_CCA1-Chol.png")

png("ForLaura_ExprPlots_HCC10-CSC.png", width=4*4, height=4*2, units="in", res=300)
par(mfcol=c(2,4))
par(mar=c(0,2,2,0))
splits <- intensity_plot(CCA1, HCC10, "EPCAM", name1="Epcam", ylab1="CCA1", ylab2="HCC10")
my_legend(splits)
splits <- intensity_plot(CCA1, HCC10, "EZR", name1="EZR", ylab1="", ylab2="")
my_legend(splits)
splits <- intensity_plot(CCA1, HCC10, "PSMD14", name1="PSMD14", ylab1="", ylab2="")
my_legend(splits)
splits <- intensity_plot(CCA1, HCC10, "KRT18", name1="KRT18", ylab1="", ylab2="")
my_legend(splits)
dev.off()

png("ForLaura_ExprPlots_HCC10-HEP.png", width=4*6, height=4*2, units="in", res=300)
par(mfcol=c(2,6))
par(mar=c(0,2,2,0))
splits <- intensity_plot(CCA1, HCC10, "CYP2E1", name1="CYP2E1", ylab1="CCA1", ylab2="HCC10")
my_legend(splits)
splits <- intensity_plot(CCA1, HCC10, "AFM", name1="AFM", ylab1="", ylab2="")
my_legend(splits)
splits <- intensity_plot(CCA1, HCC10, "ALDH6A1", name1="ALDH6A1", ylab1="", ylab2="")
my_legend(splits)
splits <- intensity_plot(CCA1, HCC10, "G6PC", name1="G6PC", ylab1="", ylab2="")
my_legend(splits)
splits <- intensity_plot(CCA1, HCC10, "ADH6", name1="ADH6", ylab1="", ylab2="")
my_legend(splits)
splits <- intensity_plot(CCA1, HCC10, "SLC51A", name1="SLC51A", ylab1="", ylab2="")
my_legend(splits)
dev.off()



### Novel Markers GO enrichments ###
require("gProfileR")
hep_GO <- gprofiler(rownames(Mature_HCC10_unsup))
chol_GO <- gprofiler(rownames(Mature_CCA1_unsup))
stem_GO <- gprofiler(rownames(Common_Stem_unsup))
### 

### Other Lines ###
plot_line <- function(file) {

	require("fastICA")
	SCE <- readRDS(file)
	marker_genes <- as.character(rowData(SCE)[rowData(SCE)$fine_marker_is.Feature, "Symbol"])
	mat <- assays(SCE)[[ expr_type ]]; mat <- mat[rowData(SCE)$Symbol %in% marker_genes,]
	clust_cols <- get_group_cols(SCE)[SCE$clusters_clean]
	CC_cols <- CC_col[SCE$CC_state_new]

	get_gene_colours <- function(SCE, genes, splits) {
	        vals <- exprs(SCE)[rowData(SCE)$Symbol %in% genes,];
		if (!is.null(dim(vals))) {
			vals <- colMeans(vals)
		}
		splits[1] <- 0; splits[length(splits)] <- max(vals, splits);
	        bins <- cut(vals, breaks=splits, include.lowest=TRUE);
	        my_cols <- colorRampPalette(brewer.pal(8, "Blues"))(length(levels(bins)));	
		return(my_cols[bins])
	}

	# Fancy Plotting:
	require("igraph")
	require("dbscan")
	set.seed(12530)
	icas <- fastICA(mat, n.comp = 50, method="C")$A
	this_sNN <- sNN(t(icas), 5)
	adj_mat <- matrix(0, nrow=ncol(icas), ncol=ncol(icas))
	for (i in 1:ncol(icas)) {
	        adj_mat[i, this_sNN$id[i,]] <- this_sNN$shared[i,]+1
	}
	colnames(adj_mat) = rownames(adj_mat) = colnames(SCE)
	graph <- graph_from_adjacency_matrix(adj_mat, mode="directed", weighted=TRUE)

	mat_to_graph <- match(V(graph)$name, colnames(SCE))
	v_colour_cc <- CC_cols[mat_to_graph]
	v_colour_g <- clust_cols[mat_to_graph]

	set.seed(103)
	par(mar=c(0,0,0,0))
	L <- layout_with_fr(graph)
	plot(graph, layout=L, vertex.color=get_gene_colours(SCE, rownames(Common_Stem), stem_splits)[mat_to_graph], 
		vertex.size=3, vertex.label=NA, edge.arrow.size=0, edge.arrow.width=0, edge.lty=3, edge.color="white", 
		vertex.frame.color=NA)
	plot(graph, layout=L, vertex.color=get_gene_colours(SCE, rownames(Mature_Hep_HCC10), hepa_splits)[mat_to_graph], 
		vertex.size=3, vertex.label=NA, edge.arrow.size=0, edge.arrow.width=0, edge.lty=3, edge.color="white",
		vertex.frame.color=NA)
	plot(graph, layout=L, vertex.color=get_gene_colours(SCE, rownames(Mature_Chol_CCA1), chol_splits)[mat_to_graph], 
		vertex.size=3, vertex.label=NA, edge.arrow.size=0, edge.arrow.width=0, edge.lty=3, edge.color="white",
		vertex.frame.color=NA)
	return(list(layout=L, mat_to_graph = mat_to_graph))
}

png("Cheating_Lineage_plot_others.png", width=4*3, height=4*3, units="in", res=300)
par(mfrow=c(3,3))
hcc6_graph <- plot_line("HCC6_noCC_SC3.rds")
d3dm_graph <- plot_line("D3DM_noCC_SC3.rds")
d3em_graph <- plot_line("D3EM_noCC_SC3.rds")
dev.off()

q("n")

# ICAs of other lines:

make_stuff <- function(file) {
	require("fastICA")
        SCE <- readRDS(file)
        marker_genes <- as.character(rowData(SCE)[rowData(SCE)$fine_marker_is.Feature, "Symbol"])
	marker_genes <- marker_genes[marker_genes != ""]
	SCE <- SCE[rowData(SCE)$Symbol %in% marker_genes,]
        mat <- assays(SCE)[[expr_type]];

        clust_cols <- get_group_cols(SCE)[SCE$clusters_clean]
        CC_cols <- CC_col[SCE$CC_state_new]
	return(list(M=mat, CCcol=CC_cols, Clustcol=clust_cols, Symbol=rowData(SCE)$Symbol))
}

get_icas <- function(mat, n=2) {
        # Fancy Plotting:
        require("igraph")
        require("dbscan")
        set.seed(12530)
        icas <- fastICA(mat, n.comp = n, method="C")
	return(icas)
}

my_colmeans <- function(mat, rows) {
	if (sum(rows) > 1) {
		return(colMeans(mat[rows,]))
	} else if (sum(rows) == 1) {
		return(mat[rows,])
	} else {
		return(rep(0, times=ncol(mat)))
	}
}

make_table <- function(icas, symbol) {
	a <- symbol %in% Chol_lineage[ Chol_lineage[,2] == "Chol-Prog" ,1]
	b <- symbol %in% Chol_lineage[ Chol_lineage[,2] == "Chol-Mature" ,1]
	c <- symbol %in% Hep_lineage[ Hep_lineage[,2] == "Hep-Prog" ,1]
	d <- symbol %in% Hep_lineage[ Hep_lineage[,2] == "Hep-Mature" ,1]

	a <- my_colmeans(icas$S, a)
	b <- my_colmeans(icas$S, b)
	c <- my_colmeans(icas$S, c)
	d <- my_colmeans(icas$S, d)

	tab <- rbind(a,b,c,d)
	rownames(tab) <- c("C-Prog", "C", "H-Prog", "H")
	return(tab)
}


hcc6_stuff <- make_stuff("HCC6_noCC_SC3.rds")
hcc6_icas <- get_icas(hcc6_stuff$M, 3)
tab <- make_table(hcc6_icas, hcc6_stuff$Symbol)


dm_stuff <- make_stuff("D3DM_noCC_SC3.rds")
dm_icas <- get_icas(dm_stuff$M, 4)
tab <- make_table(dm_icas, dm_stuff$Symbol)


em_stuff <- make_stuff("D3EM_noCC_SC3.rds")
em_icas <- get_icas(em_stuff$M, 4)
tab <- make_table(em_icas, em_stuff$Symbol)






#### Other Stuff ####

# merged
tab <- rbind(
colMeans(icas_merged$S[gene_symbols %in% Chol_lineage[ Chol_lineage[,2] == "Prog" ,1],]),
colMeans(icas_merged$S[gene_symbols %in% Chol_lineage[ Chol_lineage[,2] == "Chol" ,1],]),
colMeans(icas_merged$S[gene_symbols %in% Hep_lineage[ Hep_lineage[,2] == "Prog" ,1],]),
colMeans(icas_merged$S[gene_symbols %in% Hep_lineage[ Hep_lineage[,2] == "Hep" ,1],]))
rownames(tab) <- c("C-Prog", "C", "H-Prog", "H")

plot(icas_merged$A[2,], icas_merged$A[1,], col=c(cca1_cols, hcc10_cols), pch=16, xlab="Chol -------- Hep", ylab="Diff ------- Stem")
plot(icas_merged$A[2,], icas_merged$A[1,], col=c(cca1_CC, hcc10_CC), pch=16, xlab="Chol -------- Hep", ylab="Diff ------- Stem")
legend("bottomright", levels(HCC10$CC_state_new), col=CC_col, pch=16, bty="n")
plot(icas_merged$A[2,], icas_merged$A[1,], col=rep(type_col, times=c(ncol(cca1_mat), ncol(hcc10_mat))), pch=16, xlab="Chol -------- Hep", ylab="Diff ------- Stem")
legend("bottomright", c("CCA1", "HCC10"), col=type_col, pch=16, bty="n")

thing <- t(icas_merged$A) %*% t(tab)
thing_scaled <- apply(thing, 2, function(x){x<-x-min(x);x <- x/max(x)})



# Fancy Plotting:
require("igraph")
require("dbscan")
set.seed(20391)
icas <- fastICA(cca1_mat, n.comp = 50, method="C")$A
cca1_sNN <- sNN(t(icas), 5)
#cca1_sNN_smoothed <- sapply(1:ncol(icas), function(i) {colSums(t(icas[,cca1_sNN$id[i,]])*cca1_sNN$shared[i,])/(max(sum(cca1_sNN$shared[i,]), 1))})
adj_mat <- matrix(0, nrow=ncol(icas), ncol=ncol(icas))
for (i in 1:ncol(icas)) {
	adj_mat[i, cca1_sNN$id[i,]] <- cca1_sNN$shared[i,]+1
}
colnames(adj_mat) = rownames(adj_mat) = colnames(CCA1)
graph <- graph_from_adjacency_matrix(adj_mat, mode="directed", weighted=TRUE)

v_colour_cc <- cca1_CC[match(V(graph)$name, colnames(CCA1))]
v_colour_g <- cca1_cols[match(V(graph)$name, colnames(CCA1))]
plot(graph, layout=layout_with_fr(graph), vertex.color=v_colour_g, vertex.size=3, vertex.label=NA, edge.arrow.size=0, edge.arrow.width=0, edge.lty=3, edge.color="white")

require("Rtsne")
set.seed(20391)
tsne <- Rtsne(t(assays(CCA1)[[ expr_type]][rowData(CCA1)$noCC_fine_marker_is.Feature,]), perplexity=20, initial_dims=10)
plot(tsne$Y[,1], tsne$Y[,2], col=cca1_cols, pch=16)
