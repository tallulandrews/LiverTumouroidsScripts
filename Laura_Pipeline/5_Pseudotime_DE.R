# R3.4
args = commandArgs(trailingOnly=TRUE)

source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")

set.seed(371938)

SCE <- readRDS(args[1])
outprefix=args[2]
expr_type <- "lognorm"

require("scater")
require("SingleCellExperiment")
require("matrixStats")

#SCE <- readRDS(args[1])
cluster_col_set <- get_group_cols(SCE)
cluster_names_for_legend <- names(cluster_col_set); cluster_names_for_legend[cluster_names_for_legend == "Outliers"] <- "O";

coef_lim = 2
fdr_lim = 0.05

### Get pseudotime-DE
# MAST

do_MAST <- function(SCE, pseudo_col, branch_col) {
	require("MAST")
	require("data.table")
	set.seed(101);

	mat <- assays(SCE)[["lognorm"]]
	fdata <- rowData(SCE)
	cdata <- colData(SCE)
	scaRAW <- FromMatrix(mat, cdata, fdata)
	scaRAW <- scaRAW[rowSums(mat > 0) > 0.05*ncol(mat),]
	colData(scaRAW)$geneson <- colSums(assay(scaRAW) > 0)
	colData(scaRAW)$Time <- as.numeric(colData(scaRAW)[,pseudo_col])
	colData(scaRAW)$Time[is.na(colData(scaRAW)$Time)] <- mean(colData(scaRAW)$Time, na.rm=T); #fix possible NAs
	colData(scaRAW)$Branch <- as.numeric(colData(scaRAW)[,branch_col])
	colData(scaRAW)$Branch[is.na(colData(scaRAW)$Branch)] <- max(as.numeric(colData(scaRAW)$Branch), na.rm=T)+1 # Fix possible NAs
	colData(scaRAW)$Branch <- factor(colData(scaRAW)$Branch)

	zlmOut <- zlm(~Time+Branch+geneson, scaRAW)
	contrast_out <- summary(zlmOut, doLRT='Time')

	summaryDt <- contrast_out$datatable
	outTABLE <- as.data.frame(summaryDt[summaryDt$contrast=="Time" & summaryDt$component=="H",])
	outTABLE2 <- as.data.frame(summaryDt[summaryDt$contrast=="Time" & summaryDt$component=="logFC",])
	outTABLE[,5:8] <- outTABLE2[5:8]
	outTABLE$fdr <- p.adjust(outTABLE[,"Pr(>Chisq)"], method="fdr")
	outTABLE$Symbol <- rowData(scaRAW)$Symbol

	return(outTABLE);
}

do_MAST_with_clusters <- function(SCE, pseudo_col, branch_col, cluster_col) {
        require("MAST")
        require("data.table")
        set.seed(101);
	SCE <- SCE[,SCE$clusters_clean != "Outliers"]

        mat <- assays(SCE)[["lognorm"]]
        fdata <- rowData(SCE)
        cdata <- colData(SCE)
        scaRAW <- FromMatrix(mat, cdata, fdata)
        scaRAW <- scaRAW[rowSums(mat > 0) > 0.05*ncol(mat),]
        colData(scaRAW)$geneson <- colSums(assay(scaRAW) > 0)
        colData(scaRAW)$Time <- as.numeric(colData(scaRAW)[,pseudo_col])
        colData(scaRAW)$Time[is.na(colData(scaRAW)$Time)] <- mean(colData(scaRAW)$Time, na.rm=T); #fix possible NAs
        colData(scaRAW)$Branch <- as.numeric(colData(scaRAW)[,branch_col])
        colData(scaRAW)$Branch[is.na(colData(scaRAW)$Branch)] <- max(as.numeric(colData(scaRAW)$Branch), na.rm=T)+1 # Fix possible NAs
        colData(scaRAW)$Branch <- factor(colData(scaRAW)$Branch)
        colData(scaRAW)$Cluster <- as.character(colData(scaRAW)[,cluster_col])
        colData(scaRAW)$Cluster <- factor(colData(scaRAW)$Cluster, levels=levels(colData(scaRAW)[,cluster_col]))
	colData(scaRAW)$Cluster <- factor(colData(scaRAW)$Cluster);

	outTABLES <- list()

        #zlmOut <- zlm(~Time+Branch+Cluster+geneson, scaRAW)
        zlmOut <- zlm(~Time+Cluster+geneson, scaRAW)

	# time
        contrast_out <- summary(zlmOut, doLRT='Time')
        summaryDt <- contrast_out$datatable
        outTABLE <- as.data.frame(summaryDt[summaryDt$contrast=="Time" & summaryDt$component=="H",])
        outTABLE2 <- as.data.frame(summaryDt[summaryDt$contrast=="Time" & summaryDt$component=="logFC",])
        outTABLE[,5:8] <- outTABLE2[5:8]
        outTABLE$fdr <- p.adjust(outTABLE[,"Pr(>Chisq)"], method="fdr")
	outTABLE$Symbol <- rowData(scaRAW)$Symbol

	outTABLES[["Time"]] <- outTABLE;

	# clusters
	possible_clusters <- levels(colData(scaRAW)$Cluster)
	for (C in possible_clusters[2:length(possible_clusters)]) {
		# probably get an error for reference cluster ?
		this_name <- paste("Cluster", C, sep="");
		
		print(this_name)
		out <- summary(zlmOut, doLRT=this_name)
		
		summaryDt <- out$datatable
	        outTABLE <- as.data.frame(summaryDt[summaryDt$contrast==this_name & summaryDt$component=="H",])
	        outTABLE2 <- as.data.frame(summaryDt[summaryDt$contrast==this_name & summaryDt$component=="logFC",])
	        outTABLE[,5:8] <- outTABLE2[5:8]
	        outTABLE$fdr <- p.adjust(outTABLE[,"Pr(>Chisq)"], method="fdr")
		outTABLE$Symbol <- rowData(scaRAW)$Symbol

		outTABLES[[C]] <- outTABLE
	}

        return(outTABLES);
}

Monocle_Mega_DE <- do_MAST_with_clusters(SCE, which(colnames(colData(SCE)) == "Monocle_time_all"), which(colnames(colData(SCE)) == "Monocle_branch_all"), which(colnames(colData(SCE))=="clusters_clean"))

Mega_table <- data.frame(Gene=Monocle_Mega_DE[[1]][,1], Symbol=Monocle_Mega_DE[[1]]$Symbol)
for (t in names(Monocle_Mega_DE)) {
	this_table <- Monocle_Mega_DE[[t]]
	t_prefix <- as.character(this_table$contrast[1])
	this_table <- this_table[match(Mega_table$Gene ,this_table$primerid),]
	colnames(this_table) <- paste(t_prefix, colnames(this_table), sep="_")

	Mega_table <- cbind(Mega_table, this_table[,c(7,9)]) 
}

Mega_table$is.DE_Time <-Mega_table[,4] < 0.05
clust_cols <- seq(from=6, to=ncol(Mega_table), by=2)
Mega_table$is.DE_Cluster <- rowSums(Mega_table[,clust_cols] < 0.05/length(clust_cols)) > 0

print(paste("Monocle time DE:",sum(Mega_table$is.DE_Time)))
print(paste("Cluster DE:",sum(Mega_table$is.DE_Cluster)))


write.table(Mega_table, file=paste(outprefix, "Monocle_MegaDE.txt", sep="_"))

DPT_Mega_DE <- do_MAST_with_clusters(SCE, which(colnames(colData(SCE)) == "DPT_time"), which(colnames(colData(SCE)) == "DPT_branch"), which(colnames(colData(SCE))=="clusters_clean"))

Mega_table <- data.frame(Gene=DPT_Mega_DE[[1]][,1], Symbol=DPT_Mega_DE[[1]]$Symbol)
for (t in names(DPT_Mega_DE)) {
        this_table <- DPT_Mega_DE[[t]]
        t_prefix <- as.character(this_table$contrast[1])
        this_table <- this_table[match(Mega_table$Gene ,this_table$primerid),]
        colnames(this_table) <- paste(t_prefix, colnames(this_table), sep="_")

        Mega_table <- cbind(Mega_table, this_table[,c(7,9)]) 
}
Mega_table$is.DE_Time <-Mega_table[,4] < 0.05
clust_cols <- seq(from=6, to=ncol(Mega_table), by=2)
Mega_table$is.DE_Cluster <- rowSums(Mega_table[,clust_cols] < 0.05/length(clust_cols)) > 0

print(paste("Dpt time DE:",sum(Mega_table$is.DE_Time)))
print(paste("Cluster DE:",sum(Mega_table$is.DE_Cluster)))

write.table(Mega_table, file=paste(outprefix, "DPT_MegaDE.txt", sep="_"))


#monocle_DE <- do_MAST(SCE, which(colnames(colData(SCE)) == "Monocle_time_all"), which(colnames(colData(SCE)) == "Monocle_branch_all"))
monocle_DE <- Monocle_Mega_DE[[1]]
toadd <- monocle_DE[match(rownames(SCE),monocle_DE$primerid),]
rowData(SCE)$Monocle_time_DE_coef <- toadd[,"coef"]
rowData(SCE)$Monocle_time_DE_q.value <- toadd[,"fdr"]


#DPT_DE <- do_MAST(SCE, which(colnames(colData(SCE)) == "DPT_time"), which(colnames(colData(SCE)) == "DPT_branch"))
DPT_DE <- DPT_Mega_DE[[1]]
toadd <- DPT_DE[match(rownames(SCE),DPT_DE$primerid),]
rowData(SCE)$DPT_time_DE_coef <- toadd[,"coef"]
rowData(SCE)$DPT_time_DE_q.value <- toadd[,"fdr"]

saveRDS(SCE, file=paste(outprefix, "PseudoDE.rds", sep="_"))

OUT_mono <- monocle_DE[abs(monocle_DE[,"coef"]) > coef_lim & monocle_DE[,"fdr"] < fdr_lim & sign(monocle_DE[,"ci.hi"]) == sign(monocle_DE[,"ci.lo"]),]
reorder <- apply(OUT_mono[ , c("ci.hi", "ci.lo")], 1, function(x) {min(abs(x)) * sign(min(x))})
OUT_mono <- OUT_mono[order(-reorder),]
#write.table(OUT_mono, file=paste(outprefix, "Monocle_PseudoDE.png", sep="_"))


OUT_dpt <- DPT_DE[abs(DPT_DE[,"coef"]) > coef_lim & DPT_DE[,"fdr"] < fdr_lim & sign(DPT_DE[,"ci.hi"]) == sign(DPT_DE[,"ci.lo"]),]
reorder <- apply(OUT_dpt[ , c("ci.hi", "ci.lo")], 1, function(x) {min(abs(x)) * sign(min(x))})
OUT_dpt <- OUT_dpt[order(-reorder),]
#write.table(OUT_dpt, file=paste(outprefix, "DPT_PseudoDE.png", sep="_"))


# pseudotime-heatmap
pseudo_time_heatmap <- function(SCE, pseudotime, genes) {
        require("CellTypeProfiles")
        require("RColorBrewer")
	SCE <- SCE[,!is.na(pseudotime)]
	pseudotime <- pseudotime[!is.na(pseudotime)]

	pseudo_order <- order(pseudotime);
	type_order <- sapply(split(seq(length(pseudotime)), SCE$clusters_clean), function(a){median(pseudotime[a])})
	type_order <- type_order[as.numeric(SCE$clusters_clean)]

	cell_order <- order(type_order, pseudotime)

        heat_data <- assays(SCE)[[expr_type]]
	heat_data <- heat_data[,cell_order]
	gene_names <- rowData(SCE)$Symbol
	gene_names <- gene_names[genes]
	heat_data <- heat_data[genes,];
        columnLab <- SCE$clusters_clean[cell_order]
	pseudotime <- pseudotime[cell_order]
        set.seed(666)

	# Order genes #

	gene_2_pseudo <- function(gene) {
		gene <- gene/sum(gene) # rescale so total == 1
		pseudo <- sum(gene * rank(pseudotime));
		return(pseudo);
	}

	#gene_peaks <- apply(heat_data, 1, gene_2_pseudo);

	#gene_order <- order(gene_peaks);

        #heat_data <- heat_data[gene_order,]

        #heat_data <- heat_data/apply(heat_data,1,max)
        heat_data <- t(apply(heat_data,1,function(a){(a-min(a))/max(a)}))
        heat_colours <- c("#f7fcfd","#e0ecf4","#bfd3e6","#9ebcda","#8c96c6","#8c6bb1","#88419d","#810f7c","#4d004b")
	
	pseudo_colours <- colorRampPalette(c("white", "black"))(10)
	Pseudotime <- pseudo_colours[cut(pseudotime, breaks=10)]


        Cluster <-cluster_col_set[columnLab]
        plout <- heatmap.3(heat_data, trace="none", col=heat_colours, Colv=FALSE, 
			Rowv=TRUE, dendrogram="row",
                        ColSideColors=cbind(Cluster, Pseudotime),ColSideColorsSize=2, labRow=gene_names,
                        labCol="")
        return(plout)
}

if (sum(rownames(SCE) %in% OUT_mono$primerid) > 3) {
	png(paste(outprefix, "Monocle_PseudotimeHeatmap.png", sep="_"), width=8, height=8, units="in", res=300)
	junk <- pseudo_time_heatmap(SCE, SCE$Monocle_time_all, rownames(SCE) %in% OUT_mono$primerid)
	dev.off()
}

if (sum(rownames(SCE) %in% OUT_dpt$primerid) > 3) {
	png(paste(outprefix, "DPT_PseudotimeHeatmap.png", sep="_"), width=8, height=8, units="in", res=300)
	junk<-pseudo_time_heatmap(SCE, SCE$DPT_time, rownames(SCE) %in% OUT_dpt$primerid)
	dev.off()
}

if (sum(rowData(SCE)$fine_marker_is.Feature) > 3) {
	png(paste(outprefix, "Marker_PseudotimeHeatmap.png", sep="_"), width=8, height=8, units="in", res=300)
	junk<-pseudo_time_heatmap(SCE, SCE$DPT_time, rowData(SCE)$fine_marker_is.Feature)
	dev.off()
}
