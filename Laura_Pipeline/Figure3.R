# Lineage/Pseudotime

require("slingshot")
require("fastICA")
require("M3Drop")
require("CellTypeProfiles")
require("SingleCellExperiment")
require("scater")
require("destiny")
source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")



# New Plan : tSNE coloured by lineage scores from global. - too large scale.


# Read in each manual clustering
plotting_rds <- c("CCA1_PlottingObj.rds", "CCA5_PlottingObj.rds", "HCC6_PlottingObj.rds", 
		"HCC23_PlottingObj.rds", "HCC10_PlottingObj.rds", "HCC24_PlottingObj.rds", 
		"D3DM_PlottingObj.rds", "D3EM_PlottingObj.rds", "D9DM_PlottingObj.rds", 
		"D9EM_PlottingObj.rds");
rds_names <- c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24", "D3DM", "D3EM", "D9DM", "D9EM")
expr_type <- "norm_exprs";

Global <- readRDS("Global_SCE.rds")

#### Global DE genes - all the things - most of these genes are garbage?
#require("MAST")
#mat <- assays(Global)[["counts"]]
#fdata <- rowData(Global)
#cdata <- colData(Global)
#scaRAW <- FromMatrix(mat, cdata, fdata)
#scaRAW <- scaRAW[rowSums(mat > 0) > 0.05*ncol(mat),]
#colData(scaRAW)$geneson <- colSums(assay(scaRAW) > 0)
#colData(scaRAW)$totcount <- colSums(assay(scaRAW))


#zlmOut <- zlm(~Hep_score + Chol_score + Stem_score + Proliferating + Donor + geneson + totcount, scaRAW)
#
#Hep_de <- summary(zlmOut, doLRT='Hep_score')$datatable
#Hep_de_p.values <- as.data.frame(Hep_de[Hep_de$contrast=="Hep_score" & Hep_de$component=="H",])
#Hep_de_logFC <- as.data.frame(Hep_de[Hep_de$contrast=="Hep_score" & Hep_de$component=="logFC",])
#
#Chol_de <- summary(zlmOut, doLRT='Chol_score')$datatable
#Chol_de_p.values <- as.data.frame(Chol_de[Chol_de$contrast=="Chol_score" & Chol_de$component=="H",])
#Chol_de_logFC <- as.data.frame(Chol_de[Chol_de$contrast=="Chol_score" & Chol_de$component=="logFC",])

#Stem_de <- summary(zlmOut, doLRT='Stem_score')$datatable
#Stem_de_p.values <- as.data.frame(Stem_de[Stem_de$contrast=="Stem_score" & Stem_de$component=="H",])
#Stem_de_logFC <- as.data.frame(Stem_de[Stem_de$contrast=="Stem_score" & Stem_de$component=="logFC",])

#Prolif_de <- summary(zlmOut, doLRT='ProliferatingTRUE')$datatable
#Prolif_de_p.values <- as.data.frame(Prolif_de[Prolif_de$contrast=="ProliferatingTRUE" & Prolif_de$component=="H",])
#Prolif_de_logFC <- as.data.frame(Prolif_de[Prolif_de$contrast=="ProliferatingTRUE" & Prolif_de$component=="logFC",])

#Global <- Global[rownames(Global) %in% rownames(scaRAW),]

#rowData(Global)$HepDE.q.value <- p.adjust(Hep_de_p.values[,4], method="fdr")
#rowData(Global)$HepDE.logFC <- Hep_de_logFC$coef
#rowData(Global)$CholDE.q.value <- p.adjust(Chol_de_p.values[,4], method="fdr")
#rowData(Global)$CholDE.logFC <- Chol_de_logFC$coef
#rowData(Global)$StemDE.q.value <- p.adjust(Stem_de_p.values[,4], method="fdr")
#rowData(Global)$StemDE.logFC <- Stem_de_logFC$coef
#rowData(Global)$Prolif.q.value <- p.adjust(Prolif_de_p.values[,4], method="fdr")
#rowData(Global)$Prolif.logFC <- Prolif_de_logFC$coef

#saveRDS(Global, file="Global_DE.rds")

GO <- read.delim("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/GO/hsapiens_80_GO_Annoations_Emsembl.out", sep="\t", header=F)
GO_cc <- GO[GO[,3] == "cell cycle",1]

#within_stemy <- (rowData(Global)[(
#	rowData(Global)$Prolif.q.value < 0.05 & 
#	rowData(Global)$Prolif.logFC > 60 & 
#	!rownames(Global) %in% GO_cc ),"Symbol"])

#require("gProfileR")
#richments <- gprofiler(within_stemy, custom_bg=rowData(Global)$Symbol)

### DM + slingshot + MAST
Lineage <- read.table("Cleaned_Lineage.txt", header=F)


my_row_mean_aggregate <- function (mat, groups) {
    MAT <- as.matrix(mat)
    x <- split(seq(ncol(MAT)), groups)
    result <- sapply(x, function(a) {
		if (length(a) > 1) {rowMeans(MAT[, a])} 
		else {	MAT[,a] } })
    return(result)
}


do_slingshot_DM <- function(SCE) {
	require("CellTypeProfiles")
	# remove small clusters
	exclude <- table(SCE$Manual_Clusters)
	SCE <- SCE[,! SCE$Manual_Clusters %in% names(exclude)[exclude < 5]]

	# extract expre mat
	dm_mat <- assays(SCE)[["norm_exprs"]]
	g_keep <- rowData(SCE)$Marker.q.value < 0.05 & rowData(SCE)$Marker.q.value > -1
	dm_mat <- dm_mat[g_keep,]

	# set root & leaves
	scores <- my_row_mean_aggregate(dm_mat, SCE$Manual_Clusters)
	chol <- rowData(SCE)$Symbol[g_keep] %in% Lineage[Lineage[,2] == "Chol-Mature",1] 
	hep <- rowData(SCE)$Symbol[g_keep] %in% Lineage[Lineage[,2] == "Hep-Mature",1] 
	stem <- rowData(SCE)$Symbol[g_keep] %in% Lineage[grep("Prog", Lineage[,2]) ,1] 
	chol <- colMeans(scores[chol,])
	hep <- colMeans(scores[hep,])
	stem <- colMeans(scores[stem,])
	root <- names(stem[stem==max(stem)])
	leaf_hep <- (hep-mean(hep))/mean(hep) - (stem-mean(stem))/mean(stem);
	leaf_chol <- (chol-mean(chol))/mean(chol) - (stem-mean(stem))/mean(stem);
	if (max(leaf_chol) > max(leaf_hep)) {
		leaf <- names(leaf_chol[leaf_chol==max(leaf_chol)])
	} else {
		leaf <- names(leaf_hep[leaf_hep==max(leaf_hep)])
	}


	# Dim Reduce with DM
	DM <- DiffusionMap(t(dm_mat))
	DM_dims <- eigenvectors(DM)
	reducedDims(SCE)[["DM"]]=DM_dims[,1:2]
	
	dm_scores <- colMeans(DM_dims[SCE$Manual_Clusters==root,])- colMeans(DM_dims[SCE$Manual_Clusters==leaf,])
	time_dim <- which(abs(dm_scores) == max(abs(dm_scores)))
	if(dm_scores[time_dim] > 0) {
		time <- -DM_dims[,time_dim]
	} else {
		 time <- DM_dims[,time_dim]
	}
	require("Hmisc")
	gene_scores <- rcorr(t(dm_mat), time, type="spearman")
	gene_out <- data.frame(r=gene_scores$r[,ncol(gene_scores$r)], p=gene_scores$P[,ncol(gene_scores$P)])
	gene_out <- gene_out[-nrow(gene_out),]
	gene_out$q <- p.adjust(gene_out$p, method="fdr")
	gene_out$CC <- rownames(gene_out) %in% GO_cc
	gene_out$Symbol <- rowData(SCE)$Symbol[g_keep]
	gene_out <- cbind(gene_out, scores)

	# Slingshot - fails utterly on HCC10
	#SCE <- slingshot(SCE, clusterLabels='Manual_Clusters', reducedDim='DM', start.clus = root, end.clus= leaf)
        #fitted <- SlingshotDataSet(SCE)

	colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
	plot(reducedDims(SCE)$DM, col = colors[cut(time,breaks=100)], pch=16, asp = 1)
	#lines(SlingshotDataSet(SCE), lwd=3, col="black") #type="lineages"

	#plot(reducedDims(SCE)$DM, col = SCE@metadata$palette[SCE$Manual_Clusters], pch=16, asp = 1)
	#lines(SlingshotDataSet(SCE), lwd=3, col="black", type="lineages")

	return(gene_out)
}

HEP_genes <- do_slingshot_DM("HCC10_PlottingObj.rds")
CHOL_genes <- do_slingshot_DM("CCA1_PlottingObj.rds")

a <- HEP_genes$Symbol[HEP_genes$r < -0.5]
b <- CHOL_genes$Symbol[CHOL_genes$r < -0.5 & !CHOL_genes$CC]
stem <- a[a %in% b]

a <- HEP_genes$Symbol[HEP_genes$r > 0.5]
b <- CHOL_genes$Symbol[CHOL_genes$r > 0.5]
chol <- b[ ! b %in% a]
hep <-  a[ ! a %in% b]



do_MAST <- function(SCE, pseudo_col, class_col, CC_col) {
        require("MAST")
        require("data.table")
        set.seed(101);

        mat <- assays(SCE)[["counts"]]
        fdata <- rowData(SCE)
        cdata <- colData(SCE)
        scaRAW <- FromMatrix(mat, cdata, fdata)
        scaRAW <- scaRAW[rowSums(mat > 0) > 0.05*ncol(mat),]
        colData(scaRAW)$geneson <- colSums(assay(scaRAW) > 0)
        colData(scaRAW)$Time <- as.numeric(colData(scaRAW)[,pseudo_col])
        colData(scaRAW)$Time[is.na(colData(scaRAW)$Time)] <- mean(colData(scaRAW)$Time, na.rm=T); #fix possible NAs
        colData(scaRAW)$Class <- as.numeric(colData(scaRAW)[,class_col])
        colData(scaRAW)$Class[is.na(colData(scaRAW)$Class)] <- max(as.numeric(colData(scaRAW)$Class), na.rm=T)+1 # Fix possible NAs
        colData(scaRAW)$Class <- factor(colData(scaRAW)$Class)
        colData(scaRAW)$CC <- as.numeric(colData(scaRAW)[,CC_col])
        colData(scaRAW)$CC[is.na(colData(scaRAW)$CC)] <- max(as.numeric(colData(scaRAW)$CC), na.rm=T)+1 # Fix possible NAs

        zlmOut <- zlm(~Time+Class+geneson+CC, scaRAW)
        contrast_out <- summary(zlmOut, doLRT='Time')

        summaryDt <- contrast_out$datatable
        outTABLE <- as.data.frame(summaryDt[summaryDt$contrast=="Time" & summaryDt$component=="H",])
        outTABLE2 <- as.data.frame(summaryDt[summaryDt$contrast=="Time" & summaryDt$component=="logFC",])
        outTABLE[,5:8] <- outTABLE2[5:8]
        outTABLE$fdr <- p.adjust(outTABLE[,"Pr(>Chisq)"], method="fdr")
        outTABLE$Symbol <- rowData(scaRAW)$Symbol

        return(outTABLE);
}



#plotting_stuff <- readRDS("All_Plotting_stuff.rds")


# what was this?
#plot(reducedDims(SCE)$tSNE[,1:2], col=cell_colours, pch=16)
#globaldat <- colData(Global)[colnames(Global) %in% colnames(SCE),]
#lineage <- globaldat[,grep("score", colnames(globaldat))]
#gpca <- globaldat[,grep("PCA", colnames(globaldat))]

# Trio of intensities
# Orange = Hep
oranges <- rbind(c(255,245,235), c(254,230,206), c(253,208,162), c(253,174,107), 
		 c(253,141,60), c(241,105,19), c(217,72,1), c(140,45,4))
oranges <- apply(oranges, 1, function(x) {rgb(x[1]/255, x[2]/255, x[3]/255)})
# Blue = Chol
blues <- rbind(c(247,251,255), c(222,235,247), c(198,219,239), c(158,202,225), 
	       c(107,174,214), c(66,146,198), c(33,113,181), c(8,69,148))
blues <- apply(blues, 1, function(x) {rgb(x[1]/255, x[2]/255, x[3]/255)})
# Green = Stem
greens <- rbind(c(247,252,245), c(229,245,224), c(199,233,192), c(161,217,155), 
		c(116,196,118), c(65,171,93), c(35,139,69), c(0,90,50))
greens <- apply(greens, 1, function(x) {rgb(x[1]/255, x[2]/255, x[3]/255)})


Laura_Markers <- list(Chol=c("NDRG1", "DAPK1", "CA9"), Hep=c("AFM", "HMGCS2"), Stem=c("EXOSC9", "HMGB2", "HMMR"))

Lineage <- read.table("Cleaned_Lineage.txt", header=F)
lin_expr_mat <-assays(Global)[["lognorm"]]
rownames(lin_expr_mat) <- rowData(Global)$Symbol
#lin_expr_mat <- lin_expr_mat[rownames(lin_expr_mat) %in% Lineage[,1],]
require("proxy")
cors_out <- simil(lin_expr_mat, t(as.matrix(colData(Global)[,grep("score", colnames(colData(Global)))])))
c_out <- cbind(cors_out[,1], cors_out[,2], cors_out[,3])
colnames(c_out) <- colnames(cors_out)
c_out <- as.data.frame(c_out)
c_out$Lin <- Lineage[match(rownames(cors_out), Lineage[,1]),2]

SCEs <- list()
features <- vector()
for (i in 1:length(plotting_rds)) {
	sce <- readRDS(plotting_rds[i])
	sce$Global.Cycle <- colData(Global)$CC_state[colnames(Global) %in% colnames(sce)]
	SCEs[[plotting_rds[i]]] <- sce
	sce <- SCEs[[i]]
	features <- c(features,as.character(rowData(sce)$Symbol[rowData(sce)$Marker.q.value < 0.05]))
}
sce <- SCEs[[5]]
features_hep <- as.character(rowData(sce)$Symbol[rowData(sce)$Marker.q.value < 0.05])
sce <- SCEs[[1]]
features_chol <- as.character(rowData(sce)$Symbol[rowData(sce)$Marker.q.value < 0.05])

freq_feat <- table(features)
c_out$is.Hep.Feature <- rownames(c_out) %in% features_hep
c_out$is.Chol.Feature <- rownames(c_out) %in% features_chol
c_out$freq_feat <- freq_feat[match(rownames(c_out), names(freq_feat))]
c_out$freq_feat[is.na(c_out$freq_feat)] <- 0

head(c_out[order(c_out[,1]*as.numeric(c_out$is.Hep.Feature), decreasing=T),], 50)
head(c_out[order(c_out[,2]*as.numeric(c_out$freq_feat), decreasing=T),], 50)
head(c_out[order(c_out[,3]*as.numeric(c_out$freq_feat), decreasing=T),], 50)

hep_genes <- rownames(c_out)[ c_out$Hep_score > 0.2 & c_out$Chol_score < -0.1 & c_out$Stem_score < -0.1 & c_out$is.Hep.Feature ]
chol_genes <- rownames(c_out)[ c_out$Chol_score > 0.2 & c_out$Hep_score < -0.1 & c_out$Stem_score < -0.1 & c_out$is.Chol.Feature ]
stem_genes <- rownames(c_out)[ c_out$Stem_score > 0.2 & c_out$Hep_score < -0.1 & c_out$Chol_score < -0.1 & (c_out$is.Chol.Feature | c_out$is.Hep.Feature)]

mat <- matrix("", ncol=3, nrow=max(length(hep_genes), length(chol_genes), length(stem_genes)))
mat[1:length(hep_genes),1] <- hep_genes
mat[1:length(chol_genes),2] <- chol_genes
mat[1:length(stem_genes),3] <- stem_genes
colnames(mat) <- c("Hep", "Chol", "Stem")
write.table(mat, file="Figure3_geneslists.csv", sep=",", row.names=F, col.names=T)

# No just going to do plot for each score.
#stem_scaling <- mean(colData(Global)$Proliferating
#hep_scaling <- mean(Global$Donor == "HCC10")
#chol_scaling <- mean(Global$Donor %in% c("CCA1", "CCA5", "D3DM", "D3EM", "D9DM", "D9EM"))

get_colour_break_pts <- function(g, colours) {
        rows <- which(rowData(Global)$Symbol %in% g)
	if (length(rows) == 1) {
	        vals <- scale(assays(Global)[["norm_exprs"]][rows,])
	} else {
		#vals <- rowMeans(apply(assays(Global)[["norm_exprs"]][rows,], 1, scale)) # no-scale because over-inflates hepatocyte signature.
		vals <- colMeans(assays(Global)[["norm_exprs"]][rows,])
	}
        gaps <- seq(from=min(vals), to=max(vals), by=(max(vals)-min(vals))/length(colours))
        gaps[length(gaps)] <- gaps[length(gaps)]+1;
	return(list(vals=vals, breaks=gaps))
}


line_plots <- function(line_names, g) {
	intens_cols <- rev(brewer.pal(8, "RdBu"))
	intensities <- get_colour_break_pts(g, intens_cols)
	
	d1 <- ceiling(sqrt(length(line_names)))
	d2 <- ceiling(length(line_names)/d1)
	par(mfrow=c(d1, d2))
	par(mar=c(0,0,2,0))	

	for (l in line_names) {
		sce <- SCEs[[paste(l, "PlottingObj.rds", sep="_")]]
		xes <- reducedDims(sce)$tSNE[,1]
		yes <- reducedDims(sce)$tSNE[,2]
		pt_cols <- intensities$vals[colnames(Global) %in% colnames(sce)]
		pt_cols <- cut(pt_cols, breaks=intensities$breaks, include.lowest=T)
		pt_cols <- intens_cols[pt_cols]
		plot(xes, yes, pch=16, col=pt_cols, xaxt="n", yaxt="n", main=l, bty="n")
	}
}

multi_plots <- function(line_names, gene_sets) {
	d1 <- ceiling(length(line_names))
        d2 <- ceiling(length(gene_sets))
        par(mfcol=c(d1, d2))
	for(gi in 1:length(gene_sets)) {
		g <- gene_sets[[gi]]
		intens_cols <- rev(brewer.pal(8, "RdBu"))
	        intensities <- get_colour_break_pts(g, intens_cols)

		for (li in 1:length(line_names)) {
			l <- line_names[li]
			par(mar=c(0,1,1,0))
			sce <- SCEs[[paste(l, "PlottingObj.rds", sep="_")]]
			xes <- reducedDims(sce)$tSNE[,1]
			yes <- reducedDims(sce)$tSNE[,2]
			pt_cols <- intensities$vals[colnames(Global) %in% colnames(sce)]
			pt_cols <- cut(pt_cols, breaks=intensities$breaks, include.lowest=T)
			pt_cols <- intens_cols[pt_cols]
			plot(xes, yes, pch=16, col=pt_cols, xaxt="n", yaxt="n", bty="n")
			if (gi == 1) {
				title(ylab=l, line=0, font=2, cex=1.5)	
			}
			if (li == 1) {
				title(main=names(gene_sets)[gi], line=0, cex=1.5)	
			}
		}
	}
}


####fancy_trio_plot 
# colour by all three lineages at once - colour = which lineage is highest? + intensity = absolute value of expression.
# Need to create palettes for chol-stem and hep-stem
stem_chol_palette <- colorRampPalette(c(rev(greens), oranges))
stem_hep_palette <- colorRampPalette(c(rev(greens), blues))
gene_sets=list(Hep=hep_genes, Chol=chol_genes, Stem=stem_genes)
col_sets=list(hep=oranges, chol=blues, stem=greens)

line_plots(c("CCA1", "HCC10"), gene_sets$Hep)
line_plots(c("CCA1", "HCC10"), gene_sets$Chol)
line_plots(c("CCA1", "HCC10"), gene_sets$Stem)


png("Figure3.png", width=4*3, height=2*4, units="in", res=300) ########## HERE is the actual Figure!!!
multi_plots(c("CCA1", "HCC10"), gene_sets)
dev.off()

### VV Other stuff VV ###

hep <- get_colour_break_pts(gene_sets$hep, col_sets$hep)
chol <- get_colour_break_pts(gene_sets$chol, col_sets$chol)
stem <- get_colour_break_pts(gene_sets$stem, col_sets$stem)

# HCC10 
sce <- SCEs[[5]]
hep_local <- hep$vals[colnames(Global) %in% colnames(sce)]
chol_local <- chol$vals[colnames(Global) %in% colnames(sce)]
stem_local <- stem$vals[colnames(Global) %in% colnames(sce)]
xes <- reducedDims(sce)$tSNE[,1]
yes <- reducedDims(sce)$tSNE[,2]



# slingshot requires reduced dims before hand.
# First do ICA? PCA?
# Global-scaled lineage colouring
#### Do Slingshot ####
	SCE <- readRDS(plotting_rds[i])
	palette <- SCE@metadata$palette
        cell_colours <- palette[SCE$Manual_Clusters]
        nCs <- factor_counts(SCE$Manual_Clusters)

        # Feature Selection
        genes <- rownames(SCE)[rowData(SCE)$is.Figure.Feature]
	# PCA
	pca <- SCE@metadata$pca
	n_dims <- ncol(reducedDims(SCE)[["PCA"]])
        if (n_dims < 3) {n_dims <- 3}
	pcs <- pca$rotation[,1:n_dims]

	# Slingshot
	globaldat <- colData(Global)[colnames(Global) %in% colnames(SCE),]
	lineage <- globaldat[,grep("score", colnames(globaldat))]
	gpca <- globaldat[,grep("PCA", colnames(globaldat))]
	#reducedDims(SCE) <- SimpleList(PCA=pcs, Lin=as.matrix(lineage), gPCA=as.matrix(gpca))
	SCE <- slingshot(SCE, clusterLabels='Manual_Clusters', reducedDim='PCA')

	fitted <- SlingshotDataSet(SCE)
	#thing <- apply(fitted, 2, scale)
	thing1 <- fitted@curves$curve1$s
	chol_score_curve <- (-thing1[,1]-thing1[,4])/2
	hep_score_curve <- (thing1[,1]-thing1[,2]-thing1[,4])/3
	stem_score_curve <- (-thing1[,1]-thing1[,2]+thing1[,4]+thing1[,3])/4

	thing2 <- gpca;
	chol_score_pts <- (-thing2[,1]-thing2[,4])/2
        hep_score_pts <- (thing2[,1]-thing2[,2]-thing2[,4])/3
        stem_score_pts <- (-thing2[,1]-thing2[,2]+thing2[,4]+thing2[,3])/4

	if (median(hep_score_pts) > median(chol_score_pts)) {
		plot(hep_score_pts, stem_score_pts, col=cell_colours, pch=16, xlab="Cholangiocyte ---------- Hepatocyte", ylab="Stemness")
		reorder <- fitted@curves$curve1$ord
		lines(hep_score_curve[reorder], stem_score_curve[reorder], lwd=2)
	} else {
		plot(-chol_score_pts, stem_score_pts, col=cell_colours, pch=16, xlab="Cholangiocyte ---------- Hepatocyte", ylab="Stemness")
		reorder <- fitted@curves$curve1$ord
		lines(- chol_score_curve[reorder], stem_score_curve[reorder], lwd=2)
	}
	#boxplot(SCE$slingPseudotime_1~globaldat$CC_state)
	#plot(hep_score_pts-chol_score_pts, stem_score_pts, col=c("grey10", "grey30", "grey50", "grey70", "grey90")[cut(SCE$slingPseudotime_1, breaks=5)], pch=16)

	require(gam)
	t <- SCE$slingPseudotime_1

	Y <- log1p(assays(SCE)$norm_exprs)
	# fit a GAM with a loess term for pseudotime
	gam.pval <- apply(Y,1,function(z){
	  d <- data.frame(z=z, t=t)
	  tmp <- gam(z ~ lo(t), data=d)
	  p <- summary(tmp)[4][[1]][1,5]
	  p
	})







# Common scale for lineage markers

require("RColorBrewer")
colours <- brewer.pal(9, "GnBu")
Laura_Breaks <- list(Chol=list(), Hep=list(), Stem=list())


for (type in c("Chol", "Hep", "Stem")) {
        for (g in Laura_Markers[[type]]) {
                row <- which(rowData(Global)$Symbol == g)
                vals <- assays(Global)[["norm_exprs"]][row,]
                #quants <- quantile(vals, probs=seq(from=0, to=1, length=length(colours)+1))
                gaps <- seq(from=0, to=max(vals), by=max(vals)/length(colours))
                gaps[length(gaps)] <- gaps[length(gaps)]+1;
                Laura_Breaks[[type]][[g]] <- gaps;
        }
}

Mega_plots <- function(line_name) {
        i <- which(rds_names == line_name);
        cluGraph <- nn_graph[[i]]$graph
        cluGraph_layout <- nn_graph[[i]]$layout
        palette <- V(cluGraph)$col
        par(mar=c(0,2,0,0))
        plot(cluGraph, vertex.color = V(cluGraph)$col, edge.width=E(cluGraph)$weight,
                edge.lty=E(cluGraph)$type, edge.color="black", vertex.label.cex=3, vertex.label=NA,
                vertex.size=V(cluGraph)$size, vertex.label.color="black", rescale=FALSE,
                xlim=c(min(cluGraph_layout[,1]-0.1),max(cluGraph_layout[,1])+0.1), ylim=c(0,1), layout=cluGraph_layout)
        axis(side=2, at=c(0,0.33, 0.66, 1), labels=c("", "", "", ""))
        mtext("Diff", side=2, at=0.33/2, line=0.5)
        mtext("Prog", side=2, at=(0.33+0.66)/2, line=0.5)
        mtext("Stem", side=2, at=(0.66+1)/2, line=0.5)

        # Legend
        #par(mar=c(0,0,0,0))
        #blank_plot()
        #legend("topright", line_specific_groups[[line_name]], fill=palette, bty="n")

        par(mar=c(4,4,2,1))
        xes <- cell_coords[[i]][,1]
        yes <- cell_coords[[i]][,2]
        plot(xes, yes, xlab="Dim 1", ylab="Dim 2",
                main=line_name, col=cell_colors[[i]], pch=16)

       for (type in c("Chol", "Hep", "Stem")) {
               for (g_i in 1:2) {
                       g <- Laura_Markers[[type]][g_i]
                       g_breaks <- Laura_Breaks[[type]][[g]]
                       par(mar=c(0,0,1.5,0))
                       vals <- expr_mats[[i]][rownames( expr_mats[[i]] ) == g, ]
                       
                       plot(xes, yes, pch=16, col=colours[cut(vals, breaks=g_breaks, include.lowest=T)], 
                               xaxt="n", yaxt="n", main=g, xlab="", ylab="", bty="n")  
               }
        }

}

