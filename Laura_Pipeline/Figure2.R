# Internal variability.


# SC3 clustering
# 

source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")



Laura_Markers <- list(Chol=c("NDRG1", "DAPK1", "CA9"), Hep=c("AFM", "AFM"), Stem=c("EXOSC9", "HMGB2", "HMMR"))

# Lineage markers
#Chol_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Chol.txt", header=TRUE)
#Hep_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Hep.txt", header=TRUE)

#Hep_both <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & grepl("Hep", Hep_lineage[,2]), 1]
#Chol_both <- Chol_lineage[ grepl("Prog", Chol_lineage[,2]) & grepl("Chol", Chol_lineage[,2]), 1]
#Prog_both <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Prog",1], 1]

#Conflict1 <- Hep_lineage[ grepl("Hep", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Chol",1], 1]
#Conflict2 <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Chol",1], 1]
#Conflict3 <- Hep_lineage[ grepl("Hep", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Prog",1], 1]

#Conflicts <- c(as.character(Conflict1), as.character(Conflict2), as.character(Conflict3))

#Chol_lineage <- Chol_lineage[Chol_lineage[,1] %in% marker_genes & Chol_lineage[,1] %in% keep_genes,]
#Hep_lineage <- Hep_lineage[Hep_lineage[,1] %in% marker_genes & Hep_lineage[,1] %in% keep_genes,]

#Chol_lineage[,2] <- as.character(Chol_lineage[,2])
#Chol_lineage[Chol_lineage[,2] == "Prog",2] <- "Chol-Prog"
#Chol_lineage[Chol_lineage[,2] == "Chol",2] <- "Chol-Mature"
#Chol_lineage <- Chol_lineage[!(Chol_lineage[,1] %in% Conflicts) & !(Chol_lineage[,1] %in% Chol_both),]

#Hep_lineage[,2] <- as.character(Hep_lineage[,2])
#Hep_lineage[Hep_lineage[,2] == "Prog",2] <- "Hep-Prog"
#Hep_lineage[Hep_lineage[,2] == "Hep",2] <- "Hep-Mature"
#Hep_lineage <- Hep_lineage[!(Hep_lineage[,1] %in% Conflicts),]

#Lineage <- rbind(Chol_lineage,Hep_lineage)
#Lineage[Lineage[,1] %in% Prog_both,2] <- "Common-Prog"
#Lineage <- Lineage[!duplicated(Lineage[,1]),]
#Lineage[,1] <- as.character(Lineage[,1])
#Lineage[Lineage[,1] == "05-Mar",1] <- "MARCH5"


# Read in each manual clustering
clustered_rds <- c("CCA1_manual_SC3.rds", "CCA5_manual_SC3.rds", "HCC6_manual_SC3.rds", "HCC23_manual_SC3.rds", "HCC10_manual_SC3.rds", "HCC24_manual_SC3.rds", "D3DM_manual_SC3.rds", "D3EM_manual_SC3.rds", "D9DM_manual_SC3.rds", "D9EM_manual_SC3.rds");
rds_names <- c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24", "D3DM", "D3EM", "D9DM", "D9EM")
expr_type <- "norm_exprs";
expr_mats <- list()
cell_colors <- list();
nn_graph <- list();
cell_coords <- list();
params <- list(CCA1=cbind(c(0.5, 0.5+0.1, 0.5, 0.5-0.1),c(0.5,0,1,0)), 
		CCA5=cbind(c(0.5, 0.5+0.25, 0.5, 0.5-0.3, 0.5-0.0, 0.5+0.1), c(0.5, 0.3, 1, 0, 0.25, 0.25)), 
		HCC6=cbind(c(0.5-0.15, 0.5, 0.5+0.15, 0.5), c(0.5, 0, 0.5, 1)), 
		HCC23=cbind(c(0.5, 0.5, 0.5), c(1, 0, 0.5)), 
		HCC10=cbind(c(0.5-0.07, 0.5+0.2, 0.5, 0.5+0.2, 0.5-0.07, 0.5-0.3, 0.5-0.3), c(0.5, 0, 1, 0.5, 0, 0.4, 0.25)), 
		HCC24=cbind(c(0.5+0.15, 0.5-0.15, 0.5, 0.5, 0.5+0.3, 0.5), c(0.20, 0.20, 0.6, 1, 0, 0)),
		D3DM=cbind(c(0.5+0.2, 0.5, 0.5-0.2, 0.5), c(0, 0.5, 0, 1)),
		D3EM=cbind(c(0.5, 0.5+0.17, 0.5-0.075, 0.5), c(0, 0.5, 0.5, 1)),
		D9DM=cbind(c(0.5, 0.5+0.1, 0.5-0.1, 0.5+0.2, 0.5, 0.5), c(0.5, 0.25, 0.75, 0, 0, 1)),
		D9EM=cbind(c(0.5, 0.5, 0.5), c(0.5, 1, 0))
		); # graph layouts
line_specific_genes <- list(CCA1=c("ZWINT", "NDRG1", "EIF5A", "HERPUD1"),
			    CCA5=c("HMGB2", "LDHA", "CD81", "ANXA4"),
			    HCC6=c("MCM7", "GAL", "CST3", "IGFBP3"),
			    HCC23=c("HMGB2", "F10", "SERPINA1", "CBS"),
			    HCC10=c("HMGB2","CEBPA","ABCA1","AFP"),
			    HCC24=c("ZWINT","ALB","KLK1","TUSC3"),
			    D3DM=c("CDK2AP1","CD24","TESC","P4HA1"),
			    D3EM=c("YWHAB","ENO2","BTG1","MIF"),
			    D9DM=c("CEACAM6","SERPINA1","NDRG1","KCNT1"),
			    D9EM=c("EIF5A","HERPUD1","CA9","EFNA1")
			) # genes for violin plots
line_specific_groups <- list(CCA1=c("Prog", "Chol", "CSC", "Hypoxic"),
			    CCA5=c("Chol", "Hypoxic", "CSC", "Unk", "Hep"),
			    HCC6=c("Prog1", "Stress", "Prog2", "CSC"),
			    HCC23=c("CSC", "Clot-Hep", "Prog"),
			    HCC10=c("Prog1", "Hep", "CSC", "iHep", "Hypoxic", "Prog2"),
			    HCC24=c("Clot-Hep", "Prog1", "Prog2", "CSC", "Hep1", "Hep2"),
			    D3DM=c("Chol", "Prog", "Hypoxic", "Cycling"),
			    D3EM=c("Hypoxic", "Chol1", "Chol2", "Prog"),
			    D9DM=c("Prog1", "Prog2", "Prog3", "Clot-Hep", "Chol1", "Chol2"),
			    D9EM=c("Chol1", "Cycling", "Chol2")
			) # cluster names
			
pca_plot_infos <- list();
for (i in 1:length(clustered_rds)) {
	# save cluster IDs & markers
	# add lineage annotations to markers
	# do tSNE & sil_nn graph
	# save coords
	require("SingleCellExperiment")
	require("scater")
	require("Rtsne")
	require("M3Drop")
	require("CellTypeProfiles")
	perp <- 30
	# Set up
	SCE <- readRDS(clustered_rds[i])
	SCE <- SCE[!is.na(rowData(SCE)$biotype),]
	SCE <- SCE[rowData(SCE)$feature_symbol != "",]
	SCE <- SCE[rowData(SCE)$biotype == "protein_coding",]
	palette <- cluster_col(max(SCE$Manual_Clusters))
	cell_colours <- palette[SCE$Manual_Clusters]
	nCs <- factor_counts(SCE$Manual_Clusters)

	# save colours
	names(palette) <- 1:length(palette)
	SCE@metadata$palette <- palette
	SCE@metadata$C_keep <- nCs > 10
	SCE@metadata$C_names <- line_specific_groups[[rds_names[i]]]

	# Feature Selection
	m3d_mat <- M3DropConvertData(assays(SCE)[["norm_exprs"]], is.counts=FALSE, is.log=TRUE)
	FS <- M3DropFeatureSelection(m3d_mat, suppress.plot=TRUE)
	mark_cells <- SCE$Manual_Clusters %in% names(nCs)[nCs > 10]
	mark_mat <- assays(SCE)[["norm_exprs"]][,mark_cells]
	marks <- complex_markers(mark_mat, factor(SCE$Manual_Clusters[mark_cells]))
	m_FS <- marks[marks$AUC >= 0.8 & marks$q.value < 0.05,]
	if (i == 5) {m_FS <- marks[marks$AUC > 0.75 & marks$q.value < 0.05,]}
	genes <- union(FS$Gene[1:min(500, nrow(FS))], rownames(m_FS))
	#M3DropExpressionHeatmap(genes, m3d_mat, cell_labels=SCE$Manual_Clusters)
	
	marks <- marks[match(rownames(marks), rownames(SCE)),]
	colnames(marks) <- paste("Marker", colnames(marks), sep=".")
	rowData(SCE) <- rowData(SCE)[,-grep("marker", colnames(rowData(SCE)))]
	rowData(SCE) <- cbind(rowData(SCE), marks);	
	rowData(SCE)$is.Figure.Feature <- rownames(SCE) %in% genes;

	# PCA
	#X11()
	pca <- prcomp(assays(SCE)[["norm_exprs"]][rownames(SCE) %in% genes,])
	#plot(pca$sdev[1:20]^2/sum(pca$sdev^2)*100, type="b", xlab="Component", ylab="Variance (%)")
	p <- apply(pca$rotation[,1:50], 2, function(x) {summary(lm(x~SCE$Manual_Clusters))$coefficients[2,4]})
	n_dims <- which(p > 0.05);
	n_dims <- min(n_dims[n_dims > min(which(p < 0.05))]) -1
	if (n_dims < 3) {n_dims <- 3}
	#abline(v=n_dims+0.5, lty=3, col="grey65")
	pca_plot_infos[[i]] <- list(dat=pca$sdev[1:20]^2/sum(pca$sdev^2)*100, line=n_dims+0.5);


	# tSNE
	set.seed(28198)
	tSNE <- Rtsne(t(assays(SCE)[["norm_exprs"]][rownames(SCE) %in% genes,]), initial_dims=n_dims, dims=2, perplexity=perp)
	#plot(tSNE$Y[,1], tSNE$Y[,2], col=cell_colours, pch=16, xlab="Dim 1", ylab="Dim 2")
	coords <- tSNE$Y

	reducedDims(SCE) <- SimpleList(PCA=pca$rotation[,1:n_dims], tSNE=tSNE$Y)
	SCE@metadata$pca <- pca

	# Silhouette
	mat <- t(assays(SCE)[["norm_exprs"]][rownames(SCE) %in% genes,])
	require("proxy")
	distances <- proxy::dist(mat, method="Euclidean")

	sil <- cluster::silhouette(SCE$Manual_Clusters, distances)
	neighbours <- table(factor(sil[,1], levels=names(nCs)), 
			    factor(sil[,2], levels=names(nCs)))/nCs

	# Graph
	trimm <- function(x) {
	        thresh <- max(x)*0.5
	        x[x < thresh] <- 0
	        return(x)
	}
	nn_sil <- t(apply(neighbours, 1, trimm))
	recip <- which(nn_sil > 0, arr.ind=T)
	recip <- cbind(apply(recip,1,min), apply(recip,1,max))
	recip <- recip[order(recip[,1], recip[,2]),]
	e_type <- 1-as.numeric(duplicated(recip))
	e_type[which(e_type == 0)-1] <- 2
	e_type <- e_type[e_type != 0];

	require("igraph")
	g_sil <- igraph::graph_from_adjacency_matrix(nn_sil+t(nn_sil), mode="undirected", weighted=TRUE)
	V(g_sil)$col <- palette
	size_max=20
	V_size <- nCs/sum(nCs)*size_max #nCs-min(nCs); V_size<-(V_size/max(V_size)+0.5)*30
	V_size[V_size > size_max*1.5/length(nCs)] <- size_max*1.5/length(nCs)
	V(g_sil)$size <- V_size*size_max/max(V_size)
	E(g_sil)$type <- 3-e_type
	E(g_sil)$weight <- round(E(g_sil)$weight*3, digits=1)
	lay=params[[i]]*0.8+0.1
	set.seed(1053)
	#par(mar=c(0,0,0,0))
	#plot(g_sil, vertex.color = V(g_sil)$col, edge.width=E(g_sil)$weight,
        #	edge.lty=E(g_sil)$type,
	#       	edge.color="black", vertex.size=V(g_sil)$size, label.color="black"
	#	, rescale=FALSE, xlim=c(0,1), ylim=c(0,1), layout=lay
	#	)
	table(SCE$Cycle, SCE$Manual_Clusters)

	SCE@metadata$sil_nn <- nn_sil

	# Save Results
	nn_graph[[i]] <- list(graph=g_sil, layout=lay)
	cell_coords[[i]] <- coords
	mat <- assays(SCE)[["norm_exprs"]]
	rownames(mat) <- rowData(SCE)$feature_symbol
	expr_mats[[i]] <- mat
	cell_colors[[i]] <- cell_colours

	saveRDS(SCE, file=paste(rds_names[i], "PlottingObj.rds", sep="_"))
}

### Supplementary
names(pca_plot_infos) <- rds_names
suppl_plot_dims <- ceiling(sqrt(length(pca_plot_infos)))
suppl_plot_dims[2] <- ceiling(length(pca_plot_infos)/suppl_plot_dims)
par_defaults <- par();
png("Figure2_Supplementary_SampleSpecific_Choosing_PCs.png", width=2*suppl_plot_dims[1], height=2*suppl_plot_dims[2], units="in", res=300)
par(mfrow=c(suppl_plot_dims[2], suppl_plot_dims[1]))
par(mar=c(4,4,2,1))
for (i in 1:length(pca_plot_infos)) {
	plot(pca_plot_infos[[i]]$dat, type="b", xlab="Component", ylab="Variance (%)", main=names(pca_plot_infos)[i])
	abline(v=pca_plot_infos[[i]]$line, lty=3, col="grey65")
	legend("topright", paste("dims =", floor(pca_plot_infos[[i]]$line)), col="white", bty="n")
}
dev.off()
	

### Mega plots with violins of markers.
STUFF <- list(nn_graphs=nn_graph, cell_coords=cell_coords, cell_colors=cell_colors, layouts=params)
saveRDS(STUFF , file="All_Plotting_stuff.rds")


# Cluster Plots
# sil_nn graph + tSNE + 4 marker gene violin plots	
source("~/R-Scripts/violin_plot.R")
source("~/R-Scripts/Blank_plot.R")
# for markers being used by Laura/Meri Group?
	# Read in Global Obj
	# Create common colour scale for lineage markers across all lines
Global <- readRDS("Global_SCE.rds")

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

	nCs <- table(factor(cell_colors[[i]], levels=palette))
	outliers <- which(nCs < 5)

	if (length(outliers) > 0) {
		cluGraph <- nn_graph[[i]]$graph
		keep <- V(cluGraph) != outliers
		cluGraph <- igraph::induced_subgraph(cluGraph, which(keep))
		cluGraph_layout <- cluGraph_layout[keep,]
		palette <- palette[keep]

	}
		
	keep_cells <- cell_colors[[i]] %in% palette

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
	xes <- cell_coords[[i]][keep_cells,1]
	yes <- cell_coords[[i]][keep_cells,2]
	plot(xes, yes, xlab="Dim 1", ylab="Dim 2", 
		main=line_name, col=cell_colors[[i]][keep_cells], pch=16)

	for (g in line_specific_genes[[i]]) {
		# Need to add size-threshold for clusters!!!!!
		vio_data <- list();
		for (ele in 1:length(palette)) {
			cells <- cell_colors[[i]] == palette[ele]
			if (length(cells) < 5) {vio_data[[i]] <- NA; next;}
			grow <- which(rownames( expr_mats[[i]] ) == g)
			if (length(grow) != 1) {
				a <- rowMeans(expr_mats[[i]][grow,])
				grow <- grow[a == max(a)]
				grow <- grow[1];
			}
			vio_data[[ele]] <- expr_mats[[i]][grow, cells]
		}
		par(mar=c(4,4,2,1))
		vioplot(vio_data, col=palette, drawRect=FALSE, names=line_specific_groups[[line_name]], las=2)
		title(ylab="Expression", main=g)
	#for (type in c("Chol", "Hep", "Stem")) {
	#	for (g_i in 1:2) {
	#		g <- Laura_Markers[[type]][g_i]
	#		g_breaks <- Laura_Breaks[[type]][[g]]
	#		par(mar=c(0,0,1.5,0))
	#		vals <- expr_mats[[i]][rownames( expr_mats[[i]] ) == g, ]
	#		
	#		plot(xes, yes, pch=16, col=colours[cut(vals, breaks=g_breaks, include.lowest=T)], 
	#			xaxt="n", yaxt="n", main=g, xlab="", ylab="", bty="n")	
	#	}
	#}
	}

}


png("Figure2_CCA1_MegaViolin.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))
Mega_plots("CCA1")
dev.off()
png("Figure2_HCC10_MegaViolin.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))
Mega_plots("HCC10")
dev.off()

png("Figure2S_CCA5_MegaViolin.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))
Mega_plots("CCA5")
dev.off()
png("Figure2S_HCC6_MegaViolin.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))
Mega_plots("HCC6")
dev.off()
png("Figure2S_HCC23_MegaViolin.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))
Mega_plots("HCC23")
dev.off()
png("Figure2S_HCC24_MegaViolin.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))
Mega_plots("HCC24")
dev.off()
png("Figure2S_D3DM_MegaViolin.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))
Mega_plots("D3DM")
dev.off()
png("Figure2S_D3EM_MegaViolin.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))
Mega_plots("D3EM")
dev.off()
png("Figure2S_D9DM_MegaViolin.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))
Mega_plots("D9DM")
dev.off()
png("Figure2S_D9EM_MegaViolin.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))
Mega_plots("D9EM")
dev.off()


graphics::layout(rbind(c(1,2,3,5), c(1,2,4,6)), 
		widths=c(1, 3, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5 ))

Mega_plots("CCA5")
Mega_plots("HCC6")
Mega_plots("HCC23")
Mega_plots("HCC10")
Mega_plots("HCC24")
Mega_plots("D3DM")
Mega_plots("D3EM")
Mega_plots("D9DM")
Mega_plots("D9EM")


#graphics::layout(rbind(c(1,2,3,5,7), c(1,2,4,6,8)), 
#		widths=c(1, 3, 1.5, 1.5, 1.5), 
#		heights=c(1.5, 1.5, 1.5, 1.5))

png("Figure2_CCA1_Mega.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5,7), c(1,2,4,6,8)), 
		widths=c(1, 3, 1.5, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5, 1.5))
Mega_plots("CCA1")
dev.off()

png("Figure2_HCC10_Mega.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5,7), c(1,2,4,6,8)), 
		widths=c(1, 3, 1.5, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5, 1.5))
Mega_plots("HCC10")
dev.off()

png("Figure2_HCC6_Mega.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5,7), c(1,2,4,6,8)), 
		widths=c(1, 3, 1.5, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5, 1.5))
Mega_plots("HCC6")
dev.off()

png("Figure2_HCC23_Mega.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5,7), c(1,2,4,6,8)), 
		widths=c(1, 3, 1.5, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5, 1.5))
Mega_plots("HCC23")
dev.off()

png("Figure2_HCC24_Mega.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5,7), c(1,2,4,6,8)), 
		widths=c(1, 3, 1.5, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5, 1.5))
Mega_plots("HCC24")
dev.off()

png("Figure2_CCA5_Mega.png", width=8.5*1.7, height=6, units="in", res=300)
graphics::layout(rbind(c(1,2,3,5,7), c(1,2,4,6,8)), 
		widths=c(1, 3, 1.5, 1.5, 1.5), 
		heights=c(1.5, 1.5, 1.5, 1.5))
Mega_plots("CCA5")
dev.off()


####### Mega Marker Tables #######
source("~/Collaborations/LiverOrganoids/Laura_Pipeline/99_Functions.R")
require("CellTypeProfiles")
require("SingleCellExperiment")
require("scater")

CC_genes <- load_CC(set="cycling")
CC_genes <- c(as.character(CC_genes$Whitfield[,2]), as.character(CC_genes$Tirosh[,1]))

GO <- read.delim("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/GO/hsapiens_80_GO_Annoations_Emsembl.out", sep="\t", header=F)
GO_cc <- GO[GO[,3] == "cell cycle",1]

schwalie_StemCell <- c("MFAP2", "HSPB6", "GBP3", "MYC", "CCND2", "RAC2", "ZFP36L2", "VGLL4", "ADGRG1", "IDH2", "CDK6", "SH3BGRL", "IFITM3", "LCP1", "ETV6", "RPGRIP1", "ERI3", "SDCBP", "DAP", "JTB", "GSTO1", "DAPP1", "CTSZ")
GO_telo <- unique(GO[GO[,3] == "telomere maintenance",1])
LigRec <- read.delim("~/Data/LigandReceptorPairs.csv", sep=",", header=T)
GO_tf <- unique(GO[grep("transcription factor activity", GO[,3]),1])

Lineage <- read.table("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/Cleaned_Lineage.txt")

add_anno <- function(marker_tab) {
        marker_tab$is.cycle <- rownames(marker_tab) %in% GO_cc
        marker_tab$is.telo <- rownames(marker_tab) %in% GO_telo
        marker_tab$is.tf <- rownames(marker_tab) %in% GO_tf
        marker_tab$ReceptorOf <- LigRec[match(marker_tab$Symbol,LigRec$Receptor.ApprovedSymbol),"Ligand.ApprovedSymbol"]
        marker_tab$LigandTo <- LigRec[match(marker_tab$Symbol,LigRec$Ligand.ApprovedSymbol),"Receptor.ApprovedSymbol"]
        marker_tab$Lineage <- Lineage[match(marker_tab$Symbol,Lineage[,1]),2]
        return(marker_tab)
}

# Each line:
#	GeneID	GeneSymbol means, cell-type marker-ness, annotations, biotype
# Cluster table:
#	ClusterID ClusterName nCells, CCStage, nSigMarkers, most similar clusters, 
clustered_rds <- c("CCA1_manual_SC3.rds", "CCA5_manual_SC3.rds", "HCC6_manual_SC3.rds", "HCC23_manual_SC3.rds", "HCC10_manual_SC3.rds", "HCC24_manual_SC3.rds", "D3DM_manual_SC3.rds", "D3EM_manual_SC3.rds", "D9DM_manual_SC3.rds", "D9EM_manual_SC3.rds");
rds_names <- c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24", "D3DM", "D3EM", "D9DM", "D9EM")
expr_type <- "norm_exprs";

Global <- readRDS("Global_SCE.rds")

make_cluster_marker_tables <- function(line_name) {
        # Set up
	i <- which(rds_names == line_name);
	file = clustered_rds[i]

        SCE <- readRDS(file)
	# Cluster stats
        nCs <- factor_counts(SCE$Manual_Clusters)
	Cs_names <- line_specific_groups[[line_name]]
	names(Cs_names) <- names(nCs)[nCs > 10]
	Cs_ids <- names(nCs);
	if (!identical(colnames(Global)[Global$Donor == line_name], colnames(SCE))) {stop("Global-local mismatch")}  # checking
	Cs_cc <- table(Global$CC_state[Global$Donor == line_name], SCE$Manual_Clusters)

        # Calculate Markers
        mark_cells <- SCE$Manual_Clusters %in% names(nCs)[nCs > 10]
        mark_mat <- assays(SCE)[[expr_type]][,mark_cells]
        marks <- complex_markers(mark_mat, factor(SCE$Manual_Clusters[mark_cells]))
	marks$is.GoodMarker <- marks$AUC > 0.75 & marks$q.value < 0.05

	# Mean Expression
	expr_table <- my_row_mean_aggregate(assays(SCE)[["logcounts"]][,mark_cells], SCE$Manual_Clusters[mark_cells])

	# Assemble Marker Table
	MARKERS <- cbind(marks, expr_table, rowData(SCE)$feature_symbol, rowData(SCE)$biotype)
	colnames(MARKERS)[grep("symbol", colnames(MARKERS))] <- "Symbol"
	colnames(MARKERS)[grep("biotype", colnames(MARKERS))] <- "GeneType"
	MARKERS <- add_anno(MARKERS); 
	MARKERS <- MARKERS[MARKERS$AUC != -1 ,]

	write.table(MARKERS, file=paste(line_name, "ManualClustering_MarkerTable.csv", sep="_"), sep=",")

	Clust_marks <- MARKERS[MARKERS$is.GoodMarker, ]
	unique_marks <- Clust_marks[,2:( ncol(marks) -3)]
	unique_marks <- unique_marks[rowSums(unique_marks) == 1,]
	Cs_marks <- colSums(unique_marks)

	# Lineage markers by cluster
	C_assigned <- apply(unique_marks, 1, function(x){which(x==1)})
	C_assigned <- factor(as.character(C_assigned), levels=colnames(Cs_cc))
	Cs_marks_lin <- table(MARKERS$Lineage[rowSums(MARKERS[,2:( ncol(marks) -3)]) == 1 & MARKERS$is.GoodMarker], C_assigned)

	# Cluster-cluster similarity
	c_expr_table <- expr_table[marks$is.GoodMarker,]
	Cs_cSim <- cor(c_expr_table)
	rownames(Cs_cSim) <- paste("Cor_with", rownames(Cs_cSim), sep="_")
		
	# Assemble Cluster Table
	Cs_names <- Cs_names[match(names(nCs), names(Cs_marks))]
	Cs_marks_lin <- Cs_marks_lin[,match(names(nCs), names(Cs_marks))]
	Cs_cSim <- Cs_cSim[,match(names(nCs), names(Cs_marks))]
	Cs_marks <- Cs_marks[match(names(nCs), names(Cs_marks))]

	CLUSTERS <- rbind(Cs_ids, Cs_names, nCs, Cs_cc, Cs_marks, Cs_marks_lin, round(Cs_cSim, digits=2))
	rownames(CLUSTERS)[rownames(CLUSTERS) == "nCs"] <- "n_cells"
	rownames(CLUSTERS)[rownames(CLUSTERS) == "Cs_ids"] <- "ID_Num"
	rownames(CLUSTERS)[rownames(CLUSTERS) == "Cs_names"] <- "Suggest_Name"
	rownames(CLUSTERS)[rownames(CLUSTERS) == "Cs_marks"] <- "Unique_Markers"
	CLUSTERS <- CLUSTERS[rownames(CLUSTERS) != "G0",]
	
	write.table(CLUSTERS, file=paste(line_name,"ManualClustering_ClusterTable.csv", sep="_"), sep=",")
}

tmp_names <- rds_names[6:10]
for (line in tmp_names) {
	print(line)
	make_cluster_marker_tables(line)
}
