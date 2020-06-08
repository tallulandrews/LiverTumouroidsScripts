library(SingleCellExperiment)
require("scater")
set.seed(1932)
# Hepatocyte Ref
Camp <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/camp_human.rds")
Camp <- toSingleCellExperiment(Camp);

# Chol Ref
D3Diff <- readRDS("D3EM_manual_SC3.rds")
D9Diff <- readRDS("D9DM_manual_SC3.rds")
D3Exp <- readRDS("D3DM_manual_SC3.rds")
D9Exp <- readRDS("D9EM_manual_SC3.rds")

con_symbol <- rowData(D3Diff)$Symbol[rowData(D3Diff)$Symbol %in% rowData(D9Diff)$Symbol]
D3Diff <- D3Diff[rowData(D3Diff)$Symbol %in% con_symbol,]
D9Diff <- D9Diff[rowData(D9Diff)$Symbol %in% con_symbol,]

Camp <- Camp[rownames(Camp) %in% con_symbol,]
D3Diff <- D3Diff[match(rownames(Camp), rowData(D3Diff)$Symbol),]
D9Diff <- D9Diff[match(rownames(Camp), rowData(D9Diff)$Symbol),]
Camp <- Camp[rownames(Camp) %in% rowData(D3Diff)$Symbol,]


# Merge Refs
Camp_logcounts <- assays(Camp)[["logcounts"]]
Chol1_logcounts <- assays(D3Diff)[["logcounts"]]
Chol2_logcounts <- assays(D9Diff)[["logcounts"]]
cell_type <- c(as.character(Camp$cell_type1), rep("Chol", times=ncol(Chol1_logcounts)+ncol(Chol2_logcounts)));
mat <- cbind(Camp_logcounts, Chol1_logcounts, Chol2_logcounts)
sf <- colSums(mat)
mat <- t(t(mat)/sf*median(sf));


# Make Ref
library("scmap")
RefSCE<-SingleCellExperiment(assays=list(logcounts=as.matrix(mat)), colData=data.frame(cell_type1=cell_type), rowData=data.frame(feature_symbol=sub(" ", "-", rownames(mat))))
RefSCE <- RefSCE[,!RefSCE$cell_type1 %in% c("erythroblasts", "lymphoblasts", "hepatic", "Kupffer")]
assays(RefSCE)[["counts"]] <- assays(RefSCE)[["logcounts"]]
#png("All_Liver_Ref_FS.png", width=5, height=5, units="in", res=300)
RefSCE <- scmap::selectFeatures(RefSCE, suppress_plot=FALSE, n_features=1000)
#dev.off()
RefSCE <- scmap::indexCluster(RefSCE)

require("gplots")
#png("All_Liver_Ref_Profile.png", width=5, height=5, units="in", res=300)
heatmap.2(as.matrix(log2(metadata(RefSCE)$scmap_cluster_index+1)), trace="none")
#dev.off()

# Read in each line and map them

#Visualization Stuff
#STUFF <- list(nn_graphs=nn_graph, cell_coords=cell_coords, cell_colors=cell_colors, layouts=params)
require("igraph")
plotting_stuff <- readRDS(file="All_Plotting_stuff_Alt.rds")


clustered_rds <- c("CCA1_manual_SC3.rds", "CCA5_manual_SC3.rds", "HCC6_manual_SC3.rds", "HCC23_manual_SC3.rds", "HCC10_manual_SC3.rds", "HCC24_manual_SC3.rds", "D3DM_manual_SC3.rds", "D3EM_manual_SC3.rds", "D9DM_manual_SC3.rds", "D9EM_manual_SC3.rds");
names(clustered_rds) <- c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24", "D3DM", "D3EM", "D9DM", "D9EM");

for (i in 1:6) {

#i = 1
#i = 5
#i = 6
palette <- V(plotting_stuff$nn_graph[[i]]$graph)$col
this_SCE <- readRDS(clustered_rds[i])
keep_cells <- plotting_stuff$cell_colors[[i]] %in% palette
this_SCE <- this_SCE[,keep_cells]
tmpRef <-RefSCE[rowData(RefSCE)$feature_symbol %in% rowData(this_SCE)$Symbol,]
this_SCE <- this_SCE[match(rowData(tmpRef)$feature_symbol, rowData(this_SCE)$Symbol),]
set.seed(817)
out <- scmapCluster(this_SCE, list(ref=metadata(tmpRef)$scmap_cluster_index), threshold=-1)

# Visualize Results

require("RColorBrewer")
colours <- c("white", brewer.pal(6,"Oranges")[1:3], brewer.pal(6,"Reds")[4:6], "forestgreen", brewer.pal(3,"Blues")[2:3], "darkorchid", "grey65")
names(colours) <- c(
		"ipsc",
		"definitive endoderm",
		"hepatic endoderm", 
		"immature hepatoblast", 
		"fetal hepatocytes", 
		"mature hepatocyte", 
		"adult hepatocytes",
		"Chol",
		"mesenchymal stem cell",
		"stellate", "endothelial", "unassigned")

stuff <- table(out$scmap_cluster_labs, this_SCE$Manual_Clusters)
stuff_col <- colours[names(colours) %in% rownames(stuff)]
stuff <- stuff[match(names(stuff_col), rownames(stuff)),]

leg_names <- names(stuff_col)
leg_names[leg_names=="definitive endoderm"] <- "d-endo"
leg_names[leg_names=="hepatic endoderm"] <- "h-endo"
leg_names[leg_names=="immature hepatoblast"] <- "i-hepato"
leg_names[leg_names=="mature hepatocyte"] <- "m-hepato"
leg_names[leg_names=="mesenchymal stem cell"] <- "msc"
leg_names[leg_names=="adult hepatocytes"] <- "a-hepato"
leg_names[leg_names=="fetal hepatocytes"] <- "f-hepato"
leg_names[leg_names=="endothelial"] <- "endoth"

#png(paste(names(clustered_rds)[i], "scmap_barplot_Alt_munich.png", sep="_"), width=4, height=4, units="in", res=300)
par(mar=c(4,4,1,1))
barplot(stuff, col=stuff_col)
#dev.off();
if (!is.null(nrow(stuff))) {
thingsums <- colSums(stuff)
} else {
thingsums<- stuff
stuff <- table(out$scmap_cluster_labs, this_SCE$Manual_Clusters)
}
if (thingsums[1] > thingsums[ncol(stuff)]) {
#legend("topright", leg_names, fill=stuff_col, bty="n")
} else {
#legend("topleft", leg_names, fill=stuff_col, bty="n")
}

table(out$scmap_cluster_labs, this_SCE$Proliferating)
table(out$scmap_cluster_labs, this_SCE$Manual_Clusters)

#png(paste(names(clustered_rds)[i], "scmap_scatter_Alt_munich.png", sep="_"), width=6, height=6, units="in", res=300)
	par(mar=c(4,4,2,1))
	xes <- plotting_stuff$cell_coords[[i]][keep_cells,1]
	yes <- plotting_stuff$cell_coords[[i]][keep_cells,2]
	plot(xes, yes, 
		xlab="Dim 1", ylab="Dim 2", 
		bg=stuff_col[factor(out$scmap_cluster_labs, 
		levels=names(stuff_col))], pch=21, cex=2)
	#legend("topleft", leg_names, fill=stuff_col, bty="n")
#dev.off()	

out$scmap_cell_Cols <- stuff_col[factor(out$scmap_cluster_labs, levels=names(stuff_col))]
names(out$scmap_cell_Cols) <- colnames(this_SCE)[keep_cells];
saveRDS(out, paste(names(clustered_rds)[i], "scmap_output.rds", sep="_"))
}


#png("scmap_Alt_munich_legend.png", width=3, height=4, units="in", res=300)
plot(1,1, xaxt="n", yaxt="n", col="white", xlab="", ylab="", main="", bty="none")
legend("left", fill=stuff_col, leg_names)
#dev.off()


