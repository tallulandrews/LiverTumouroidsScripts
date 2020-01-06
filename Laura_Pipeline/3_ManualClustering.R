# R3.4
# After Talking with Laura & Meri Decided:
# 1) Manually picking k seems to be better than algorithmic approach.
# 2) Regressed Cell-cycle looks like the best approach.


args <- commandArgs(trailingOnly=TRUE) # rds file for data, prefix for output, number of cores

source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("~/Collaborations/LiverOrganoids/Laura_Pipeline/99_Functions.R")

# Input = Merged SCE.
# And Set-Up
n.cores = 4;
k_min = 4
expr_type="norm_exprs" # CC-regressed

require("scater")
SCE <- readRDS("D3EM_merged_SC3.rds")
outprefix = "D3EM";
if (class(SCE)[1] == "SCESet") {
	SCE <- scater::toSingleCellExperiment(SCE)
}

set.seed(123)

max_ks <- ceiling(dim(SCE)[2]/25)

require("SC3")
source("~/R-Scripts/Ensembl_Stuff.R")

rowData(SCE)$feature_symbol <- General_Map(rownames(SCE), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")

SCE_orig <- SCE

## Feature Selection ##
SCE <- SCE[rowData(SCE)$feature_symbol != "",]
SCE <- SCE[!is.na(rowData(SCE)$biotype),]
SCE <- SCE[rowData(SCE)$biotype == "protein_coding",]
SCE <- SCE[rowSums(assays(SCE)[["norm_exprs"]]) > 0,]

#thing <- M3DropConvertData(assays(SCE)[["norm_exprs"]], is.log=TRUE, is.counts=FALSE)
#fs <- M3DropFeatureSelection(thing, mt_method="fdr", mt_threshold=0.05, suppress.plot=FALSE)
#fs2 <- rownames(SCE)[order(rowVars(assays(SCE)[["norm_exprs"]]), decreasing=TRUE)]
#fs2 <- fs2[1:500]
#Features <- union(fs$Gene, fs2)
Features <- rownames(SCE);

### Run SCE Clustering ###
assays(SCE)[["logcounts"]] <- assays(SCE)[[expr_type]]
SCE <- SCE[rownames(SCE) %in% Features,]
SCE <- sc3_prepare(SCE)
SCE <- sc3_estimate_k(SCE)
estimated_k <- SCE@metadata$sc3$k_estimation # Does this work?
if( estimated_k > max_ks ) {
	max_ks <- estimated_k
}
set.seed(3810)
SCE <- sc3(SCE, ks=2:max_ks, n_cores=as.numeric(n.cores), biology=FALSE)

print("SC3 - finished")

### Manual Selecting K & Justification ###

op <- par(ask=TRUE)
stats <- list()
scores <- vector(length=max_ks)
sils <- vector(length=max_ks)
for (i in 2:max_ks) {
	stats[[i]] <- my_clustering(SCE, as.character(i), suppress.plot=TRUE)
	sils[i] <- stats[[i]]$overallSil
	scores[i] <- stats[[i]]$overallScore
}
par(ask=FALSE)

#print(force_crash)

set.seed(123)
K <- 4 # hcc10=7, cca1 = 4, hcc6=4, cca5=6, hcc23=3, hcc24=6, d3dm=4, d3em=4, d9dm=6, d9em=3
sils <- sils[sils>0]
scores <- scores[scores>0]
png(paste(outprefix, "manual_clustering_scores.png", sep="_"), width=7, height=7, units="in", res=300)
plot((1:length(sils))+1, sils, ylim=c(0,1), type="b", xlab="K", ylab="Score", main=outprefix)
lines((1:length(sils))+1, scores, lty=2)
points((1:length(sils))+1, scores)
# chosen K
pos <- max(sils[K-1], scores[K-1])
arrows(K,pos+0.1, K, pos+0.01, col="purple", length=0.1, lwd=2)
text(K, pos+0.11, paste("k =",K), pos=3, col="purple", font=2)
# estimated K
pos <- max(sils[estimated_k-1], scores[estimated_k-1])
arrows(estimated_k,pos+0.1, estimated_k, pos+0.01, col="grey50", length=0.1, lwd=1.5)
text(estimated_k, pos+0.11, paste("est.",estimated_k), pos=3, col="grey50")
legend("topright", c("Silhouette", "Tightness"), pch=1, lty=c(1,2), bty="n")
dev.off()

png(paste(outprefix, "manual_k", K, "SC3_consensus.png", sep="_"), width=7, height=7, units="in", res=300)
out <- my_clustering(SCE, as.character(K), suppress.plot=FALSE)
dev.off()
palette <- cluster_col(max(out$Cs))
source("~/R-Scripts/Blank_plot.R")

png(paste(outprefix, "manual_k", K, "legend.png", sep="_"), width=2.5, height=5, units="in", res=300)
blank_plot()
legend("left", as.character(1:max(out$Cs)), fill=palette, bty="n")
dev.off()



Best_Clusters <- stats[[K]]$Cs
nCs <- factor_counts(Best_Clusters)
cell_colours <- palette[Best_Clusters]

## Calculate Markers ##
require("CellTypeProfiles")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/My_R_packages/CellTypeProfiles/R/Markers.R")

#colData(SCE_orig)$Manual_Clusters <- factor(Best_Clusters);
dat_for_markers <- assays(SCE)[[ expr_type ]]
dat_for_markers <- dat_for_markers[,!Best_Clusters %in% names(nCs)[nCs < 10] ]
groups_for_clusters <- Best_Clusters[!Best_Clusters %in% names(nCs)[nCs < 10]]
best_markers <- complex_markers(dat_for_markers, groups_for_clusters)
best_markers$is.Feature <- best_markers$q.value < 0.05 & best_markers$q.value >= 0 & best_markers$AUC > 0.7

write.table(best_markers, file=paste(outprefix, "manual_Markers.csv", sep="_"), sep="\t")

good <- best_markers[best_markers$is.Feature,]
good <- good[order(good$AUC, decreasing=TRUE),]
good$Symbol <- General_Map(rownames(good), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")

assigned <- good[,-c(1, ncol(good), ncol(good)-1, ncol(good)-2)]
a_names <- colnames(assigned)
#assigned <- t(apply(assigned, 1, function(a){
#                    a <- as.numeric(a)
#                    if (mean(a) > 0.5){return (1-a)}
#                    else {return(a)}
#                    }))
colnames(assigned) <- a_names


ons <- rowSums(assigned)
keep <- (ons/ncol(assigned) <= 0.5 | ons <= 2) & !grepl("-", good$Symbol)
assigned <- assigned[keep,]
good <- good[keep,]
classes <- apply(assigned, 1, paste, collapse="")
to_test <- factor_counts(factor(classes));
to_test <- names(to_test)[to_test > 20]

background <- unique(rowData(SCE)$feature_symbol)
richment_table <- vector();
require("gProfileR")
#for( s in to_test ) {
#	gene_set <- unique(good$Symbol[classes == s])
#	exclude <- c(grep("RPS", gene_set), grep("RPL", gene_set), grep("NDUF", gene_set), grep("PSM", gene_set))
#	if (length(exclude) > 1) {
#		gene_set <- gene_set[-exclude]
#	}
#	gene_set <- gene_set[1:min(length(gene_set), 250)]
#	richments <- gprofiler(gene_set, organism="hsapiens", ordered_query=T, significant=T, underrep=F, 
#			min_set_size=15, max_set_size=3000, correction_method="fdr", custom_bg=background, 
#			src_filter=c("GO:BP", "GO:MF", "KEGG", "REAC"), hier_filtering="moderate")
#	richments <- richments[order(richments$p.value),]
#	richments$OnOff <- s;
#	richment_table <- rbind(richment_table, richments[1:20,])
#	head(richments$term.name, 20)
#	key_genes <- as.vector(unlist(strsplit(head(richments$intersection, 20), ",")))
#	key_genes <- factor_counts(key_genes)
#	key_genes <- rev(sort(key_genes))
#	key_genes <- key_genes[key_genes > max(key_genes)/2]
#	gene_set[gene_set %in% names(key_genes)]
#}

marks <- c()
for(i in 1:ncol(assigned)) {
	if (sum(assigned[,i]==1 & rowSums(assigned)==1) > 10) {
		tmp <- assigned[assigned[,i]==1 & rowSums(assigned)==1,]
	} else {
		tmp <- assigned[assigned[,i]==1 & !rownames(assigned) %in% marks ,]
	}
	marks <- c(marks, rownames(tmp)[1:round(40/ncol(assigned))])
}

heat_data <- assays(SCE)[[expr_type]][rownames(SCE) %in% marks,]
heat_data <- t(apply(heat_data, 1, scale))
h <- heatmap.2(heat_data, trace="none", Colv=TRUE, dendrogram="none", key.title="", key.xlab="Z-score",
		hclustfun=function(x){hclust(x, method="ward.D")})
reorder <- order.dendrogram(h$colDendrogram)
tmp1 <- heat_data[,reorder]
tmp2 <- cell_colours[reorder]
reorder <- order(tmp2)
rownames(tmp1) <- General_Map(rownames(tmp1), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")

png(paste(outprefix, "manual_markers.png", sep="_"), width=6, height=6, units="in", res=300)
heatmap.2(tmp1[,reorder], col=rev(c(rev(brewer.pal(5, "Reds")),"white", brewer.pal(5,"Blues"))), 
		breaks=c(-10, -5, -2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2, 5, 10),
		ColSideColors=tmp2[reorder], trace="none", Colv=FALSE, dendrogram="none", key.title="", key.xlab="Z-score",
		hclustfun=function(x){hclust(x, method="ward.D")})
dev.off();

#Phangorn tree?
#library(phangorn)
#pharngorn_data <- as.matrix(t(assigned))
#rownames(pharngorn_data) <- paste("c", rownames(pharngorn_data), sep="")
#pharngorn_data <- phyDat(pharngorn_data, type="USER", levels=c("0", "1"))
#init_tree <- rtree(ncol(assigned), rooted=FALSE, tip.label=names(pharngorn_data))
#fit <- pml(init_tree, pharngorn_data)
#fit <- optim.pml(fit, optNni=TRUE, optRooted=FALSE, subs=c(0)) #weird error that is difficult to trace!



# Save results
colData(SCE_orig)$Manual_Clusters <- Best_Clusters


add_markers <- best_markers
colnames(add_markers) <- paste("best_marker", colnames(best_markers), sep="_");
if (identical(rownames(best_markers), rownames(rowData(SCE_orig))) ) {
	rowData(SCE_orig) <- cbind(rowData(SCE_orig), add_markers)
}

saveRDS(SCE_orig, file=paste(outprefix,"manual_SC3.rds", sep="_"))
