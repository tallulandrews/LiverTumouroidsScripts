sils <- sils[sils>0]
scores <- scores[scores>0]
png(paste(outprefix, "manual_clustering_scores.png", sep="_"), width=7, height=7, units="in", res=300)
plot((1:length(sils))+1, sils, ylim=c(0,1), type="b", xlab="K", ylab="Score", main=outprefix)
lines((1:length(sils))+1, scores, lty=2)
points((1:length(sils))+1, scores)
pos <- max(sils[K-1], scores[K-1])
arrows(K,pos+0.1, K, pos+0.01, col="red", length=0.1, lwd=1.5)
text(K, pos+0.11, paste("k =",K), pos=3, col="red")
legend("topright", c("Silhouette", "Tightness"), pch=1, lty=c(1,2), bty="n")
dev.off()

png(paste(outprefix, "manual_k", K, "SC3_consensus.png", sep="_"), width=7, height=7, units="in", res=300)
out <- my_clustering(SCE, as.character(K), suppress.plot=FALSE)
dev.off()
palette <- cluster_col(max(out$Cs))
source("~/R-Scripts/Blank_plot.R")

png(paste(outprefix, "manual_k", K, "legend.png", sep="_"), width=4, height=5, units="in", res=300)
blank_plot()
legend("left", as.character(1:max(out$Cs)), fill=palette, bty="n")
dev.off()



Best_Clusters <- stats[[K]]$Cs

## Calculate Markers ##
require("CellTypeProfiles")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/My_R_packages/CellTypeProfiles/R/Markers.R")

#colData(SCE_orig)$Manual_Clusters <- factor(Best_Clusters);
best_markers <- complex_markers(assays(SCE_orig)[[ expr_type ]], Best_Clusters)
best_markers$is.Feature <- best_markers$q.value < 0.05 & best_markers$q.value >= 0 & best_markers$AUC > 0.7

## Refine Clusters ##
good <- best_markers[best_markers$is.Feature,]
assigned <- good[,-c(1, ncol(good), ncol(good)-1, ncol(good)-2)]
a_names <- colnames(assigned)
assigned <- t(apply(assigned, 1, function(a){
                    a <- as.numeric(a)
                    if (mean(a) > 0.5){return (1-a)}
                    else {return(a)}
                    }))
colnames(assigned) <- a_names

key_markers <- assigned[rowSums(assigned) == 1,]
shared_markers <- assigned[rowSums(assigned) > 1,]

keepC <- colnames(assigned)[colSums(key_markers) > median(colSums(key_markers))]
nCs <- factor_counts(Best_Clusters)

newCs <- Best_Clusters;
nn_sil <- out$sil_nn
nn_mark <- matrix(0,nrow=ncol(assigned), ncol=ncol(assigned))
rownames(nn_mark) <- colnames(assigned);
colnames(nn_mark) <- colnames(assigned);
for (c in colnames(assigned)) {
	this_C_markers <- shared_markers[,colnames(shared_markers) == c] == 1
        m_score <- colSums(shared_markers[this_C_markers,])/sum(this_C_markers)
	m_score[c] <- 0;
	nn_mark[c,] <- m_score
}
# Mannually check merge list does it make sense to merge any? 
# Also graph of relationships!
trimm <- function(x) {
	thresh <- max(x)*0.75
	x[x < thresh] <- 0
	return(x)
}

nn_mark <- t(apply(nn_mark, 1, trimm))
nn_sil <- t(apply(nn_sil, 1, trimm))

g_mark <- igraph::graph_from_adjacency_matrix(nn_mark, mode="undirected", weighted=TRUE)
V(g_mark)$col <- palette[match(V(g_mark), names(palette))]
png(paste(outprefix, "manual_markNet.png", sep="_"), width=6, height=6, units="in", res=300)
plot(g_mark, vertex.col = V(g_mark)$col, edge.width=E(g_mark)$weight*3)
dev.off()

g_sil <- igraph::graph_from_adjacency_matrix(nn_sil, mode="undirected", weighted=TRUE)
V(g_sil)$col <- palette[match(V(g_sil), names(palette))]

png(paste(outprefix, "manual_silNet.png", sep="_"), width=6, height=6, units="in", res=300)
plot(g_sil, vertex.col = V(g_sil)$col, edge.width=E(g_sil)$weight*3)
dev.off()




# Save results
colData(SCE_orig)$Manual_Clusters <- newCs


add_markers <- complex_markers(assays(SCE_orig)[[ expr_type ]], newCs)
add_markers$is.Feature <- add_markers$q.value < 0.05 & add_markers$q.value >= 0 & add_markers$AUC > 0.7
colnames(add_markers) <- paste("best_marker", colnames(best_markers), sep="_");
if (identical(rownames(best_markers), rownames(rowData(SCE_orig))) ) {
	rowData(SCE_orig) <- cbind(rowData(SCE_orig), add_markers)
}

saveRDS(SCE_orig, file=paste(outprefix,"manual_SC3.rds", sep="_"))
