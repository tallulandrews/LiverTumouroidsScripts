require("CellTypeProfiles")
require("RColorBrewer")
get_pseudotime <- function(gene_list, SCE) {
	set.seed(101)
        toplot <- exprs(SCE)[rownames(SCE) %in% gene_list,]
        require ("destiny")
        dm <- DiffusionMap(t(toplot))
        pseudotime = eigenvectors(dm)[,1]
        pseudotime_order = (1:ncol(toplot))[order(pseudotime)]
        if (mean(pseudotime[SCE$Proliferating]) > mean(pseudotime[!SCE$Proliferating])) {
                pseudotime_order <- rev(pseudotime_order);
        }
        return(list(time=pseudotime, order=pseudotime_order))
}
get_pseudotime_2D <- function(gene_list, SCE) {
	set.seed(101)
        toplot <- exprs(SCE)[rownames(SCE) %in% gene_list,]
        require ("destiny")
        dm <- DiffusionMap(t(toplot))
	p1 <- eigenvectors(dm)[,1]
	p2 <-  eigenvectors(dm)[,2]
	diff1 <- mean(p1[SCE$Proliferating]) - mean(p1[!SCE$Proliferating])
	diff2 <- mean(p2[SCE$Proliferating]) - mean(p2[!SCE$Proliferating])
	pseudotime <- -1*sign(diff1)*p1 + -1*sign(diff2)*p2

        pseudotime_order = (1:ncol(toplot))[order(pseudotime)]
        if (mean(pseudotime[SCE$Proliferating]) > mean(pseudotime[!SCE$Proliferating])) {
                pseudotime_order <- rev(pseudotime_order);
        }
        return(list(time=pseudotime, order=pseudotime_order))
}

pseudo_time_heatmap <- function(heat_data, pseudotime_out, column_labs, markers, genes_to_label=rownames(heat_data)) {
        require("CellTypeProfiles")
        require("RColorBrewer")
        heat_data <- heat_data[, pseudotime_out$order]
        columnLab <- column_labs[pseudotime_out$order]
	sig <- markers$q.value < 0.05 & markers$q.value > -1 & markers$AUC > 0.7
	heat_data<-heat_data[sig,]
	markers<-markers[sig,]

        set.seed(666)
#        markers <- complex_markers(as.matrix(heat_data), columnLab)

        markers_assigned <- markers[,-c(1, ncol(markers), ncol(markers)-1)]
	C_order <- aggregate(1:length(columnLab), list(factor(columnLab)), mean)
	refactor_levels <- C_order[order(C_order[,2]),1]
	markers_assigned <- markers_assigned[,match(refactor_levels, colnames(markers_assigned))]
        #markers_assigned[rowSums(markers_assigned) > ncol(markers_assigned)/2,] <- 1-markers_assigned[rowSums(markers_assigned) > ncol(markers_assigned)/2,] # Flip for negative markers
        pattern <- apply(markers_assigned, 1, paste, collapse="")



        ordering_fxn <- function(x) {
                out <- which(x==1)
                out <- c(out, rep(0, times=length(x)-length(out)))

                return(out);
        }

        gene_order <- do.call(order, c(decreasing=FALSE, as.data.frame(t(apply(markers_assigned, 1, ordering_fxn)))))

        heat_data <- heat_data[gene_order,]
        row_groups <- pattern[gene_order]
        fcs <- factor_counts(factor(row_groups))
        row_groups[row_groups %in% names(fcs)[fcs < 5]] = "Other"
        row_groups <- factor(row_groups)

        new_rownames <- rownames(heat_data); new_rownames[!new_rownames %in% genes_to_label] <- ""

        #heat_data <- t(scale(t(heat_data), center=FALSE))
        heat_data <- heat_data/apply(heat_data,1,max)
        heat_colours <- c("#f7fcfd","#e0ecf4","#bfd3e6","#9ebcda","#8c96c6","#8c6bb1","#88419d","#810f7c","#4d004b")
#heat_colours <- c("#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026")

        rowColVec <- c(rainbow(length(levels(row_groups))-1), "grey85")[row_groups]
        colColVec <- brewer.pal("Set2", n=length(unique(columnLab)))[columnLab]
        heatmap.3(heat_data, trace="none", col=heat_colours, Colv=FALSE, Rowv=FALSE, dendrogram="none",
                        RowSideColors=matrix(rowColVec, nrow=1), ColSideColors=matrix(colColVec, ncol=1),
                        labRow = new_rownames, labCol="")
	names(rowColVec) <- rownames(heat_data)
	#return(cbind(row_groups, rowColVec))
	return(rowColVec)

}

require("CellTypeProfiles")
CCA1 <- readRDS("CCA1_SC3_Prolif.rds")
CCA1_markers <- complex_markers(exprs(CCA1), CCA1$sc3_6_clusters)
HCC6 <- readRDS("HCC6_SC3_Prolif.rds")
HCC6_markers <- complex_markers(exprs(HCC6), HCC6$sc3_6_clusters)
HCC10 <- readRDS("HCC10_SC3_Prolif.rds")
HCC10_markers <- complex_markers(exprs(HCC10), HCC10$sc3_6_clusters)

CCA1_pseudo <- get_pseudotime(rownames(CCA1_markers)[CCA1_markers$q.value < 0.05], CCA1)
HCC6_pseudo <- get_pseudotime_2D(rownames(HCC6_markers)[HCC6_markers$q.value < 0.05], HCC6)
HCC10_pseudo <- get_pseudotime(rownames(HCC10_markers)[HCC10_markers$q.value < 0.05], HCC10)

png("CCA1_pseudotime_marker_heatmap.png", width=8, height=8, units="in", res=300)
CCA1_gene_groups <- pseudo_time_heatmap(exprs(CCA1), CCA1_pseudo, CCA1$sc3_6_clusters, CCA1_markers) 
dev.off()
png("HCC6_pseudotime_marker_heatmap.png", width=8, height=8, units="in", res=300)
HCC6_gene_groups <- pseudo_time_heatmap(exprs(HCC6), HCC6_pseudo, HCC6$sc3_6_clusters, HCC6_markers) 
dev.off()
png("HCC10_pseudotime_marker_heatmap.png", width=8, height=8, units="in", res=300)
HCC10_gene_groups <- pseudo_time_heatmap(exprs(HCC10), HCC10_pseudo, HCC10$sc3_6_clusters, HCC10_markers) 
dev.off()

get_gprofiler_terms <- function(gene_groups, markers) {
	require("gProfileR")
	groups <- factor(gene_groups);
	
	get_stuff <- function(lvl) {
		terms <- gprofiler(rownames(markers[rownames(markers) %in% names(gene_groups)[gene_groups == lvl],]))
		terms <- terms[terms$domain %in% c("BP", "MF", "hpa"),]
		on_in <- colnames(markers)[which(markers[rownames(markers) == (names(gene_groups)[gene_groups == lvl])[1] ,] ==1)]
		stuff <- c(paste(on_in, collapse=","), terms$term.name)
		return(stuff);
	}
	out <- lapply(levels(groups), get_stuff)
	return(out)		
}
CCA1_gene_group_char <- get_gprofiler_terms(CCA1_gene_groups, CCA1_markers)
HCC6_gene_group_char <- get_gprofiler_terms(HCC6_gene_groups, HCC6_markers)
HCC10_gene_group_char <- get_gprofiler_terms(HCC10_gene_groups, HCC10_markers)

as.numeric(unique(factor(CCA1_gene_groups)))
tmp <- factor(CCA1_gene_groups); plot(1:length(levels(tmp)), col=levels(tmp), pch=16, cex=5)
