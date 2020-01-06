

my_get_clusters <- function(consensus_mat, k) {
#	consensus_mat <- SCE@sc3$consensus[[as.character(k)]]
#        consensus_mat <- consensus_mat[[1]]
        distances <- as.dist(1-consensus_mat)
        Cs <- cutree(hclust(distances, method="average"), k=as.numeric(k))
        Cs2 <- cutree(hclust(distances, method="average"), h=1-10^-10)
        if (max(Cs2) > max(Cs)) {Cs <-Cs2}
	return(Cs)
}

my_clustering <- function(SCE, k, suppress.plot=TRUE) {
	require("CellTypeProfiles")
	
	
	if (class(SCE)[1]=="SCESet") {
		consensus_mat <- SCE@sc3$consensus[[k]]
	        consensus_mat <- consensus_mat[[1]]
	} else if (class(SCE)[1] =="SingleCellExperiment") {
		consensus_mat <- SCE@metadata$sc3$consensus[[k]]
	        consensus_mat <- consensus_mat[[1]]
	}
	distances <- as.dist(1-consensus_mat)
	Cs <- my_get_clusters(consensus_mat, k);

        #distances <- as.dist(1-consensus_mat)
        #Cs <- cutree(hclust(distances, method="average"), k=as.numeric(k))
        #Cs2 <- cutree(hclust(distances, method="average"), h=1-10^-10)
	#if (max(Cs2) > max(Cs)) {Cs <-Cs2}

	Cns <- CellTypeProfiles::factor_counts(factor(Cs))

	palette <- cluster_col(max(Cs))
	if (!suppress.plot) {
		require("gplots")
		heatmap.2(consensus_mat, distfun=function(x){as.dist(1-x)}, hclustfun=function(x){hclust(x, method="average")},
			trace="none", ColSideColors=palette[Cs])
	}

	require("cluster")
	sil <- cluster::silhouette(Cs, distances)
	#m <- mean(sil[,3])
	#s <- sd(sil[,3])
	#p <- p.adjust(pnorm(abs(sil[,3]-m)/s, lower.tail=FALSE), method="fdr")
	silhouettes <- mean(sil[,3])
	neighbours <- table(sil[,1], sil[,2])/Cns

	sets <- split(seq(length(Cs)), factor(Cs))
        score <- sapply(sets, function(a) {mean(consensus_mat[a,a]) - mean(consensus_mat[a,-a])})
        scores <- sum(score*Cns)/ncol(consensus_mat)

	return(list(sil_nn=neighbours, c_scores=score, overallSil=silhouettes, overallScore=scores, Cs=Cs, cell_colours=palette[Cs]));
}

get_optimal_k <- function(SCE) {
	silhouettes <- vector()
	scores <- vector()

	if (class(SCE)[1] == "SCESet") {
		consensus <- SCE@sc3$consensus
	} else if (class(SCE)[1] == "SingleCellExperiment") {
		consensus <- SCE@metadata$sc3$consensus
	}

	for (i in names(consensus)) {
		clust_name <- paste("sc3", i, "clusters", sep="_");
		stuff <- my_clustering(SCE, i)

		if (class(SCE)[1] == "SCESet") {
			pData(SCE)[,clust_name] <- stuff$Cs;
		} else if (class(SCE)[1] == "SingleCellExperiment") {
			colData(SCE)[,clust_name] <- stuff$Cs;
		}

		silhouettes <- c(silhouettes, stuff$overallSil)
		scores <- c(scores, stuff$overallScore)
	}
	stability <- vector()
	for (i in names(consensus)) {
		stability <- c(stability, calculate_stability(SCE, i)$overall)
	}

	t1 <- silhouettes-min(silhouettes); t1 <- t1/max(t1)
	t2 <- scores-min(scores); t2 <- t2/max(t2)
	t3 <- stability-min(stability); t3 <- t3/max(t3)
	composite_score <- apply(cbind(t1, t2), 1, min);
	composite_score2 <- apply(cbind(t1, t2, t3), 1, min);
	optimal_k <- which(composite_score == max(composite_score))+1 # Select optimal K
	fine_k <- optimal_k
	#if (optimal_k < k_min) {
		#composite_score[optimal_k-1] = 0
		fine_k <- which(composite_score2 == max(composite_score2))+1
	#}
	return(list(optim = optimal_k, fine = fine_k))
}


refine_clusters <- function(SCE, expr_type, clusters, markers, lim.AUC=0.7, lim.size_pct=5, lim.expect_bernoulli=0.05, lim.sil=0.75, lim.markers=0.4) {
	require("CellTypeProfiles")
	if (length(unique(clusters)) < 3) {
		return(list(newCs=clusters, markers=markers))
	}

	#markers <- complex_markers(get_exprs(SCE, expr_type), pData(SCE_orig)$clusters_fine)
	sig <- markers[markers$q.value < 0.05 & markers$q.value > 0,]
	good <- sig[sig$AUC > lim.AUC,]

	markers$is.Feature <- rownames(markers) %in% rownames(good);	

	assigned <- good[,-c(1, ncol(good), ncol(good)-1, ncol(good)-2)]
	a_names <- colnames(assigned)
	assigned <- t(apply(assigned, 1, function(a){
			a <- as.numeric(a)
			if (mean(a) > 0.5){return (1-a)}
			else {return(a)}
			}))
	colnames(assigned) <- a_names

	# Unique vs Shared Markers
	key_markers <- assigned[rowSums(assigned) == 1,]
	shared_markers <- assigned[rowSums(assigned) > 1,]
	# How evenly are they spread?
	expected <- qbinom(lim.expect_bernoulli/ncol(assigned), size=nrow(key_markers), prob=1/ncol(key_markers))
	# How big are the clusters?
	nCs <- factor_counts(clusters)
	min_C_size <- ncol(SCE)*lim.size_pct/100

	keepC <- colnames(assigned)[colSums(key_markers) > expected] # lots of unique markers == real

	# Do the Refining
	allC <- as.character(levels(clusters))
	rawCs <- as.character(clusters)
	newCs <- rawCs
	for(i in allC) {
		new <- i
		if (nCs[as.numeric(i)] == 1) {
			new <- "Outliers"
			newCs[rawCs==i] <- new[1];
			next;
		}
		# Small & cohesive = outliers
		if (min_C_size > nCs[as.numeric(i)] & out$c_scores[as.numeric(i)] >= min(out$c_scores[as.numeric(keepC)])) {
			new <- "Outliers"
		} 
		if ( !(i %in% keepC) ) {
			# if more cohesive than a keptC keep it too
			#if (out$c_scores[as.numeric(i)] > min(out$c_scores[as.numeric(keepC)])) {
	                #        new <- i
	                #} else {
				# Proportion of shared markers that are shared with cluster X
				this_C_markers <- shared_markers[,colnames(shared_markers) == i] == 1
				m_score <- colSums(shared_markers[this_C_markers,])/sum(this_C_markers)
				m_score[colnames(shared_markers) == i] <- 0
				# Proportion of cells where the next closest cluster is cluster X
				sil_score <- out$sil_nn[as.numeric(i),]

				m_closest <- names(m_score)[m_score == max(m_score)]
				sil_closest <- names(sil_score)[sil_score == max(sil_score)]

				if ( length(intersect(m_closest, sil_closest)) == 1) {
					# Aggreement on closest cluster
					if (max(m_score) > lim.markers & max(sil_score) > lim.sil) {
					# similar enough to be equivalent
						new <- unique(newCs[rawCs==intersect(m_closest, sil_closest)])
					}
				} else {
					if (min_C_size > nCs[as.numeric(i)]) {
						new <- "Outliers"
					}
				}
			#}
		}
		newCs[rawCs==i] <- new[1];
	}
	return(list(newCs=newCs, markers=markers))
}

load_CC <- function(set="all") {
	cellcycle <- read.table("~/Data/Whitfield_CC.txt")
	cellcycle_simple <- as.matrix(cellcycle[cellcycle[,1] != "CC",])
	cellcycle_simple[cellcycle_simple[,1] == "G2",1] = "G2M";
	cellcycle_simple[cellcycle_simple[,1] == "S",1] = "G1S";
	cellcycle_simple = cellcycle_simple[cellcycle_simple[,1] != "MG1",];
	G0_genes <- read.table("~/Data/Reactome_G0.txt", header=F) # Update this with genes from Laura/MiSigDB?
	Quiescence <- read.table("~/Data/Quiescence.txt")
	new_cellcycle <- read.table("~/Collaborations/LiverOrganoids/New_CC_171117.txt", header=FALSE) #Tirsoh2016Nature
	if (set == "all") {
		return(list(Whitfield=cellcycle, Simple=cellcycle_simple, G0=G0_genes, Quiescence=Quiescence, Tirosh=new_cellcycle))
	} 
	if (set == "cycling") {
		return(list(Whitfield=cellcycle, Simple=cellcycle_simple, Tirosh=new_cellcycle))
	}
}


calculate_stability <- function(SCE, k) {
    #hc <- consensus[[as.character(k)]]$hc
    #labs <- reindex_clusters(hc, k)
    
    clust_name <- paste("sc3", k, "clusters", sep="_");
	if (class(SCE)[1] == "SCESet") {
    		labs <- as.numeric(pData(SCE)[,clust_name])
	} else if (class(SCE)[1] == "SingleCellExperiment") {
    		labs <- as.numeric(colData(SCE)[,clust_name])
	}
	if (class(SCE)[1] == "SCESet") {
		consensus <- SCE@sc3$consensus
	} else if (class(SCE)[1] == "SingleCellExperiment") {
		consensus <- SCE@metadata$sc3$consensus
	}
    names(labs) <- colnames(SCE)
    Cns <- CellTypeProfiles::factor_counts(factor(labs))
    
    ks <- as.numeric(names(consensus))
    kRange <- length(ks)
    
    stability <- rep(0, max(labs))
    
    for (i in 1:max(labs)) {
        inds <- names(labs[labs == i])
        # sum over k range
        for (k2 in ks) {
            if (k2 != k) {
    		clust_name2 <- paste("sc3", k2, "clusters", sep="_");
		if (class(SCE)[1] == "SCESet") {
    			labs2 <- as.numeric(pData(SCE)[,clust_name2])
		} else if (class(SCE)[1] == "SingleCellExperiment") {
    			labs2 <- as.numeric(colData(SCE)[,clust_name2])
		}
		names(labs2) <- colnames(SCE)
                clusts <- as.numeric(names(table(labs2[names(labs2) %in% inds])))
                N <- length(clusts)
                # sum over new clusters, taking into account new cells from other clusters
                for (j in clusts) {
                  inds2 <- names(labs2[labs2 == j])
                  s <- length(inds[inds %in% inds2])/length(inds2)/N/N/kRange
                  stability[i] <- stability[i] + s
                }
            }
        }
    }
    return(list(per_cluster=stability, overall=sum(Cns*stability)/sum(Cns)))
}

remove_cell_cycle <- function(SCE) {
	CC_genes <- load_CC("cycling")
	CC_genes <- c(as.character(CC_genes$Whitfield[,2]), as.character(CC_genes$Tirosh[,1]))
	SCE <- SCE[!( fData(SCE)$Symbol %in% CC_genes ),]
	GO_cc <- read.delim("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/GO/hsapiens_80_GO_Annoations_Emsembl.out", sep="\t", header=F)
	GO_cc <- GO_cc[GO_cc[,3] == "cell cycle",1]
	olap <- rownames(SCE) %in% GO_cc
	if (sum(olap) > 0) {
		SCE <- SCE[!( olap ),]
	}
	return(SCE)
}

