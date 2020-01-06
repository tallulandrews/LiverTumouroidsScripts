library("scater")
library("RColorBrewer")
group_cols <- brewer.pal(n=6, "Set2")
map = read.table("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/Hsap_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
        new = as.character(map[match(x, map[,1]),2])
        new[is.na(new)] = as.character(x[is.na(new)])
        new[duplicated(new)] = x[duplicated(new)]
        return(new)
}



#### Epigenetic regulators ####
CCA1 <- readRDS("CCA1_SC3_Prolif.rds")
HCC6 <- readRDS("HCC6_SC3_Prolif.rds")
HCC10 <- readRDS("HCC10_SC3_Prolif.rds")

get_log_norm <- function(SCE) {
        dat <- counts(SCE);
        fac <- colSums(dat);
        norm <- t(t(dat)/fac*median(fac))
        return(log(norm +1)/log(2));
}

exprs(CCA1) <- get_log_norm(CCA1)
exprs(HCC6) <- get_log_norm(HCC6)
exprs(HCC10) <- get_log_norm(HCC10)

fData(CCA1)$Symbol = ensg2symbol(rownames(fData(CCA1)))
fData(HCC6)$Symbol = ensg2symbol(rownames(fData(HCC6)))
fData(HCC10)$Symbol = ensg2symbol(rownames(fData(HCC10)))


emt_genes1 <- read.table("/nfs/users/nfs_t/ta6/Data/EMT_Anastassiou_GSEA_Hsap.txt", header=F)
emt_genes2 <- read.table("/nfs/users/nfs_t/ta6/Data/EMT_TGFb_Reactome_Hsap.txt", header=F)


make_PathScore_plot <- function(SCE, path_genes, plot_type=c("PCA","hist","box")) {
	path_genes <- as.character(path_genes);

	require(matrixStats)
	require(mixtools)
	loggednorm <- exprs(SCE)
	detected <- rowSums(loggednorm > 0) > 5
	
	sigma <- rowVars(loggednorm); sigma[sigma == 0] = 1; sigma = sqrt(sigma)
	Zed <- (loggednorm-rowMeans(loggednorm))/sigma

	score <- colSums(Zed[fData(SCE)$Symbol %in% path_genes & detected,]) 
	score2 <- apply(Zed[fData(SCE)$Symbol %in% path_genes & detected,], 2, max) 
	print(sum(fData(SCE)$Symbol %in% path_genes & detected));

	intensity_col <- colorRampPalette(c("white", "black"))(100)
	bins <- c(seq(from=0, to=quantile(score, prob=0.99), length=100), max(score)+1)

	require(TeachingDemos)
	source("/nfs/users/nfs_t/ta6/R-Scripts/Colour_bar.R")
	myplot <- function() {
	
		min=0; max=length(bins)-1; ticks.at=seq(from=0, to=length(bins)-1, by=10); ticks.lab=round(bins[(0:10)*10+1]); title='';
		scale = (length(intensity_col))/(max-min)
		#par(mar=c(0,4,1,0))
	        plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
	        axis(2, at=ticks.at, labels=ticks.lab, las=1)
	        for (i in 1:(length(intensity_col))) {
	             y = (i-1)/scale + min
	             rect(0,y,10,y+1/scale, col=intensity_col[i], border=NA)
	        }
	        mtext(side=1, line=2.5, title, font=2, cex=1.1)
	}

	if (plot_type[1] == "PCA") {
		set.seed(1)
		PCA <- prcomp(loggednorm, scale=TRUE)
		par(mar=c(4,4,1,1))
		par(cex=1)
		plot(PCA$rotation[,1], PCA$rotation[,2], col=intensity_col[cut(score, breaks=bins)], pch=16, xlab=paste("PC1 (",round(PCA$sdev[1]/sum(PCA$sdev)*100, digits=2)," %)",sep=""), ylab=paste("PC2 (",round(PCA$sdev[2]/sum(PCA$sdev)*100, digits=2)," %)",sep=""))
		par(cex=0.75)
		subplot(myplot(),x=grconvertX(c(0.1,0.15), from="npc"), y=grconvertY(c(0,0.4),from="npc"))
		par(cex=1)
		return(list(score=score, assign=intensity_col[cut(score, breaks=bins)], pca=PCA));
	} 
	if (plot_type[1] == "hist") {
		mix <- normalmixEM(score)
		plot(mix, which=2, main2="", xlab2="Score")
		
		p1 <- dnorm(score, mean = mix$mu[1], sd=mix$sigma[1])
	        p2 <- dnorm(score, mean = mix$mu[2], sd=mix$sigma[2])
		m1 = mix$mu[1]
	        m2 = mix$mu[2]

		means_test <- abs(m1-m2) > min(mix$sigma)


	        if (m1 < m2 & means_test) {
	                assign <- p2 > p1
	        } else if (m2 < m1 & means_test) {
	                assign <- p1 > p2
	        } else {
	                assign <- rep(FALSE, times=length(p1))
		}
		print(cbind(aggregate(assign, by=list(pData(SCE)$sc3_6_clusters), sum),summary(factor(pData(SCE)$sc3_6_clusters))))
		return(list(score=score, assign=assign, pca=NULL));


	}
	if (plot_type[1] == "box") {
		boxplot(score~pData(SCE)$sc3_6_clusters, notch=T, col = group_cols, xlab="Group", ylab="Score")
		return(list(score=score, assign=NULL, pca=NULL))
	}
}

CCA1_hout  <- make_PathScore_plot(CCA1,  emt_genes1[,1], "hist")
HCC6_hout  <- make_PathScore_plot(HCC6,  emt_genes1[,1], "hist")
HCC10_hout <- make_PathScore_plot(HCC10, emt_genes1[,1], "hist")

png("CCA1_EMT_Score_PCA_plot.png", width=6, height=6, units="in", res=300)
CCA1_pout  <- make_PathScore_plot(CCA1,  emt_genes1[,1], "PCA")
dev.off()
png("HCC6_EMT_Score_PCA_plot.png", width=6, height=6, units="in", res=300)
HCC6_pout  <- make_PathScore_plot(HCC6,  emt_genes1[,1], "PCA")
dev.off()
png("HCC10_EMT_Score_PCA_plot.png", width=6, height=6, units="in", res=300)
HCC10_pout <- make_PathScore_plot(HCC10, emt_genes1[,1], "PCA")
dev.off()

png("CCA1_EMT_Score_Box_plot.png", width=6, height=6, units="in", res=300)
CCA1_pout  <- make_PathScore_plot(CCA1,  emt_genes1[,1], "box")
dev.off()
png("HCC6_EMT_Score_Box_plot.png", width=6, height=6, units="in", res=300)
HCC6_pout  <- make_PathScore_plot(HCC6,  emt_genes1[,1], "box")
dev.off()
png("HCC10_EMT_Score_Box_plot.png", width=6, height=6, units="in", res=300)
HCC10_pout <- make_PathScore_plot(HCC10, emt_genes1[,1], "box")
dev.off()

# Heatmap
fancy_heatmap <- function(SCE, genes) {
	cosine_D <- function(x,y) {x %*% y / sqrt(x%*%x * y%*%y)}
	cosine_dist <- function(x) {
		nrow = length(x[,1])
		out = outer(1:nrow,1:nrow, FUN = Vectorize( function(i,j) {cosine_D(x[i,],x[j,])} ) )
		out[out==0] <- min(out[out>0]/2) # prevent infinite distances
		return(as.dist(out))
	}

        require("gplots");
	require("scater")
        group_cols <- group_cols

        toplot = genes
	loggednorm <- exprs(SCE)
	rownames(loggednorm) <- fData(SCE)$Symbol
	detected <- rowSums(loggednorm > 0) > 5
	loggednorm <- loggednorm[detected,]

	# Variance-Adjusted bins
	nbins_approx <- 100;
	M <- loggednorm[rownames(loggednorm) %in% toplot,]
	require("matrixStats")
	x = rowMeans(M); y = rowVars(M)
	model = lm(y~poly(x,3))
	blockbreaks = 0:ceiling(max(M))
	blocks = cut(x, breaks=blockbreaks)
	density = as.vector(unlist(by(fitted(model), blocks, mean)))
	density[is.na(density)] = min(density, na.rm=T)
	density <- ceiling(density*(nbins_approx/( sum(density) )))
	bins = unique(unlist(sapply(1:length(density),function(i) { seq(from=i-1, to=i, length=density[i]+1) })))
	


	# Colour sets
        col_col <- group_cols[pData(SCE)$sc3_6_clusters]
#        heatcols <- colorRampPalette(c("black","blue","cyan","green","yellow"))(length(bins)-1)
        heatcols <- colorRampPalette(c("black","#081d58","#253494","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#edf8b1"))(length(bins)-1)

#	if (length(genes) > 100) {
#        heatmap.2(loggednorm[rownames(loggednorm) %in% toplot,], trace="n", ColSideColors=col_col, col=heatcols, breaks=bins, key.title="", key.xlab="Log2 Expression", hclustfun=function(x){hclust(x, method="ward.D")}, distfun=function(x){as.dist(1-cor(t(x)))})
#	} else {
        heatmap.2(loggednorm[rownames(loggednorm) %in% toplot,], trace="n", ColSideColors=col_col, col=heatcols, key.title="", key.xlab="Log2 Expression", hclustfun=function(x){hclust(x, method="ward.D")}, distfun=function(x){as.dist(1-cor(t(x)))})
#	}
#        heatmap.2(M, trace="n", ColSideColors=col_col, col=heatcols, breaks=bins, key.title="", key.xlab="Log2 Expression", hclustfun=function(x){hclust(x, method="average")}, distfun=function(x){as.dist(1/cosine_dist(x))})
	return(rownames(loggednorm)[rownames(loggednorm) %in% toplot])
}

png("CCA1_EMT_Heatmap.png", width=6, height=6, units="in", res=300)
cca1_emt_genes = fancy_heatmap(CCA1, emt_genes1[,1])
dev.off()
png("HCC6_EMT_Heatmap.png", width=6, height=6, units="in", res=300)
hcc6_emt_genes = fancy_heatmap(HCC6, emt_genes1[,1])
dev.off()
png("HCC10_EMT_Heatmap.png", width=6, height=6, units="in", res=300)
hcc10_emt_genes = fancy_heatmap(HCC10, emt_genes1[,1])
dev.off()

L_marker = read.delim("~/Collaborations/LiverOrganoids/Laura_Markers.txt", header=F, sep=" ")

png("CCA1_LMark_Heatmap.png", width=6, height=10, units="in", res=300)
cca1_lmark_genes = fancy_heatmap(CCA1, L_marker[,1])
dev.off()
png("HCC6_LMark_Heatmap.png", width=6, height=10, units="in", res=300)
hcc6_lmark_genes = fancy_heatmap(HCC6, L_marker[,1])
dev.off()
png("HCC10_LMark_Heatmap.png", width=6, height=10, units="in", res=300)
hcc10_lmark_genes = fancy_heatmap(HCC10, L_marker[,1])
dev.off()

