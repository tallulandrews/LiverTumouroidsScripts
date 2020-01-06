require("matrixStats")
require("destiny")
obj <- readRDS("ForRevewPaper_stuff.rds")

norm <- log(obj$norm+1)/log(2)
keep = rowSums(norm) > 0 & rowVars(norm) > 0
norm<-norm[keep,]

dm <- DiffusionMap(t(norm))
dm_dims <- eigenvectors(dm)

require("Rtsne")
set.seed(123)
tsne2 <- Rtsne(t(obj$norm), perplexity=50)

require("scater")
CCA1 <- readRDS("CCA1_SC3_Prolif.rds")
CCA1 <- rownames(pData(CCA1))
HCC6 <- readRDS("HCC6_SC3_Prolif.rds")
HCC6 <- rownames(pData(HCC6))
HCC10 <- readRDS("HCC10_SC3_Prolif.rds")
HCC10 <- rownames(pData(HCC10))

cell_lab <- as.character(colnames(norm))
cell_lab[cell_lab %in% CCA1] = 1
cell_lab[cell_lab %in% HCC6] = 2
cell_lab[cell_lab %in% HCC10] = 3
cell_lab <- as.numeric(cell_lab)

my_colours <- obj$colours[cell_lab]

png("ForReview_DR.png", width=3.5*2,height=3.5*2, units="in", res=300)
par(mfrow=c(2,2))
par(mar=c(3,3.5,2,0.5))
plot(obj$PCA$rotation[,1], obj$PCA$rotation[,2], col=my_colours, 
	xlab="", ylab="", pch=16, xaxt="n", yaxt="n")
title(xlab="Component 1", ylab="Component 2", line=1)
title(main="PCA",line=0.5)
mtext("A", font=2, las=1, at=0.04, line=2, side=2, cex=1.5)

plot(dm_dims[,1], dm_dims[,2], col=my_colours, 
	xlab="", ylab="", pch=16, xaxt="n", yaxt="n")
title(xlab="Dimension 1", ylab="Dimension 2", line=1)
title(main="DM",line=0.5)
mtext("B", font=2, las=2, at=0.85, line=2, side=2, cex=1.5)

plot(obj$tsne$Y[,1], obj$tsne$Y[,2], col=my_colours, 
	xlab="", ylab="", pch=16, xaxt="n", yaxt="n")
title(xlab="Dimension 1", ylab="Dimension 2", line=1)
title(main="tSNE (10)",line=0.5)
mtext("C", font=2, las=2, at=45, line=2, side=2, cex=1.5)

plot(tsne2$Y[,1], tsne2$Y[,2], col=my_colours, 
	xlab="", ylab="", pch=16, xaxt="n", yaxt="n")
title(xlab="Dimension 1", ylab="Dimension 2", line=1)
title(main="tSNE (50)",line=0.5)
mtext("D", font=2, las=2, at=15, line=2, side=2, cex=1.5)

dev.off()
