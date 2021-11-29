# Human Liver - published
# Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations - MacParland et al. 2018
#GSE115469
#https://github.com/BaderLab/singleLiverCells
#https://github.com/BaderLab/HumanLiver

#devtools::install_github("BaderLab/HumanLiver")
require("HumanLiver")
viewHumanLiver()
require("fpc")
out <- dbscan(ds$dr_viz, eps=2)
summary(factor(out$cluster))
cols=c("grey50","red","pink","green","blue","orange","purple","navy","orchid","forestgreen","firebrick","goldenrod1","cornflowerblue","salmon","purple4", "turquoise", "magenta", "darkgreen", "darkorange")
plot(ds$dr_viz[,1], ds$dr_viz[,2], col=cols[factor(out$cluster)])
labels <- out$cluster
expr_mat <- ds$nge
plot(ds$dr_viz[,1], ds$dr_viz[,2], col=cols[factor(out$cluster)])
legend("bottomleft", levels(factor(out$cluster)), pch=16, col=cols, bty="n")
labels <- as.character(labels)
labels[labels==0] <- "outliers"
labels[labels==1] <- "LSECs"
labels[labels==2] <- "Cholangiocytes"
labels[labels==3] <- "Macrophages"
labels[labels==4] <- "abTcells"
labels[labels==5] <- "NKcells"
labels[labels==6] <- "gdTcells"
labels[labels==7] <- "Hepatocytes5"
labels[labels==8] <- "Endothelial"
labels[labels==9] <- "gdTcells2"
labels[labels==10] <- "Hepatocytes"
labels[labels==11] <- "Bcells"
labels[labels==12] <- "Stellate"
labels[labels==13] <- "PlasmaCells"
labels[labels==14] <- "Erythroid"
labels[labels==15] <- "Hepatocyptes3"
labels[labels==16] <- "unknown"
labels[labels==17] <- "unknown2"
labels[labels==18] <- "unknown3"
require("SingleCellExperiment")
cc <- ds$md$Phase
tsne <- ds$dr_viz
cdata <- data.frame(cbind(tsne, labels, cc))
colnames(cdata) <- c("tSNE_1", "tSNE_2", "cell_type1", "CC_Phase")
rdata <- data.frame(rownames(expr_mat))
colnames(rdata) <- c("feature_symbol")
obj <- SingleCellExperiment(assay=list(logcounts=expr_mat), colData=cdata, rowData=rdata)
saveRDS(obj, file="HLiver_MacParland.rds")
