# R3.4
args <- commandArgs(trailingOnly=TRUE) # rds file for data, prefix for output

source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/99_Functions.R")

set.seed(1973)

expr_type <- "lognorm"


require("scater")
require("SingleCellExperiment")
SCE <- readRDS(args[1]);
SCE <- toSingleCellExperiment(SCE);

require("CycleMix")
require("mclust")
source("~/NetworkInferencePipeline/Dropouts/My_R_packages/CycleMix/R/Analysis.R")
source("~/NetworkInferencePipeline/Dropouts/My_R_packages/CycleMix/R/Plotting.R")

CC_genes <- rbind(HGeneSets$Tirosh, HGeneSets$Quiesc)

out <- classifyCells(SCE, CC_genes, expr_name="lognorm", do.scale=FALSE, symbol_column="Symbol", allow.multi=FALSE)

colData(SCE)$CC_state <- factor(out$phase, levels=c("None", "G0", "G1S", "G2M"))
colData(SCE)$Proliferating <- ! (out$phase %in% c("G0", "None"))

# Make Plot
png(paste(args[2], "_ProlifMix.png", sep=""), width=7, height=7, units="in", res=300)
par(mfrow=c(2,2))
for(i in out$fits) {
	plotMixture(i, BIC=FALSE)
}
dev.off()

if (sum(out$phase == "G1S") > 0 & sum(out$phase == "G2M") > 0) { 
	SCE <- regressCycle_partial(SCE, out, expr_name="lognorm", method="phase", phases=c("G1S", "G2M"))
} else {
	assays(SCE)[["norm_exprs"]] <- assays(SCE)[["lognorm"]]
}

# Check Regression
if (sum( colData(SCE)$CC_state == "G1S" ) > 2 &
    sum( colData(SCE)$CC_state == "G2M" ) > 2) {
	de <- apply(assays(SCE)[["norm_exprs"]], 1, function(x) {
		wilcox.test(x[colData(SCE)$CC_state == "G1S"], 
			    x[colData(SCE)$CC_state == "G2M"])$p.value})
	de_orig <- apply(assays(SCE)[["lognorm"]], 1, function(x) {
		wilcox.test(x[colData(SCE)$CC_state == "G1S"], 
			    x[colData(SCE)$CC_state == "G2M"])$p.value})
} else {
	de <- NA
	de_orig <- NA
}

# Save some output
sink(file=paste(args[2], "Prolif_regDE.txt", sep="_"))
summary(colData(SCE)$CC_state)
paste("Raw DE:", sum(p.adjust(de_orig, method="fdr") < 0.05, na.rm=T))
paste("Corrected DE:", sum(p.adjust(de, method="fdr") < 0.05, na.rm=T))
sink();


#names(assays(SCE))[names(assays(SCE)) == "norm_exprs"] <- "noCC_exprs"


saveRDS(SCE, file=paste(args[2], "Prolif_v2.rds", sep="_"))
