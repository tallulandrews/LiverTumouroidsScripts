# array data for hiPSC-derived cholangiocyte-like cells 
# Sampaziotis et al. (2015). Cholangiocytes derived from human induced pluripotent stem cells for disease modeling and drug validation. Nature Biotechnology. 33:845-852
# PMID: 26167629
# ArrayExpress: E-MTAB-2965
#require("lumi") # I think this is already QCed & quantile normalized

probe_lvls <- read.delim("fotis_sampaziotis_2014-09-15_groutset_quantile_nobkgd_ArrayExpress_DataFile.txt", header=TRUE, stringsAsFactors=FALSE)

probe_lvls <- probe_lvls[-c(1:6),]
colnames(probe_lvls) <- probe_lvls[1,]
probe_lvls <- probe_lvls[-1,]

sample_anno <- read.delim("E-MTAB-2965.sdrf.txt", header=TRUE, stringsAsFactors=FALSE)
probe_anno <- read.delim("A-GEOD-10558.adf.txt", header=TRUE, stringsAsFactors=FALSE)


library("illuminaHumanv4.db")
thing <- mget(x=probe_lvls[,1], envir= illuminaHumanv4SYMBOL)
probe_lvls$Gene <- unlist(thing)

require("limma")
limma_obj <- read.ilmn("fotis_sampaziotis_2014-09-15_groutset_quantile_nobkgd_ArrayExpress_DataFile.txt")
conditions <- sample_anno[match(colnames(limma_obj), sample_anno[,1]), c("Source.Name", "Description")]
conditions <- matrix(unlist(strsplit(conditions[,2], " ")), ncol=2,byrow=TRUE)

probe_filter <- rowSums(limma_obj$other[[1]] < 0.05) >= 3
limma_obj<- limma_obj[probe_filter,] 
cond <- factor(conditions[,1], levels=c("hPSCs", "HBs", "CLCs", "BIL"))
mod <- model.matrix(~cond)

fit <- lmFit(limma_obj$E, mod)

#contrast.matrix <- makeContrasts(condCLCs-condBIL, levels=mod)
#fit_immature_vs_mature_Chol <- contrasts.fit(fit, contrast.matrix)
#fit_immature_vs_mature_Chol <- eBayes(fit_immature_vs_mature_Chol)
#results_imm_mat <- decideTests(fit_immature_vs_mature_Chol)

# List of genes I want: 
#	CLCs & BIL vs hPSCs & HBs = Chol
#	CLCs vs hPSCs & HBs OR BIL vs hPSCs & HBs = Chol
#	CLCs vs BIL & in one of the above = immature/mature Chol

contrast.matrix <- makeContrasts(condBIL-condCLCs, condCLCs-Intercept, condBIL-Intercept, condHBs-Intercept, levels=mod)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)
fit_contrasts <- eBayes(fit_contrasts)
results <- decideTests(fit_contrasts)

CholUp <- results[,1] > 0 & results[,2] > 0 & results[,3] > 0 & results[,4] < 1


cond <- factor(conditions[,1], levels=c("BIL", "HBs", "CLCs", "hPSCs"))
mod <- model.matrix(~cond)
fit <- lmFit(limma_obj$E, mod)
contrast.matrix <- makeContrasts(condhPSCs-Intercept, condhPSCs-condCLCs, condhPSCs-condHBs, levels=mod)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)
fit_contrasts <- eBayes(fit_contrasts)
results <- decideTests(fit_contrasts)


hPSCsUp <- results[,1] > 0 & results[,2] > 0 & results[,3] > 0
averages <- my_row_mean_aggregate(limma_obj$E, factor(conditions[,1]))


library("illuminaHumanv4.db")
symbol <-unlist(mget(x=rownames(averages), envir=illuminaHumanv4SYMBOL))

keep <- !is.na(symbol)
CholUp <- CholUp[keep]
hPSCsUp <- hPSCsUp[keep]
averages <- averages[keep,]
symbol <- symbol[keep]
averages <- as.data.frame(averages)
averages$symbol <- symbol

averages$CholScore <- apply(averages[,1:2], 1, min)/apply(averages[,3:4], 1, max)
averages$StemScore <- averages[,4]/apply(averages[,1:3], 1, max)
averages$isChol <- CholUp
averages$ishPSC <- hPSCsUp
score_threshold <- 2
keep_score <- averages[,6] > score_threshold | averages[,7] > score_threshold
averages<- averages[keep_score,]
Chol_score <- aggregate(averages$CholScore, list(averages$symbol), mean)
PSC_score <- aggregate(averages$StemScore, list(averages$symbol), mean)
PSC_score <- PSC_score[PSC_score[,2] > 2,]
Chol_score <- Chol_score[Chol_score[,2] > 2,]


OUT <- list(CholUp=Chol_score, hPSCsUp=PSC_score)

saveRDS(OUT, file="Sampaziotis_CholDimension.rds")
