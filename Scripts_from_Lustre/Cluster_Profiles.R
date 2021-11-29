require("M3Drop")
require("scater")
require("matrixStats")
require("RColorBrewer")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/DiffExpr/DE_functions.R")
source("Cluster_Profiles_Functions.R")

CCA1 <- readRDS("CCA1_SC3.rds")
HCC6 <- readRDS("HCC6_SC3.rds")
HCC10 <- readRDS("HCC10_SC3.rds")

markers <- readRDS("Laura_SC3_k6_Markers.rds")

exprs(CCA1) <- get_log_norm(CCA1)
exprs(HCC6) <- get_log_norm(HCC6)
exprs(HCC10) <- get_log_norm(HCC10)

## KW test
CCA1_kw <- apply(exprs(CCA1), 1, function(x) {kruskal.test(x ~ pData(CCA1)$sc3_6_clusters)$p.value})
HCC6_kw <- apply(exprs(HCC6), 1, function(x) {kruskal.test(x ~ pData(HCC6)$sc3_6_clusters)$p.value})
HCC10_kw <- apply(exprs(HCC10), 1, function(x) {kruskal.test(x ~ pData(HCC10)$sc3_6_clusters)$p.value})

png("CCA1_Cluster_KWDE_heatmap.png", width=8, height=8, units="in", res=300)
fancy_heatmap(CCA1, CCA1_kw, markers$CCA1)
dev.off()
png("HCC6_Cluster_KWDE_heatmap.png", width=8, height=8, units="in", res=300)
fancy_heatmap(HCC6, HCC6_kw, markers$HCC6)
dev.off()
png("HCC10_Cluster_KWDE_heatmap.png", width=8, height=8, units="in", res=300)
fancy_heatmap(HCC10, HCC10_kw, markers$HCC10)
dev.off()

CCA1_genes <- c(names(CCA1_kw)[p.adjust(CCA1_kw, method="bon") < 0.05],markers$CCA1$Gene[markers$CCA1$pvalue < 3.266693e-07])
HCC6_genes <- c(names(HCC6_kw)[p.adjust(HCC6_kw, method="bon") < 0.05],markers$HCC6$Gene[markers$HCC6$pvalue < 3.266693e-07])
HCC10_genes <- c(names(HCC10_kw)[p.adjust(HCC10_kw, method="bon") < 0.05],markers$HCC10$Gene[markers$HCC10$pvalue < 3.266693e-07])


###### Removing Line Specific Effects #######

### Combine datasets ###
require("scater")
CCA1_pData <- pData(CCA1)[,1:4]
HCC6_pData <- pData(HCC6)[,1:4]
HCC10_pData <- pData(HCC10)[,1:4]

CCA1_pData$Final_Clusters <- paste("CCA1", pData(CCA1)$sc3_6_clusters, sep="_")
HCC6_pData$Final_Clusters <- paste("HCC6", pData(HCC6)$sc3_6_clusters, sep="_")
HCC10_pData$Final_Clusters <- paste("HCC10", pData(HCC10)$sc3_6_clusters, sep="_")

Combined_counts <- cbind(counts(CCA1), counts(HCC6), counts(HCC10))
Combined_pData <- rbind(CCA1_pData, HCC6_pData, HCC10_pData)
Combined_fData <- fData(CCA1)[,1:2]

gd <- new("AnnotatedDataFrame", data=Combined_fData)
pd <- new("AnnotatedDataFrame", data=Combined_pData)
CombinedSCE <- newSCESet(countData = Combined_counts, phenoData = pd, featureData = gd)

CombinedSCE <- calculateQCMetrics(CombinedSCE)
batch <- factor(pData(CombinedSCE)$Type)

# Match clusters
CCA1_HCC6 <- match_clusters(CCA1, HCC6)
CCA1_HCC10 <- match_clusters(CCA1, HCC10)
HCC6_HCC10 <- match_clusters(HCC6, HCC10)

recip_table <- cbind(c( paste("CCA1", CCA1_HCC6$recip[,1], sep="_"),
			paste("CCA1", CCA1_HCC10$recip[,1], sep="_"),
			paste("HCC6", HCC6_HCC10$recip[,1], sep="_")), 
		     c( paste("HCC6", CCA1_HCC6$recip[,2], sep="_"),
			paste("HCC10", CCA1_HCC10$recip[,2], sep="_"),
			paste("HCC10", HCC6_HCC10$recip[,2], sep="_")))

new_clusters <- pData(CombinedSCE)$Final_Clusters

new_clusters[new_clusters == "HCC10_1"] = "Stem"
new_clusters[new_clusters == "HCC6_4"] = "Stem"
new_clusters[new_clusters == "CCA1_4"] = "Stem"
new_clusters[new_clusters == "CCA1_6"] = "Chol"
new_clusters[new_clusters == "HCC6_5"] = "Chol"
new_clusters[new_clusters == "HCC10_3"] = "Chol"
new_clusters[new_clusters == "CCA1_3"] = "Prog"
new_clusters[new_clusters == "HCC6_2"] = "Prog"
new_clusters[new_clusters == "HCC10_6"] = "Prog"
new_clusters[new_clusters == "HCC10_5"] = "Prog"
new_clusters[new_clusters == "HCC6_1"] = "Prog2"
new_clusters[new_clusters == "HCC10_4"] = "Prog2"

### Combat ###
require("sva")
#CombinedSCE <- newSCESet(countData = Combined_counts, phenoData = pd, featureData = gd)
#Combined_data <- get_log_norm(CombinedSCE)
Combined_data <- exprs(CombinedSCE)
mod = model.matrix(~as.factor(new_clusters))
mod0 = model.matrix(~1, data=pData(CombinedSCE))
# latent variables
#n.sv = num.sv(Combined_data,mod,method="leek")
#svobj = sva(Combined_data,mod,mod0,n.sv=n.sv)
# DE test - accounting for latent variables
#modSv = cbind(mod,svobj$sv)
#mod0Sv = cbind(mod0,svobj$sv)
#pValues = f.pvalue(Combined_data,modSv,mod0Sv)
#qValues = p.adjust(pValues,method="BH")
combat_data = ComBat(dat=Combined_data, batch=factor(batch), mod=mod0, par.prior=TRUE, prior.plots=FALSE)

png("ComBat_Profiles_Heatmap.png", width=6, height=6, units="in", res=300)
cluster_relative_heatmap( combat_data, pData(CombinedSCE)$Final_Clusters, npermute=100)
dev.off()

### Tung et al/Blischak model - this is very very slow ###
require("limma")
#require("edgeR")
#require("ggplot2")
#theme_set(theme_bw(base_size = 12))
#source("functions.R")
require("Humanzee")


##design <- model.matrix(~ 1 + Final_Clusters, data = pData(CombinedSCE)) # This didn't work - maybe because of confounded clusters?
#design <- model.matrix(~1 , data = pData(CombinedSCE))
#block <- batch
#dup_corrs <- duplicateCorrelation(Combined_data, design = design, block = block)
#gls_fit <- Humanzee::ruv_mixed_model(Combined_data, ndups = 1, design = design, block = block, correlation = dup_corrs$cons)

#residual_expr <- t( design %*% t(gls_fit$coef) ) + gls_fit$resid
#colnames(molecules_final) <- colnames(molecules_cpm_trans)

#cluster_relative_heatmap(residual_expr, pData(CombinedSCE)$Final_Clusters, DEonly=FALSE)

### RUVs ###
# Runs fine does not correct sufficiently!
require("RUVSeq")

scIdx <- matrix(-1, ncol = max(table(batch)), nrow = 3)
batch_labs = levels(batch)
for( b in 1:length(batch_labs)) {
	tmp <- which(batch == batch_labs[b])
	scIdx[b, 1:length(tmp)] <- tmp
}
cIdx <- rownames(Combined_data)
ruvs <- RUVs(Combined_data, cIdx, k = 10, scIdx = scIdx, isLog = TRUE)
png("RUVs_Profiles_Heatmap.png", width=6, height=6, units="in", res=300)
cluster_relative_heatmap(ruvs$normalizedCounts, pData(CombinedSCE)$Final_Clusters, npermute=100)
dev.off()

### Relative Expression Profiles ###

# Naive
batch_labs = levels(batch)
corrected <- Combined_data
for( b in 1:length(batch_labs)) {
	corrected[,batch==batch_labs[b]] <- corrected[,batch==batch_labs[b]]-rowMeans(corrected[,batch==batch_labs[b]])
}

png("NaiveMean_Profiles_Heatmap.png", width=6, height=6, units="in", res=300)
cluster_relative_heatmap(corrected, pData(CombinedSCE)$Final_Clusters, npermute=100)
dev.off()

# Known Clusters
corrected2 <- Combined_data
for( b in 1:length(batch_labs)) {
#	row_means <- t(aggregate(t(corrected2[,batch==batch_labs[b]]), by=list(pData(CombinedSCE)$Final_Clusters[batch==batch_labs[b]]), mean))
	cluster_means <- my_row_mean_aggregate(corrected2[,batch==batch_labs[b]], pData(CombinedSCE)$Final_Clusters[batch==batch_labs[b]])

        corrected2[,batch==batch_labs[b]] <- corrected2[,batch==batch_labs[b]]-rowMeans(cluster_means)
}

png("WeightedMean_Profiles_Heatmap.png", width=6, height=6, units="in", res=300)
cluster_relative_heatmap(corrected2, pData(CombinedSCE)$Final_Clusters, npermute=100)
dev.off()

my_profiles <- cluster_means <- my_row_mean_aggregate(corrected2,  pData(CombinedSCE)$Final_Clusters)

write.table(my_profiles, file="WeightedMean_Profiles.csv", sep=",")
