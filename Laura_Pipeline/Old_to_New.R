# CCA1, HCC6, HCC10 old vs new clusters

require("scater")
new_dir <- "/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output"
old_dir <- "/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/"


match_SCEs <- function(SCE1, SCE2) {
	cols <- colnames(SCE1)[colnames(SCE1) %in% colnames(SCE2)]
	SCE1 <- SCE1[,match(cols, colnames(SCE1))]
	SCE2 <- SCE2[,match(cols, colnames(SCE2))]
	rows <- rownames(SCE1)[rownames(SCE1) %in% rownames(SCE2)]
	SCE1 <- SCE1[match(rows, rownames(SCE1)),]
	SCE2 <- SCE2[match(rows, rownames(SCE2)),]
	return(list(SCE1, SCE2))
}


CCA1_new <- readRDS(paste(new_dir, "CCA1_SC3.rds", sep="/"))
CCA1_old <- readRDS(paste(old_dir, "CCA1_SC3.rds", sep="/"))

fixed <- match_SCEs(CCA1_old, CCA1_new)
Assignment <- table(fixed[[1]]$sc3_6_clusters, fixed[[2]]$clusters_clean)
rownames(Assignment) <- paste("old", rownames(Assignment))
colnames(Assignment) <- paste("new", colnames(Assignment))

write.table(Assignment, file="CCA1_SC3_old_vs_new.txt", col.names=T, row.names=T)

HCC6_new <- readRDS(paste(new_dir, "HCC6_SC3.rds", sep="/"))
HCC6_old <- readRDS(paste(old_dir, "HCC6_SC3.rds", sep="/"))

fixed <- match_SCEs(HCC6_old, HCC6_new)
Assignment <- table(fixed[[1]]$sc3_6_clusters, fixed[[2]]$clusters_clean)
rownames(Assignment) <- paste("old", rownames(Assignment))
colnames(Assignment) <- paste("new", colnames(Assignment))

write.table(Assignment, file="HCC6_SC3_old_vs_new.txt", col.names=TRUE, row.names=TRUE)

HCC10_new <- readRDS(paste(new_dir, "HCC10_SC3.rds", sep="/"))
HCC10_old <- readRDS(paste(old_dir, "HCC10_SC3.rds", sep="/"))

fixed <- match_SCEs(HCC10_old, HCC10_new)
Assignment <- table(fixed[[1]]$sc3_6_clusters, fixed[[2]]$clusters_clean)
rownames(Assignment) <- paste("old", rownames(Assignment))
colnames(Assignment) <- paste("new", colnames(Assignment))

write.table(Assignment, file="HCC10_SC3_old_vs_new.txt", col.names=TRUE, row.names=TRUE)

