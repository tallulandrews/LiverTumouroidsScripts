#source("Cluster_Profiles_Functions.R")
get_log_norm <- function(SCE) {
        dat <- counts(SCE);
        fac <- colSums(dat);
        norm <- t(t(dat)/fac*median(fac))
        return(log(norm +1)/log(2));
}


library("CellTypeProfiles")
CCA1 <- readRDS("CCA1_SC3.rds")
HCC6 <- readRDS("HCC6_SC3.rds")
HCC10 <- readRDS("HCC10_SC3.rds")

exprs(CCA1) <- get_log_norm(CCA1)
exprs(HCC6) <- get_log_norm(HCC6)
exprs(HCC10) <- get_log_norm(HCC10)

CCA1_markers <- complex_markers(exprs(CCA1), pData(CCA1)$sc3_6_clusters, n_max=5, strict_only=FALSE)
HCC6_markers <- complex_markers(exprs(HCC6), pData(HCC6)$sc3_6_clusters, n_max=5, strict_only=FALSE)
HCC10_markers <- complex_markers(exprs(HCC10), pData(HCC10)$sc3_6_clusters, n_max=5, strict_only=FALSE)

source("~/R-Scripts/Ensembl_Stuff.R")
CCA1_markers$symbol <- map_symbol_ensg(rownames(CCA1_markers), "Hsap", "ensg")
HCC6_markers$symbol <- map_symbol_ensg(rownames(HCC6_markers), "Hsap", "ensg")
HCC10_markers$symbol <- map_symbol_ensg(rownames(HCC10_markers), "Hsap", "ensg")

obj <- list(CCA1=CCA1_markers, HCC6=HCC6_markers, HCC10=HCC10_markers)
saveRDS(obj, "Laura_Complex_Markers.rds")

# Strict only
CCA1_markers <- complex_markers(exprs(CCA1), pData(CCA1)$sc3_6_clusters, n_max=5, strict_only=TRUE)
HCC6_markers <- complex_markers(exprs(HCC6), pData(HCC6)$sc3_6_clusters, n_max=5, strict_only=TRUE)
HCC10_markers <- complex_markers(exprs(HCC10), pData(HCC10)$sc3_6_clusters, n_max=5, strict_only=TRUE)

source("~/R-Scripts/Ensembl_Stuff.R")
CCA1_markers$symbol <- map_symbol_ensg(rownames(CCA1_markers), "Hsap", "ensg")
HCC6_markers$symbol <- map_symbol_ensg(rownames(HCC6_markers), "Hsap", "ensg")
HCC10_markers$symbol <- map_symbol_ensg(rownames(HCC10_markers), "Hsap", "ensg")

obj <- list(CCA1=CCA1_markers, HCC6=HCC6_markers, HCC10=HCC10_markers)
saveRDS(obj, "Laura_Strict_Complex_Markers.rds")
