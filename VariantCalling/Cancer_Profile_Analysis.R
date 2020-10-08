#require("infercnv")
profiled_genes <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/VariantCalling/CancerProfileGenes.csv", stringsAsFactors=FALSE)
source("~/R-Scripts/Ensembl_Stuff.R")
profiled_genes$ensg <- General_Map(profiled_genes[,1], in.org="Hsap", in.name="symbol", out.org="Hsap", out.name="ensg")

load("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/CNVAnalysis/inferCNVRes/run.final.infercnv_obj")
G <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/Global_SCE.rds")


### Plot variants by XXXX ###
All_Calls <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/CNVAnalysis/inferCNV_allcall_CtrlNullBon.rds")
call_type <- "CtrlNullBon";

gene_order <- read.table("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/CNVAnalysis/inferCNV_input_geneorder.txt", header=F)

sex_genes <- gene_order[gene_order[,2] %in% c("chrX", "chrY"),1]

dir <- "/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output"
SCE_files <- c("CCA1_manual_SC3.rds", "HCC24_manual_SC32.rds", "CCA5_manual_SC3.rds", "HCC6_manual_SC3.rds", "HCC10_manual_SC32.rds", "HCC23_manual_SC3.rds")
plot_objs <- c("CCA1_PlottingObj_Alt.rds", "HCC24_PlottingObj_Alt.rds", "CCA5_PlottingObj_Alt.rds", "HCC6_PlottingObj_Alt.rds", "HCC10_PlottingObj_Alt.rds", "HCC23_PlottingObj_Alt.rds")

names(SCE_files) <- c("CCA1", "HCC24", "CCA5", "HCC6", "HCC10", "HCC23")
names(plot_objs) <- c("CCA1", "HCC24", "CCA5", "HCC6", "HCC10", "HCC23")
source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")

# Rename clusters

line_specific_groups <- list(CCA1=c("Progenitor", "Differentiated1", "TICs", "Differentiated2"),
                            CCA5=c("Chol", "Stress", "CSC", "Unk", "Hep"),
                            HCC6=c("Prog1", "Stress", "Prog2", "CSC"),
                            HCC23=c("CSC", "Clot-Hep", "Prog"),
                            #HCC10=c("Prog1", "Hep", "CSC", "iHep", "Stress", "Prog2"),
                            HCC10=c("TICs", "Progenitor", "Differentiated"),
                            #HCC24=c("Clot-Hep", "Prog1", "Prog2", "CSC", "Hep1", "Hep2"),
			    HCC24=c("Differentiated", "Progenitor1", "Progenitor2", "TICs"),
                            D3DM=c("Chol", "Prog", "Stress", "Cycling"),
                            D3EM=c("Stress", "Chol1", "Chol2", "Prog"),
                            D9DM=c("Prog1", "Prog2", "Prog3", "Clot-Hep", "Chol1", "Chol2"),
                            D9EM=c("Chol1", "Cycling", "Chol2")
                        ) # cluster names



del_col="cornflowerblue"
dup_col="brown1"

# May need to mod below here

for (t in names(SCE_files)) {
	print(t)
	call_mat <- All_Calls[[t]]
	call_mat <- call_mat[!rownames(call_mat) %in% as.character(sex_genes),]
	up <- colSums(call_mat > 0);
	down <- colSums(call_mat < 0);

	profiled_call_mat <- call_mat[match(profiled_genes$ensg, rownames(call_mat)),]
	profiled_call_mat <- profiled_call_mat[!is.na(rowSums(profiled_call_mat)) & rowSums(abs(profiled_call_mat)) > 5,]

	

	# by cluster #
	require("scater")
	sce <- readRDS(paste(dir, SCE_files[names(SCE_files) == t], sep="/"))
	groups <- sce$Manual_Clusters

	sce@metadata$C_names <- line_specific_groups[[t]]
        c_names <- sce@metadata$C_names	
	names(c_names) <- 1:length(c_names)

	groups_char <- as.character(groups);
	for (i in 1:length(c_names)) {
		id <- names(c_names)[i]
		groups_char[groups_char == id] <- c_names[i];
	}
	groups <- factor(groups_char);

	exclude <- !(groups %in% c_names) | ! (colnames(sce) %in% colnames(profiled_call_mat))
	groups <- factor(groups_char[!exclude])
	up <- up[!exclude]
	down <- down[!exclude]
	sce <- sce[,!exclude]
	profiled_call_mat <- profiled_call_mat[,match(colnames(sce), colnames(profiled_call_mat))]


	require("CellTypeProfiles")
	require("gplots")
	require("RColorBrewer")
	hot <- brewer.pal(4, "Reds")
	cold <- brewer.pal(4, "Blues")
	heatcol <-  c(rev(cold), "white", hot)

	profiled <- my_row_mean_aggregate(profiled_call_mat, groups)
	profiled_var <- my_row_var_aggregate(profiled_call_mat, groups)
	sqrt_n <- as.vector(sqrt(table(groups)))
	profiled_se <- t(t(profiled_var)/sqrt_n)
	mt_confidence_95 <- qnorm(1-0.025/nrow(profiled))
	profiled_up <- profiled + t(t(profiled_se)*mt_confidence_95)
	profiled_dw <- profiled - t(t(profiled_se)*mt_confidence_95)
	rownames(profiled) <- profiled_genes[match(rownames(profiled), profiled_genes$ensg),1]
	write.table(cbind(profiled, profiled_up, profiled_dw), paste(t, call_type, "profiled_genes.txt", sep="_"), row.names=T, col.names=T)
	# Heatmap
	pdf(paste(t, call_type, "profiled_genes_heatmap.pdf", sep="_"), width=8, height=8)
	par(mar=c(8, 8, 3, 1))
	heatmap.2(profiled, trace="none", scale="none", col=heatcol, symbreaks=TRUE)
	dev.off()
	print("Done")

}

exit();
