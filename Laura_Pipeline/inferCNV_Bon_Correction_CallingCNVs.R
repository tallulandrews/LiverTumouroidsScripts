require("infercnv")
load("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/CNVAnalysis/inferCNVRes/run.final.infercnv_obj")
G <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/Global_SCE.rds")

vals <- infercnv_obj@expr.data
ref_cells <- unlist(infercnv_obj@reference_grouped_cell_indices)
clean_vals <- function(g, ref_cells) {
	g <- g-1;
	t <- quantile(abs(g[ref_cells]), p=0.99)
	g[abs(g) < t] <- 0;
	return(g)
}

cleaned_vals <- apply(vals, 2, clean_vals, ref_cells)

up <- colSums(cleaned_vals > 0.1)
down <- colSums(cleaned_vals < -0.1)
boxplot(up+down ~ G$Proliferating+G$Donor, las=2, col=c("blue", "red"), notch=TRUE)

# Attempt 2 (Sliding window call & clean)

gene_order <- infercnv_obj@gene_order
identical(rownames(vals), rownames(gene_order));


# Sliding Window to call CNVs
# Bonferroni correction across all cells for each window
G <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/Global_SCE.rds")
tumour <- c("CCA1",  "CCA5",  "HCC6",  "HCC23", "HCC10", "HCC24")
ref_cells <- unlist(infercnv_obj@reference_grouped_cell_indices)
All_Calls <- list();
for (t in tumour) {

	ngenes_window = 20
	min_var <- 0.05
	calls <- matrix(0, nrow=nrow(vals), ncol=sum(G$Donor == t));
	colnames(calls) <- colnames(vals)[G$Donor == t]
	rownames(calls) <- rownames(vals);
	sig_threshold <- 0.05;
	for (i in 1:(nrow(vals)-ngenes_window)) {
		res_t<- apply(vals[i:(i+ngenes_window), G$Donor == t], 2, function(x) {c(mean(x), sd(x))})
		res_t[2,res_t[2,] < min_var] <- min_var;

		res_r<- apply(vals[i:(i+ngenes_window), ref_cells], 2, function(x) {c(mean(x), sd(x))})
		res_r[2,res_r[2,] < min_var] <- min_var;

		null_hypothesis <- 1
		std_err <- res_t[2,]

		p.vals <- pnorm(abs(res_t[1,]-null_hypothesis)/std_err, lower.tail=FALSE)
		sig <- p.vals < sig_threshold/ncol(res_t)
		dir <- sign(res_t[1,]-null_hypothesis)
		calls[i:(i+ngenes_window), sig & dir==1] <- 1
		calls[i:(i+ngenes_window), sig & dir==-1] <- -1
	}
	All_Calls[[t]] <- calls;
}
saveRDS(All_Calls, "inferCNV_allcall_SimpleNullBon.rds")
