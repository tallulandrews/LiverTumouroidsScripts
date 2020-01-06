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

		null_hypothesis <- mean(res_r[1,])
		std_err <- max(res_t[2,], sd(res_r[1,]))
		#std_err <- res_t[2,]

		p.vals <- pnorm(abs(res_t[1,]-null_hypothesis)/std_err, lower.tail=FALSE)
		sig <- p.vals < sig_threshold/ncol(res_t)
		dir <- sign(res_t[1,]-null_hypothesis)
		calls[i:(i+ngenes_window), sig & dir==1] <- 1
		calls[i:(i+ngenes_window), sig & dir==-1] <- -1
	}
	All_Calls[[t]] <- calls;
}
saveRDS(All_Calls, "inferCNV_allcall_CtrlNullBon.rds")

### Plot variants by XXXX ###
All_Calls <- readRDS("inferCNV_allcall_CtrlNullBon.rds")
call_type <- "CtrlNullBon";

gene_order <- read.table("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/CNVAnalysis/inferCNV_input_geneorder.txt", header=F)

sex_genes <- gene_order[gene_order[,2] %in% c("chrX", "chrY"),1]

dir <- "/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output"
SCE_files <- c("CCA1_manual_SC3.rds", "HCC24_manual_SC3.rds", "CCA5_manual_SC3.rds", "HCC6_manual_SC3.rds", "HCC10_manual_SC3.rds", "HCC23_manual_SC3.rds")
names(SCE_files) <- c("CCA1", "HCC24", "CCA5", "HCC6", "HCC10", "HCC23")
source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")


del_col="cornflowerblue"
dup_col="brown1"

for (t in names(SCE_files)) {
	call_mat <- All_Calls[[t]]
	call_mat[!rownames(call_mat) %in% as.character(sex_genes),]
	up <- colSums(call_mat > 0);
	down <- colSums(call_mat < 0);

	# by cluster #
	require("scater")
	sce <- readRDS(paste(dir, SCE_files[names(SCE_files) == t], sep="/"))
	groups <- sce$Manual_Clusters

	png(paste(t, call_type, "CNV_by_celltype.png", sep="_"), width=9, height=6, units="in", res=300)
	par(mfrow=c(1,2))
	tot <- data.frame(nCNV=c(up, down),
			 dir=c(rep("up", length(up)), rep("down", length(down))), 
			 group=c(groups, groups), 
			cc=c(as.character(sce$CC_state_new), as.character(sce$CC_state_new)));
	outliers <- pnorm(tot[,1], mean=mean(tot[,1]), sd=sd(tot[,1]), lower.tail=FALSE) < 10^-5
	tot <- tot[!outliers,]
	loc <- 1:(length(unique(groups))*3);
	loc <- loc[loc %% 3 != 0]
	loc_names <- rep(levels(factor(groups)), each=2)
	loc_names[seq(from=1, to=length(loc_names), by=2)] <- "";

	group_colours <- cluster_col(max(sce$Manual_Clusters))
	boxplot(tot[,1]~tot[,2]+tot[,3], col=c(del_col, dup_col), notch=TRUE, at=loc, names=loc_names, xaxt="n", ylab="Total Variant Size", xlab="Cluster");
	axis(1, at=(loc[seq(from=1, to=length(loc), by=2)]+loc[seq(from=2, to=length(loc), by=2)])/2,
		levels(factor(groups)))
	legend("top", fill=c(del_col, dup_col), bty="n", c("Del", "Dup"), horiz=TRUE)
	title(main=t)

	# by proliferation #
	tot[,4] <- as.character(tot[,4])
	tot[tot[,4]=="None",4] <- "G0";
	loc <- 1:(length(unique(tot[,4]))*3)
	loc <- loc[loc %% 3 != 0]
	loc_names <- rep(levels(factor(tot[,4])), each=2)
	loc_names[seq(from=1, to=length(loc_names), by=2)] <- "";
	tot[,4] <- factor(tot[,4])
	
	boxplot(tot[,1]~tot[,2]+tot[,4], col=c(del_col, dup_col), at=loc, names=loc_names, notch=T, xaxt="n", ylab="Total Variant Size", xlab="Cell Cycle Stage");
	axis(1, at=(loc[seq(from=1, to=length(loc), by=2)]+loc[seq(from=2, to=length(loc), by=2)])/2,
		levels(factor(tot[,4])))
	legend("top", fill=c(del_col, dup_col), bty="n", c("Del", "Dup"), horiz=TRUE)
	title(main=t)
	dev.off()

	png(paste(t, call_type, "CNV_by_cellcycle_line.png", sep="_"), width=6, height=6, units="in", res=300)
	# by proliferation - lines #
	l_up <- aggregate(tot[tot[,2]=="up",1], by=list(tot[tot[,2]=="up",4]), mean)
	l2_up <- aggregate(tot[tot[,2]=="up",1], by=list(tot[tot[,2]=="up",4]), sd)
	l_down <- aggregate(tot[tot[,2]=="down",1], by=list(tot[tot[,2]=="down",4]), mean)
	l2_down <- aggregate(tot[tot[,2]=="down",1], by=list(tot[tot[,2]=="down",4]), sd)
	l2_up[,2] <- l2_up[,2]/sqrt(sum(tot[,2]=="up"))
	l2_down[,2] <- l2_down[,2]/sqrt(sum(tot[,2]=="down"))

	plot(1:3, l_up[,2], type="both", col=dup_col, xaxt="n", xlab="", ylab="Total Variant Size",
	ylim=c(min(l_up[,2]-2*l2_up[,2], l_down[,2]-2*l2_down[,2]), 
		max(l_up[,2]+2*l2_up[,2], l_down[,2]+2*l2_down[,2])))
	lines(1:3, l_up[,2]-2*l2_up[,2], lty=2, col=dup_col)
	lines(1:3, l_up[,2]+2*l2_up[,2], lty=2, col=dup_col)
	lines(1:3, l_down[,2]+2*l2_down[,2], lty=2, col=del_col)
	lines(1:3, l_down[,2]-2*l2_down[,2], lty=2, col=del_col)
	lines(1:3, l_down[,2], type="both", col=del_col)
	dev.off()

}

### Karyoplot ###
gene_order <- read.table("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/CNVAnalysis/inferCNV_input_geneorder.txt", header=F)

for (t in names(All_Calls)) {
	mat <- All_Calls[[t]];
	if (is.null(rownames(mat))) {rownames(mat) <- rownames(infercnv_obj@expr.data)}
	mat1 <- mat[,!G[,G$Donor==t]$Proliferating]
	dup_freq <- rowMeans(mat1 > 0);
	del_freq <- rowMeans(mat1 < 0);

	mat2 <- mat[,G[,G$Donor==t]$Proliferating]
	dup_freq_prolif <- rowMeans(mat2 > 0);
	del_freq_prolif <- rowMeans(mat2 < 0);

	g_order <- gene_order[gene_order[,1] %in% rownames(mat),]
	reorder <- match(g_order[,1], rownames(mat))

	stats <- data.frame(chr=g_order[,2], pos=g_order[,3], dup=dup_freq[reorder], del=del_freq[reorder], dup2=dup_freq_prolif[reorder], del2=del_freq_prolif[reorder]);

	require("karyoploteR")
	scale_fac <- max(stats[,3], stats[,4])
	png(paste(t, "part1", call_type ,"karyo.png", sep="_"), width=10, height=3, units="in", res=300)
	kp <- plotKaryotype(genome="hg19", plot.type=4, chromosomes=paste("chr", 1:10, sep=""))
        kpAxis(kp, ymin=0, ymax=scale_fac, r0=0, r1=1, side=1, data.panel=1)
        kpLines(kp, chr=stats$chr, x=stats$pos, y=stats$dup/scale_fac, col=dup_col, lwd=1.75)
        kpLines(kp, chr=stats$chr, x=stats$pos, y=stats$del/scale_fac, col=del_col, lwd=1.75)
        kpLines(kp, chr=stats$chr, x=stats$pos, y=stats$dup2/scale_fac, col=dup_col, lwd=1.75, lty=2)
        kpLines(kp, chr=stats$chr, x=stats$pos, y=stats$del2/scale_fac, col=del_col, lwd=1.75, lty=2)
	dev.off()

	png(paste(t, "part2", call_type, "karyo.png", sep="_"), width=10, height=3, units="in", res=300)
	kp <- plotKaryotype(genome="hg19", plot.type=4, chromosomes=paste("chr", 11:22, sep=""))
        kpAxis(kp, ymin=0, ymax=scale_fac, r0=0, r1=1, side=1, data.panel=1)
        kpLines(kp, chr=stats$chr, x=stats$pos, y=stats$dup/scale_fac, col=dup_col, lwd=1.75)
        kpLines(kp, chr=stats$chr, x=stats$pos, y=stats$del/scale_fac, col=del_col, lwd=1.75)
        kpLines(kp, chr=stats$chr, x=stats$pos, y=stats$dup2/scale_fac, col=dup_col, lwd=1.75, lty=2)
        kpLines(kp, chr=stats$chr, x=stats$pos, y=stats$del2/scale_fac, col=del_col, lwd=1.75, lty=2)
	dev.off()

}


