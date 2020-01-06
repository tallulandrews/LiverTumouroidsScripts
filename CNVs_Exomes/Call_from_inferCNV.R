require("infercnv")
load("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/inferCNVRes/run.final.infercnv_obj")
G <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/Global_SCE.rds")
Cells <- which(G$Donor %in% unique(G$Donor));
#Cells <- infercnv_obj@observation_grouped_cell_indices$malignant_HCC10
vals <- infercnv_obj@expr.data[,Cells]
counts <- infercnv_obj@count.data[,Cells]
NAME <- "All"


# Order genes by location - already in order
#vals <- round(vals, digits=2)
gene_order <- infercnv_obj@gene_order
identical(rownames(vals), rownames(gene_order));

# Sliding Window to call CNVs
# Bonferroni correction across all cells for each window
ngenes_window = 20
null_hypothesis = 1;
min_var <- 0.05
calls <- matrix(0, nrow=nrow(vals), ncol=ncol(vals));
colnames(calls) <- colnames(vals)
sig_threshold <- 0.05;
for (i in 1:(nrow(vals)-ngenes_window)) {
	res <- apply(vals[i:(i+ngenes_window), ], 2, function(x) {c(mean(x), sd(x))})
	res[2,res[2,] < min_var] <- min_var;
	p.vals <- pnorm(abs(res[1,]-null_hypothesis)/res[2,], lower.tail=FALSE)
	sig <- p.vals < sig_threshold/ncol(res)
	dir <- sign(res[1,]-null_hypothesis)
	calls[i:(i+ngenes_window), sig & dir==1] <- 1
	calls[i:(i+ngenes_window), sig & dir==-1] <- -1
}

# Require calls to have been seen in at least 5% of cells
rep_threshold=0.05
png("infercnv_calls.png", width=8, height=4, units="in", res=300)
plot(rowMeans(calls>0), ylim=c(-1,1), type="l", col="red", lwd=2, ylab="del --- Proportion --- dup", xlab="Position")
lines(1:nrow(calls),-rowMeans(calls<0), ylim=c(-1,1), col="blue", lwd=2)
chr_lines <- vector()
for(c in unique(gene_order[,1])) {
	chr_lines <- c(chr_lines, min(which(gene_order[,1]==c)))
}
abline(v=chr_lines, col="grey65", lwd=2, lty=2)
abline(h=c(rep_threshold, -rep_threshold), lwd=1, lty=1, col="black")
dev.off();

calls_tidy <- calls;
calls_tidy[rowMeans(calls>0) < rep_threshold & rowMeans(calls<0) < rep_threshold,] <- 0

# Require calls to be at least one full window long 
# Reformat calls to be bedr style
REGIONS <- list();
REGIONS_dir <- list();

for(c in 1:ncol(calls_tidy)) {
	chr_curr = 0;
	dir_curr = 0;
	st_curr=-1;
	end_curr=-1;
	regions <- vector();
	regions_dir <- vector();
	regions_n <- vector()
	n = 0
	for (g in 1:nrow(calls_tidy)) {
		dir <- calls_tidy[g,c]
		chr <- gene_order[g,1]
		if (chr != chr_curr) {
#			region <- paste(chr_curr, ":", st_curr, "-", end_curr, sep="");
			region <- c(chr_curr, st_curr, end_curr);
			if (chr_curr != 0 & st_curr != -1) {
				regions <- c(regions, region);
				regions_dir <- c(regions_dir, dir_curr);
				regions_n <- c(regions_n, n)
				n = 0;
			}
			st_curr <- -1;
			dir_curr = dir;
			chr_curr = chr;
		}
		if (dir_curr != 0 & dir != dir_curr) {
#			region <- paste(chr_curr, ":", st_curr, "-", end_curr, sep="");
			region <- c(chr_curr, st_curr, end_curr);
			regions <- c(regions, region)
			regions_dir <- c(regions_dir, dir_curr);
			regions_n <- c(regions_n, n)
			st_curr <- -1;
			n = 0;
		}
		if (dir != 0) {
			if (st_curr < 0) {
				st_curr <- gene_order[g,2]
				end_curr <- gene_order[g,3]
				dir_curr <- dir;
				n = 1;
			} else {
				end_curr <- gene_order[g,3]
				n = n+1;
			}
		}
		dir_curr <- dir;

	}
	t <- data.frame(region=regions, dir=regions_dir, n=regions_n)
	REGIONS[[ colnames(calls_tidy)[c] ]] <- t[t[,3] > ngenes_window,];
}
saveRDS(REGIONS, paste(NAME, "inferCNV_called_CNVs.rds", sep="_"))
