NB.fit <- function(x, lib.size=rep(1, length(x)), max_r=10^10) {
	# Defunct
	# Fit
	mus <- mean(x)*lib.size
	obs_err <- sum( (x - mus)^2 )
	rg <- sum( mus^2 )/(obs_err - sum(mus))
	if (rg <= 0) {rg <- max_r}
	return(list(mu=mean(x), size=rg, d=0));
}

ZINB.fit <- function(x, lib.size=rep(1, length(x)), max_r=10^10, e=0.00001) {
	# Defunct
	# Fit
	d.obs <- mean(x==0)
	if (d.obs == 0) {return(NB.fit(x, lib.size))}
	d_curr <- d.obs;
        d_prev <- -100;
	nc <- length(x);
	while( abs(d_curr-d_prev) > e ) {
		mus <- sum(x)/(nc-d_curr*nc)*lib.size
		weights <- rep(1, times=length(x))
		weights[x == 0] <- (1-d_curr/d.obs)
		obs_err <- sum( (x - mus)^2*weights )
		rg <- sum( mus^2*weights )/(obs_err - sum(mus*weights))
                if (rg <= 0) {rg <- max_r}
		pds <- (1 + mus/rg)^(-rg)
                d.exp <- mean(pds)
		d_prev <- d_curr
                d_curr <- (d.obs - d.exp)
		if (d_curr <= 0) {d_curr <- d_prev}
	}
	return(list(mu=mus[1]/lib.size[1], size=rg, d=d_prev))
}

scale_mat_to_ref <- function(gene_expr, ref_cells, lib.size) { # Tested!
	# Defunct
	lib.size <- lib.size/mean(lib.size)
	zinb.params <- ZINB.fit(gene_expr[ref_cells], lib.size[ref_cells]);
	mus <- as.numeric(zinb.params$mu*lib.size)
	p0s <- zinb.params$d+(1-zinb.params$d)*sapply(mus, function(m) {pnbinom(0, mu=m, size=zinb.params$size)})
	ps <- sapply(1:length(gene_expr), function(i) {pnbinom(gene_expr[i], mu=mus[i], size=zinb.params$size)})
	ps[gene_expr==0] <- p0s[gene_expr==0];
	z <- qnorm(ps)
	z[z>10] <- 10
	return(z);
}

scale_mat_to_ref_simple <- function(gene_expr, ref_cells, lib.size) { # Tested!
	# Defunct
	lib.size <- lib.size/mean(lib.size)
	gene_expr <- log2(gene_expr/lib.size+1)
	gene_expr <- ( gene_expr-median(gene_expr[ref_cells]))/ sd(gene_expr[ref_cells])
	return(gene_expr)
}

quick_normalize <- function(raw_counts, lib.size) {
	lib.size <- lib.size/mean(lib.size)
	norm <- log2(t(t(raw_counts)/lib.size)+1)
	return(norm);
}


get_null_distr <- function(gene_expr, ref_cells, lib.size) {
	gene_expr_null <- gene_expr[ref_cells]
	return( c(mean(gene_expr_null[gene_expr_null >0]), sd(gene_expr_null[gene_expr_null > 0]), mean(gene_expr_null == 0)) )
}

id_cnv_total <- function(norm_expr, gene_order, ref_cells, 
			window_size=100, window_step=5, 
			exclude_sex=TRUE, sig_threshold=0.01){
	# Order scaled matrix by gene order
	gene_order <- gene_order[gene_order$Gene %in% rownames(norm_expr),]
	norm_expr <- norm_expr[match(rownames(gene_order), rownames(norm_expr)),]
	# Exclude sex chromosomes if requested
	if (exclude_sex) {
		sex <- gene_order$chr %in% c("chrX", "chrY");
		gene_order <- gene_order[!sex,]
		norm_expr <- norm_expr[!sex,]
	}
	OUT <- list()
	
	# Sliding window along each chromosome.
	for (chr in as.character(unique(gene_order$chr))) {
		print(chr);
		chr_vals <- norm_expr[gene_order$chr == chr,]
		out <- matrix(0, nrow=nrow(chr_vals), ncol= ncol(chr_vals));
		if (nrow(chr_vals) < window_size+5) {next;}
		# Sliding window
		slide_window <- seq(from=1, to=(nrow(chr_vals)-window_size), by=window_step)
		for (i in 1:slide_window) {
			# Is the proportion of total expression in this window higher/lower than
			# in controls?
			V <- colSums(chr_vals[i:(i+window_size),])/colSums(norm_expr)
			N <- V[ref_cells]
			effect <- (V-mean(N))/mean(N)
			p <- pnorm(abs((V-mean(N))/sd(N)), lower.tail=FALSE);
			# Ignore statistics and just compare to distribution of controls:
			sig <- p < quantile(p[ref_cells], probs=sig_threshold);
			# Check if more than expected in tumor cells
			check <- prop.test(sum(sig[!ref_cells]), sum(!ref_cells), 
					 p=sum(sig[ref_cells])/sum(ref_cells), alternative="greater")
			effect[!sig] <- 0;
			if (check$p.value < 0.05) {
				# There is evidence of CNVs!
				# What genes contribute?
				gene_ref <- rowMeans(chr_vals[i:(i+window_size),ref_cells])
				replace <- sign(chr_vals[i:(i+window_size),]-gene_ref)
				replace <- t( apply(replace,1, function(row){
						row <- row*sign(effect); 
						row[row < 0] <- 0;
						row <- row*sign(effect); 
						return(row);
						}))
				
				# Save results
				un_set <- which(rowSums(out) == 0);
				replace_start <- max(min(un_set[un_set >= i]), i)
				out[replace_start:(i+window_size),] <- replace[ (replace_start-i+1):nrow(replace), ]
			}
		}
		colnames(out) <- colnames(chr_vals);
		rownames(out) <- rownames(chr_vals);
		OUT[[chr]] <- out;
	}
	return(OUT)
}

id_cnv_detect_rate <- function(norm_expr, null_distributions, gene_order, ref_cells,
			fc_threshold=0.5, sig_threshold=0.05, 
			window_size=100, window_step=5, exclude_sex=TRUE) {
	# UNFINISHED

	# Order scaled matrix by gene order
	gene_order <- gene_order[gene_order$Gene %in% rownames(norm_expr),]
	norm_expr <- norm_expr[match(rownames(gene_order), rownames(norm_expr)),]
	null_distributions <- null_distributions[match(rownames(gene_order), rownames(null_distributions)),]
	# Exclude sex chromosomes if requested
	if (exclude_sex) {
		sex <- gene_order$chr %in% c("chrX", "chrY");
		gene_order <- gene_order[!sex,]
		norm_expr <- norm_expr[!sex,]
		null_distributions <- null_distributions[!sex,]
	}

	
	zeds <- (norm_expr - null_distributions[,1])/null_distributions[,2];
	null_prop_up <- sum(zeds[norm_expr>0] > fc_threshold)/sum(norm_expr>0)
	null_prop_down <- sum(zeds[norm_expr>0] < -fc_threshold)/sum(norm_expr>0)
	
	OUT <- list()
	
	# Sliding window along each chromosome.
	for (chr in as.character(unique(gene_order$chr))) {
		print(chr);
		chr_vals <- norm_expr[gene_order$chr == chr,]
		chr_zeds <- zeds[gene_order$chr == chr,]
		out <- matrix(0, nrow=nrow(chr_vals), ncol= ncol(chr_vals));
		if (nrow(chr_vals) < window_size+5) {next;}
		# Sliding window
		slide_window <- seq(from=1, to=(nrow(chr_vals)-window_size), by=window_step)
		for (i in 1:slide_window) {
			# DE non-zeros
			Z <- chr_zeds[i:(i+window_size),]
			V <- chr_vals[i:(i+window_size),]
			N <- null_distributions[i:(i+window_size),]
			Z[V==0] <- 0;
			n_up <- colSums(Z > fc_threshold)
			n_dw <- colSums(Z < -1*fc_threshold)
			pos <- colSums(V>0);
			p_pos_up <- apply(cbind(n_up, pos), 1, function(x){prop.test(x[1],x[2], p=null_prop_up)$p.value})
			p_pos_dw <- apply(cbind(n_dw, pos), 1, function(x){prop.test(x[1],x[2], p=null_prop_down)$p.value})
			effect1 <- (n_up/pos - null_prop_up)
			effect12 <- (null_prop_down - n_dw/pos)
			effect1[p_pos_up > quantile(p_pos_up[ref_cells], probs=sig_threshold)] <- 0
			effect12[p_pos_dw > quantile(p_pos_dw[ref_cells], probs=sig_threshold)] <- 0
			effect1[effect1*effect12 < 0] <- 0
			effect12[effect1*effect12 < 0] <- 0
			effect1 <- effect1+effect12;

			# DE zeros
			zeros <- colSums(V==0);
			zero_exp <- sum(N[,3]);
			zero_sd <- sqrt(sum(N[,3]*(1-N[,3])));
			p_zero <- pnorm(abs(zeros-zero_exp)/zero_sd, lower.tail=FALSE)

			sig <- apply(cbind(p_pos, p_zero), 1, min)
			is.var <- sig < quantile(sig[ref_cells], probs=sig_threshold)
			effect2 <- (zero_exp-zeros)/nrow(V)
			effect2[p_zero > quantile(p_zero[ref_cells], probs=sig_threshold)] <- 0

			# Direction:
			sign(effect1)
			
			# Contributing Genes
		}
	}
}


				






test_cnv_detect_rate <- function(gene_expr, null_distr, prop.null=0.5, fc_threshold=2, sig_threshold=0.05) {
	out <- rep(0, length(gene_expr))
	# non-zeros
	pos <- which(gene_expr > 0);
	zeds <- (gene_expr[pos] - null_distr[pos,1])/null_distr[pos,2]
	n_up <- sum(zeds > fc_threshold)
	n_down <- sum(zeds < -1*fc_threshold)
	if (n_up + n_down ==0) {p_pos=1} else {
		p_pos <- prop.test(n_up, length(pos), p=prop.null)$p.value
	}
	# zeros
	n_0_obs <- sum(gene_expr == 0);
	n_0_exp <- sum(null_distr[,3]);
	n_0_sd <- sqrt(sum(null_distr[,3]*(1-null_distr[,3])));
	z <- (n_0_obs-n_0_exp)/n_0_sd
	p_zero <- pnorm(abs(z), lower.tail=FALSE)
	# Sig?
	if (min(p_pos, p_zero) > sig_threshold) {return(out)}
	# Which direction?
	dir = 0;
	effect <- 0;
	effect1 <- (mean(gene_expr)-mean(null_distr[,1]))/mean(null_distr[,1]);
	effect2 <- (mean(null_distr[,3])-mean(gene_expr==0))/mean(null_distr[,3]);
	if (n_0_obs < n_0_exp & n_up > n_down) {
		dir <- 1;
		effect <- max(effect1, effect2);
	} else if (n_0_obs > n_0_exp & n_down > n_up) {
		dir <- -1;
		effect <- min(effect1, effect2);
	} else {
		if (p_zero < p_pos) {
			dir <- sign(n_0_exp-n_0_obs)
			effect <- effect2
		} else {
			dir <- sign(n_up-n_down)
			effect <- effect1;
		}
	}
	# contributing genes:
	if (dir == -1) {
		top_n <- round(n_0_obs-n_0_exp)
		g <- null_distr[gene_expr==0,];
		g <- g[order(g[,3], decreasing=FALSE),]
		g <- rownames(g)[1:top_n]
		contrib <- rownames(null_distr) %in% g | zeds < -1*fc_threshold;
		out[contrib] <- -1*abs(effect)
	} else if (dir == 1) {
		top_n <- round(n_0_exp-n_0_obs)
		g <- null_distr[gene_expr>0,];
		g <- g[order(g[,3], decreasing=TRUE),]
		g <- rownames(g)[1:top_n]
		contrib <- rownames(null_distr) %in% g | zeds > 1*fc_threshold;
		out[contrib] <- 1*abs(effect)
	}
	return(out)
}

test_cnv <- function(z_scores, null_up=1, null_down=null_up, fc_threshold=1, sig_threshold=0.05) {
	out <- rep(0, length(z_scores))
	n_up <- sum(z_scores > fc_threshold)
	n_down <- sum(z_scores < -1*fc_threshold)
	if (n_up == 0 & n_down == 0){return(out)} 
	p <- prop.test(n_up, n_up+n_down, p=null_up/(null_up+null_down))$p.value
	# Geometric mean of fc
	if (p < sig_threshold & n_up > n_down) {
		contrib <- z_scores > 0
		fc <- exp( sum(log(abs(z_scores[contrib])))*(1/sum(contrib)))
		out[contrib] <- fc
	}
	if (p < sig_threshold & n_down > n_up) {
		contrib <- z_scores < 0
		fc <- exp( sum(log(abs(z_scores[contrib])))*(1/sum(contrib)))
		out[contrib] <- -1*fc
	}
	return(out);
}

id_cnvs<- function(scaled_count_mat, gene_order, null_distr=NULL, window_size=100, window_step=5, fc_threshold=2, tailor_null=TRUE, exclude_sex=TRUE) {
	# Order scaled matrix by gene order
	gene_order <- gene_order[gene_order$Gene %in% rownames(scaled_count_mat),]
	scaled_count_mat <- scaled_count_mat[match(rownames(gene_order), rownames(scaled_count_mat)),]
	# Exclude sex chromosomes if requested
	if (exclude_sex) {
		sex <- gene_order$chr %in% c("chrX", "chrY");
		gene_order <- gene_order[!sex,]
		scaled_count_mat <- scaled_count_mat[!sex,]
	}

	OUT <- list()	
	# Expected ratio of up/down: either average across all dataset or 1:1
	if (tailor_null) {
		null_up <- sum(scaled_count_mat > fc_threshold);
		null_down <- sum(scaled_count_mat < -fc_threshold);
	} else {
		null_up <- 1
		null_down <- 1
	}
	# Sliding window along each chromosome.
	for (chr in as.character(unique(gene_order$chr))) {
		print(chr);
		chr_vals <- scaled_count_mat[gene_order$chr == chr,]
		out <- matrix(0, nrow=nrow(chr_vals), ncol= ncol(chr_vals));
		if (nrow(chr_vals) < window_size+5) {next;}
		# Sliding window
		slide_window <- seq(from=1, to=(nrow(chr_vals)-window_size), by=window_step)
		for (i in 1:slide_window) {
			if (!is.null(null_distr)) {
				res <- apply(chr_vals[i:(i+window_size),], 2, test_cnv_detect_rate, null_distr=null_distr[i:(i+window_size),], fc_threshold=fc_threshold, sig_threshold=0.05/ncol(chr_vals))
			} else {
				res <- apply(chr_vals[i:(i+window_size),], 2, test_cnv, fc_threshold=fc_threshold, sig_threshold=0.05/ncol(chr_vals))
			}
			un_set <- which(rowSums(out) == 0);
			replace_start <- max(min(un_set[un_set >= i]), i)
			out[replace_start:(i+window_size),] <- res[ (replace_start-i+1):nrow(res), ]

		}
		rownames(out) <- rownames(chr_vals);
		colnames(out) <- colnames(chr_vals);
		OUT[[chr]] <- out;
	}
	return(OUT);
}

call_cnv_regions <- function(id_cnv_output, gene_order, max_gap_g=100, max_gap_kb=20000) {
	all_regions <- list();
	for (chr in names(id_cnv_output)) {
		vals <- id_cnv_output[[chr]]

		for (c in 1:ncol(vals)) {
			OUT <- vector();

			current_region <- list(chr=chr, st=-1, end=-1, ngene=0, dir=10, last_gene=-1)
			for (i in which(vals[,c] != 0)) {
				g_info <- gene_order[gene_order[,1] == rownames(vals)[i], ]
				this_dir <- sign(vals[i,c])
				new_region=FALSE
				# Check gene gap
				if (i - current_region$last_gene > max_gap_g) {
					new_region=TRUE
				}
				# Check kb gap
				if (g_info[3] - current_region$end > max_gap_kb*1000) {
					new_region=TRUE
				}
				# Check direction
				if (this_dir != current_region$dir) {
					new_region=TRUE
				}

				if (new_region) {
					if (current_region$st != -1) {
						# Save current
						OUT <- rbind(OUT, unlist(current_region))
					}
					# Reset
					current_region$st <- g_info[3]
					current_region$end <- g_info[4]
					current_region$ngene <- 1
					current_region$dir <- this_dir
					current_region$last_gene <- i;
				} else {
					# Extend current
					current_region$end <- max(g_info[4], current_region$end)
					current_region$ngene <- current_region$ngene+1
					current_region$last_gene <- i;
				}

			}
			OUT <- rbind(OUT, unlist(current_region))
			if (length(all_regions) >= c) {
				all_regions[[c]] <- rbind(all_regions[[c]], OUT);
			} else {
				 all_regions[[c]] <- OUT;
			}
		}
	}
	return(all_regions);
}

exit()

### Notes ###

require("SingleCellExperiment")
G <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/Global_SCE.rds")

raw_counts <- G@assays[["counts"]]
control_cells <- G$Donor %in% c("D3DM", "D3EM", "D9DM", "D9EM")
raw_counts <- raw_counts[rowMeans(raw_counts[,control_cells] > 0) > 0.01,]
raw_counts <- raw_counts[rowMeans(raw_counts[,!control_cells] > 0) > 0.01,]

# Scale it
lib.size <- colSums(raw_counts)
scaled_mat <- t(apply(raw_counts, 1, scale_mat_to_ref_simple, control_cells, lib.size))
scaled_mat_zinb <- t(apply(raw_counts, 1, scale_mat_to_ref, control_cells, lib.size))

norm_mat <- quick_normalize(raw_counts, lib.size);
null_distributions <- t(apply(norm_mat, 1, get_null_distr, control_cells, lib.size))
colnames(null_distributions) <- c("mean", "sd", "detect")

# Call CNVs
gene_order <- read.table("inferCNV_input_geneorder.txt")
colnames(gene_order) <- c("Gene", "chr", "st", "end")
rownames(gene_order) <- gene_order[,1]

# Test method:
test_out <- id_cnv_total(norm_mat, gene_order, control_cells, sig_threshold=0.05)

# Exclude lineage markers
line_genes <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Cleaned_PatLineage.txt")
source("~/R-Scripts/Ensembl_Stuff.R")
line_genes <- General_Map(as.character(line_genes[,1]), in.org="Hsap", out.org="Hsap", in.name="symbol", out.name="ensg")
scaled_mat <- scaled_mat[! rownames(scaled_mat) %in% line_genes,]
scaled_mat_zinb <- scaled_mat_zinb[! rownames(scaled_mat_zinb) %in% line_genes,]
null_distributions <- null_distributions[! rownames(null_distributions) %in% line_genes,]
norm_mat <- norm_mat[! rownames(norm_mat) %in% line_genes,]

RES <- id_cnvs(scaled_mat, gene_order, fc_threshold=2, window_size=200)
REGIONS <- call_cnv_regions(RES, gene_order)
saveRDS(RES, "MyCNVs_Simplenorm.rds")
saveRDS(REGIONS, "MyCNVs_regions_Simplenorm.rds")

RES <- id_cnvs(scaled_mat_zinb, gene_order, window_size=300)
REGIONS <- call_cnv_regions(RES, gene_order)
saveRDS(RES, "MyCNVs_NBnorm.rds")
saveRDS(REGIONS, "MyCNVs_regions_NBnorm.rds")

### Line specific frequencies ###
plot_cnv_freq_karyo <- function(cells, REGIONS, name) {
	require("bedr")
	require("karyoploteR")
	dups <- list()
	dels <- list()
	for (c in cells) {
		regions <- REGIONS[[c]]
		if (sum(regions[,"dir"] == 1) > 0) {
			#dups[[c]] <- regions[regions[,"dir"] == 1,]
			tmp <- paste(regions[regions[,"dir"] == 1,1], ":", regions[regions[,"dir"] == 1,2],"-", regions[regions[,"dir"] == 1,3], sep="")
			dups[[length(dups)+1]] <- tmp#regions[regions[,"dir"] == 1,1:3]
		} 
		if (sum(regions[,"dir"] == -1) > 0) {
			#dels[[c]] <- regions[regions[,"dir"] == -1,]
			tmp <- paste(regions[regions[,"dir"] == -1,1], ":", regions[regions[,"dir"] == -1,2],"-", regions[regions[,"dir"] == -1,3], sep="")
			dels[[length(dels)+1]] <-tmp
		}
	}

	consistency_dups <- bedr.join.multiple.region(dups)
	consistency_dels <- bedr.join.multiple.region(dels)

	# need to add splitting rownames of consistency to get positions back out.

	dup_stats <- data.frame(chr=c(consistency_dups[,1], consistency_dups[,1]), pos=c(consistency_dups[,2], consistency_dups[,3]), val=as.numeric(c(consistency_dups$n.overlaps, rep(0, nrow(consistency_dups)))))
	dup_stats <- dup_stats[!duplicated(dup_stats[,c(1,2)]),]
	del_stats <- data.frame(chr=c(consistency_dels[,1], consistency_dels[,1]), pos=c(consistency_dels[,2], consistency_dels[,3]), val=as.numeric(c(consistency_dels$n.overlaps, rep(0, nrow(consistency_dels)))))
	del_stats <- del_stats[!duplicated(del_stats[,c(1,2)]),]

	#scale_fac<-max(dup_stats$val, del_stats$val)
	scale_fac <- length(cells);

	#png(paste(name, "Karyotype_myCNVs_part1.png", sep="_"), width=10, height=2, units="in", res=300)
	X11()
	kp <- plotKaryotype(genome="hg19", plot.type=4, chromosomes=paste("chr", 1:10, sep=""))
	kpAxis(kp, ymin=0, ymax=scale_fac, r0=0, r1=1, side=1, data.panel=1)
	kpLines(kp, chr=dup_stats$chr, x=dup_stats$pos, y=dup_stats$val/scale_fac, col="red", lwd=1.75)
	kpLines(kp, chr=del_stats$chr, x=del_stats$pos, y=del_stats$val/scale_fac, col="blue", lwd=1.75)
#	dev.off()

#	png(paste(name, "Karyotype_myCNVs_part2.png", sep="_"), width=10, height=2, units="in", res=300)
	X11()
	kp <- plotKaryotype(genome="hg19", plot.type=4, chromosomes=paste("chr", 11:22, sep=""))
	kpAxis(kp, ymin=0, ymax=scale_fac, r0=0, r1=1, side=1, data.panel=1)
	kpLines(kp, chr=dup_stats$chr, x=dup_stats$pos, y=dup_stats$val/scale_fac, col="red", lwd=1.75)
	kpLines(kp, chr=del_stats$chr, x=del_stats$pos, y=del_stats$val/scale_fac, col="blue", lwd=1.75)
#	dev.off()
}

REGIONS<-readRDS("MyCNVs_regions_Simplenorm.rds")
for (donor in unique(G$Donor)) {
	cells <- which(G$Donor == donor);
	name <- paste(donor, "Simplenorm", sep="_")
	plot_cnv_freq_karyo(cells, REGIONS, name)
	cells <- which(G$Donor == donor & G$Proliferating);
	name <- paste(donor, "Simplenorm_Prolif", sep="_")
	plot_cnv_freq_karyo(cells, REGIONS, name)
}

REGIONS<-readRDS("MyCNVs_regions_NBnorm.rds")
for (donor in unique(G$Donor)) {
	cells <- which(G$Donor == donor);
	name <- paste(donor, "NBnorm", sep="_")
	plot_cnv_freq_karyo(cells, REGIONS, name)
	cells <- which(G$Donor == donor & G$Proliferating);
	name <- paste(donor, "NBnorm_Prolif", sep="_")
	plot_cnv_freq_karyo(cells, REGIONS, name)
}

