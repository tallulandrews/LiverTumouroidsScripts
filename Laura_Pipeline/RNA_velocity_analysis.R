args <- commandArgs(trailingOnly=TRUE) # velocity counts file, SCE object

# Testing
args <- c("FirstExp_velocity_counts.rds", "CCA1_SC3.rds", "HCC10_SC3.rds", "HCC6_SC3.rds")
#args <- c("Donor3_velocity_counts.rds", "D3DM_SC3.rds", "D3EM_SC3.rds")
#args <- c("Donor9_velocity_counts.rds", "D9DM_SC3.rds", "D9EM_SC3.rds")

dat <- readRDS(args[1])

require("scater")
source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")

for (f in 2:length(args)) {
	file <- args[f]
	prefix <- strsplit(file, "_")[[1]][1]

	cell.color.sets <- list()
	SCE <- readRDS(file)
	c.col <- get_group_cols(SCE)
	cell.colors <- c.col[factor(SCE$clusters_clean, levels=names(c.col))] # cluster colours

	# check cellIDs
	tmp <- strsplit(colnames(SCE)[1], "_")[[1]]
	if (length(tmp) == 2) {
		# add runID
		runID <- strsplit(colnames(dat$emat), "_")[[1]][1]
		colnames(SCE) <- sub("X", paste(runID, "_", sep=""), colnames(SCE))
	} else {
		colnames(SCE) <- sub("X", "", colnames(SCE))
	}
	names(cell.colors) <- colnames(SCE)
	dat_local <- dat
	local_cells <- colnames(SCE)
	dat_local$emat <- dat_local$emat[,local_cells]
	dat_local$iomat <- dat_local$iomat[,local_cells]
	dat_local$smat <- dat_local$smat[,local_cells]

	require("velocyto.R")

	#emb <- #dimensionality reduction coordinates

	# How does RNA velocity deal with names? - fixed

	# Check plots?
	#hist(log10(Matrix::rowSums(dat_local$emat)+1),col='wheat',
	#	xlab='log10[ number of reads + 1]',main='number of reads per gene')

	# Filter genes
	# exonic read (spliced) expression matrix
	emat <- dat_local$emat;
	# intronic read (unspliced) expression matrix
	nmat <- dat_local$iomat;
	# spanning read (intron+exon) expression matrix
	smat <- dat_local$smat;
	# filter expression matrices based on some minimum max-cluster averages
	emat <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 5, do.preview=TRUE)
	nmat <- filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 1, do.preview=TRUE)
	smat <- filter.genes.by.cluster.expression(smat,cell.colors,min.max.cluster.average = 0.5, do.preview=TRUE)
# look at the resulting gene set
	str(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))

	# Robust estimate?
	rvel.qf <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 10,fit.quantile = 0.02)

	saveRDS(rvel.qf, file=paste(prefix, "rnaVelo_output.rds", sep="_"))

	# Plotting
	png(paste(prefix, "rnaVelo_PCA.png", sep="_"), width=8, height=8, units="in", res=300)
	pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cell.colors=ac(cell.colors,alpha=0.9),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1))
	dev.off()
	# Plot gene
	#gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = 0.02,old.fit=rvel.qf,show.gene='Chga',cell.emb=emb,cell.colors=cell.colors)
}
