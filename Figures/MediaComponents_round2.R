CustomPalette <- function (low = "white", high = "red", mid = NULL, k = 50) 
{
    low <- col2rgb(col = low)/255
    high <- col2rgb(col = high)/255
    if (is.null(x = mid)) {
        r <- seq(from = low[1], to = high[1], len = k)
        g <- seq(from = low[2], to = high[2], len = k)
        b <- seq(from = low[3], to = high[3], len = k)
    }
    else {
        k2 <- round(x = k/2)
        mid <- col2rgb(col = mid)/255
        r <- c(seq(from = low[1], to = mid[1], len = k2), seq(from = mid[1], 
            to = high[1], len = k2))
        g <- c(seq(from = low[2], to = mid[2], len = k2), seq(from = mid[2], 
            to = high[2], len = k2))
        b <- c(seq(from = low[3], to = mid[3], len = k2), seq(from = mid[3], 
            to = high[3], len = k2))
    }
    return(rgb(red = r, green = g, blue = b))
}


my_DotPlot_raw <- function(size.val=size_matrix, colour.val=colour_matrix, 
	colours = CustomPalette(low = "magenta", high = "yellow", mid = "black", k = 21),
	min.cex=0, max.cex=3, low.size.threshold = 0, low.col.threshold=-Inf) {
	if (! identical(dim(size.val), dim(colour.val)) ) {
		stop("Size and colour matrices are not the same dimension.")
	}

	exclude <- size.val < low.size.threshold | colour.val < low.col.threshold;
	size.val[exclude] <- NA;

	if (min(colour.val) < 0 & max(colour.val) > 0) {
		bins <- seq(from=0, to=max(abs(colour.val)), length=ceiling((length(colours)+1)/2));
		bins[length(bins)] <- bins[length(bins)]+1
		bins <- unique(c(-1*rev(bins), bins))
	} else {
		bins <- seq(from=min(colour.val), to=max(colour.val), length=length(colours+1));
		bins[length(bins)] <- bins[length(bins)]+1
		bins[1] <- bins[1]-1
	}
	x_coords <- matrix(rep(1:nrow(size.val), each=ncol(size.val)), ncol=ncol(size.val), byrow=T)
	y_coords <- matrix(rep(1:ncol(size.val), each=nrow(size.val)), ncol=ncol(size.val), byrow=F)
	dot.colour <- colours[cut(as.vector(colour.val), breaks=bins)]
	size.cex = size.val-min(size.val, na.rm=T); size.cex = size.cex/max(size.cex, na.rm=T)*max.cex;
	size.cex <- as.vector(as.matrix(round(size.cex, digits=1)))
	plot(as.vector(x_coords), as.vector(y_coords), col=dot.colour, cex=size.cex, pch=16, 
		xlab="", ylab="", xaxt="n", yaxt="n", bty="L", xaxs="i", xlim=c(0,max(x_coords)+1))
	axis(1, at=1:max(x_coords), label=rownames(size.val), las=2)
	axis(2, at=1:max(y_coords), label=colnames(size.val), las=2)
	invisible(list(xes=as.vector(x_coords), yes=as.vector(y_coords), col=dot.colour, cex=size.cex))

}

files <- Sys.glob("*Marker*Nov*.csv")

pathways <- read.table("B27_componnents_genesets_Harmonizome database_aggregate_gene_list.csv", sep=",", header=T)

categories <- unique(pathways[,1])

max_genes_dotplot = 50

for (f in files) {
	ID <- unlist(strsplit(f, "_"))[1]
	tab <- read.table(f, sep=",", header=T)
	for (c in categories) {
		# Add column to tab
		tab[,c] <- tab$Symbol %in% pathways[ pathways[,1] == c ,2]

		# Select data for dotplot
		this_tab <- tab[tab[,c],]
		this_tab <- this_tab[order(this_tab$p.values, decreasing=FALSE),]
		if (nrow(this_tab) > max_genes_dotplot) {
			dat_tab <- this_tab[1:max_genes_dotplot,]
		} else {	
			dat_tab <- this_tab
		}
		# Make dotplot
		manual_clusters <- 3:which(colnames(this_tab) == "p.values")-1
		manual_clusters_name <- colnames(this_tab)[manual_clusters]

		size_mat <- dat_tab[,manual_clusters+max(manual_clusters)+2] # Expression level
		on_off_mat <- dat_tab[,manual_clusters]
		col_mat <- matrix( log10(rep(dat_tab$p.values, times=ncol(size_mat))) , ncol=ncol(size_mat) )
		col_mat[on_off_mat == 1] = abs(col_mat[on_off_mat == 1])
		colnames(size_mat) <- manual_clusters_name;
		rownames(size_mat) <- dat_tab$Symbol

	
		png(paste(ID, sub(" ", "_", c), "GeneDotplot.png", sep="_"), width=12, height=8, units="in", res=200)
		par(mar=c(6, 8, 2, 1))
		my_colours <- CustomPalette(low = "magenta", high = "yellow", mid = "black", k = 21)
		my_colours <- my_colours[c(1:5, 10, 11, 16:20)]
		my_DotPlot_raw(size.val=size_mat, colour.val=col_mat, colours = my_colours)
		significance <- 0.05/nrow(this_tab)
		abline(v = max(which(dat_tab$p.values < significance))+0.5, lty=2)
		dev.off()

	}
}

#Add whether each gene is associated with eac pathway from the media component to the giant marker tables

#Make dotplots where dotsize = expression level, colour = pvalue for those where it is "on"
