x = read.table("E-MTAB-3929_counts.txt", header=T)
labels = read.table("E-MTAB-3929_labels.txt", header=F)
require(M3Drop)
data_list = M3Drop_Clean_Data(x, labels=labels, is.counts=T)
png("E-MTAB-3929_counts_M3Drop_DE.png", width=6, height=6, units="in", res=300)
DE = M3Drop_Differential_Expression(data_list$data)
dev.off()
png("E-MTAB-3929_counts_M3Drop_heat.png", width=7, height=7, units="in", res=300)
heatout = M3Drop_Expression_Heatmap(DE[,1], data_list$data, cell_labels=unlist(data_list$labels))
dev.off()
C = M3Drop_Get_Heatmap_Cell_Clusters(heatout, k=7)
markers = M3Drop_getmarkers(data_list$data, unlist(C))
sig_markers = markers[p.adjust(markers$pval, method="bonferroni")<0.05,]

my_box_plot<-function(data,labels,groups,gene_name,colours, group_names=groups, plottitle=gene_name) {
        box_data <-list();
        N = length(groups);
        thresh = 0.05
        pvalues <- matrix(nrow=N, ncol=N);
        for (i in 1:length(groups)){
                box_data[[i]] = data[rownames(data) == gene_name, labels == groups[i]]
                if (i > 1) {
                        for(j in (1:(i-1))) {
                                p = suppressWarnings(wilcox.test(box_data[[i]],box_data[[j]])$p.value)
                                p = p*(N)*(N-1)/2
                                pvalues[i,j] = p
                                pvalues[j,i] = p
                        }
                }
        }
        pos = boxplot(box_data, col=colours, main="", ylab="", xlab="", notch=TRUE, names=group_names)
        title(ylab="Expression", line=2)
        title(main=plottitle, line=0.5)
        for (i in 1:length(groups)) {
                worstp = max(pvalues[i,],na.rm=T)
                if (worstp < thresh) {
                        level = abs(round(log(worstp)/log(10)))
                        stars = paste(rep("*",times=min(level,5)),collapse="")
                        text(i,pos$stats[4,i]*1.1,stars)
                }
        }
}

my_box_plot(data_list$data, C, unique(C), "SOX2", c("grey50","grey90","grey50","grey25","grey35","grey75","black"))
my_box_plot(data_list$data, C, unique(C), "PDGFRA", c("grey50","grey90","grey50","grey25","grey35","grey75","black"))
my_box_plot(data_list$data, C, unique(C), "TDGF1", c("grey50","grey90","grey50","grey25","grey35","grey75","black"))
my_box_plot(data_list$data, C, unique(C), "NANOG", c("grey50","grey90","grey50","grey25","grey35","grey75","black"))
my_box_plot(data_list$data, C, unique(C), "POU5F1", c("grey50","grey90","grey50","grey25","grey35","grey75","black"))
my_box_plot(data_list$data, C, unique(C), "GATA2", c("grey50","grey90","grey50","grey25","grey35","grey75","black"))
