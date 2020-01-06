source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/99_Functions.R")

set.seed(1973)

args <- commandArgs(trailingOnly=TRUE) # rds file for data, prefix for output
expr_type <- "lognorm"


require("scater")
SCE <- readRDS(args[1]);

gene_lists <- load_CC("all")
cellcycle <- gene_lists$Whitfield
cellcycle_simple <- gene_lists$Simple
GO_genes <- gene_lists$G0
Quiescence <- gene_lists$Quiescence
new_cellcycle <- gene_lists$Tirosh
#cellcycle <- read.table("~/Data/Whitfield_CC.txt")
#cellcycle_simple <- as.matrix(cellcycle[cellcycle[,1] != "CC",])
#cellcycle_simple[cellcycle_simple[,1] == "G2",1] = "G2M";
#cellcycle_simple[cellcycle_simple[,1] == "S",1] = "G1S";
#cellcycle_simple = cellcycle_simple[cellcycle_simple[,1] != "MG1",];
#G0_genes <- read.table("~/Data/Reactome_G0.txt", header=F) # Update this with genes from Laura/MiSigDB?
#Quiescence <- read.table("~/Data/Quiescence.txt")
#new_cellcycle <- read.table("~/Collaborations/LiverOrganoids/New_CC_171117.txt", header=FALSE) #Tirsoh2016Nature

get_prolif <- function(CC, SCE, suppress.plot=TRUE) {
	require(mixtools)
	cc_g <- fData(SCE)$Symbol %in% CC
	score <- colSums(get_exprs(SCE, expr_type)[cc_g,])
	mix <- normalmixEM(score)
	if (!suppress.plot) {
		plot(mix, which=2)
	}
	p1 <- dnorm(score, mean = mix$mu[1], sd=mix$sigma[1]) 
	p2 <- dnorm(score, mean = mix$mu[2], sd=mix$sigma[2])
	if (mix$mu[1] < mix$mu[2]) {
		assign <- p2 > p1
	} else {
		assign <- p1 > p2
	}
	return(list(score=score, assign=assign))
}

get_prolif2 <- function(CC, SCE, name="cycle genes", suppress.plot=TRUE, min_detected = 0.05) {
	require(mixtools)
        require(matrixStats)
        loggednorm <- get_exprs(SCE, expr_type)
        cc_g <- fData(SCE)$Symbol %in% CC
	loggednorm <- loggednorm[cc_g,]
	loggednorm <- loggednorm[rowVars(loggednorm) > 0,]
	loggednorm <- loggednorm[rowMeans(loggednorm > 0) >= min_detected,]
        loggednorm <- (loggednorm-rowMeans(loggednorm))/rowVars(loggednorm)
        score <- colMeans(loggednorm)
        mix <- normalmixEM(score)
        p1 <- dnorm(score, mean = mix$mu[1], sd=mix$sigma[1])
        p2 <- dnorm(score, mean = mix$mu[2], sd=mix$sigma[2])
	
	
	m1 = mix$mu[1]
	m2 = mix$mu[2]
	means_test <- abs(m1-m2) > min(mix$sigma)
	

        if (m1 < m2 & means_test) {
                assign <- p2 > p1
        } else if (m2 < m1 & means_test) {
                assign <- p1 > p2
        } else {
		assign <- rep(FALSE, times=length(p1))
		name=paste(name, "- No Difference")
	}
        if (!suppress.plot) {
                plot(mix, which=2, main2=name, xlab2="Score")
        }
        return(list(score=score, assign=assign))

}

get_prolif3 <- function(CC, SCE, name="cycle genes", suppress.plot=TRUE, min_detected = 0.05) {
        require(mixtools)
        require(matrixStats)
        loggednorm <- get_exprs(SCE, expr_type)
        cc_g <- fData(SCE)$Symbol %in% CC[,1]
        loggednorm <- loggednorm[cc_g,]
	rownames(loggednorm) <-fData(SCE)$Symbol[cc_g]; 
        loggednorm <- loggednorm[rowVars(loggednorm) > 0,]
        loggednorm <- loggednorm[rowMeans(loggednorm > 0) >= min_detected,]
        loggednorm <- (loggednorm-rowMeans(loggednorm))/rowVars(loggednorm)
	
	rownames(CC) <- CC[,1]
	weights <- CC[match(rownames(loggednorm), CC[,1]),]
	loggednorm <- loggednorm*weights[,2]
	
        score <- colMeans(loggednorm)
        mix <- normalmixEM(score)
        p1 <- dnorm(score, mean = mix$mu[1], sd=mix$sigma[1])
        p2 <- dnorm(score, mean = mix$mu[2], sd=mix$sigma[2])
        
        
        m1 = mix$mu[1]
        m2 = mix$mu[2]
        means_test <- abs(m1-m2) > min(mix$sigma)
        

        if (m1 < m2 & means_test) {
                assign <- p2 > p1
        } else if (m2 < m1 & means_test) {
                assign <- p1 > p2
        } else {
                assign <- rep(FALSE, times=length(p1))
                name=paste(name, "- No Difference")
        }
        if (!suppress.plot) {
                plot(mix, which=2, main2=name, xlab2="Score")
        }
        return(list(score=score, assign=assign))

}


prolif_analysis<- function(SCE, name="Data") {
	
	# Get Prolif
	png(paste(name, "Proliferation_MixtureModels.png", sep="_"), width=8, height=8*2/3, units="in", res=300)
	par(mfrow=c(2,3))
	sets <- levels(factor(cellcycle[,1]))
	out <- matrix(0, nrow=length(pData(SCE)[,1]), ncol=length(sets))
	for( i in 1:length(sets) ) {
		s = sets[i]
		assign<-get_prolif2(cellcycle[cellcycle[,1]==s,2], SCE, name=s, suppress.plot=F)$assign
		out[,i] = assign
	}
	dev.off()

	summary(factor(rowSums(out)))
	prolif = rowSums(out) >= 3

	# Output
	return(prolif)
}

prolif_analysis2 <- function(SCE, name="Data") {
	
#	g0 <- get_prolif2(G0_genes[,1], SCE, name="G0", suppress.plot=F)$assign
	png(paste(name, "MixtureModels2_Figure.png", sep="_"), width=7*3/2, height=7, units="in", res=300)
	par(mfrow=c(2,3))
	par(mar=c(4,4,1,1))
	g1 <- get_prolif2(cellcycle_simple[cellcycle_simple[,1] == "G1S",2], SCE, name="G1S", suppress.plot=F)
	g1_new <- get_prolif2(new_cellcycle[new_cellcycle[,2] == "G1S",1], SCE, name="(new) G1S", suppress.plot=F)
	g2 <- get_prolif2(cellcycle_simple[cellcycle_simple[,1] == "G2M",2], SCE, name="G2M", suppress.plot=F)
	g2_new <- get_prolif2(new_cellcycle[new_cellcycle[,2] == "G2M",1], SCE, name="(new) G2M", suppress.plot=F)
	g0 <- get_prolif3(Quiescence, SCE, name="G0", suppress.plot=F)
	dev.off()

	Assign <- cbind(g1$assign, g2$assign, g0$assign)
	colnames(Assign) = c("G1S", "G2M", "G0")
	Score <- cbind(g1$score, g2$score, g0$score)
	colnames(Score) = c("G1S", "G2M", "G0")
	
	state <- apply(Score, 1, function(x) {colnames(Score)[which(x == max(x))]})
	state[state=="G0" & !Assign[,3]] <- "None"
	state[state=="G1S" & !Assign[,1]] <- "None"
	state[state=="G2M" & !Assign[,2]] <- "None"

	state1 <- factor(state, levels=c("None", "G0", "G1S", "G2M"))

	Assign <- cbind(g1_new$assign, g2_new$assign, g0$assign)
	colnames(Assign) = c("G1S", "G2M", "G0")
	Score <- cbind(g1_new$score, g2_new$score, g0$score)
	colnames(Score) = c("G1S", "G2M", "G0")
	
	state <- apply(Score, 1, function(x) {colnames(Score)[which(x == max(x))]})
	state[state=="G0" & !Assign[,3]] <- "None"
	state[state=="G1S" & !Assign[,1]] <- "None"
	state[state=="G2M" & !Assign[,2]] <- "None"

	state2 <- factor(state, levels=c("None", "G0", "G1S", "G2M"))

	# Output
	return(list(state1, state2))
}

prolif_out <- prolif_analysis(SCE, args[2])
pData(SCE)$Proliferating <- prolif_out
prolif2_out <- prolif_analysis2(SCE, args[2])
pData(SCE)$CC_state <- prolif2_out[[1]]
pData(SCE)$CC_state_new <- prolif2_out[[2]]
print(table(pData(SCE)$CC_state, pData(SCE)$Proliferating))
print(table(pData(SCE)$CC_state, pData(SCE)$CC_state_new))
saveRDS(SCE, file=paste(args[2], "Prolif.rds", sep="_"))
