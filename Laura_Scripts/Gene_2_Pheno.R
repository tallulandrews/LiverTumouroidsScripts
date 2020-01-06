obj <- readRDS("/warehouse/team218_wh01/tallulah/TCGA/All_TCGA_Transcriptome.rds")
my_stage <- as.character(obj$ann[,3])
my_stage[my_stage == "stage:stage i, "] = 1
my_stage[my_stage == "stage:stage ii, "] = 2
my_stage[my_stage == "stage:stage iii, "] = 3
my_stage[my_stage == "stage:stage iiia, "] = 3
my_stage[my_stage == "stage:stage iiib, "] = 3
my_stage[my_stage == "stage:stage iiic, "] = 3
my_stage[my_stage == "stage:stage iv, "] = 4
my_stage[my_stage == "stage:stage iva, "] = 4
my_stage[my_stage == "stage:stage ivb, "] = 4
my_stage[my_stage == "stage:not reported, "] = NA

my_death_time <- as.character(obj$ann[,5])
my_death_time<-matrix(unlist(strsplit(my_death_time, ":")), ncol=2, byrow=T)
my_death_time<-sub(", ","", my_death_time[,2])
my_death_time[my_death_time =="null"] = NA
my_death_time<-as.numeric(my_death_time)

my_followup_time <- as.character(obj$ann[,7])
my_followup_time<-matrix(unlist(strsplit(my_followup_time, ":")), ncol=2, byrow=T)
my_followup_time<-sub(", ","", my_followup_time[,2])
my_followup_time[my_followup_time =="null"] = NA
my_followup_time<-as.numeric(my_followup_time)


norm<- colSums(obj$data)
normed <- t(t(obj$data)/norm*1000000)
normed <- log(normed+1)/log(2)

profiles <- read.delim("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/WeightedMean_Profiles.csv", sep=",", header=T)

profile_score <- (profiles$CCA1_4+profiles$CCA1_1+profiles$HCC6_4+profiles$HCC10_1)-rowSums(profiles[,!(colnames(profiles) %in% c("CCA1_4","CCA1_1","HCC6_4","HCC10_1"))])

# Laura-esque
normal <- apply(profiles[,!(colnames(profiles) %in% c("CCA1_4","CCA1_1","HCC6_4","HCC10_1"))], 1, mean)
stemy <- apply(profiles[,(colnames(profiles) %in% c("CCA1_4","CCA1_1","HCC6_4","HCC10_1"))], 1, min)
#profile_score2 <- (stemy-normal)/abs(normal) # no work b/c of vals close to 0.
profile_score2<- stemy-normal

profile_score = profile_score2

profiles <- profiles[order(-1*profile_score),]

# Consistent gene set #
normed <- normed[match(rownames(profiles),rownames(normed)),]
profiles<-profiles[rownames(profiles) %in% rownames(normed),]
normed<-normed[rownames(normed) %in% rownames(profiles),]

# General score #
stem_score_weights <- (profiles$CCA1_4+profiles$CCA1_1+profiles$HCC6_4+profiles$HCC10_1)/4
pat_stem_score <- t(normed) %*% stem_score_weights
png("Stage_vs_Score.png", width=6, height=6, units="in", res=300)
boxplot(pat_stem_score~my_stage, notch=TRUE, xlab="Stage", ylab="Stem-ness Score", col=c("#4575b4","#fee090","#d73027","#a50026"))
p = wilcox.test(pat_stem_score[my_stage=="1"], pat_stem_score[my_stage=="3"])
lines(c(1-0.15,3+0.15), c(11000,11000), lwd=2)
text(2, 11000, paste("p =", round(p$p.value, digits=2)), font=2, pos=3)
dev.off()

# Exteme score #
pat_extreme_score <- colSums(normed[1:100,]) - colSums(normed[length(normed[,1])-1:100,]) # tried 20, 50, 100


hi <- pat_stem_score > median(pat_stem_score)
png("Score_vs_Death.png", width=6, height=6, units="in", res=300)
boxplot(my_death_time~hi, notch=TRUE, names=c("Low", "High"), xlab="Stem-ness Score", ylab="Time to Death (days)", col=c("#4575b4","#d73027"))
p <- wilcox.test(my_death_time[hi], my_death_time[!hi])$p.value
lines(c(1-0.15,2+0.15), c(1600,1600), lwd=2)
text(1.5, 1600, paste("p =",round(p, digits=4)), font=2, pos=3)
dev.off()

hi <- pat_extreme_score > median(pat_extreme_score)
png("Extreme_Score_vs_Death.png", width=6, height=6, units="in", res=300)
boxplot(my_death_time~hi, notch=TRUE, names=c("Low", "High"), xlab="Stem-ness Score", ylab="Time to Death (days)", col=c("#4575b4","#d73027"))
p <- wilcox.test(my_death_time[hi], my_death_time[!hi])$p.value
lines(c(1-0.15,2+0.15), c(1600,1600), lwd=2)
text(1.5, 1600, paste("p =",round(p, digits=4)), font=2, pos=3)
dev.off()


# Survival Test
require("survival")
dead = !is.na(my_death_time);
times <- my_followup_time;
times[dead] <- my_death_time[dead]

predictor <- as.numeric(pat_stem_score > median(pat_stem_score))
#predictor <- as.numeric(pat_extreme_score > median(pat_extreme_score))
predictor <- factor(predictor)

keep <- (dead | !is.na(my_followup_time)) & !is.na(predictor);
predictor<-predictor[keep];
times<-times[keep]
dead <- dead[keep]

surv_obj <- Surv(times, dead)
surv_test <- survdiff(surv_obj~predictor) # Discrete predictor
curves <- survfit(surv_obj ~ predictor)
png("Stem-ness_Survival_Curve.png", width=6, height=6, units="in", res=300)
plot(curves, col=c("blue", "red"), mark.time=TRUE, conf.int=TRUE, xlab="Days", ylab="Survival")
legend("bottomleft", c("High", "Low"), col=c("red","blue"), title="Stem-ness Score", lty=1, bty="n")
legend("topright", bty="n", paste("p =", signif(1-pchisq(surv_test$chisq, df=1), digits=1)))
dev.off()


predictor<-pat_stem_score
predictor<-predictor[keep];
surv_test1 <- coxph(surv_obj~predictor) # Continuous predictor?
surv_test2 <- survreg(surv_obj~predictor, dist="logistic") # Continuous predictor

sng_gene_test <- function(x) {
	require("survival")
	dead = !is.na(my_death_time);
	times <- my_followup_time;
	times[dead] <- my_death_time[dead]
	predictor<-x
	keep <- (dead | !is.na(my_followup_time)) & !is.na(predictor);
	predictor<-predictor[keep];
	times<-times[keep]
	dead <- dead[keep]

	surv_obj <- Surv(times, dead)
	surv_test1 <- coxph(surv_obj~predictor) # Cox proportional hazards model 
#	surv_test2 <- survreg(surv_obj~predictor, dist="logistic") # Continuous predictor
	p = pchisq(surv_test1$wald.test, df=1, lower.tail=FALSE)
	eff <- surv_test1$coefficients[1]
	return(c(eff,p))
}

out <- t(apply(normed, 1, sng_gene_test))
q.vals <- p.adjust(out[,2], method="fdr")

source("~/R-Scripts/Ensembl_Stuff.R")
ForLaura <- profiles[,c(6,17,9,5,18,8,2,15,12,11,9,3,14,13,1,16,4,7)]
ForLaura <- data.frame(ForLaura)
ForLaura$profile_score <- profile_score
ForLaura$effect_size =out[,1] 
ForLaura$p.value=out[,2]
ForLaura$q.value=q.vals
ForLaura$Symbol <- map_symbol_ensg(rownames(ForLaura), is.org="Hsap", is.name="ensg")

write.table(ForLaura, file="WeightedMean_Profile_PlusTCGA.csv", sep=",")






# gene-wise #

#sng_gene_test <- function(x) {
#	tmp = summary(glm(my_death_time~x+my_stage))$coeff
#	return(c(tmp[2,1], tmp[2,4]))
#}

#out <- apply(normed, 1, sng_gene_test)
#q.vals <- p.adjust(t(out)[,2], method="fdr")
#glm_tab <- cbind(t(out), q.vals)
#colnames(glm_tab) <- c("eff_size", "p.value", "q.value")
#require("gplots")
#heatmap.2(as.matrix(profiles[glm_tab[,3] < 0.05 & glm_tab[,1] < -100,]), col=rev(colorRampPalette(c("#d73027","white","#4575b4"))(20)), symbreaks=TRUE, trace="none")

ForLaura <- profiles[,c(6,17,9,5,18,8,2,15,12,11,9,3,14,13,1,16,4,7)]
ForLaura <- cbind(ForLaura, glm_tab)

source("~/R-Scripts/Ensembl_Stuff.R")
ForLaura <- data.frame(ForLaura)
ForLaura$Symbol <- map_symbol_ensg(rownames(ForLaura), is.org="Hsap", is.name="ensg")


#sum(ForLaura$q.value[1:100] < 0.05 & ForLaura$eff_size[1:100] < 0)


my_row_mean_aggregate <- function(mat, groups) {
        # Much faster version
        MAT <- as.matrix(mat)
        x <- split(seq(ncol(MAT)), groups)
        result <- sapply(x, function(a) rowMeans(MAT[,a]))

        return(result);
}

survival_analysis_discrete_predictor <- function(predictor, death_time, last_followup){
#	https://www.openintro.org/download.php?file=survival_analysis_in_R&referrer=/stat/surv.php
	require("survival")
	Surv # Create response variable
	survdiff # Test for difference in survival curves
	survfit # Create survival curve

	dead = !is.na(death_time);
	times <- last_followup;
	times[dead] <- death_time[dead]

	keep <- dead | !is.na(last_followup);
	times<-times[keep]
	dead <- dead[keep]
	predictor<-predictor[keep];
	
	surv_obj <- Surv(times, dead)

	surv_test <- survdiff(surv_obj~predictor) # Discrete predictor
	surv_test <- coxph(surv_obj~predictor) # Continuous predictor?
	surv_test <- survreg(surv_obj~predictor, dist="weibull") # Continuous predictor

	curves <- survfit(surv_obj ~ predictor)
	plot(curves, main="Survival Plot", xlab="Time (days)", ylab="Survival", col=c("blue", "red"))
	#curves$strata

	tabSurv<-survdiff(ttime  ~ c(rep("top", nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN)) ))
	plot(survfit(ttime ~ c(rep("low", nrow(cfu_onlyTOP)), rep("high", nrow(cfu_onlyDOWN)))), col = c("green", "red"),main= titlePlot,xlab="Days",ylab="Survival")

}
