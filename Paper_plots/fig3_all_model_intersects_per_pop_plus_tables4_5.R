#Make another Figure 3 Boxplot. Where the gene intersects per MESA training subpopulation is used per algorithm.
#e.g In AFA, all gene intersects of EN, RF, SVR, KNN

library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
"%&%" = function(a,b) paste (a,b,sep="")

#df <- NULL

algs <- c("en", "knn", "rf", "svr")
#pops <- c("AFA", "AFHI", "ALL", "CAU", "HIS")


#AFA (for common genes in AFA)
en_AFA <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/en_AFA_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
rf_AFA <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/rf_AFA_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
svr_AFA <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/svr_AFA_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
knn_AFA <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/knn_AFA_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))

AFA <- inner_join(en_AFA, rf_AFA, by = "gene")
AFA <- inner_join(AFA, svr_AFA, by = "gene")
AFA <- inner_join(AFA, knn_AFA, by = "gene")
colnames(AFA) <- c("gene", "en", "rf", "svr", "knn")


#CAU (for common genes in CAU)
en_CAU <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/en_CAU_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
rf_CAU <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/rf_CAU_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
svr_CAU <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/svr_CAU_2_METS_corr_filt.txt",
                 header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
knn_CAU <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/knn_CAU_2_METS_corr_filt.txt",
                 header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))

CAU <- inner_join(en_CAU, rf_CAU, by = "gene")
CAU <- inner_join(CAU, svr_CAU, by = "gene")
CAU <- inner_join(CAU, knn_CAU, by = "gene")
colnames(CAU) <- c("gene", "en", "rf", "svr", "knn")


#HIS (for common genes in HIS)
en_HIS <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/en_HIS_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
rf_HIS <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/rf_HIS_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
svr_HIS <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/svr_HIS_2_METS_corr_filt.txt",
                 header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
knn_HIS <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/knn_HIS_2_METS_corr_filt.txt",
                 header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))

HIS <- inner_join(en_HIS, rf_HIS, by = "gene")
HIS <- inner_join(HIS, svr_HIS, by = "gene")
HIS <- inner_join(HIS, knn_HIS, by = "gene")
colnames(HIS) <- c("gene", "en", "rf", "svr", "knn")


#AFHI (for common genes in AFHI)
en_AFHI <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/en_AFHI_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
rf_AFHI <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/rf_AFHI_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
svr_AFHI <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/svr_AFHI_2_METS_corr_filt.txt",
                 header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
knn_AFHI <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/knn_AFHI_2_METS_corr_filt.txt",
                 header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))

AFHI <- inner_join(en_AFHI, rf_AFHI, by = "gene")
AFHI <- inner_join(AFHI, svr_AFHI, by = "gene")
AFHI <- inner_join(AFHI, knn_AFHI, by = "gene")
colnames(AFHI) <- c("gene", "en", "rf", "svr", "knn")


#ALL (for common genes in ALL)
en_ALL <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/en_ALL_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
rf_ALL <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/rf_ALL_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
svr_ALL <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/svr_ALL_2_METS_corr_filt.txt",
                 header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
knn_ALL <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/knn_ALL_2_METS_corr_filt.txt",
                 header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))

ALL <- inner_join(en_ALL, rf_ALL, by = "gene")
ALL <- inner_join(ALL, svr_ALL, by = "gene")
ALL <- inner_join(ALL, knn_ALL, by = "gene")
colnames(ALL) <- c("gene", "en", "rf", "svr", "knn")


df <- NULL
for (alg in algs){
  
  pop1 <- data.frame(spearman=AFA[[alg]], Model=toupper(alg), mesa="AFA")
  pop2 <- data.frame(spearman=CAU[[alg]], Model=toupper(alg), mesa="CAU")
  pop3 <- data.frame(spearman=HIS[[alg]], Model=toupper(alg), mesa="HIS")
  pop4 <- data.frame(spearman=AFHI[[alg]], Model=toupper(alg), mesa="AFHI")
  pop5 <- data.frame(spearman=ALL[[alg]], Model=toupper(alg), mesa="ALL")
  df <- rbind(df, pop1, pop2, pop3, pop4, pop5)
}

mesa2mets <- mutate(df,mesa=factor(mesa,levels=c("AFA","HIS","CAU","AFHI", "ALL")), 
                    Model=factor(Model,levels=c("EN","RF","SVR","KNN")))
fig <- ggplot(mesa2mets, aes(x=mesa, y=spearman, fill=Model)) + geom_boxplot() + theme_classic(18) + 
  xlab("Population") + scale_y_continuous(breaks=seq(-1.0, 1.0, 0.25), limits=c(-1.0, 1.0)) + 
  ylab(expression(paste("Spearman Correlation ", rho)))
#print(fig)
tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/used_pcair_10pc_10peer_obs_exp/Fig3_all_model_intersects_per_pop.tiff", 
     width = 18, height = 14, units = 'cm', res = 300, compression = 'lzw')
fig
dev.off()

# means <- aggregate(spearman ~ Model, mesa2mets, mean)
# 
# library(ggpubr)
# ggboxplot(mesa2mets, "mesa", "spearman", fill="Model", bxp.errorbar = T, repel=T) + theme_classic(18) + 
#   xlab("Population") + scale_y_continuous(breaks=seq(-1.0, 1.0, 0.25), limits=c(-1.0, 1.0)) + 
#   ylab(expression(paste("Spearman Correlation ", rho))) + stat_summary(fun.y="mean") +
#   geom_text(data=means, aes(y=spearman))

#Make a table for the mean prediction performance of each algorithm per population. Note, this is on their intersect genes

algs <- c("en", "rf", "svr", "knn")
pops <- c("ALL","AFHI", "AFA", "CAU", "HIS")

tabut <- matrix(nrow=5, ncol=3) #ttest pvalue # Table 5 in Paper
colnames(tabut) <- c("RF", "SVR", "KNN")
rownames(tabut) <- pops

tabum <- matrix(nrow=5, ncol=4) #Mean # Table 4 in paper
colnames(tabum) <- algs
rownames(tabum) <- pops

ALL <- mutate(ALL, pop="ALL")
AFHI <- mutate(AFHI, pop="AFHI")
AFA <- mutate(AFA, pop="AFA")
CAU <- mutate(CAU, pop="CAU")
HIS <- mutate(HIS, pop="HIS")

popbox <- rbind(ALL, AFHI, AFA, CAU, HIS)


for (i in 1:length(pops)){
  for (j in 2:length(algs)){
    df <- subset(popbox, pop==pops[i])
    tt <- t.test(df[["en"]], df[[algs[j]]], paired=TRUE)
    tabut[i,j-1] <- tt$p.value
  }
  df[,c(1,6)] <- NULL
  tabum[i,] <- colMeans(df)
}

