
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(viridis)
"%&%" = function(a,b) paste (a,b,sep="")

algs <- c("en", "knn", "rf", "svr")
pops <- c("AFA", "CAU", "HIS", "AFHI", "ALL")

model <- c("EN", "EN", "EN", "RF", "RF", "RF", "SVR", "SVR", "SVR", "KNN", "KNN", "KNN")
#threshold <- c("all","p>0", "p>0.1")

#tabcol <- c("population", "Model", "Threshold", "Count")

df <- matrix(nrow=5, ncol=12)
colnames(df) <- model
rownames(df) <- pops

for (i in 1:length(pops)){
  
  ml <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/en_" %&% pops[i] %&% "_2_METS_corr_filt.txt",
              header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene", "spearman"))
  df[i,1] <- nrow(ml)
  ml <- subset(ml, spearman>0)
  df[i,2] <- nrow(ml)
  ml <- subset(ml, spearman>0.1)
  df[i,3] <- nrow(ml)
}


for (i in 1:length(pops)){
  
  ml <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/rf_" %&% pops[i] %&% "_2_METS_corr_filt.txt",
              header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene", "spearman"))
  df[i,4] <- nrow(ml)
  ml <- subset(ml, spearman>0)
  df[i,5] <- nrow(ml)
  ml <- subset(ml, spearman>0.1)
  df[i,6] <- nrow(ml)
}


for (i in 1:length(pops)){
  
  ml <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/svr_" %&% pops[i] %&% "_2_METS_corr_filt.txt",
              header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene", "spearman"))
  df[i,7] <- nrow(ml)
  ml <- subset(ml, spearman>0)
  df[i,8] <- nrow(ml)
  ml <- subset(ml, spearman>0.1)
  df[i,9] <- nrow(ml)
}


for (i in 1:length(pops)){
  
  ml <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/knn_" %&% pops[i] %&% "_2_METS_corr_filt.txt",
              header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene", "spearman"))
  df[i,10] <- nrow(ml)
  ml <- subset(ml, spearman>0)
  df[i,11] <- nrow(ml)
  ml <- subset(ml, spearman>0.1)
  df[i,12] <- nrow(ml)
}


afa <- data.frame(population=rep("AFA", 12), Model=model, Threshold=rep(threshold, 4), Count=df[1,])

cau <- data.frame(population=rep("CAU", 12), Model=model, Threshold=rep(threshold, 4), Count=df[2,])

his <- data.frame(population=rep("HIS", 12), Model=model, Threshold=rep(threshold, 4), Count=df[3,])

afhi <- data.frame(population=rep("AFHI", 12), Model=model, Threshold=rep(threshold, 4), Count=df[4,])

all <- data.frame(population=rep("ALL", 12), Model=model, Threshold=rep(threshold, 4), Count=df[5,])

table2 <- rbind(afa, cau, his, afhi, all)

table2 <- mutate(table2,mesa=factor(population,levels=c("AFA","HIS","CAU","AFHI", "ALL")), 
                    Model=factor(Model,levels=c("EN","RF","SVR","KNN")))

fig <- ggplot(data=table2, aes(x=mesa, y=Count, fill=Threshold)) + 
  geom_bar(stat="identity", position=position_dodge()) + scale_fill_brewer(palette="Paired", aesthetics = "colour") + 
  theme_minimal(20)  + xlab("Population") + ylab("Genes") + scale_fill_viridis(discrete=T, option="A")+ facet_wrap(~Model) + 
  scale_fill_discrete(name =expression(rho), labels = c("all", ">0", ">0.1"))
#print(fig)
tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/Fig2.tiff", width = 22, height = 18, 
     units = 'cm', res = 300, compression = 'lzw')
fig
dev.off()
