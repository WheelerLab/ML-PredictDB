#Remake figure 4 where we compared EN and other ML models performance on test data set (METS)

library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

"%&%" = function(a,b) paste (a,b,sep="")

df <- NULL

algs <- c("knn", "rf", "svr")
pops <- c("AFA", "AFHI", "ALL", "CAU", "HIS")

for (pop in pops){
  
  en <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/en_" %&% pop %&% "_2_METS_corr_filt.txt",
              header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
  colnames(en) <- c("gene", "ENrho")
  
  for (alg in algs){
    
    ml <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/" %&% alg %&% "_" %&% pop %&% "_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
    colnames(ml) <- c("gene", "rho")
    
    ml <- inner_join(en, ml, by = "gene")
    ml <- mutate(ml, MLmodel=toupper(alg),pop=pop)
    df <- rbind(df, ml)
    
  }
}

fwrite(df, file="Z:/data/ml_paper_reviewers_corrections/paper_figs/fig4_rhodf.txt", 
       row.names=F, quote=F, sep="\t")

#Plot
data <- mutate(df,pop=factor(pop,levels=c("ALL","AFHI","AFA","HIS","CAU")),MLmodel=factor(MLmodel,levels=c("RF","SVR","KNN")))

# fig <- ggplot(data,aes(x=ENrho,y=rho)) + geom_point(shape=".") + geom_abline(intercept=0, slope=1, color="blue") +
#   geom_smooth(method="lm", color="red", lwd=0.5) + xlim(-0.5,1) + ylim(-0.5,1) + 
#   xlab(expression(paste("Elastic Net ", rho))) + ylab(expression(paste("Other Machine Learning Model ", rho))) +
#   theme_bw(16) + facet_grid(pop~MLmodel) #xlim(0.1,1) + ylim(0.1,1) + old
# 
# #print(fig)

tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/used_pcair_10pc_10peer_obs_exp/Fig4.tiff", 
     width = 16, height = 20, units = 'cm', res = 300, compression = 'lzw')

ggscatter(data, x="ENrho", y="rho", size=0.8, add = "reg.line", add.params = list(color="red", size=1), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 2.5) + geom_abline(intercept=0, slope=1, color="blue", lwd=1) + 
  xlim(-0.5,1) + ylim(-0.5,1) + xlab(expression(paste("Elastic Net ", rho))) + 
  ylab(expression(paste("Other Machine Learning Model ", rho))) + theme_bw(16) + facet_grid(pop~MLmodel)

#fig
dev.off()




##############
#Table 6
##Make a table with the p-values of the t-test between EN vs ML across each sub population

tb <- matrix(nrow=5,ncol=3)

pops <- c("ALL", "AFHI", "AFA", "CAU", "HIS")
algs <- c("RF", "SVR", "KNN")

colnames(tb) <- algs
rownames(tb) <- pops


for (i in 1:length(pops)){
  
  for (j in 1:length(algs)){
    
    ml <- subset(df, df$MLmodel==algs[j] & df$pop==pops[i])
    tt <- t.test(ml[["ENrho"]], ml[["rho"]], paired=TRUE)
    tb[i,j] <- tt$p.value
  }
}

#Take the mean as well for table 6
# tm <- matrix(nrow=5,ncol=6)
# colnames(tm) <- c("EN", "RF", "EN", "SVR", "EN", "KNN")

ten <- matrix(nrow=5,ncol=3)
colnames(ten) <- algs
rownames(ten) <- pops

tml <- matrix(nrow=5,ncol=3)
colnames(tml) <- algs
rownames(tml) <- pops


for (i in 1:length(pops)){
  
  for (j in 1:length(algs)){
    
    ml <- subset(df, df$MLmodel==algs[j] & df$pop==pops[i])
    ten[i,j] <- mean(ml[["ENrho"]]) 
    tml[i,j] <- mean(ml[["rho"]])
  }
}


#####Table 7

###
en <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/en_ALL_2_METS_corr_filt.txt",
            header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
en <- subset(en, spearman>0.1)

rf <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/rf_ALL_2_METS_corr_filt.txt",
            header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
rf <- subset(rf, spearman>0.1)

svr <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/svr_ALL_2_METS_corr_filt.txt",
             header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
svr <- subset(svr, spearman>0.1)

knn <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/knn_ALL_2_METS_corr_filt.txt",
             header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
knn <- subset(knn, spearman>0.1)

#Overlap
enrf <- inner_join(en,rf,by="gene")
ensvr <- inner_join(en,svr,by="gene")
enknn <- inner_join(en,knn,by="gene")

#Unique
en_uni_rf <- anti_join(en, rf, by="gene")#compared to RF
en_uni_rf_genes <- en_uni_rf$gene
rf_uni <- anti_join(rf, en, by="gene")
rf_genes <- rf_uni$gene

en_uni_svr <- anti_join(en, svr, by="gene") #compared to SVR
svr_uni <- anti_join(svr, en, by="gene")

en_uni_knn <- anti_join(en, knn, by="gene") #compared to KNN
knn_uni <- anti_join(knn, en, by="gene")

#Uniques of each ML model not in EN combined
uni_of_rfsvr <- inner_join(rf_uni, svr_uni, by="gene")
uni_of_rfsvrknn <- inner_join(uni_of_rfsvr, knn_uni, by="gene")
ml_genes <- uni_of_rfsvrknn$gene #
#write.table(ml_genes, "/Users/okoro/OneDrive/Desktop/mlgenes.txt", quote=F, row.names=F, col.names=F)

#Uniques of EN model combined
en_uni_of_rfsvr <- inner_join(en_uni_rf, en_uni_svr, by="gene")
en_uni_of_rfsvrknn <- inner_join(en_uni_of_rfsvr, en_uni_knn, by="gene")
en_genes <- en_uni_of_rfsvrknn$gene

#Check expression levels of these genes found only in ML models

mets_obs <- fread("Z:/data/ml_paper_reviewers_corrections/mets_observed_expressions/METS_age_F_only_king_10PCS_QN_RN_expression10_peer_factor_adjusted.txt",
                  header=T)
obs_genes <- mets_obs$gene_id
#remove decimal places from the gene_id
for (i in 1:length(obs_genes)){
  obs_genes[i] <- gsub('\\.[0-9]+','',obs_genes[i])
} #just to remove the decimal places in the gene_id
mets_obs$gene_id <- obs_genes

#Get ML genes expression levels
ml_exp_levels <- subset(mets_obs, gene_id %in% ml_genes) #Keep only those ML unique genes
ml_exp_levels <- as.data.frame(t(ml_exp_levels))
colnames(ml_exp_levels) <- ml_exp_levels[1,]
ml_exp_levels <- ml_exp_levels[-1,]
ml_exp_levels <- as.data.frame(sapply(ml_exp_levels, as.numeric))
dfml <- NULL
for (i in 1:length(ml_genes)){
  ml <- data.frame(exp=ml_exp_levels[[ml_genes[i]]], gene=ml_genes[i])
  dfml <- rbind(dfml,ml)
}
dfml <- mutate(dfml, Type="ML Unique")
ggplot(dfml,aes(x=exp, colour=gene)) + geom_density() + theme(legend.position = "none")

#Get en genes expression levels
en_exp_levels <- subset(mets_obs, gene_id %in% en_genes) #Keep only those en unique genes
en_exp_levels <- as.data.frame(t(en_exp_levels))
colnames(en_exp_levels) <- en_exp_levels[1,]
en_exp_levels <- en_exp_levels[-1,]
en_exp_levels <- as.data.frame(sapply(en_exp_levels, as.numeric))
dfen <- NULL
for (i in 1:length(en_genes)){
  ml <- data.frame(exp=en_exp_levels[[en_genes[i]]], gene=en_genes[i])
  dfen <- rbind(dfen,ml)
}
dfen <- mutate(dfen, Type="EN Unique")
ggplot(dfen,aes(x=exp, colour=gene)) + geom_density() + theme(legend.position = "none")

df <- rbind(dfml,dfen)
p <- ggplot(df,aes(x=exp, colour=gene)) + geom_density() + xlab("Gene Expression Levels") +
  theme_classic(16) + facet_wrap(~Type)
p + theme(legend.position = "none")

##Compare the EN and ML expression levels

ml_mean <- colMeans(ml_exp_levels)
en_mean <- colMeans(en_exp_levels)
t.test(ml_mean, en_mean)

# #or test for homogeneity of variances with fligner test
# grp1 <- data.frame(mean=ml_mean, group="grp1")
# grp2 <- data.frame(mean=en_mean, group="grp2")
# grp_df <- rbind(grp1,grp2)
# fligner.test(mean ~ group, data=grp_df) #Assumes non normality
# bartlett.test(mean ~ group, data=grp_df) #Assumes normality

ml_var <- sapply(ml_exp_levels, var)
en_var <- sapply(en_exp_levels, var)
t.test(ml_var, en_var)

mldf <- data.frame(variance=ml_var, Model="ML")
endf <- data.frame(variance=en_var, Model="EN")
df <- rbind(mldf,endf)

ggplot(df, aes(x=Model, y=variance)) + geom_violin()
ggplot(df, aes(x=Model, y=variance)) + geom_boxplot()



##Get rf genes expression levels
rf_exp_levels <- subset(mets_obs, gene_id %in% rf_genes) #Keep only those rf unique genes
rf_exp_levels <- as.data.frame(t(rf_exp_levels))
colnames(rf_exp_levels) <- rf_exp_levels[1,]
rf_exp_levels <- rf_exp_levels[-1,]
rf_exp_levels <- as.data.frame(sapply(rf_exp_levels, as.numeric))
dfrf <- NULL
for (i in 1:length(rf_genes)){
  ml <- data.frame(exp=rf_exp_levels[[rf_genes[i]]], gene=rf_genes[i])
  dfrf <- rbind(dfrf,ml)
}
dfrf <- mutate(dfrf, Type="RF Unique")
ggplot(df,aes(x=exp, colour=gene)) + geom_density() + theme(legend.position = "none")

df <- rbind(dfrf,endf)

#Get en_uni_rf gen_uni_rfes expression levels
en_uni_rf_exp_levels <- subset(mets_obs, gene_id %in% en_uni_rf_genes) #Keep only those en_uni_rf unique gen_uni_rfes
en_uni_rf_exp_levels <- as.data.frame(t(en_uni_rf_exp_levels))
colnames(en_uni_rf_exp_levels) <- en_uni_rf_exp_levels[1,]
en_uni_rf_exp_levels <- en_uni_rf_exp_levels[-1,]
en_uni_rf_exp_levels <- as.data.frame(sapply(en_uni_rf_exp_levels, as.numeric))
endf <- NULL
for (i in 1:length(en_uni_rf_genes)){
  ml <- data.frame(exp=en_uni_rf_exp_levels[[en_uni_rf_genes[i]]], gene=en_uni_rf_genes[i])
  endf <- rbind(endf,ml)
}
endf <- mutate(endf, Type="EN Unique")
ggplot(df,aes(x=exp, colour=gene)) + geom_density() + theme(legend.position = "none")

p <- ggplot(df,aes(x=exp, colour=gene)) + geom_density() + xlab("Gene Expression Levels") +
  theme_classic(16) + facet_wrap(~Type)
p + theme(legend.position = "none")


##Compare the EN and RF expression levels

rf_mean <- colMeans(rf_exp_levels)
enrf_mean <- colMeans(en_uni_rf_exp_levels)
t.test(rf_mean, enrf_mean)

rf_var <- sapply(rf_exp_levels, var)
enrf_var <- sapply(en_uni_rf_exp_levels, var)
t.test(rf_var, enrf_var)

rfdf <- data.frame(variance=rf_mean, Model="ML")
enrfdf <- data.frame(variance=enrf_mean, Model="EN")
df <- rbind(rfdf,enrfdf)

ggplot(df, aes(x=Model, y=variance)) + geom_violin() + geom_boxplot()
ggplot(df, aes(x=Model, y=variance)) + geom_boxplot()






##Plots of the RHO vs TPM

mets_obs1 <- fread("Z:/data/ml_paper_reviewers_corrections/mets_observed_expressions/METS_age_F_only_king_10PCS_QN_RN_expression10_peer_factor_adjusted.txt",
                   header=T)
IID <- colnames(mets_obs1)

mets_obs <- fread("Z:/data/ml_paper_reviewers_corrections/mets_observed_expressions/expression_sal_WG.txt.gz",header=T)
mets_obs <- mets_obs %>% select(all_of(IID))
obs_genes <- mets_obs$gene_id
#remove decimal places from the gene_id
for (i in 1:length(obs_genes)){
  obs_genes[i] <- gsub('\\.[0-9]+','',obs_genes[i])
} #just to remove the decimal places in the gene_id
mets_obs$gene_id <- obs_genes


algs <- c("en", "knn", "rf", "svr")
pops <- c("AFA", "AFHI", "ALL", "CAU", "HIS")


gex <- mets_obs[,c(2:77)]
gex <- data.frame(gene=obs_genes, meanTPM=rowMeans(gex))

df <- NULL

for (pop in pops){
  
  for (alg in algs){
    
    ml <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/" %&% alg %&% "_" %&% pop %&% "_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
    #colnames(ml) <- c("gene", "rho")
    
    ml <- inner_join(ml, gex, by = "gene")
    ml <- mutate(ml, MLmodel=toupper(alg),pop=pop)
    df <- rbind(df, ml)
    
  }
}


tiff("/Users/okoro/OneDrive/Desktop/rhobytmp_raw.tiff", width = 24, height = 20, units = 'cm', res = 300, compression = 'lzw')
ggscatter(df, x="meanTPM", y="spearman", size=0.8, add = "reg.line", add.params = list(color="red", size=1), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 2.5) + geom_abline(intercept=0, slope=1, color="blue", lwd=1) +
  xlab("Expression Levels") + ylab("Spearman") + xlim(0,500) + theme_bw(16) + facet_grid(pop~MLmodel)
dev.off()