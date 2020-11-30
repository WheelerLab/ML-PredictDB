library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

"%&%" = function(a,b) paste (a,b,sep="")

mets_obs1 <- fread("Z:/data/ml_paper_reviewers_corrections/mets_observed_expressions/METS_age_F_only_king_10PCS_QN_RN_expression10_peer_factor_adjusted.txt",
                   header=T)
IID <- colnames(mets_obs1)

mets_obs <- fread("Z:/data/ml_paper_reviewers_corrections/mets_observed_expressions/expression_sal_WG.txt.gz",header=T) #RAW TPM
mets_obs <- mets_obs %>% select(all_of(IID))
obs_genes <- mets_obs$gene_id
#remove decimal places from the gene_id
for (i in 1:length(obs_genes)){
  obs_genes[i] <- gsub('\\.[0-9]+','',obs_genes[i])
} #just to remove the decimal places in the gene_id
mets_obs$gene_id <- obs_genes

data <- NULL

algs <- c("en", "knn", "rf", "svr")
pops <- c("AFA", "AFHI", "ALL", "CAU", "HIS")

# for (pop in pops){
# 
#   for (alg in algs){
# 
#     ml <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/" %&% alg %&% "_" %&% pop %&% "_2_METS_corr.txt",
#                 header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
#     ml <- mutate(ml, model=toupper(alg),pop=pop)
#     data <- rbind(data, ml)
# 
#   }
# }

#en <- subset(data, model=="EN" & pop=="ALL")
gex <- mets_obs[,c(2:77)]
gex <- data.frame(gene=obs_genes, meanTPM=rowMeans(gex))
#en <- inner_join(en, gex, by="gene")

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


tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/FigS7_rho_vs_raw-tpm.tiff", width = 24, height = 20, 
     units = 'cm', res = 300, compression = 'lzw')
ggscatter(df, x="meanTPM", y="spearman", size=0.8, add = "reg.line", add.params = list(color="red", size=1), conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 2.5) +
  xlab("Mean Raw Expression Levels") + ylab("Model Performance") + xlim(0,400) + theme_bw(16) + facet_grid(pop~MLmodel)
dev.off()





# #### WIth Normalized TPM
# mets_obs1 <- mets_obs1 %>% select(all_of(IID))
# # obs_genes <- mets_obs1$gene_id
# # #remove decimal places from the gene_id
# # for (i in 1:length(obs_genes)){
# #   obs_genes[i] <- gsub('\\.[0-9]+','',obs_genes[i])
# # } #just to remove the decimal places in the gene_id
# mets_obs1$gene_id <- obs_genes
# 
# data <- NULL
# 
# algs <- c("en", "knn", "rf", "svr")
# pops <- c("AFA", "AFHI", "ALL", "CAU", "HIS")
# 
# gex <- mets_obs1[,c(2:77)]
# gex <- data.frame(gene=obs_genes, meanTPM=rowMeans(gex))
# #en <- inner_join(en, gex, by="gene")
# 
# data <- NULL
# 
# for (pop in pops){
#   
#   for (alg in algs){
#     
#     ml <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/" %&% alg %&% "_" %&% pop %&% "_2_METS_corr_filt.txt",
#                 header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
#     #colnames(ml) <- c("gene", "rho")
#     
#     ml <- inner_join(ml, gex, by = "gene")
#     ml <- mutate(ml, MLmodel=toupper(alg),pop=pop)
#     data <- rbind(data, ml)
#     
#   }
# }
# 
# tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/FigS8_rho_vs_normalized-tpm.tiff", width = 24, height = 20, 
#      units = 'cm', res = 300, compression = 'lzw')
# ggscatter(data, x="meanTPM", y="spearman", size=0.8, add = "reg.line", add.params = list(color="red", size=1), conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 2.5) + geom_abline(intercept=0, slope=1, color="blue", lwd=1) +
#   xlab("Normalized Expression Levels") + ylab("Model Performance") + theme_bw(16) + facet_grid(pop~MLmodel)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##check for R2 vs exp levels of MESA
# #But first adjust the mesa exp levels
# r2 <- fread(file="Z:/ml_paper_figs/fig1df.txt",header=T,stringsAsFactors = F)
# #colnames(r2)[1] <- "gene"
# 
# #For adjusting mesa expression
# #this is the function
# adjust_for_covariates <- function(expression_vec, cov_df) {
#   combined_df <- cbind(expression_vec, cov_df)
#   expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
#   expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
#   expr_resid
# }
# 
# gexr <- fread("Z:/data/mesa_models/mesa_train_data/HIS_PF10.txt.gz", header=T)
# gexr_IID <- colnames(gexr)[2:ncol(gexr)]
# gexr <- as.data.frame(t(gexr))
# colnames(gexr) <- gexr[1,]
# gexr <- gexr[-1,]
# gexr <- as.data.frame(sapply(gexr, as.numeric))
# gexr_genes <- colnames(gexr) #take the genes names in the expression dataframe, so I can use it to iteratively adjust the expression with pc
# gexr <- cbind(gexr_IID, gexr)
# gexr$gexr_IID <- as.character(gexr$gexr_IID)
# 
# pc <- fread("Z:/data/mesa_models/mesa_train_data/HIS_3_PCs.txt", header=T) %>% select(c("IID","PC1","PC2","PC3"))
# pc_IID <- data.frame(IID=pc$IID)
# pc_IID$IID <- as.character(pc_IID$IID)
# 
# gexr <- inner_join(pc_IID, gexr, by = c("IID"="gexr_IID"))
# 
# cov_df <- pc[,c(2:4)]#take only the 3 pcs since the pc sample IID is in same order with the expression
# 
# 
# #Adjust the predicted expression
# mesa_adj <- NULL
# for (gene in gexr_genes){
#   gene_expr <- gexr[[gene]]
#   mesa_adj <- cbind(mesa_adj,adjust_for_covariates(gene_expr, cov_df))
# }
# mesa_adj <- as.data.frame(mesa_adj)
# names(mesa_adj) <- gexr_genes #put b
# mesa_tpm <- data.frame(gene=gexr_genes, meanTPM=colMeans(mesa_adj))
# 
# ml <- subset(r2, r2$MLmodel=="RF" & r2$pop=="HIS")
# en <- subset(r2, r2$MLmodel=="RF" & r2$pop=="HIS")
# 
# hismesa <- inner_join(mesa_tpm, his, by="gene")
# ggscatter(hismesa, x="meanTPM", y="ENcvR2", size=0.8, add = "reg.line", add.params = list(color="red", size=1), conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 2.5) + geom_abline(intercept=0, slope=1, color="blue", lwd=1)
# 
# 
# knn <- subset(r2, r2$MLmodel=="KNN" & r2$pop=="ALL")
# 
# #knn_twas <- fread("Z:/data/twas_mesa/ml_pred/transformed_ALL_allchrom_knn.txt", header=T)
