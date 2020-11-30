#Do a cluster analysis to see how the genes cluster around themselves
#Also color the genes with if they are shared by all ml (common) and found only in each algorithm (unique)

###
en <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/en_ALL_2_METS_corr_filt.txt",
            header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
en <- subset(en, spearman>0.1)
en <- mutate(en, model="EN")

rf <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/rf_ALL_2_METS_corr_filt.txt",
            header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
rf <- subset(rf, spearman>0.1)
rf <- mutate(rf, model="RF")

svr <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/svr_ALL_2_METS_corr_filt.txt",
             header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
svr <- subset(svr, spearman>0.1)
svr <- mutate(svr, model="SVR")

knn <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/knn_ALL_2_METS_corr_filt.txt",
             header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
knn <- subset(knn, spearman>0.1)
knn <- mutate(knn, model="KNN")


#Overlap
enrf <- inner_join(en,rf,by="gene")
ensvr <- inner_join(en,svr,by="gene")
enknn <- inner_join(en,knn,by="gene")
enrfsvr <- inner_join(enrf,svr,by="gene")
enrfsvrknn <- inner_join(enrfsvr,knn,by="gene")

enall <- data.frame(gene=enrfsvrknn[,1], spearman=enrfsvrknn$spearman.x, model="Common")
rfall <- data.frame(gene=enrfsvrknn[,1], spearman=enrfsvrknn$spearman.y, model="Common")
svrall <- data.frame(gene=enrfsvrknn[,1], spearman=enrfsvrknn$spearman.x.x, model="Common")
knnall <- data.frame(gene=enrfsvrknn[,1], spearman=enrfsvrknn$spearman.y.y, model="Common")
allgenes <- data.frame(gene=enrfsvrknn[,1], model="Common")

#rfsvrknn
rfsvr <- full_join(rf,svr,by="gene")
rfsvrknn <- full_join(rfsvr,knn,by="gene")

#enrfsvr
enrf <- full_join(en,rf,by="gene")
enrfsvr <- full_join(enrf,svr,by="gene")

#enrfknn
enrf <- full_join(en,rf,by="gene")
enrfknn <- full_join(enrf,knn,by="gene")

#ensvrknn
ensvr <- full_join(en,svr,by="gene")
ensvrknn <- full_join(ensvr, knn, by="gene")

#Uniques
enonly <- anti_join(en, rfsvrknn, by="gene")
enonlyg <- enonly[,c(1,3)]
rfonly <- anti_join(rf, ensvrknn, by="gene")
rfonlyg <- rfonly[,c(1,3)]
svronly <- anti_join(svr, enrfknn, by="gene")
svronlyg <- svronly[,c(1,3)]
knnonly <- anti_join(knn, enrfsvr, by="gene")
knnonlyg <- knnonly[,c(1,3)]

gene_df <- rbind(allgenes, enonlyg, rfonlyg, svronlyg, knnonlyg)


#read in measured expressions - Normalized TPM
mets_obs1 <- fread("Z:/data/ml_paper_reviewers_corrections/mets_observed_expressions/METS_age_F_only_king_10PCS_QN_RN_expression10_peer_factor_adjusted.txt",
                   header=T)
IID <- colnames(mets_obs1)
obs_genes <- mets_obs1$gene_id
#remove decimal places from the gene_id
for (i in 1:length(obs_genes)){
  obs_genes[i] <- gsub('\\.[0-9]+','',obs_genes[i])
} #just to remove the decimal places in the gene_id
mets_obs1$gene_id <- obs_genes

#select only gene_df genes and keep in order on the normalized and peer factor corrected TPM
obs_gene_df <- inner_join(gene_df, mets_obs1, by=c("gene"="gene_id"))
#remove gene and model colums
rownames(obs_gene_df) <- gene_df$gene
obs_gene_df[,c(1,2)] <- NULL
#Do PCA on the genes expression levels to see how the genes are related. Thus, rows are genes, while sample ids are column
data1 <- as.matrix(obs_gene_df)
pca <- prcomp((data1), scale.=T)
plot(pca$x[,1], pca$x[,2])
#plot
ggdata1 <- data.frame(gene=rownames(pca$x), PC1=pca$x[,1], PC2=pca$x[,2])
ggdata1 <- inner_join(ggdata1, gene_df, by="gene")
tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/FigS9_normalized_gene_tpm_pca_clusters.tiff", width = 16, height = 16, 
     units = 'cm', res = 300, compression = 'lzw')
ggplot(data=ggdata1, aes(x=PC1, y=PC2, color=model)) + geom_point() + theme_bw(15) +labs(color="Genes")# + xlim(-0.2,50)
dev.off()


# #USe the raw TPM
# mets_obs <- fread("Z:/data/ml_paper_reviewers_corrections/mets_observed_expressions/expression_sal_WG.txt.gz",header=T) #RAW TPM
# mets_obs <- mets_obs %>% select(all_of(IID))
# obs_genes <- mets_obs$gene_id
# #remove decimal places from the gene_id
# for (i in 1:length(obs_genes)){
#   obs_genes[i] <- gsub('\\.[0-9]+','',obs_genes[i])
# } #just to remove the decimal places in the gene_id
# mets_obs$gene_id <- obs_genes
# 
# 
# #select only gene_df genes and keep in order on the normalized and peer factor corrected TPM
# obs_gene_df <- inner_join(gene_df, mets_obs, by=c("gene"="gene_id"))
# #remove gene and model colums
# rownames(obs_gene_df) <- gene_df$gene
# obs_gene_df[,c(1,2)] <- NULL
# #Do PCA on the genes expression levels to see how the genes are related. Thus, rows are genes, while sample ids are column
# data <- as.matrix(obs_gene_df)
# pca <- prcomp((data), scale.=T)
# plot(pca$x[,1], pca$x[,2])
# #plot
# ggdata <- data.frame(gene=rownames(pca$x), PC1=pca$x[,1], PC2=pca$x[,2])
# ggdata <- inner_join(ggdata, gene_df, by="gene")
# tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/FigS13_raw_gene_tpm_pca_clusters.tiff", width = 16, height = 16, 
#      units = 'cm', res = 300, compression = 'lzw')
# ggplot(data=ggdata, aes(x=PC1, y=PC2, color=model)) + geom_point() + theme_bw(15) +labs(color="Genes") + xlim(-0.2,40)
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
# # #Unique
# # en_uni_rf <- anti_join(en, rf, by="gene")#compared to RF
# # en_uni_rf_genes <- en_uni_rf$gene
# # rf_uni <- anti_join(rf, en, by="gene")
# # rf_genes <- rf_uni$gene
# # 
# # en_uni_svr <- anti_join(en, svr, by="gene") #compared to SVR
# # svr_uni <- anti_join(svr, en, by="gene")
# # 
# # en_uni_knn <- anti_join(en, knn, by="gene") #compared to KNN
# # knn_uni <- anti_join(knn, en, by="gene")