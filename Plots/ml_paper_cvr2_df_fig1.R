library(dplyr)
en_afa <- read.table(file = "/home/pokoro/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)
en_afa$gene_id <- as.character(en_afa$gene_id)
en_afa <- en_afa[,c(1,10)]
names(en_afa) <- c("gene", "cvR2")
en_afa <- subset(en_afa, cvR2 > -1)
#start df for dot plots ML v. EN
en_afa <- mutate(en_afa, model="EN",pop="AFA")

rf_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- rf_afa[,c(1,3)]
names(rf_afa) <- c("gene", "cvR2")
rf_afa <- subset(rf_afa, cvR2 > -1)
rf_afa <- mutate(rf_afa, model="RF",pop="AFA")

svr_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_svr_all_chrom.txt", header = T)
svr_afa$Gene_ID <- as.character(svr_afa$Gene_ID)
svr_afa <- svr_afa[,c(1,3)]
names(svr_afa) <- c("gene", "cvR2")
svr_afa <- subset(svr_afa, cvR2 > -1)
svr_afa <- mutate(svr_afa, model="SVR",pop="AFA")

knn_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_knn_all_chrom.txt", header = T)
knn_afa$Gene_ID <- as.character(knn_afa$Gene_ID)
knn_afa <- knn_afa[,c(1,3)]
names(knn_afa) <- c("gene", "cvR2")
knn_afa <- subset(knn_afa, cvR2 > -1)
knn_afa <- mutate(knn_afa, model="KNN",pop="AFA")

#data.frame for	boxplots
bpdf <- rbind(en_afa, rf_afa, svr_afa, knn_afa)

#data.frame for dot plots EN v. ML model
enrf <- inner_join(en_afa,rf_afa,by=c("gene","pop"))
ensvr <- inner_join(en_afa,svr_afa,by=c("gene","pop"))
enknn <- inner_join(en_afa,knn_afa,by=c("gene","pop"))

#try facet_grid(pop ~ model)
fig1df <- rbind(enrf, ensvr, enknn)

#do same for HIS
en_afa <- read.table(file = "/home/pokoro/data/mesa_models/split_mesa/results/all_chr_HIS_model_summaries.txt", header = TRUE)
en_afa$gene_id <- as.character(en_afa$gene_id)
en_afa <- en_afa[,c(1,10)]
names(en_afa) <- c("gene", "cvR2")
en_afa <- subset(en_afa, cvR2 > -1)
#start df for dot plots ML v. EN
en_afa <- mutate(en_afa, model="EN",pop="HIS")

rf_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- rf_afa[,c(1,3)]
names(rf_afa) <- c("gene", "cvR2")
rf_afa <- subset(rf_afa, cvR2 > -1)
rf_afa <- mutate(rf_afa, model="RF",pop="HIS")

svr_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_svr_all_chrom.txt", header = T)
svr_afa$Gene_ID <- as.character(svr_afa$Gene_ID)
svr_afa <- svr_afa[,c(1,3)]
names(svr_afa) <- c("gene", "cvR2")
svr_afa <- subset(svr_afa, cvR2 > -1)
svr_afa <- mutate(svr_afa, model="SVR",pop="HIS")

knn_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_knn_all_chrom.txt", header = T)
knn_afa$Gene_ID <- as.character(knn_afa$Gene_ID)
knn_afa <- knn_afa[,c(1,3)]
names(knn_afa) <- c("gene", "cvR2")
knn_afa <- subset(knn_afa, cvR2 > -1)
knn_afa <- mutate(knn_afa, model="KNN",pop="HIS")

#data.frame for boxplots
bpdf <- rbind(bpdf, en_afa, rf_afa, svr_afa, knn_afa)

#data.frame for dot plots EN v. ML model
enrf <- inner_join(en_afa,rf_afa,by=c("gene","pop"))
ensvr <- inner_join(en_afa,svr_afa,by=c("gene","pop"))
enknn <- inner_join(en_afa,knn_afa,by=c("gene","pop"))
fig1df <- rbind(fig1df, enrf, ensvr, enknn)

#do same for CAU
en_afa <- read.table(file = "/home/pokoro/data/mesa_models/split_mesa/results/all_chr_CAU_model_summaries.txt", header = TRUE)
en_afa$gene_id <- as.character(en_afa$gene_id)
en_afa <- en_afa[,c(1,10)]
names(en_afa) <- c("gene", "cvR2")
en_afa <- subset(en_afa, cvR2 > -1)
#start df for dot plots ML v. EN
en_afa <- mutate(en_afa, model="EN",pop="CAU")

rf_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- rf_afa[,c(1,3)]
names(rf_afa) <- c("gene", "cvR2")
rf_afa <- subset(rf_afa, cvR2 > -1)
rf_afa <- mutate(rf_afa, model="RF",pop="CAU")

svr_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_svr_all_chrom.txt", header = T)
svr_afa$Gene_ID <- as.character(svr_afa$Gene_ID)
svr_afa <- svr_afa[,c(1,3)]
names(svr_afa) <- c("gene", "cvR2")
svr_afa <- subset(svr_afa, cvR2 > -1)
svr_afa <- mutate(svr_afa, model="SVR",pop="CAU")

knn_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_knn_all_chrom.txt", header = T)
knn_afa$Gene_ID <- as.character(knn_afa$Gene_ID)
knn_afa <- knn_afa[,c(1,3)]
names(knn_afa) <- c("gene", "cvR2")
knn_afa <- subset(knn_afa, cvR2 > -1)
knn_afa <- mutate(knn_afa, model="KNN",pop="CAU")

#data.frame for boxplots
bpdf <- rbind(bpdf, en_afa, rf_afa, svr_afa, knn_afa)

#data.frame for dot plots EN v. ML model
enrf <- inner_join(en_afa,rf_afa,by=c("gene","pop"))
ensvr <- inner_join(en_afa,svr_afa,by=c("gene","pop"))
enknn <- inner_join(en_afa,knn_afa,by=c("gene","pop"))
fig1df <- rbind(fig1df, enrf, ensvr, enknn)

#do same for ALL
en_afa <- read.table(file = "/home/pokoro/data/mesa_models/split_mesa/results/all_chr_ALL_model_summaries.txt", header = TRUE)
en_afa$gene_id <- as.character(en_afa$gene_id)
en_afa <- en_afa[,c(1,10)]
names(en_afa) <- c("gene", "cvR2")
en_afa <- subset(en_afa, cvR2 > -1)
#start df for dot plots ML v. EN
en_afa <- mutate(en_afa, model="EN",pop="ALL")

rf_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- rf_afa[,c(1,3)]
names(rf_afa) <- c("gene", "cvR2")
rf_afa <- subset(rf_afa, cvR2 > -1)
rf_afa <- mutate(rf_afa, model="RF",pop="ALL")

svr_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header = T)
svr_afa$Gene_ID <- as.character(svr_afa$Gene_ID)
svr_afa <- svr_afa[,c(1,3)]
names(svr_afa) <- c("gene", "cvR2")
svr_afa <- subset(svr_afa, cvR2 > -1)
svr_afa <- mutate(svr_afa, model="SVR",pop="ALL")

knn_afa <- read.table(file = "/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_knn_all_chrom.txt", header = T)
knn_afa$Gene_ID <- as.character(knn_afa$Gene_ID)
knn_afa <- knn_afa[,c(1,3)]
names(knn_afa) <- c("gene", "cvR2")
knn_afa <- subset(knn_afa, cvR2 > -1)
knn_afa <- mutate(knn_afa, model="KNN",pop="ALL")

#data.frame for boxplots
bpdf <- rbind(bpdf, en_afa, rf_afa, svr_afa, knn_afa)

#data.frame for dot plots EN v. ML model
enrf <- inner_join(en_afa,rf_afa,by=c("gene","pop"))
ensvr <- inner_join(en_afa,svr_afa,by=c("gene","pop"))
enknn <- inner_join(en_afa,knn_afa,by=c("gene","pop"))
fig1df <- rbind(fig1df, enrf, ensvr, enknn)
colnames(fig1df) <- c("gene", "ENcvR2", "EN", "pop", "cvR2", "MLmodel")

fig1df <- select(fig1df, gene, ENcvR2, MLmodel, cvR2, pop)

write.table(fig1df, "fig1df.txt", quote=F, row.names=F)
write.table(bpdf, "boxplotsdf.txt", quote=F, row.names=F)
