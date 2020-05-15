#join the best grid results for all pops
#make a pop column
#remove the time column, and combine the table
library(data.table)
library(dplyr)


#RF
afa_rf <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_rf_all_chrom.txt", header=T)
afa_rf$time.s. <- NULL
afa_rf <- mutate(afa_rf, pop="AFA")
afa_rf <- afa_rf[order(afa_rf$CV_R2),]

cau_rf <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_rf_all_chrom.txt", header=T)
cau_rf$time.s. <- NULL
cau_rf <- mutate(cau_rf, pop="CAU")
cau_rf <- cau_rf[order(cau_rf$CV_R2),]

his_rf <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_rf_all_chrom.txt", header=T)
his_rf$time.s. <- NULL
his_rf <- mutate(his_rf, pop="HIS")
his_rf <- his_rf[order(his_rf$CV_R2),]

all_rf <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header=T)
all_rf$time.s. <- NULL
all_rf <- mutate(all_rf, pop="ALL")
all_rf <- all_rf[order(all_rf$CV_R2),]

afhi_rf <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_rf_all_chrom.txt", header=T)
afhi_rf$time.s. <- NULL
afhi_rf <- mutate(afhi_rf, pop="AFHI")
afhi_rf <- afhi_rf[order(afhi_rf$CV_R2),]

#join
rf <- rbind(afa_rf, cau_rf, his_rf, afhi_rf, all_rf)
names(rf) <- c("Gene_ID", "Gene_Name", "Cross-validated R2", "RF Trees", "Population")
fwrite(rf, file="Z:/ml_paper_figs/RF_optimum_hyperparameter.txt", row.names=F, quote=F, sep="\t")

#SVR
afa_svr <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_svr_all_chrom.txt", header=T)
afa_svr$time.s. <- NULL
afa_svr <- mutate(afa_svr, pop="AFA")
afa_svr <- afa_svr[order(afa_svr$CV_R2),]

cau_svr <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_svr_all_chrom.txt", header=T)
cau_svr$time.s. <- NULL
cau_svr <- mutate(cau_svr, pop="CAU")
cau_svr <- cau_svr[order(cau_svr$CV_R2),]

his_svr <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_svr_all_chrom.txt", header=T)
his_svr$time.s. <- NULL
his_svr <- mutate(his_svr, pop="HIS")
his_svr <- his_svr[order(his_svr$CV_R2),]

all_svr <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header=T)
all_svr$time.s. <- NULL
all_svr <- mutate(all_svr, pop="ALL")
all_svr <- all_svr[order(all_svr$CV_R2),]

afhi_svr <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_svr_all_chrom.txt", header=T)
afhi_svr$time.s. <- NULL
afhi_svr <- mutate(afhi_svr, pop="AFHI")
afhi_svr <- afhi_svr[order(afhi_svr$CV_R2),]

svr <- rbind(afa_svr, cau_svr, his_svr, afhi_svr, all_svr)
names(svr) <- c("Gene_ID", "Gene_Name", "Cross-validated R2", "Kernel", "Degree", "C", "Population")
fwrite(svr, file="Z:/ml_paper_figs/SVR_optimum_hyperparameter.txt", row.names=F, quote=F, sep="\t")

#KNN
afa_knn <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_knn_all_chrom.txt", header=T)
afa_knn$time.s. <- NULL
afa_knn <- mutate(afa_knn, pop="AFA")
afa_knn <- afa_knn[order(afa_knn$CV_R2),]

cau_knn <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_knn_all_chrom.txt", header=T)
cau_knn$time.s. <- NULL
cau_knn <- mutate(cau_knn, pop="CAU")
cau_knn <- cau_knn[order(cau_knn$CV_R2),]

his_knn <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_knn_all_chrom.txt", header=T)
his_knn$time.s. <- NULL
his_knn <- mutate(his_knn, pop="HIS")
his_knn <- his_knn[order(his_knn$CV_R2),]

all_knn <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_knn_all_chrom.txt", header=T)
all_knn$time.s. <- NULL
all_knn <- mutate(all_knn, pop="ALL")
all_knn <- all_knn[order(all_knn$CV_R2),]

afhi_knn <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_knn_all_chrom.txt", header=T)
afhi_knn$time.s. <- NULL
afhi_knn <- mutate(afhi_knn, pop="AFHI")
afhi_knn <- afhi_knn[order(afhi_knn$CV_R2),]

knn <- rbind(afa_knn, cau_knn, his_knn, afhi_knn, all_knn)
names(knn) <- c("Gene_ID", "Gene_Name", "Cross-validated R2", "Number of Neigbors", "Weights", "P", "Population")
fwrite(knn, file="Z:/ml_paper_figs/KNN_optimum_hyperparameter.txt", row.names=F, quote=F, sep="\t")

