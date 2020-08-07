#Hyperopt standardization
#Maximum evaluation = 30

library(dplyr)
library(ggplot2)
"%&%" <- function(a,b) paste(a,b, sep = "")

#Elastic Net
en <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    en <- rbind(en, read.table(file="Z:/data/paper_hyperopt/max_evals_30/afa_en_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
en$X <- NULL
en <- en[,c(1,34)]
en <- subset(en, X29 > -0.5)
names(en) <- c("gene_id", "cvr2")


#Random Forest, chr22
rf <- NULL
for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:22) {
    rf <- rbind(rf, read.table(file="Z:/data/paper_hyperopt/RF/max_evals_30/afa_rf_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
rf$X <- NULL
rf <- rf[,c(1,34)]
rf <- subset(rf, X29 > -0.5)
names(rf) <- c("gene_id", "cvr2")


#SVR
svr <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    svr <- rbind(svr, read.table(file="Z:/data/paper_hyperopt/max_evals_30/afa_svr_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
svr$X <- NULL
svr <- svr[,c(1,34)]
svr <- subset(svr, X29 > -0.5)
names(svr) <- c("gene_id", "cvr2")


#KNN
knn <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    knn <- rbind(knn, read.table(file="Z:/data/paper_hyperopt/max_evals_30/afa_knn_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
knn$X <- NULL
knn <- knn[,c(1,34)]
knn <- subset(knn, X29 > -0.5)
names(knn) <- c("gene_id", "cvr2")


#compare the hyperopts of EN against RF, SVR or KNN

#RF
enrf <- inner_join(en, rf, by = c("gene_id"="gene_id"))
names(enrf) <- c("gene", "EN", "RF")
enrf <- subset(enrf, (EN>-0.5 & RF>-0.5)) #9441
ggplot(data = enrf, aes(x=EN, y=RF)) + geom_point() + geom_point() + theme_bw(24) + 
  geom_abline(intercept=0, slope=1, color="blue") + 
  geom_smooth(method="lm", color="red", lwd=0.5)

#SVR
ensvr <- inner_join(en, svr, by = c("gene_id"="gene_id"))
names(ensvr) <- c("gene", "EN", "SVR")
ensvr <- subset(ensvr, (EN>-0.5 & SVR>-0.5)) #9441
ggplot(data = ensvr, aes(x=EN, y=SVR)) + geom_point() + geom_point() + theme_bw(24) + 
  geom_abline(intercept=0, slope=1, color="blue") + 
  geom_smooth(method="lm", color="red", lwd=0.5)


#KNN
enknn <- inner_join(en, knn, by = c("gene_id"="gene_id"))
names(enknn) <- c("gene", "EN", "KNN")
enknn <- subset(enknn, (EN>-0.5 & KNN>-0.5)) #9441
ggplot(data = enknn, aes(x=EN, y=KNN)) + geom_point() + geom_point() + theme_bw(24) + 
  geom_abline(intercept=0, slope=1, color="blue") + 
  geom_smooth(method="lm", color="red", lwd=0.5)


#Make df
#AFA

endf <- mutate(en, model="EN", pop="AFA")
rfdf <- mutate(rf, model="RF", pop="AFA")
svrdf <- mutate(svr, model="SVR", pop="AFA")
knndf <- mutate(knn, model="KNN", pop="AFA")

#df <- rbind(endf, rfdf, svrdf, knndf)
enrf <- inner_join(endf,rfdf,by=c("gene_id","pop"))
ensvr <- inner_join(endf,svrdf,by=c("gene_id","pop"))
enknn <- inner_join(endf,knndf,by=c("gene_id","pop"))

df <- rbind(enrf, ensvr, enknn)
colnames(df) <- c("gene", "ENcvR2", "EN", "pop", "cvR2", "MLmodel")

data <- mutate(df,pop=factor(pop,levels=c("AFA")),MLmodel=factor(MLmodel,levels=c("RF","SVR","KNN")))

tiff("/Users/okoro/OneDrive/Desktop/FigS6.tiff", width = 16, height = 20, units = 'cm', res = 300, compression = 'lzw')

ggplot(data,aes(x=ENcvR2,y=cvR2)) + geom_point(shape=".") + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red", lwd=0.5) + xlim(-0.5,1) + ylim(-0.5,1) + 
  xlab(expression(paste("Elastic Net ", R^{2}))) + ylab(expression(paste("Machine Learning Model   ", R^{2}))) +
  theme_bw(16) + facet_grid(pop~MLmodel) 

dev.off()





######
#AFA

#Elastic Net
en_afa <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    en_afa <- rbind(en_afa, read.table(file="Z:/data/paper_hyperopt/max_evals_30/afa_en_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
en_afa$X <- NULL
en_afa <- en_afa[,c(1,34)]
en_afa <- subset(en_afa, X29 > -0.5)
names(en_afa) <- c("gene_id", "cvr2")


#Random Forest, chr22
rf_afa <- NULL
for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:22) {
    rf_afa <- rbind(rf_afa, read.table(file="Z:/data/paper_hyperopt/RF/max_evals_30/afa_rf_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
rf_afa$X <- NULL
rf_afa <- rf_afa[,c(1,34)]
rf_afa <- subset(rf_afa, X29 > -0.5)
names(rf_afa) <- c("gene_id", "cvr2")


#SVR
svr_afa <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    svr_afa <- rbind(svr_afa, read.table(file="Z:/data/paper_hyperopt/max_evals_30/afa_svr_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
svr_afa$X <- NULL
svr_afa <- svr_afa[,c(1,34)]
svr_afa <- subset(svr_afa, X29 > -0.5)
names(svr_afa) <- c("gene_id", "cvr2")


#KNN
knn_afa <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    knn_afa <- rbind(knn_afa, read.table(file="Z:/data/paper_hyperopt/max_evals_30/afa_knn_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
knn_afa$X <- NULL
knn_afa <- knn_afa[,c(1,34)]
knn_afa <- subset(knn_afa, X29 > -0.5)
names(knn_afa) <- c("gene_id", "cvr2")



##################################################################################
######
#CAU

#Elastic Net
en_cau <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    en_cau <- rbind(en_cau, read.table(file="Z:/data/paper_hyperopt/max_evals_30/cau_en_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
en_cau$X <- NULL
en_cau <- en_cau[,c(1,34)]
en_cau <- subset(en_cau, X29 > -0.5)
names(en_cau) <- c("gene_id", "cvr2")


#Random Forest, chr22
rf_cau <- NULL
for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:22) {
    rf_cau <- rbind(rf_cau, read.table(file="Z:/data/paper_hyperopt/RF/max_evals_30/cau_rf_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
rf_cau$X <- NULL
rf_cau <- rf_cau[,c(1,34)]
rf_cau <- subset(rf_cau, X29 > -0.5)
names(rf_cau) <- c("gene_id", "cvr2")


#SVR
svr_cau <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    svr_cau <- rbind(svr_cau, read.table(file="Z:/data/paper_hyperopt/max_evals_30/cau_svr_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
svr_cau$X <- NULL
svr_cau <- svr_cau[,c(1,34)]
svr_cau <- subset(svr_cau, X29 > -0.5)
names(svr_cau) <- c("gene_id", "cvr2")


#KNN
knn_cau <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    knn_cau <- rbind(knn_cau, read.table(file="Z:/data/paper_hyperopt/max_evals_30/cau_knn_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
knn_cau$X <- NULL
knn_cau <- knn_cau[,c(1,34)]
knn_cau <- subset(knn_cau, X29 > -0.5)
names(knn_cau) <- c("gene_id", "cvr2")







##################################################################################
######
#HIS

#Elastic Net
en_his <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    en_his <- rbind(en_his, read.table(file="Z:/data/paper_hyperopt/max_evals_30/his_en_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
en_his$X <- NULL
en_his <- en_his[,c(1,34)]
en_his <- subset(en_his, X29 > -0.5)
names(en_his) <- c("gene_id", "cvr2")


#Random Forest, chr22
rf_his <- NULL
for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:22) {
    rf_his <- rbind(rf_his, read.table(file="Z:/data/paper_hyperopt/RF/max_evals_30/his_rf_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
rf_his$X <- NULL
rf_his <- rf_his[,c(1,34)]
rf_his <- subset(rf_his, X29 > -0.5)
names(rf_his) <- c("gene_id", "cvr2")


#SVR
svr_his <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    svr_his <- rbind(svr_his, read.table(file="Z:/data/paper_hyperopt/max_evals_30/his_svr_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
svr_his$X <- NULL
svr_his <- svr_his[,c(1,34)]
svr_his <- subset(svr_his, X29 > -0.5)
names(svr_his) <- c("gene_id", "cvr2")


#KNN
knn_his <- NULL
for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    knn_his <- rbind(knn_his, read.table(file="Z:/data/paper_hyperopt/max_evals_30/his_knn_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
knn_his$X <- NULL
knn_his <- knn_his[,c(1,34)]
knn_his <- subset(knn_his, X29 > -0.5)
names(knn_his) <- c("gene_id", "cvr2")



#

#Make df

#AFA
endf <- mutate(en_afa, model="EN", pop="AFA")
rfdf <- mutate(rf_afa, model="RF", pop="AFA")
svrdf <- mutate(svr_afa, model="SVR", pop="AFA")
knndf <- mutate(knn_afa, model="KNN", pop="AFA")

#df <- rbind(endf, rfdf, svrdf, knndf)
enrf <- inner_join(endf,rfdf,by=c("gene_id","pop"))
ensvr <- inner_join(endf,svrdf,by=c("gene_id","pop"))
enknn <- inner_join(endf,knndf,by=c("gene_id","pop"))
df <- rbind(enrf, ensvr, enknn)

#CAU
endf <- mutate(en_cau, model="EN", pop="CAU")
rfdf <- mutate(rf_cau, model="RF", pop="CAU")
svrdf <- mutate(svr_cau, model="SVR", pop="CAU")
knndf <- mutate(knn_cau, model="KNN", pop="CAU")

#df <- rbind(endf, rfdf, svrdf, knndf)
enrf <- inner_join(endf,rfdf,by=c("gene_id","pop"))
ensvr <- inner_join(endf,svrdf,by=c("gene_id","pop"))
enknn <- inner_join(endf,knndf,by=c("gene_id","pop"))
df <- rbind(df, enrf, ensvr, enknn)


#HIS
endf <- mutate(en_his, model="EN", pop="HIS")
rfdf <- mutate(rf_his, model="RF", pop="HIS")
svrdf <- mutate(svr_his, model="SVR", pop="HIS")
knndf <- mutate(knn_his, model="KNN", pop="HIS")

#df <- rbind(endf, rfdf, svrdf, knndf)
enrf <- inner_join(endf,rfdf,by=c("gene_id","pop"))
ensvr <- inner_join(endf,svrdf,by=c("gene_id","pop"))
enknn <- inner_join(endf,knndf,by=c("gene_id","pop"))

df <- rbind(df, enrf, ensvr, enknn)
colnames(df) <- c("gene", "ENcvR2", "EN", "pop", "cvR2", "MLmodel")


#Plot

data <- select(df, gene, ENcvR2, MLmodel, cvR2, pop)
data <- mutate(data,pop=factor(pop,levels=c("AFA","HIS","CAU")),MLmodel=factor(MLmodel,levels=c("RF","SVR","KNN")))

tiff("/Users/okoro/OneDrive/Desktop/FigS2.tiff", width = 16, height = 20, units = 'cm', res = 300, compression = 'lzw')

ggplot(data,aes(x=ENcvR2,y=cvR2)) + geom_point(shape=".") + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red", lwd=0.5) + xlim(-0.5,1) + ylim(-0.5,1) + 
  xlab(expression(paste("Elastic Net ", R^{2}))) + ylab(expression(paste("Machine Learning Model ", R^{2}))) +
  theme_bw(16) + facet_grid(pop~MLmodel) 

#tiff("Fig1.tiff", width = 16, height = 20, units = 'cm', res = 300, compression = 'lzw')
#fig
dev.off()
