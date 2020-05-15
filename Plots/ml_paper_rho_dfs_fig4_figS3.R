#Make Figure 2 and 3
#Figure 2


library(dplyr)
#####AFHI
afhi_2_mets_en <- read.table(file="Z:/data/ml_paper/afhi_2_mets_en_rho_all.txt", header=T)
afhi_2_mets_en <- subset(afhi_2_mets_en, spearman > -0.5)
afhi_2_mets_en <- mutate(afhi_2_mets_en, model="EN",pop="AFHI")

afhi_2_mets_rf <- read.table(file="Z:/data/ml_paper/afhi_2_mets_rf_rho_all.txt", header=T)
afhi_2_mets_rf <- subset(afhi_2_mets_rf, spearman > -0.5)
afhi_2_mets_rf <- mutate(afhi_2_mets_rf, model="RF",pop="AFHI")

afhi_2_mets_svr <- read.table(file="Z:/data/ml_paper/afhi_2_mets_svr_rho_all.txt", header=T)
afhi_2_mets_svr <- subset(afhi_2_mets_svr, spearman >-0.5)
afhi_2_mets_svr <- mutate(afhi_2_mets_svr, model="SVR",pop="AFHI")

afhi_2_mets_knn <- read.table(file="Z:/data/ml_paper/afhi_2_mets_knn_rho_all.txt", header=T)
afhi_2_mets_knn <- subset(afhi_2_mets_knn, spearman >-0.5)
afhi_2_mets_knn <- mutate(afhi_2_mets_knn, model="KNN",pop="AFHI")

#data.frame for	boxplots
bpdf <- rbind(afhi_2_mets_en, afhi_2_mets_rf, afhi_2_mets_svr, afhi_2_mets_knn)

#data.frame for dot plots EN v. ML model
enrf <- inner_join(afhi_2_mets_en, afhi_2_mets_rf,by=c("gene","pop"))
ensvr <- inner_join(afhi_2_mets_en, afhi_2_mets_svr,by=c("gene","pop"))
enknn <- inner_join(afhi_2_mets_en, afhi_2_mets_knn,by=c("gene","pop"))
#try facet_grid(pop ~ model)
fig1df <- rbind(enrf, ensvr, enknn)

#####ALL
all_2_mets_en <- read.table(file="Z:/data/ml_paper/all_2_mets_en_rho_all.txt", header=T)
all_2_mets_en <- subset(all_2_mets_en, spearman >-0.5)
all_2_mets_en <- mutate(all_2_mets_en, model="EN",pop="ALL")

all_2_mets_rf <- read.table(file="Z:/data/ml_paper/all_2_mets_rf_rho_all.txt", header=T)
all_2_mets_rf <- subset(all_2_mets_rf, spearman >-0.5)
all_2_mets_rf <- mutate(all_2_mets_rf, model="RF",pop="ALL")

all_2_mets_svr <- read.table(file="Z:/data/ml_paper/all_2_mets_svr_rho_all.txt", header=T)
all_2_mets_svr <- subset(all_2_mets_svr, spearman >-0.5)
all_2_mets_svr <- mutate(all_2_mets_svr, model="SVR",pop="ALL")

all_2_mets_knn <- read.table(file="Z:/data/ml_paper/all_2_mets_knn_rho_all.txt", header=T)
all_2_mets_knn <- subset(all_2_mets_knn, spearman >-0.5)
all_2_mets_knn <- mutate(all_2_mets_knn, model="KNN",pop="ALL")

#data.frame for	boxplots
bpdf <- rbind(bpdf, all_2_mets_en, all_2_mets_rf, all_2_mets_svr, all_2_mets_knn)

#data.frame for dot plots EN v. ML model
enrf <- inner_join(all_2_mets_en, all_2_mets_rf,by=c("gene","pop"))
ensvr <- inner_join(all_2_mets_en, all_2_mets_svr,by=c("gene","pop"))
enknn <- inner_join(all_2_mets_en, all_2_mets_knn,by=c("gene","pop"))
#try facet_grid(pop ~ model)
fig1df <- rbind(fig1df, enrf, ensvr, enknn)

#CAU
cau_2_mets_en <- read.table(file="Z:/data/ml_paper/cau_2_mets_en_rho_all.txt", header=T)
cau_2_mets_en <- subset(cau_2_mets_en, spearman >-0.5)
cau_2_mets_en <- mutate(cau_2_mets_en, model="EN",pop="CAU")

cau_2_mets_rf <- read.table(file="Z:/data/ml_paper/cau_2_mets_rf_rho_all.txt", header=T)
cau_2_mets_rf <- subset(cau_2_mets_rf, spearman >-0.5)
cau_2_mets_rf <- mutate(cau_2_mets_rf, model="RF",pop="CAU")

cau_2_mets_svr <- read.table(file="Z:/data/ml_paper/cau_2_mets_svr_rho_all.txt", header=T)
cau_2_mets_svr <- subset(cau_2_mets_svr, spearman >-0.5)
cau_2_mets_svr <- mutate(cau_2_mets_svr, model="SVR",pop="CAU")

cau_2_mets_knn <- read.table(file="Z:/data/ml_paper/cau_2_mets_knn_rho_all.txt", header=T)
cau_2_mets_knn <- subset(cau_2_mets_knn, spearman >-0.5)
cau_2_mets_knn <- mutate(cau_2_mets_knn, model="KNN",pop="CAU")

#data.frame for	boxplots
bpdf <- rbind(bpdf, cau_2_mets_en, cau_2_mets_rf, cau_2_mets_svr, cau_2_mets_knn)

#data.frame for dot plots EN v. ML model
enrf <- inner_join(cau_2_mets_en, cau_2_mets_rf,by=c("gene","pop"))
ensvr <- inner_join(cau_2_mets_en, cau_2_mets_svr,by=c("gene","pop"))
enknn <- inner_join(cau_2_mets_en, cau_2_mets_knn,by=c("gene","pop"))
#try facet_grid(pop ~ model)
fig1df <- rbind(fig1df, enrf, ensvr, enknn)

#HIS
his_2_mets_en <- read.table(file="Z:/data/ml_paper/his_2_mets_en_rho_all.txt", header=T)
his_2_mets_en <- subset(his_2_mets_en, spearman >-0.5)
his_2_mets_en <- mutate(his_2_mets_en, model="EN",pop="HIS")

his_2_mets_rf <- read.table(file="Z:/data/ml_paper/his_2_mets_rf_rho_all.txt", header=T)
his_2_mets_rf <- subset(his_2_mets_rf, spearman >-0.5)
his_2_mets_rf <- mutate(his_2_mets_rf, model="RF",pop="HIS")

his_2_mets_svr <- read.table(file="Z:/data/ml_paper/his_2_mets_svr_rho_all.txt", header=T)
his_2_mets_svr <- subset(his_2_mets_svr, spearman >-0.5)
his_2_mets_svr <- mutate(his_2_mets_svr, model="SVR",pop="HIS")

his_2_mets_knn <- read.table(file="Z:/data/ml_paper/his_2_mets_knn_rho_all.txt", header=T)
his_2_mets_knn <- subset(his_2_mets_knn, spearman >-0.5)
his_2_mets_knn <- mutate(his_2_mets_knn, model="KNN",pop="HIS")

bpdf <- rbind(bpdf, his_2_mets_en, his_2_mets_rf, his_2_mets_svr, his_2_mets_knn)

#data.frame for dot plots EN v. ML model
enrf <- inner_join(his_2_mets_en, his_2_mets_rf,by=c("gene","pop"))
ensvr <- inner_join(his_2_mets_en, his_2_mets_svr,by=c("gene","pop"))
enknn <- inner_join(his_2_mets_en, his_2_mets_knn,by=c("gene","pop"))
#try facet_grid(pop ~ model)
fig1df <- rbind(fig1df, enrf, ensvr, enknn)

#AFA
afa_2_mets_en <- read.table(file="Z:/data/ml_paper/afa_2_mets_en_rho_all.txt", header=T)
afa_2_mets_en <- subset(afa_2_mets_en, spearman >-0.5)
afa_2_mets_en <- mutate(afa_2_mets_en, model="EN",pop="AFA")

afa_2_mets_rf <- read.table(file="Z:/data/ml_paper/afa_2_mets_rf_rho_all.txt", header=T)
afa_2_mets_rf <- subset(afa_2_mets_rf, spearman >-0.5)
afa_2_mets_rf <- mutate(afa_2_mets_rf, model="RF",pop="AFA")

afa_2_mets_svr <- read.table(file="Z:/data/ml_paper/afa_2_mets_svr_rho_all.txt", header=T)
afa_2_mets_svr <- subset(afa_2_mets_svr, spearman >-0.5)
afa_2_mets_svr <- mutate(afa_2_mets_svr, model="SVR",pop="AFA")

afa_2_mets_knn <- read.table(file="Z:/data/ml_paper/afa_2_mets_knn_rho_all.txt", header=T)
afa_2_mets_knn <- subset(afa_2_mets_knn, spearman >-0.5)
afa_2_mets_knn <- mutate(afa_2_mets_knn, model="KNN",pop="AFA")

bpdf <- rbind(bpdf, afa_2_mets_en, afa_2_mets_rf, afa_2_mets_svr, afa_2_mets_knn)

#data.frame for dot plots EN v. ML model
enrf <- inner_join(afa_2_mets_en, afa_2_mets_rf,by=c("gene","pop"))
ensvr <- inner_join(afa_2_mets_en, afa_2_mets_svr,by=c("gene","pop"))
enknn <- inner_join(afa_2_mets_en, afa_2_mets_knn,by=c("gene","pop"))
#try facet_grid(pop ~ model)
fig1df <- rbind(fig1df, enrf, ensvr, enknn)

colnames(fig1df) <- c("gene", "ENrho", "EN", "pop", "rho", "MLmodel")

fig1df <- select(fig1df, gene, ENrho, MLmodel, rho, pop)

write.table(fig1df, "Z:/ml_paper_figs/fig4df.txt", quote=F, row.names=F)
write.table(bpdf, "Z:/ml_paper_figs/figs3densityplotsdf.txt", quote=F, row.names=F)
