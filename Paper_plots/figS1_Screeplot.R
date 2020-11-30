#RAN on WHEELER LAB 3
#path<- "/home/paul/lauren_mesa/used_in_training/"

library(data.table)
library(dplyr)
library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(ggplot2)
# "%&%" = function(a,b) paste (a,b,sep="")
# 
# df <- NULL
eigenfacs <- c(1:32)
# 
# pccol <- vector(mode="character", length=32)
# for (i in 1:32){
#   pccol[i] <- paste("PC", i, sep="")
# }
# 
# #ALL
# gds <- snpgdsOpen("all/all.gds")
# snpset <- snpgdsLDpruning(gds, method="r", slide.max.bp=10e6, maf=0.01, missing.rate=0.01, ld.threshold=0.3, verbose=T)
# pruned <- unlist(snpset, use.names=F)
# snpgdsClose(gds)
# geno <- GdsGenotypeReader("all/all.gds")
# genodata <- GenotypeData(geno)
# mypcair <- pcair(genodata, snp.include=pruned)
# 
# pcs <- as.data.frame(mypcair$vectors)
# colnames(pcs) <- pccol
# pcval <- as.data.frame(mypcair$values)
# colnames(pcval) <- "Eigenvalues"
# fwrite(pcs, "all/mesa_ALL_onlytrain_pc.txt", quote=F, sep="\t")
# fwrite(pcval, "all/mesa_ALL_onlytrain_pcval.txt", quote=F, sep="\t")
# 
# eigenval <- data.frame(Components=eigenfacs, Eigenvalues=pcval$Eigenvalues, Population="ALL")
# df <- rbind(df, eigenval)
# 
# 
# #AFHI
# gds <- snpgdsOpen("afhi/afhi.gds")
# snpset <- snpgdsLDpruning(gds, method="r", slide.max.bp=10e6, maf=0.01, missing.rate=0.01, ld.threshold=0.3, verbose=T)
# pruned <- unlist(snpset, use.names=F)
# snpgdsClose(gds)
# geno <- GdsGenotypeReader("afhi/afhi.gds")
# genodata <- GenotypeData(geno)
# mypcair <- pcair(genodata, snp.include=pruned)
# 
# pcs <- as.data.frame(mypcair$vectors)
# colnames(pcs) <- pccol
# pcval <- as.data.frame(mypcair$values)
# colnames(pcval) <- "Eigenvalues"
# fwrite(pcs, "afhi/mesa_AFHI_onlytrain_pc.txt", quote=F, sep="\t")
# fwrite(pcval, "afhi/mesa_AFHI_onlytrain_pcval.txt", quote=F, sep="\t")
# 
# eigenval <- data.frame(Components=eigenfacs, Eigenvalues=pcval$Eigenvalues, Population="AFHI")
# df <- rbind(df, eigenval)
# 
# 
# 
# #AFA
# gds <- snpgdsOpen("afa/afa.gds")
# snpset <- snpgdsLDpruning(gds, method="r", slide.max.bp=10e6, maf=0.01, missing.rate=0.01, ld.threshold=0.3, verbose=T)
# pruned <- unlist(snpset, use.names=F)
# snpgdsClose(gds)
# geno <- GdsGenotypeReader("afa/afa.gds")
# genodata <- GenotypeData(geno)
# mypcair <- pcair(genodata, snp.include=pruned)
# 
# pcs <- as.data.frame(mypcair$vectors)
# colnames(pcs) <- pccol
# pcval <- as.data.frame(mypcair$values)
# colnames(pcval) <- "Eigenvalues"
# fwrite(pcs, "afa/mesa_AFA_onlytrain_pc.txt", quote=F, sep="\t")
# fwrite(pcval, "afa/mesa_AFA_onlytrain_pcval.txt", quote=F, sep="\t")
# 
# eigenval <- data.frame(Components=eigenfacs, Eigenvalues=pcval$Eigenvalues, Population="AFA")
# df <- rbind(df, eigenval)
# 
# 
# 
# #CAU
# gds <- snpgdsOpen("cau/cau.gds")
# snpset <- snpgdsLDpruning(gds, method="r", slide.max.bp=10e6, maf=0.01, missing.rate=0.01, ld.threshold=0.3, verbose=T)
# pruned <- unlist(snpset, use.names=F)
# snpgdsClose(gds)
# geno <- GdsGenotypeReader("cau/cau.gds")
# genodata <- GenotypeData(geno)
# mypcair <- pcair(genodata, snp.include=pruned)
# 
# pcs <- as.data.frame(mypcair$vectors)
# colnames(pcs) <- pccol
# pcval <- as.data.frame(mypcair$values)
# colnames(pcval) <- "Eigenvalues"
# fwrite(pcs, "cau/mesa_CAU_onlytrain_pc.txt", quote=F, sep="\t")
# fwrite(pcval, "cau/mesa_CAU_onlytrain_pcval.txt", quote=F, sep="\t")
# 
# eigenval <- data.frame(Components=eigenfacs, Eigenvalues=pcval$Eigenvalues, Population="CAU")
# df <- rbind(df, eigenval)
# 
# 
# 
# #HIS
# gds <- snpgdsOpen("his/his.gds")
# snpset <- snpgdsLDpruning(gds, method="r", slide.max.bp=10e6, maf=0.01, missing.rate=0.01, ld.threshold=0.3, verbose=T)
# pruned <- unlist(snpset, use.names=F)
# snpgdsClose(gds)
# geno <- GdsGenotypeReader("his/his.gds")
# genodata <- GenotypeData(geno)
# mypcair <- pcair(genodata, snp.include=pruned)
# 
# pcs <- as.data.frame(mypcair$vectors)
# colnames(pcs) <- pccol
# pcval <- as.data.frame(mypcair$values)
# colnames(pcval) <- "Eigenvalues"
# fwrite(pcs, "his/mesa_HIS_onlytrain_pc.txt", quote=F, sep="\t")
# fwrite(pcval, "his/mesa_HIS_onlytrain_pcval.txt", quote=F, sep="\t")
# 
# eigenval <- data.frame(Components=eigenfacs, Eigenvalues=pcval$Eigenvalues, Population="HIS")
# df <- rbind(df, eigenval)
# 
# fwrite(df, "/home/paul/lauren_mesa/used_in_training/screeplot_df.txt", row.names=F, quote=F, sep="\t")



#ON ROCKS
df <- fread("Z:/data/lauren_mesa/lauren_train_mesa_plink/screeplot_df.txt", header=T)

#Add the eignevals of mets and twas_mesa

pcval <- fread("Z:/data/ml_paper_reviewers_corrections/genotypes/mets/mets_male_removed_pcair_king_eigenval_ld_threshold_0.3.txt",
               header=F)
eigenval <- data.frame(Components=eigenfacs, Eigenvalues=pcval$V2, Population="METS")
df <- rbind(df, eigenval)

#TWAS MESA
pcval <- fread("Z:/data/ml_paper_reviewers_corrections/genotypes/twas_mesa/mesa_notrain_pcair_nokin_eigenval_ld_threshold_0.3.txt",
               header=F)
eigenval <- data.frame(Components=eigenfacs, Eigenvalues=pcval$V2, Population="TWAS_MESA")
df <- rbind(df, eigenval)


#fwrite(df, "Z:/data/ml_paper_reviewers_corrections/paper_figs/new_screeplot_df.txt", row.names=F, quote=F, sep="\t")


#Make plots
#df1 <- subset(df, Population=="AFHI")

tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/FigS1_Screeplot.tiff",
     width = 24, height = 20, units = 'cm', res = 500, compression = 'lzw')

ggplot(df,aes(x=Components, y=Eigenvalues)) + geom_point() + geom_line() + xlim(0,30) +
  theme_classic(16) + facet_wrap(~Population, scales="free_y") + scale_x_continuous(breaks=c(seq(0,30,3)))

dev.off()
