library(data.table)
library(dplyr)
"%&%" <- function(a,b) paste(a,b, sep = "")
library(ggplot2)

gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)

#cetp <- "ENSG00000087237"
#st8sia4 <- "ENSG00000113532"

pop <- "" #this means ALL
trait <- "hdl"

#for (trait in c("hdl", "ldl", "chol", "trig", "plt5")){}

#for (pop in c("AFA", "CAU", "HIS")){}

en <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "en_rank_" %&% trait %&% "_assoc.txt", header=T)
rf <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "rf_rank_" %&% trait %&% "_assoc.txt", header=T)
svr <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "svr_rank_" %&% trait %&% "_assoc.txt", header=T)
knn <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "knn_rank_" %&% trait %&% "_assoc.txt", header=T)


#prepare data for manhattan plot

gencode$gene_id <- as.character(gencode$gene_id)
gencode$gene_type <- as.character(gencode$gene_type)
gencode <- subset(gencode, gene_type=="protein_coding") #this to reduce the df size to only protein coding genes
#remove the decimal in the gene id
for (i in 1:nrow(gencode)){
  gencode$gene_id[i] <- gsub('\\.[0-9]+','',gencode$gene_id[i])
} #just to remove the decimal places in the gene_id

gencode1 <- gencode[,c(1,2,3,4)]

en_man <- inner_join(gencode1, en, by = c("gene_id" = "gene"))
en_man <- en_man[,c(1,3,4,7)]#drop gene_id and use gene_name to help annotation in the plot

rf_man <- inner_join(gencode1, rf, by = c("gene_id" = "gene"))
rf_man <- rf_man[,c(1,3,4,7)] #drop gene_id and use gene_name to help annotation in the plot

svr_man <- inner_join(gencode1, svr, by = c("gene_id" = "gene"))
svr_man <- svr_man[,c(1,3,4,7)] #drop gene_id

knn_man <- inner_join(gencode1, knn, by = c("gene_id" = "gene"))
knn_man <- knn_man[,c(1,3,4,7)] #drop gene_id


#arrange the data by chrom number and start pos

#take each chrom out, sort it by "start", and rbind them
enassocman <- NULL
rfassocman <- NULL
svrassocman <- NULL
knnassocman <- NULL

#en
for (i in 1:22){
  i <- as.character(i)
  a <- subset(en_man, chr==i)
  a <- a[order(a$start),]
  enassocman <- rbind(enassocman,a)
}
names(enassocman) <- c("CHR", "SNP", "BP", "P")
enassocman$CHR <- as.numeric(enassocman$CHR)
# manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Elastic Net HDL", suggestiveline = F, ylim=c(0,14),
#           col=c("red", "blue"))



#rf
for (i in 1:22){
  i <- as.character(i)
  a <- subset(rf_man, chr==i)
  a <- a[order(a$start),]
  rfassocman <- rbind(rfassocman,a)
}
names(rfassocman) <- c("CHR", "SNP", "BP", "P")
rfassocman$CHR <- as.numeric(rfassocman$CHR)
# manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Random Forest HDL", suggestiveline = F,
#           col=c("red", "blue"), ylim=c(0,14))



#svr
for (i in 1:22){
  i <- as.character(i)
  a <- subset(svr_man, chr==i)
  a <- a[order(a$start),]
  svrassocman <- rbind(svrassocman,a)
}
names(svrassocman) <- c("CHR", "SNP", "BP", "P")
svrassocman$CHR <- as.numeric(svrassocman$CHR)
# manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Support Vector HDL", suggestiveline = F, ylim=c(0,14),
#           col=c("red", "blue"))


#knn
for (i in 1:22){
  i <- as.character(i)
  a <- subset(knn_man, chr==i)
  a <- a[order(a$start),]
  knnassocman <- rbind(knnassocman,a)
}
names(knnassocman) <- c("CHR", "SNP", "BP", "P")
knnassocman$CHR <- as.numeric(knnassocman$CHR)
#knnassocman$P <- as.numeric(knnassocman$P)
# manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="KNN HDL", suggestiveline = F, ylim=c(0,14),
#           col=c("red", "blue"))


# Fugure 5 - Dataframe for Manhattan Plot
fwrite(enassocman, file="Z:/ml_paper_figs/en_assoc_manplot_df.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassocman, file="Z:/ml_paper_figs/rf_assoc_manplot_df.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassocman, file="Z:/ml_paper_figs/svr_assoc_manplot_df.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassocman, file="Z:/ml_paper_figs/knn_assoc_manplot_df.txt", row.names=F, quote=F, sep="\t")



# Figure S4 - Dataframe for QQ Plot

#Make qq plot with the P-values
en_p <- enassocman[,4]
en_logp<- -log10(en_p)
en_logp <- sort(en_logp, decreasing = T, na.last = NA) %>% as.data.frame()

rf_p <- rfassocman[,4]
rf_logp<- -log10(rf_p)
rf_logp <- sort(rf_logp, decreasing = T, na.last = NA) %>% as.data.frame()

svr_p <- svrassocman[,4]
svr_logp<- -log10(svr_p)
svr_logp <- sort(svr_logp, decreasing = T, na.last = NA) %>% as.data.frame()

knn_p <- knnassocman[,4]
knn_logp<- -log10(knn_p)
knn_logp <- sort(knn_logp, decreasing = T, na.last = NA) %>% as.data.frame()

create_qq_input<-function(df, limit=nrow(df),range=nrow(df)){ #takes in one column df of pvalues - assumes no NAs, can also limit to the top n hits. Default of this will be just all pvalues
  print("Size of parental set is presumed to be " %&% range)
  print("Provided number of pvalues is " %&% nrow(df))
  print("limiting output to top " %&% limit %&% " rows")
  count <- nrow(df) #count number of pvals that are not na
  ExpP <- -log10((1:count)/(range+1)) %>% as.data.frame() #generate a list of expected pvalues equal in length to the obs pvalues that are not na
  qqvals<-cbind.data.frame(df[1:limit,],ExpP[1:limit,])
  colnames(qqvals)<-c("Observed","Expected")
  return(qqvals)
}

en_qq_input<-create_qq_input(en_logp)
rf_qq_input<-create_qq_input(rf_logp)
svr_qq_input<-create_qq_input(svr_logp)
knn_qq_input<-create_qq_input(knn_logp)

en_qq_input <- mutate(en_qq_input, model="Elastic Net")
rf_qq_input <- mutate(rf_qq_input, model="Random Forest")
svr_qq_input <- mutate(svr_qq_input, model="Support Vector")
knn_qq_input <- mutate(knn_qq_input, model="K Nearest Neighbor")

qq_input <- rbind(en_qq_input, rf_qq_input, svr_qq_input, knn_qq_input)
fwrite(qq_input, file="Z:/ml_paper_figs/qq_input_ALL.txt", row.names = F, quote = F, sep = "\t")