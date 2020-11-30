#Make ML paper Figure 7

library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(data.table)

#compare the t statistic performance of EN and RF
#The predicted expression was done with ALL trained model
#read in the assoc

#color the intersection genes black, and others red
#read in the assoc
assoc_en <- fread("Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/en_twas_mesa_hdl_assoc_r2_pc3.txt", header=T)
assoc_en <- assoc_en[,c(1,3)] #cols are gene and t
names(assoc_en) <- c("gene","en")

assoc_rf <- fread("Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/rf_twas_mesa_hdl_assoc_r2_pc3.txt", header=T)
assoc_rf <- assoc_rf[,c(1,3)] #cols are gene and t
names(assoc_rf) <- c("gene","rf")

#join both and keep all rows
en_rf <- full_join(assoc_en, assoc_rf, by = c("gene"="gene"))
en_rf$type <- "" #this type column is so i can keep track of the intercept and unique genes and color it so

#name gene intercept if there is no NA vice versa
for (i in 1:nrow(en_rf)){
  ifelse(is.na(en_rf[i,2]) | is.na(en_rf[i,3]), en_rf$type[i]<-"unique", en_rf$type[i]<-"common")
}

#now replace all NAs with 0
en_rf[is.na(en_rf)] <- 0

#plot them
ggplot(en_rf, aes(x=en, y=rf, col=type)) + geom_point(lwd=3)  +
  xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(30) + geom_abline(intercept=0, slope=1, color="blue") +
  scale_color_manual(values = c("black","red")) #  + ggtitle("EN vs RF on HDL")
#geom_smooth(method="lm", color="red")


###put gene name

# #join both and keep all rows
# en_rf <- full_join(assoc_en, assoc_rf, by = c("gene"="gene"))

#merge the gene pheno assoc results with gencode annotation
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
gencode$gene_id <- as.character(gencode$gene_id)
gencode$gene_type <- as.character(gencode$gene_type)
gencode <- subset(gencode, gene_type=="protein_coding") #this to reduce the df size to only protein coding genes
#remove the decimal in the gene id
for (i in 1:nrow(gencode)){
  gencode$gene_id[i] <- gsub('\\.[0-9]+','',gencode$gene_id[i])
} #just to remove the decimal places in the gene_id

gencode <- gencode[,c(2,3)]
gencode$gene_name <- as.character(gencode$gene_name)

#annotate the df in order to tag genes cetp or st8sia4
en_rf_gene <- inner_join(gencode, en_rf, by = c("gene_id"="gene")) %>% 
  mutate(is_annotate=ifelse(gene_name == "CETP" | gene_name == "ST8SIA4", "yes", "no"))

# for (i in 1:nrow(en_rf_gene)){
#   ifelse(is.na(en_rf_gene[i,2]) | is.na(en_rf_gene[i,3]), en_rf_gene$type[i]<-"unique", en_rf_gene$type[i]<-"common")
# }
# 
# #now replace all NAs with 0
# en_rf_gene[is.na(en_rf_gene)] <- 0

fwrite(en_rf_gene, file="Z:/data/ml_paper_reviewers_corrections/paper_figs/en_rf_tstatistic_df.txt", 
       row.names=F, quote=F, sep="\t")


fig <- ggplot(en_rf_gene, aes(x=en, y=rf, colour=type)) + geom_point(lwd=1.5)  +
  xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(20) + geom_abline(intercept=0, slope=1, color="blue") +
  scale_color_manual(values = c("black","red")) + xlim(-7,5) + ylim(-7.5,5) +
  geom_label_repel(data=en_rf_gene[en_rf_gene$is_annotate=="yes",], aes(label=as.factor(gene_name)), size=5, force=1.3) +
  labs(colour="Type")
print(fig)
tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/Fig7_tstatistic.tiff", width = 18, height = 16, 
     units = 'cm', res = 300, compression = 'lzw')
fig
dev.off()
