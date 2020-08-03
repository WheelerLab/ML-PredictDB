####This script is to make plots for the different random forest trees parameter cv R2 performance
# Using only AFA

library(data.table)
library(dplyr)
library(ggplot2)

trees50_500 <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_rf_100_to_500tree_cv_full_chr.txt", header=T)
trees50_500$gene_id <- as.character(trees50_500$gene_id)
trees50_500$chr <- NULL

trees5 <- fread(file = "Z:/data/mesa_models/python_ml_models/results/AFA_rf_5tree_cv_full_chr.txt", header=T)
trees5$Gene_ID <- as.character(trees5$Gene_ID)
trees5[,c(2,4)] <- NULL
names(trees5) <- c("gene_id", "5")


trees5000 <- fread(file = "Z:/data/mesa_models/python_ml_models/results/AFA_rf_5000tree_check_cv_full_chr.txt", header=T)
trees5000$Gene_ID <- as.character(trees5000$Gene_ID)
trees5000[,c(2,4)] <- NULL
names(trees5000) <- c("gene_id", "5000")

#combine the dataframes of trees 5, 50, 100, 150 ... 5000
trees5_500 <- inner_join(trees50_500, trees5, by = c("gene_id"="gene_id"))
trees5_5000 <- inner_join(trees5_500, trees5000, by = c("gene_id"="gene_id"))


#arrange the the df to be in order starting from trees 5 to 5000
trees5_5000 <- trees5_5000[,c(1,2,13,3,4,5,6,7,8,9,10,11,12,14)]
names(trees5_5000) <- c("gene_id","gene_name","5","50","100","150","200","250","300","350","400","450","500","5000")

fwrite(trees5_5000, file = "Z:/ml_paper_figs/rf_trees_cvr2_afa_allchrom.txt", row.names=F, quote=F, sep="\t")

#Now make the plot of trees 5, 50, 500, and 5000
t5_50_500_5000 <- trees5_5000[,c(1,3:14)]

#Take out only trees 5, 50, 500 and 5000
tree5 <- data.frame(cvr2=t5_50_500_5000[,2],trees="5")
tree50 <- data.frame(cvr2=t5_50_500_5000[,3], trees="50")
tree500 <- data.frame(cvr2=t5_50_500_5000[,12], trees="500")
tree5000 <- data.frame(cvr2=t5_50_500_5000[,13], trees="5000")

treedf <- rbind(tree5,tree50,tree500,tree5000)
#filter out cvr2 < -1
treedf <- subset(treedf, cvr2 > -1)

tiff("/users/okoro/OneDrive/Desktop/mesa_scripts/ml_paper_figs/FigS5.tiff", width = 18, height = 13, 
     units = 'cm', res = 300, compression = 'lzw')

ggplot(treedf, aes(x=trees, y=cvr2, fill=trees)) + geom_violin() + theme_classic(18) + ylim(-1,1) +
  ylab(expression(paste("Cross-Validated ", R^{2}))) + xlab("Number of Trees") + 
  theme(legend.position = "none") + geom_boxplot(width=0.15) # +ggtitle("Random Forest Trees Performance") + 

dev.off()


#rf <- t5_50_500_5000[,c(2,3,12,13)]
# trees <- c(5,50,500, 5000)
# trees <- as.factor(trees)
# 
# cl <- rainbow(length(rf3$`5`))
# plot(trees, rf3[1,], type = 'o', xlab = "Trees", ylab = "CV R2", main = "RF Tree", ylim = c(0.3,1), col=cl[1])
# 
# for (i in 2:length(rftrees$`50`)){
#   lines(trees, rf3[i,], type = "o", col=cl[i+3])
# }


