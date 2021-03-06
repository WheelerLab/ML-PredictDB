library(data.table)
library(dplyr)
library(ggplot2)


trees5_5000 <- fread(file = "Z:/ml_paper_figs/rf_trees_cvr2_afa_allchrom.txt", header=T)
names(trees5_5000) <- c("gene_id","gene_name","5","50","100","150","200","250","300","350","400","450","500","5000")

#Now make the plot of trees 5, 50, 500, and 5000
t5_50_500_5000 <- trees5_5000[,c(1,3:14)]
#Take out only trees 5, 50, 500 and 5000
tree5 <- data.frame(cvr2=t5_50_500_5000[,2],trees="5")
names(tree5) <- c("cvr2","trees")
tree50 <- data.frame(cvr2=t5_50_500_5000[,3], trees="50")
names(tree50) <- c("cvr2","trees")
tree500 <- data.frame(cvr2=t5_50_500_5000[,12], trees="500")
names(tree500) <- c("cvr2","trees")
tree5000 <- data.frame(cvr2=t5_50_500_5000[,13], trees="5000")
names(tree5000) <- c("cvr2","trees")
treedf <- rbind(tree5,tree50,tree500,tree5000)

#filter out cvr2 < -1
treedf <- subset(treedf, cvr2 > -1)

fig <- ggplot(treedf, aes(x=trees, y=cvr2, fill=trees)) + geom_violin() + theme_classic(18) + ylim(-1,1) + 
  ylab(expression(paste("Cross-Validated ", R^{2}))) + xlab("Number of Trees") + 
  theme(legend.position = "none") + geom_boxplot(width=0.15) # +ggtitle("Random Forest Trees Performance") + 
print(fig)
tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/FigS2_rf_trees.tiff", width = 18, height = 13, 
     units = 'cm', res = 300, compression = 'lzw')
fig
dev.off()
