library(dplyr)
library(ggplot2)
library(viridis)

data <- read.table("Z:/ml_paper_figs/boxplotsdf.txt",header=T,stringsAsFactors = F)

#rearrange pop and ML order
data <- mutate(data,pop=factor(pop,levels=c("ALL","AFA","HIS","CAU")),model=factor(model,levels=c("EN","RF","SVR","KNN")))
fig <- ggplot(data,aes(x=pop,y=cvR2,fill=model)) + geom_boxplot(outlier.size=0.01, lwd=0.2) + 
  xlab("Population") + ylab(expression(paste("Cross-validated ", R^{2}))) +
  theme_classic(18) + labs(fill="Model") #theme_classic(12) + geom_boxplot(outlier.size = 0.01,lwd=0.1)
print(fig)
tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/FigS4_boxplots.tiff", width = 17, height = 12, units = 'cm', res = 300, compression = 'lzw') #width=12 #height=6
fig
dev.off()
