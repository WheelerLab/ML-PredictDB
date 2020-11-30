library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

"%&%" = function(a,b) paste (a,b,sep="")

data <- NULL

algs <- c("en", "knn", "rf", "svr")
pops <- c("AFA", "AFHI", "ALL", "CAU", "HIS")

for (pop in pops){

  for (alg in algs){
    
    ml <- fread("Z:/data/ml_paper_reviewers_corrections/mets_spearman/used_pcair_10pc_10peer_obs_exp/" %&% alg %&% "_" %&% pop %&% "_2_METS_corr_filt.txt",
                header=T, sep="\t", stringsAsFactors=F) %>% select(c("gene","spearman"))
    ml <- mutate(ml, model=toupper(alg),pop=pop)
    data <- rbind(data, ml)
    
  }
}


data <- mutate(data,pop=factor(pop,levels=c("ALL","AFHI","AFA","HIS","CAU")),model=factor(model,levels=c("EN","RF","SVR","KNN")))
data <- subset(data, spearman > -0.5)

fig <- ggplot(data,aes(x=spearman, colour=model)) + geom_density(lwd=0.8) + 
  xlab(expression(paste("Spearman Correlation ", rho))) +
  theme_classic(16) + facet_wrap(~pop) + labs(colour="Model")
# print(fig)
# tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/used_pcair_10pc_10peer_obs_exp/FigS6.tiff", 
#      width = 24, height = 18, units = 'cm', res = 300, compression = 'lzw')
# fig
# dev.off()


print(fig)
tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/FigS6_density_plot.tiff", 
     width = 24, height = 18, units = 'cm', res = 300, compression = 'lzw')
fig
dev.off()
