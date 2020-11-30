library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)

data <- read.table(file="Z:/ml_paper_figs/fig1df.txt",header=T,stringsAsFactors = F)
data <- mutate(data,pop=factor(pop,levels=c("ALL","AFA","HIS","CAU")),MLmodel=factor(MLmodel,levels=c("RF","SVR","KNN")))

tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/Fig1.tiff", 
     width = 16, height = 20, units = 'cm', res = 300, compression = 'lzw')

ggscatter(data, x="ENcvR2", y="cvR2", size=0.8, add = "reg.line", add.params = list(color="red", size=1), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 2.5) + geom_abline(intercept=0, slope=1, color="blue", lwd=1) + 
  xlim(-0.5,1) + ylim(-0.5,1) + xlab(expression(paste("Elastic Net ", R^{2}))) + 
  ylab(expression(paste("Other Machine Learning Model ", R^{2}))) + theme_bw(16) + facet_grid(pop~MLmodel)

dev.off()


##############
##Make a table with the p-values of the t-test between EN vs ML across each sub population

tb <- matrix(nrow=4,ncol=3)

pops <- c("ALL", "AFA", "CAU", "HIS")
algs <- c("RF", "SVR", "KNN")

colnames(tb) <- algs
rownames(tb) <- pops


for (i in 1:length(pops)){
  
  for (j in 1:length(algs)){
    
    ml <- subset(data, data$MLmodel==algs[j] & data$pop==pops[i])
    tt <- t.test(ml[["ENcvR2"]], ml[["cvR2"]], paired=TRUE)
    tb[i,j] <- tt$p.value
  }
}

#Calculate the mean prediction performance
tm <- matrix(nrow=4,ncol=4)

pops <- c("ALL", "AFA", "CAU", "HIS")
algs <- c("RF", "SVR", "KNN")

colnames(tm) <- c("EN", "RF", "SVR", "KNN")
rownames(tm) <- pops


for (i in 1:length(pops)){
  
  for (j in 1:length(algs)){
    
    ml <- subset(data, data$MLmodel==algs[j] & data$pop==pops[i])
    tm[i,1] <- mean(ml[["ENcvR2"]]) 
    tm[i,j+1] <- mean(ml[["cvR2"]])
  }
} #Note that EN vs RF in AFA, mean of EN is actually 0.05282505, but that's almost same with what we have already in the table


#I can also compare the slope as suggested by reviewer 1

# Another approach that works when you're comparing a slope with 1 (not other numbers) is to compare the linear regression model with 
# one that has an offset. An offset is a term you add to a model where the slope is not estimated, but instead fixed at 1.
# Example R code:
# r.x <- lm(y ~ x)
# r1 <- lm(y ~1 + offset(x))
# anova(r.x, r1)

ml <- subset(data, data$MLmodel=="KNN" & data$pop=="ALL")

#E.g comparison of the null hypothesis of slope = 1 versus slope != 1 in Figure 1
# Not sure if this is correct
lmdt <- cbind(ml[["ENcvR2"]], ml[["cvR2"]])
lmdt <- as.data.frame(lmdt)
reg <- lm(V2~V1, lmdt)
reg1 <- lm(V2~1 + offset(V1), lmdt)
res_aov <- anova(reg,reg1)
summary(res_aov)
