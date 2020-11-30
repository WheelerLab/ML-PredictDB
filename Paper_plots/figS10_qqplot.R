library(dplyr)
library(ggplot2)
library(viridis)
library(data.table)

#new qq for the reviewers 
newqq <- fread(file="Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/qq_input_ALL.txt", header=T, stringsAsFactors=F)
#rearrange and ML order
newqq <- mutate(newqq,model=factor(model,levels=c("Elastic Net","Random Forest","Support Vector","K Nearest Neighbor")))

tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/FigS10_qqplot.tiff", 
     width = 16, height = 16, units = 'cm', res = 300, compression = 'lzw')

ggplot(newqq,aes(x=Expected,y=Observed)) + geom_point(shape=1) + theme_classic(16) + 
  scale_colour_manual(values=c(rgb(163,0,66,maxColorValue = 255),"dark gray")) + 
  theme(legend.position = c(0.01,0.99),legend.justification = c(0,1)) + 
  geom_abline(slope=1,intercept = 0, colour="red") + 
  xlab(expression(Expected~~-log[10](italic(p)))) + 
  ylab(expression(Observed~~-log[10](italic(p)))) +
  expand_limits(x = 0, y = 0) + scale_y_continuous(breaks = seq(0,14,2), limits = c(0,14)) + facet_wrap(~model)

dev.off()
