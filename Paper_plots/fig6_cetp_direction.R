library(data.table)
library(dplyr)
library(ggplot2)


en <- fread(file="Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/en_cetp_pc3.txt", header = T)
rf <- fread(file="Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/rf_cetp_pc3.txt", header = T)
svr <- fread(file="Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/svr_cetp_pc3.txt", header=T)
knn <- fread(file="Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/knn_cetp_pc3.txt", header=T)

#use facet_wrap
en$method <- "Elastic Net"
rf$method <- "Random Forest"
svr$method <- "Support Vector"
knn$method <- "K Nearest Neighbor"

df <- rbind(en,rf,svr, knn)
fig <- ggplot(data=df, aes(x=CETP, y=Phenotype)) + geom_point(lwd=2) + xlab("Predicted gene expression") + 
  ylab("HDL (rank normalized)") + geom_smooth(method="lm", color="red", lwd=2, se=TRUE) + geom_density_2d(lwd=1.5) + 
  theme_classic(20) + scale_x_continuous(breaks=seq(-2.0, 2.0, 1)) + facet_wrap(~method)
print(fig)


tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/Fig6_cetp.tiff", width = 24, height = 19, 
     units = 'cm', res = 300, compression = 'lzw')
fig
dev.off()
