library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(data.table)

enassocman <- fread(file="Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/en_assoc_manplot_df.txt", header = T)
rfassocman <- fread(file="Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/rf_assoc_manplot_df.txt", header = T)
svrassocman <- fread(file="Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/svr_assoc_manplot_df.txt", header = T)
knnassocman <- fread(file="Z:/data/ml_paper_reviewers_corrections/twas_mesa/hdl/knn_assoc_manplot_df.txt", header = T)

#####Make the manhattan plot
#this script for mamhattan plot below was gotten from
#https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function
#I only made a few changes


en_df.tmp <- enassocman %>% group_by(CHR) %>% summarise(chr_len=as.numeric(max(BP))) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% 
  select(-chr_len) %>% 
  left_join(enassocman, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse(SNP %in% NA, "yes", "no")) %>% mutate( is_annotate=ifelse(P < 3.3e-6, "yes", "no"))
en_axisdf <- en_df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


rf_df.tmp <- rfassocman %>% group_by(CHR) %>% summarise(chr_len=as.numeric(max(BP))) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% 
  select(-chr_len) %>% 
  left_join(rfassocman, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse(SNP %in% NA, "yes", "no")) %>% mutate( is_annotate=ifelse(P < 3.3e-6, "yes", "no"))
rf_axisdf <- rf_df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


svr_df.tmp <- svrassocman %>% group_by(CHR) %>% summarise(chr_len=as.numeric(max(BP))) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% 
  select(-chr_len) %>% 
  left_join(svrassocman, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse(SNP %in% NA, "yes", "no")) %>% mutate( is_annotate=ifelse(P < 3.3e-6, "yes", "no"))
svr_axisdf <- svr_df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


knn_df.tmp <- knnassocman %>% group_by(CHR) %>% summarise(chr_len=as.numeric(max(BP))) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% 
  select(-chr_len) %>% 
  left_join(knnassocman, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse(SNP %in% NA, "yes", "no")) %>% mutate( is_annotate=ifelse(P < 3.3e-6, "yes", "no"))
knn_axisdf <- knn_df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )



#Plot the Manhattan Plots
##, 
#text=element_text(size=8)) + 

en_gg <- ggplot(en_df.tmp, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
  scale_color_manual(values = rep(c("red","blue"), 22 )) + 
  scale_x_continuous( label = en_axisdf$CHR, breaks= en_axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14), breaks = seq(0,14,2)) + ggtitle(paste0("Elastic Net")) + 
  labs(x = "Chromosome") +
  geom_hline(yintercept = -log10(5e-6), colour="red") +
  geom_label_repel(data=en_df.tmp[en_df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP)), size=5, force=1.3) + 
  theme_classic(base_size = 20) + theme(legend.position = "none", axis.text.x = element_text(hjust=0, angle=90, size=8)) +
  ylab(expression(~~-log[10](italic(p))))


rf_gg <- ggplot(rf_df.tmp, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
  scale_color_manual(values = rep(c("red","blue"), 22 )) + 
  scale_x_continuous( label = rf_axisdf$CHR, breaks= rf_axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14), breaks = seq(0,14,2)) + ggtitle(paste0("Random Forest")) + 
  labs(x = "Chromosome") +
  geom_hline(yintercept = -log10(5e-6), colour="red") +
  geom_label_repel(data=rf_df.tmp[rf_df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP)), size=5, force=1.3) + 
  theme_classic(base_size = 20) + theme(legend.position = "none", axis.text.x = element_text(hjust=0, angle=90, size=8)) + 
  ylab(expression(~~-log[10](italic(p))))


svr_gg <- ggplot(svr_df.tmp, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
  scale_color_manual(values = rep(c("red","blue"), 22 )) + 
  scale_x_continuous( label = svr_axisdf$CHR, breaks= svr_axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14), breaks = seq(0,14,2)) + ggtitle(paste0("Support Vector")) + 
  labs(x = "Chromosome") +
  geom_hline(yintercept = -log10(5e-6), colour="red") +
  geom_label_repel(data=svr_df.tmp[svr_df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP)), size=5, force=1.3) + 
  theme_classic(base_size = 20) + theme(legend.position = "none", axis.text.x = element_text(hjust=0, angle=90, size=8)) + 
  ylab(expression(~~-log[10](italic(p))))


knn_gg <- ggplot(knn_df.tmp, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
  scale_color_manual(values = rep(c("red","blue"), 22 )) + 
  scale_x_continuous( label = knn_axisdf$CHR, breaks= knn_axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14), breaks = seq(0,14,2)) + ggtitle(paste0("K Nearest Neighbor")) + 
  labs(x = "Chromosome") +
  geom_hline(yintercept = -log10(5e-6), colour="red") +
  theme_classic(base_size = 20) + theme(legend.position = "none", axis.text.x = element_text(hjust=0, angle=90, size=8)) + 
  ylab(expression(~~-log[10](italic(p))))


#plot all 4 plots in one page
library(gridExtra)
grid.arrange(en_gg,rf_gg,svr_gg,knn_gg,nrow=2)
tiff("Z:/data/ml_paper_reviewers_corrections/paper_figs/Fig5_manplot.tiff", 
     width = 34, height = 18, units = 'cm', res = 300, compression = 'lzw')
grid.arrange(en_gg,rf_gg,svr_gg,knn_gg,nrow=2)
dev.off()
