###viral tree#######
library(ggtree)
library(ggnewscale)
library(ggplot2)
V_tree<-read.tree('tree.txt')
info_V <- read.csv("info.csv")
vcluster_cat_D1<- read.csv("D1_cat.csv",row.names = 1)
vcluster_cat_D15<- read.csv("D15_cat.csv",row.names = 1)
viraltree<-ggtree(V_tree, branch.length='none',layout = 'circular') %<+% info_V + geom_tippoint(size=0.2)+
  geom_hilight(node=481, fill="gold",alpha=.3)+
  geom_hilight(node=648, fill="purple",alpha=.3)+
  geom_hilight(node=656, fill="green",alpha=.3)
viraltree_new<-open_tree(viraltree, 5)
fig <- viraltree_new + new_scale_fill()
fig_D1<-gheatmap(fig, vcluster_cat_D1, offset=0.8, width=.15, 
                #border_color='#00BFC4', 
                colnames=F) +
  scale_fill_manual(values =c("#F28127","#D3D3D3",'#FFFFFF',"#0c00e8"), name="Phage-host interaction")+new_scale_fill()
fig_D1D15<-gheatmap(fig_D1, vcluster_cat_D15, offset=10, width=.15, 
                 #color='#F8766D',
                 colnames=F) +
  scale_fill_manual(values =c("#F28127","#D3D3D3",'#FFFFFF',"#0c00e8"), name="Phage-host interaction")+
  theme(legend.position = 'left')
#######box plot####
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
df <-read.csv("abcd.csv", header = TRUE)
stat.test <- df %>%
  group_by(Variable3) %>%
  t_test(Value ~ Sample) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")%>%
  add_xy_position(x = "Sample", dodge = 0.8)
stat.test
abcd <- ggboxplot(df, x = 'Sample', y = 'Value', 
                       palette = c('#00BFC4','#F8766D'), color = 'Sample')+ 
  stat_pvalue_manual(stat.test,  label = "p.signif", size=2.5)+
  facet_wrap(.~Variable3, scale='free', ncol=4,strip.position='left')+
  theme(axis.text.x = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        axis.text=element_text(size=6))
df_FM <-read.csv("ef.csv", header = TRUE)
df_D <-read.csv("gh.csv", header = TRUE)
stat.test_FM <- df_FM %>%
  group_by(Variable3) %>%
  t_test(Value ~ Variable1) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")%>%
  add_xy_position(x = "Variable1", dodge = 0.8)
stat.test_FM
stat.test_D <- df_D %>%
  group_by(Variable3) %>%
  t_test(Value ~ Variable1) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")%>%
  add_xy_position(x = "Variable1", dodge = 0.8)
stat.test_D
ef <- ggboxplot(df_FM, x = 'Variable1', y = 'Value', 
                     palette = c('#00BFC4'), color = 'Sample')+ 
  stat_pvalue_manual(stat.test_FM,  label = "p.signif", size=2.5)+
  facet_wrap(.~Variable3, scale='free', ncol=2,strip.position='left')+
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        axis.text=element_text(size=6),
        legend.position='none')
gh <- ggboxplot(df_D, x = 'Variable1', y = 'Value', 
                     palette = c('#F8766D'), color = 'Sample')+ 
  stat_pvalue_manual(stat.test_D,  label = "p.signif", size=2.5)+
  facet_wrap(.~Variable3, scale='free', ncol=2,strip.position='left')+
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        axis.text=element_text(size=6),
        legend.position='none')
library(gridExtra)
grid.arrange(abcd, arrangeGrob(ef,gh, ncol=2),ncol=1)
#####viral copy per host###
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
d <-read.csv("D1_D15_VPH.csv", header = TRUE)
stat.test <- d %>%
  t_test(CopyPerCell ~ Sample) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test_new <- stat.test %>% add_xy_position(x = "Sample")
p.adj_new<- scientific(stat.test_new$p.adj, digits = 3)

library(RColorBrewer)
library(colorRamps)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
library(scales)
show_col(getPalette(15))
plt1<-ggplot(d, aes(x = factor(Sample,level=c('Wet','Dry')), y = CopyPerCell)) +
  geom_boxplot(aes(y = CopyPerCell),outlier.shape = NA) +
  geom_jitter(aes(y = CopyPerCell, color=BinTax), width = 0.2, size=0.8)+
  scale_color_manual(values =getPalette(15))+
  stat_pvalue_manual(stat.test_new,  label = "p.adj.signif", 
                     tip.length = 0, y.position = 1.48,size=3, vjust = 4, hjust = -1)+
  scale_y_continuous(name = "Viral copy per host (adjusted)")+
  theme_bw()+
  theme(axis.text.y = element_text(size=8, color=c('#00BFC4','#F8766D')),
        axis.text.x = element_text(size=6),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=8),
        legend.title = element_blank(), 
        legend.position = 'bottom',
        legend.text=element_text(size=rel(0.45),margin = margin(t = 0.1)),
        legend.key.size = unit(0.2, 'lines'),
        legend.spacing.x = unit(0.8, "cm"),
        plot.margin = unit(c(1, 1, 0, 1), "lines"),
        strip.text = element_text(size=8))+
  guides(color=guide_legend(ncol=2))+coord_flip()
g<- ggplot_gtable(ggplot_build(plt1))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c('#F8766D','#00BFC4')
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
library(tidyverse)
library(grid)

grid.draw(g)
df <-read.csv("D1_VPH.csv", header = TRUE)
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(data=df, aes(x=log(CopyPerCell), y=BinAbundance),fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(aes(colour =Sample), size = 2) +
    stat_smooth(method = "lm", se=TRUE, fullrange=FALSE, level=0.95) +
    scale_color_manual(values =c('#00BFC4'))+
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                       "Intercept =",signif(fit$coef[[1]],2 ),
                       " Slope =",signif(fit$coef[[2]], 2),
                       " P =",signif(summary(fit)$coef[2,4], 2)))+
    theme(legend.position = "none",plot.title = element_text(size=8), 
          axis.title=element_text(size=8), axis.text = element_text(size=6))+
    xlab('log(Viral copy per host (adjusted))')+
    ylab('Host abundance')
}
fit1 <- lm(log(CopyPerCell) ~ BinAbundance, data = df)
D1_reg<- ggplotRegression(fit1)
df <-read.csv("D15_VPH.csv", header = TRUE)
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(data=df, aes(x=log(CopyPerCell), y=BinAbundance),fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(aes(colour =Sample), size = 2) +
    stat_smooth(method = "lm", se=TRUE, fullrange=FALSE, level=0.95) +
    scale_color_manual(values =c('#F8766D'))+
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                       "Intercept =",signif(fit$coef[[1]],2 ),
                       " Slope =",signif(fit$coef[[2]], 2),
                       " P =",signif(summary(fit)$coef[2,4], 2)))+
    theme(legend.position = "none", plot.title = element_text(size=8), 
          axis.title=element_text(size=8), axis.text = element_text(size=6))+
    xlab('log(Viral copy per host (adjusted))')+
    ylab('Host abundance')
}
fit2 <- lm(log(CopyPerCell) ~ BinAbundance, data = df)
D15_reg<-ggplotRegression(fit2)
library(gridExtra)
grid.arrange(plt1, D1_reg, D15_reg,ncol=1)

