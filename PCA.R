library(vegan)
library(tidyverse)
library(FD)
library(gridExtra)
library(ggplot2)
library(ggpubr)

## PCA for environmental variables by treatment, visualized by  ##
## Followed by extraction of PCA scores, regressed against STL and ETL and SOC
## NOTE: data from LTM_manuscript_2.R
data14a <- data14 %>% 
  filter(Time=="during") %>%
  mutate(ID = paste(Block,Treatment, Yearf, sep = "_"))

env<-data14a  %>% 
  dplyr::select(moisture, pH, conductivity, SOC2, Chl_a)

row.names(env) <- data14a$ID

myrda2 <- rda(na.omit(env), scale = TRUE)

# extract values
siteout2 <- as.data.frame(scores(myrda2, choices=c(1,2), display=c("sites")))
siteout2$ID<-rownames(siteout2)
siteout2$name <- siteout2$ID

enviroout2<-as.data.frame(scores(myrda2, choices=c(1,2), display=c("species")))
enviroout2$type<-"variable"
enviroout2$name<-rownames(enviroout2)

# merge PC axes with trait data
tog <- left_join(data14a, siteout2, by="ID") 

a.pca<-ggplot(subset(tog, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W"), aes(x=PC1, y=PC2))+ 
  #ggtitle("a) PCA on environmental variables")+
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  #geom_text(aes(label = name, color = Treatment), size = 2) +
  geom_point(aes(color=Treatment, shape=Treatment), size=4)+
  scale_shape_manual(values=c(15, 16, 17, 8), 
                     guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Control", "Water", "Mannitol", "Sucrose"))+ #change labels in the legend)+
  scale_color_manual(values = c("black", "#8d4a85", "#63ba95", "#a78b47"), labels=c("Control", "Water", "Mannitol", "Sucrose"), guide = guide_legend(title = "Treatment")) +
  geom_segment(data = enviroout2,
               aes(x = 0, xend =  PC1,
                   y = 0, yend =  PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_text(data = enviroout2,
            aes(x=  PC1*1.2, y =  PC2*1.2, #we add 10% to the text to push it slightly out from arrows
                label = name), #otherwise you could use hjust and vjust
            size = 4,
            hjust = 0.5, 
            color="black") + 
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 20), legend.position="none")+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda2$CA$eig["PC1"]/myrda2$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda2$CA$eig["PC2"]/myrda2$tot.chi*100,3),"%)",sep="")) 
a.pca

#linear regression by subplot
fit1 <- lm(log(STL+1) ~ PC1, data=subset(tog, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W")) 
summary(fit1) #intercept= 6.15, slope=-1.12, r2=0.25
as.vector(fit1$coefficients)
fit2 <- lm(log(STL+1) ~ PC2, data=subset(tog, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W")) 
summary(fit2) #intercept=5.85, slope = -1.13, r2=0.31
as.vector(fit2$coefficients)
fit3 <- lm(log(ETL+1) ~ PC1, data=subset(tog, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W")) 
summary(fit3) #intercept=3.15, slope = 1.29, r2=0.17
fit4 <- lm(log(ETL+1) ~ PC2, data=subset(tog, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W")) 
summary(fit4) #intercept=3.40, slope = 0.55, r2=0.04

#equations for figures
lbl1 <- expression("y = 6.15 - 1.12 x," ~ r^2 ~ "= 0.25")
lbl2 <- expression("y = 5.85 - 1.13 x," ~ r^2 ~ "= 0.31")
lbl3 <- expression("y = 3.15 + 1.29 x," ~ r^2 ~ "= 0.17")
lbl4 <- expression("y = 3.40 + 0.55 x," ~ r^2 ~ "= 0.04")

a.1<-ggplot(data=subset(tog, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W"), 
            aes(x=PC1, y=log(STL+1), group=Treatment, color=Treatment, shape=Treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #ggtitle("a)")+
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  scale_shape_manual(values=c(15, 16, 17, 8), 
                     guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Control", "Water", "Mannitol", "Sucrose"))+ #change labels in the legend)+
  scale_color_manual(values=c( "black", "#8d4a85", "#63ba95", "#a78b47"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Control", "Water", "Mannitol", "Sucrose"))+ #change labels in the legend)+
  xlab("PC1 Scores")+
  ylab("STL")+
  #xlim(40,100)+
  #ylim(0,1000)+
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
  annotate("text", x = 0.5, y = 9, label = lbl1, size = 3, color = "black")
a.1

a.2<-ggplot(data=subset(tog, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W"), 
            aes(x=PC2, y=log(STL+1), group=Treatment, color=Treatment, shape=Treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #ggtitle("a)")+
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  scale_shape_manual(values=c(15, 16, 17, 8), 
                     guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Control", "Water", "Mannitol", "Sucrose"))+ #change labels in the legend)+
  scale_color_manual(values=c("black", "#8d4a85", "#63ba95", "#a78b47"), 
                     guide = guide_legend(title = "Treatment"), #change legend title
  labels=c("Control", "Water", "Mannitol", "Sucrose"))+ #change labels in the legend)+
  xlab("PC2 Scores")+
  ylab("STL")+
  #xlim(40,100)+
  #ylim(0,1000)+
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
  annotate("text", x = 0, y = 9, label = lbl2, size = 3, color = "black")
a.2

b.1<-ggplot(data=subset(tog, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W"), 
            aes(x=PC1, y=log(ETL+1), group=Treatment, color=Treatment, shape=Treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #ggtitle("a)")+
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  scale_shape_manual(values=c(15, 16, 17, 8), 
                     guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Control", "Water", "Mannitol", "Sucrose"))+ #change labels in the legend)+
  scale_color_manual(values=c( "black", "#8d4a85", "#63ba95", "#a78b47"), guide = guide_legend(title = "Treatment"), #change legend title
  labels=c("Control", "Water", "Mannitol", "Sucrose"))+ #change labels in the legend)+
  xlab("PC1 Scores")+
  ylab("ETL")+
  #xlim(40,100)+
  #ylim(0,1000)+
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
  annotate("text", x = -0.75, y = 5, label = lbl3, size = 3, color = "black")
b.1

b.2<-ggplot(data=subset(tog, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W"), 
            aes(x=PC2, y=log(ETL+1), group=Treatment, color=Treatment, shape=Treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #ggtitle("a)")+
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  scale_shape_manual(values=c(15, 16, 17, 8), 
                     guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Control", "Water", "Mannitol", "Sucrose"))+ #change labels in the legend)+
  scale_color_manual(values=c("black", "#8d4a85", "#63ba95", "#a78b47"), guide = guide_legend(title = "Treatment"), #change legend title
  labels=c("Control", "Water", "Mannitol", "Sucrose"))+ #change labels in the legend)+
  xlab("PC2 Scores")+
  ylab("ETL")+
  #xlim(40,100)+
  #ylim(0,1000)+
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
  annotate("text", x = -1, y = 5, label = lbl4, size = 3, color = "black")
b.2

p.1<-ggarrange(a.1, a.2, b.1, b.2, 
               ncol=2, nrow=2, common.legend=TRUE, legend="bottom",
               labels = c("b)", "c)", "d)", "e)"))
p.1

PCA_fig <- ggarrange(a.pca, p.1,
                       ncol =1, nrow =2, 
                       labels = c("a)",""))
PCA_fig


## PCA for environmental variables by treatment, visualized by  ##
## Followed by extraction of PCA scores, regressed against STL and ETL and SOC
## NOTE: data from LTM_manuscript_2.R
data14b <- data14 %>% 
  filter(Time=="recovery") %>%
  mutate(ID = paste(Block,Treatment, Yearf, sep = "_"))

env_b<-data14b  %>% 
  dplyr::select(moisture, pH, conductivity,  SOC2, Chl_a)

row.names(env_b) <- data14b$ID

myrda2_b <- rda(na.omit(env_b), scale = TRUE)

# extract values
siteout2b <- as.data.frame(scores(myrda2_b, choices=c(1,2), display=c("sites")))
siteout2b$ID<-rownames(siteout2b)
siteout2b$name <- siteout2b$ID

enviroout2b<-as.data.frame(scores(myrda2_b, choices=c(1,2), display=c("species")))
enviroout2b$type<-"variable"
enviroout2b$name<-rownames(enviroout2b)

# merge PC axes with trait data
tog_b <- left_join(data14b, siteout2b, by="ID") 

b.pca<-ggplot(subset(tog_b, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W"), aes(x=PC1, y=PC2))+ 
  ggtitle("a) PCA on environmental variables")+
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_text(aes(label = name, color = Treatment), size = 5) +
  #scale_color_manual(values = c("mediumpurple", "seagreen3","sienna1"), labels=c("Forb","Grass", "Legume"), guide = guide_legend(title = "Functional Group")) +
  geom_segment(data = enviroout2,
               aes(x = 0, xend =  PC1,
                   y = 0, yend =  PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_text(data = enviroout2,
            aes(x=  PC1*1.2, y =  PC2*1.2, #we add 10% to the text to push it slightly out from arrows
                label = name), #otherwise you could use hjust and vjust
            size = 6,
            hjust = 0.5, 
            color="black") + 
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 20), legend.position="none")+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda2$CA$eig["PC1"]/myrda2$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda2$CA$eig["PC2"]/myrda2$tot.chi*100,3),"%)",sep="")) 
b.pca


c.1<-ggplot(data=subset(tog_b, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W"), aes(x=PC1, y=log(STL+1), group=Treatment, color=Treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #ggtitle("a)")+
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
  #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("PC1 Scores")+
  ylab("STL")+
  #xlim(40,100)+
  #ylim(0,1000)+
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
c.1

c.2<-ggplot(data=subset(tog_b, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W"), aes(x=PC2, y=log(STL+1), group=Treatment, color=Treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #ggtitle("a)")+
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
  #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("PC2 Scores")+
  ylab("STL")+
  #xlim(40,100)+
  #ylim(0,1000)+
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
c.2

d.1<-ggplot(data=subset(tog_b, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W"), aes(x=PC1, y=log(ETL+1), group=Treatment, color=Treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #ggtitle("a)")+
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
  #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("PC1 Scores")+
  ylab("ETL")+
  #xlim(40,100)+
  #ylim(0,1000)+
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
d.1

d.2<-ggplot(data=subset(tog_b, Treatment=="WM"|Treatment=="WS"|Treatment=="C"|Treatment=="W"), aes(x=PC2, y=log(ETL+1), group=Treatment, color=Treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #ggtitle("a)")+
  geom_point()+
  theme_bw()+
  #theme(legend.position="none")+
  #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
  #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("PC2 Scores")+
  ylab("ETL")+
  #xlim(40,100)+
  #ylim(0,1000)+
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
d.2
