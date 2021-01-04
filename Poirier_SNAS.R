
# Poirier, Waterhouse, Dunn & Smith, SN Applied Sciences 2021

# Libraries 
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(grid)
library(cowplot)
library(car)
library(lme4)


# SNAS Datasets
# Number of compounds
dede<-read.csv(file.choose(), header=T, sep=";")
dede[,c(1,2,3,4,6)]<- lapply(dede[,c(1,2,3,4,6)], factor)
head(dede)

dede2<- melt(xtabs(~Spl+Ind+Tr+Ex, dede))
dede2<- dede2[which(dede2$value!=0),]
head(dede2)

Torder2<- c("5-extractions_short-delay","5-extractions_long-delay",
            "2-extractions_long-delay")
dede2<- transform(dede2, Tr= factor(Tr, levels=Torder2))


# Lost/Gained
LG<- read.csv(file.choose(), header=T, sep=";", dec=",")
LG[,c(1,2,3,4,5)]<- lapply(LG[,c(1,2,3,4,5)], factor)
head(LG)

Torder<- c("5 Extractions_Short delay","5 Extractions_Long delay",
           "2 Extractions_Long delay")
LG2<- transform(LG, Treatment2= factor(Treatment2, levels= Torder))
head(LG2)


# SNAS paper fig1a
Fig1a<-
ggplot(dede2, aes(x=Tr, y=value, fill=Ex))+
  geom_boxplot(outlier.shape=1, outlier.size=1)+    
  scale_fill_manual(values=c("gray35","firebrick3","red","orange","yellow1"),
                    labels=c(expression("1"^{st}), expression("2"^{nd}),
                             expression("3"^{rd}), expression("4"^{th}),
                             expression("5"^{th})))+
  xlab("Experimental condition")+ 
  ylab("Number of compounds detected in the samples")+
  labs(fill="Extraction")+
  labs(title=expression(paste(bold("a."))))+
  scale_x_discrete(labels=c("5 Extractions\nShort delay",
                            "5 Extractions\nLong delay", 
                            "2 Extractions\nMaximum delay"))+
  theme_classic(base_size=11)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="bottom")


# SNAS paper Fig1b
Fig1b<-
ggplot(LG2, aes(x=factor(Treatment2), y=Percent, fill=ELG, stat="identity"))+
  geom_bar(position=position_stack(reverse=T), colour="black", 
           stat="identity", width=0.5)+
  geom_hline(yintercept=0, colour="black")+
  scale_y_continuous(breaks=seq(-40, 10, 10))+
  scale_x_discrete(labels=c("5 Extractions\nShort delay",
                            "5 Extractions\nLong delay", 
                            "2 Extractions\nMaximum delay"))+
  scale_fill_manual(values=c("darkblue","firebrick3","royalblue","red","deepskyblue","orange",
                             "lightcyan","yellow1"),
                    labels=c(expression("2"^{nd}), expression("2"^{nd}),
                             expression("3"^{rd}), expression("3"^{rd}),
                             expression("4"^{th}), expression("4"^{th}),
                             expression("5"^{th}), expression("5"^{th})))+
  xlab("Experimental condition")+ 
  ylab("Cumulative proportion of compounds 
  gained and lost at successive extractions (in %)")+
  labs(fill="Extraction   Gained \n                    Lost")+
  labs(title=expression(paste(bold("b."))))+
  theme_classic(base_size=11)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="bottom")


abba <- plot_grid(Fig1a, Fig1b, rel_heights=c(1,1), ncol=2, align="h")
ggsave(plot=abba, filename="D:/ALICE_2019-20/Alice-Rprojects/July2019/Paper3_Stability_results/Fig1_Oct20_600dpi.tiff", 
       device="tiff", scale=1, dpi=600, units="cm", width=31, height=14)



# SNAS model
momo3<- glmer(value~Tr+Ex+(1|Ind/Spl), data=dede2, family=poisson(link="log"))
AIC(momo3)
summary(momo3)
vif(momo3)
qqnorm(residuals(momo3))
plot(residuals(momo3))


# DHARMa residuals check
library(DHARMa)

mymodel=momo3
simulationOutput<- simulateResiduals(fittedModel=mymodel, n=10000, refit=F, use.u=F)
plot(simulationOutput)



