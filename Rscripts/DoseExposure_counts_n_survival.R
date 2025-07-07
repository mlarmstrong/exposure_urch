##DOSE EXPOSURE STUDY SCRIPT ###
#code for assessing fertilization success, developmental success and survival of S. purpuratus larvae
#Madison Armstrong 
#Last Modified 7/7/2025

setwd("~/Desktop/purp dev data")

library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(lme4)
library(emmeans)
library(wesanderson)
library(patchwork)

#figure theme to keep things consistent
theme_box <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 20), 
          panel.border = element_rect(color = "black", fill = NA, size = 1)) 
}

#Fertilization## 
  #Fertilization#----
fert<-read.csv("Dose_Exposure/DE_Fertilized.csv")

fert<-fert %>% dplyr::select (1:7) %>% 
  filter(!is.na(Treatment), !is.na(Beaker)) %>% 
  rename(Two_cell=X2.cell)

fert <- fert %>% 
  mutate(sum.fert = rowSums(across(6:7), na.rm = T)) %>% 
  mutate(total= rowSums(across(5:7), na.rm=T)) %>% 
  mutate (prop.fert=(sum.fert/total)) 
 

fert.summary<- fert %>% 
  group_by(MatePair,Beaker, Treatment) %>% 
  summarize(mean.prop=mean(prop.fert), se.prop=(sd(prop.fert)/sqrt(length(prop.fert)))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "100 ppb", "500 ppb", "1000 ppb")))


fert.summary$MatePair<- factor(x=fert.summary$MatePair, labels=c(
  "1" ="Pair 1",
  "2" = "Pair 2",
  "3" = "Pair 3",
  "4" = "Pair 4"
))

#Ff<-
  ggplot(fert.summary, aes(x=Treatment, y=mean.prop, color=MatePair)) +
    geom_point(size=2)+ ylim(0,1)+
    geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Fertilization Success', x = 'Nonylphenol Concentration', y = 'Proportion of Fertilized Embryos')+
    theme_box()+
  scale_color_brewer(palette="Dark2")


#Stats: are these groups statically different?
#ANOVA https://statsandr.com/blog/anova-in-r/
Fe<-lm(mean.prop~Treatment+MatePair, data=fert.summary)
anova(Fe) # no significant differences in proportion fertilized due to treatment or MP

  #2-cell Dev Check#----
dev<-read.csv("Dose_Exposure/DE_Developed.csv")
#View(dev)

dev <- dev %>% 
  mutate(total= rowSums(across(6:7), na.rm=T)) %>% #that's why only 6:7 rather than column 4
  mutate (prop.dev=(Developed/total))

dev.summary<- dev %>% 
  group_by(MatePair,Beaker, Treatment) %>% 
  summarize(mean.prop=mean(prop.dev), se.prop=(sd(prop.dev)/sqrt(length(prop.dev)))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "100 ppb", "500 ppb", "1000 ppb")))
#View(dev.summary)

dev.summary$MatePair<- factor(x=dev.summary$MatePair, labels=c(
  "1" ="Pair 1",
  "2" = "Pair 2",
  "3" = "Pair 3",
  "4" = "Pair 4"
))

#DD<-
  ggplot(dev.summary, aes(x=Treatment, y=mean.prop, color=MatePair))+
    geom_point(size=2)+ ylim(0,1)+
    geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Developmental Success', x = 'Nonylphenol Concentration', y = 'Proportion of Developed Embryos')+
    theme_box()+
  scale_color_brewer(palette="Dark2")


fullcount <- Ff + DD + plot_layout(guides="collect")
ggsave("fullcount.png",fullcount, width=30, height=18, units = "cm")
#save figure in large format--> higher quality too!
ggsave("DE_Dev.png",DD, width=16, height=20, units = "cm") 

#Stats: are these groups statically different?
#ANOVA
De=aov(mean.prop~Treatment+MatePair, data=dev.summary)
summary(De) #no variation in development success


#Survival ----
#  looking at overall shift: blastula to 7dpf Pluteus

#read in blastula file
blast<-read.csv("Dose_Exposure/DE_survival_blastula.csv")

#get mean count data per mL and make summary datasheet
blast.summary<- blast %>% 
  group_by(Stage, Beaker, Treatment) %>% 
  summarize(mean.blast.per.mL= mean((Count)*2), blast.se.prop=(sd(Count)/sqrt(length(Count))))

#read in pluteus file
plut<-read.csv("Dose_Exposure/DE_survival_plut.csv")

#get mean count data per mL and make summary datasheet
plut.summary<- plut %>% 
  group_by(MatePair,Stage, Beaker, Treatment) %>% 
  summarize(mean.plut.per.mL =mean((Count)*2), plut.se.prop=(sd(Count)/sqrt(length(Count))))

#merge two datasets together
BP<-merge(blast.summary,plut.summary, by="Beaker")
#View(BP) #looks like it merged successfully 

#remove unwanted columns, like duplicate columns
BP<-BP %>% select(-c(Treatment.y))

#rename Treatment column 
BP<-BP %>% dplyr::select (1:9) %>% 
  rename(Treatment=Treatment.x)


#add new column for survival differences between blast and plut
BP <- BP %>% 
  mutate (prop_survival=((mean.plut.per.mL-mean.blast.per.mL)/mean.blast.per.mL)) %>% 
  mutate(Treatment = fct_relevel(Treatment.x, c("Control", "100 ppb", "500 ppb", "1000 ppb")))

BP$MatePair<- factor(x=BP$MatePair, labels=c(
  "1" ="Pair 1",
  "2" = "Pair 2",
  "3" = "Pair 3",
  "4" = "Pair 4"
))
write.csv(BP, "meancounts_survival.csv")


survivalblast_plut<-
ggplot(BP, aes(x=Treatment, y=prop_survival, color=MatePair)) +
  #geom_line(aes(color= MatePair, group=MatePair))+
  geom_point(size=4)+
  labs(title = '', x = 'Nonylphenol Concentration', 
       y = 'Proportion Survival')+
  theme_box()+
  scale_color_brewer(palette="Dark2")


ggsave("blast_plut_survival-nolines.png",survivalblast_plut, width=30, height=20, units = "cm") 

BP_stats<-lm(prop_survival~Treatment+MatePair, data=BP)
summary(BP_stats) 
anova(BP_stats) 

#write csv and let's get the comparison of survival from the control to treatments for each dose in excel
write.csv(BP, "Dose_exposure/BlastPlut_avgsurvival-summary.csv")
#in excel= look at difference in treatment - control for each dose at each matepair == prop_survival_2control

survival2C<-read.csv("Dose_Exposure/BlastPlut_avgsurvival-summary2control.csv")

survival2C <- survival2C %>% 
  mutate(Treatment = fct_relevel(Treatment, c("Control", "100 ppb", "500 ppb", "1000 ppb")))

#View(survival2C)
#relsurvival<-
  ggplot(survival2C, aes(x=Treatment, y=prop_2control, color=MatePair)) +
  geom_line(aes(color= MatePair, group=MatePair))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, color="gray")+ #no change in survival from control to treatment within mate pair
  labs(title = '', x = 'Nonylphenol Concentration', 
       y = 'Proportion Survival (Relative)')+
  theme_box()+
  scale_color_brewer(palette="Dark2")
ggsave("blast_plut_relsurvival.png",relsurvival, width=18, height=16, units = "cm") 


#boxplot
#survival_boxplot<-
ggplot(survival2C, aes(x=Treatment, y=prop_2control, fill=Treatment)) +
  geom_boxplot()+
  geom_hline(yintercept = 0, color="gray")+
  labs(title = '', x = '',
       y = 'Proportion Survival (relative)') + theme_box()+
  scale_fill_brewer(palette="Dark")
#ggsave("boxplot_survival.png",survival_boxplot, width=20, height=13, units = "cm") 


relcontrol_stats<-glm(prop_2control~Treatment+MatePair, data=survival2C)
summary(relcontrol_stats) 
anova(relcontrol_stats) 

#t.tests with rel2C
sub <- survival2C%>%filter(Treatment=="100 ppb")
t.test(sub$prop_2control)


#subset out mate pair 2
sub<-survival2C[-c (5:8),]
#View(sub)

substatsBP<-glm(prop_2control~Treatment+(MatePair), data=sub)
summary(substatsBP) 
anova(substatsBP)
subanovaBPM <-aov(prop_survival~Treatment, data=sub)
AIC(substatsBP, subanovaBPM)

TukeyHSD(subanovaBPM)
