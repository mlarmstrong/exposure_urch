#Dose Exposure---Abnormalities
#code for assessing developmental abnormalities in the three developmental stages of interest: Blastula, Gastrula and Pluteus
#Madison Armstrong 
#Last Modified 7/7/2025

setwd("~/Desktop/purp development")
#if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("patchwork"))
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(lme4)
library(emmeans)
library(RColorBrewer)
library(patchwork)

#figure theme to keep things consistent
theme_box <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 20), 
          panel.border = element_rect(color = "black", fill = NA, size = 1)) 
}
#Blastula####
abn_countB<-read.csv("Dose_exposure/Blastula_abn_counts.csv") #count data

abncountB<- abn_countB %>% 
  group_by(MatePair,BeakerID, Treatment) %>% 
  mutate (prop.abn=(abn_count/total)) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("Control", "100ppb", "500ppb", "1000ppb")))

abncountB$MatePair<- factor(x=abncountB$MatePair, labels=c(
  "1" ="Pair 1",
  "2" = "Pair 2",
  "3" = "Pair 3",
  "4" = "Pair 4"
))

blast_abncounts2<-
  ggplot(abncountB, aes(x=Treatment, y=prop.abn, color=MatePair)) +
    geom_point(size=2)+
  geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Blastula', x = 'Nonylphenol Concentration',
       y = 'Proportion of Abnormal Blastula') + theme_box()+
    scale_color_brewer(palette="Dark2")

ggsave("blastsbn2.png",blast_abncounts2, width=20, height=13, units = "cm") 

#since its a proportion you can do anova
Blm<-glm(prop.abn~Treatment+MatePair, data=abncountB)
summary(Blm) # sig diffs by mate pair 2
anova(Blm)

#write csv and let's get the comparison of blast abn from the control to treatments for each dose in excel
write.csv(abncountB, "Dose_exposure/abncountB-summary.csv")
#in excel= look at difference in treatment - control for each dose at each matepair == prop.abn.2C

abnB2C<-read.csv("Dose_Exposure/abncountB-summary2C.csv")

abnB2C<- abnB2C %>%
  group_by(MatePair,BeakerID, Treatment) %>%
  mutate(Treatment = fct_relevel(Treatment, "Control", "100ppb", "500ppb", "1000ppb"))


#blast_abncounts2control<-
  ggplot(abnB2C, aes(x=Treatment, y=prop.abn.2C, color=MatePair)) +
  geom_point(size=2)+
  geom_line(aes(color= MatePair, group=MatePair))+
    geom_hline(yintercept = 0, color="gray")+
  labs(title = 'Proportion Abnormal Blastula Relative to Control Across Treatments', x = 'Nonylphenol Concentration',
       y = 'Proportion of Abnormal Blastula') + theme_box()+
  scale_color_brewer(palette="Dark2")

ggsave("blastsbn2control.png",blast_abncounts2control, width=20, height=13, units = "cm") 

BlastC<-glm(prop.abn.2C~Treatment+MatePair, data=abnB2C)
summary(BlastC) # slight diff at MP 2 (0.0819)
anova(BlastC)
#take out mate pair 2 (LOTS of Abn)
sub_abnB2C<-abnB2C[-c (5:8),]
View(sub_abnB2C)
subblastC<-glm(prop.abn.2C~Treatment+MatePair, data=sub_abnB2C)
summary(subblastC)
anova(subblastC)

#Gastrula####
abn_countG <-read.csv("Dose_exposure/Gastrula_abn_counts.csv") #count data

abncountG<- abn_countG %>% 
  group_by(MatePair,BeakerID, Treatment) %>% 
  mutate (prop.abn=(abn_count/total)) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("Control", "100 ppb", "500 ppb", "1000 ppb")))

abncountG$MatePair<- factor(x=abncountG$MatePair, labels=c(
  "1" ="Pair 1",
  "2" = "Pair 2",
  "3" = "Pair 3",
  "4" = "Pair 4"
))

gastLM<-glm(prop.abn~Treatment+MatePair, data=abncountG)
summary(gastLM)
anova(gastLM)

#write csv and let's get the comparison of blast abn from the control to treatments for each dose in excel
write.csv(abncountG, "Dose_exposure/abncountG-summary.csv")
#in excel= look at difference in treatment - control for each dose at each matepair == prop.abn.2C

abnG2C<-read.csv("Dose_Exposure/abncountG-summary2C.csv")

abnG2C<- abnG2C %>%
  group_by(MatePair,BeakerID, Treatment) %>%
  mutate(Treatment = fct_relevel(Treatment, "Control", "100ppb", "500ppb", "1000ppb"))


#gast_abncounts2control<-
ggplot(abnG2C, aes(x=Treatment, y=prop.abn.2C, color=MatePair)) +
  geom_point(size=2)+
  geom_line(aes(color= MatePair, group=MatePair))+
  geom_hline(yintercept = 0, color="gray")+
  labs(title = 'Proportion Abnormal Gastrula Relative to Control Across Treatments', x = 'Nonylphenol Concentration',
       y = 'Proportion of Abnormal Gastrula') + theme_box()+
  scale_color_brewer(palette="Dark2")

#ggsave("gastsbn2control.png",gast_abncounts2control, width=20, height=13, units = "cm") 

GastC<-glm(prop.abn.2C~Treatment+MatePair, data=abnG2C)
summary(GastC) # slight diff at MP 2 (0.0819)
anova(GastC)
sub_abnG2C<-abnG2C[-c (5:8),]
View(sub_abnG2C)
subgastC<-glm(prop.abn.2C~Treatment+MatePair, data=sub_abnG2C)
summary(subgastC)
anova(subgastC)
#Pluteus####
abn_countP<-read.csv("Dose_exposure/Pluteus7dpf_abn-counts.csv") #count data

#count data first
abncountP<- abn_countP %>% 
  group_by(MatePair,BeakerID, Treatment) %>% 
  mutate (prop.abn=(abn_count/total)) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("Control", "100 ppb", "500 ppb", "1000 ppb")))

abncountP$MatePair<- factor(x=abncountP$MatePair, labels=c(
  "1" ="Pair 1",
  "2" = "Pair 2",
  "3" = "Pair 3",
  "4" = "Pair 4"
))

#plut_abncounts<-
ggplot(abncountP, aes(x=Treatment, y=prop.abn, color=MatePair)) +
  geom_line(aes(color= MatePair, group=MatePair))+
  geom_point(size=2)+
  labs(title = 'Proportion Abnormal Pluteus Across Treatments', x = 'Nonylphenol Concentration',
       y = 'Proportion of Abnormal Pluteus') +
  theme_box() + scale_color_manual(values = wes_palette(n=4, name="BottleRocket2"))

ggsave("plut_abncounts.png",plut_abncounts, width=20, height=13, units = "cm") 

#since its a proportion you can do anova
Plm<-glm(prop.abn~Treatment+MatePair, data=abncountP)
summary(Plm) 
anova(Plm)

#write csv and let's get the comparison of blast abn from the control to treatments for each dose in excel
write.csv(abncountP, "Dose_exposure/abncountP-summary.csv")
#in excel= look at difference in treatment - control for each dose at each matepair == prop.abn.2C

###relative to control###
abnP2C<-read.csv("Dose_Exposure/abncountP-summary2C.csv")

abnP2C<- abnP2C %>%
  group_by(MatePair,BeakerID, Treatment) %>%
  mutate(Treatment = fct_relevel(Treatment, "Control", "100ppb", "500ppb", "1000ppb"))

#plut_abncounts2control<-
ggplot(abnP2C, aes(x=Treatment, y=prop.abn.2C, color=MatePair)) +
  geom_point(size=2)+
  geom_line(aes(color= MatePair, group=MatePair))+
  geom_hline(yintercept = 0, color="gray")+
  labs(title = 'Proportion Abnormal Pluteus Relative to Control Across Treatments', x = 'Nonylphenol Concentration',
       y = 'Proportion of Abnormal Pluteus') + theme_box()+
  scale_color_brewer(palette="Dark2")

ggsave("plutsbn2control.png",plut_abncounts2control, width=20, height=13, units = "cm") 

#boxplot
plutabn_boxplot<-ggplot(abnP2C, aes(x=Treatment, y=prop.abn.2C, fill=Treatment)) +
         geom_boxplot()+
  geom_hline(yintercept = 0, color="gray")+
  labs(title = 'Proportion Abnormal Pluteus Relative to Control Across Treatments', x = 'Nonylphenol Concentration',
       y = 'Proportion of Abnormal Pluteus') + theme_box()+
  scale_fill_brewer(palette="Dark")

ggsave("boxplot_plutsbn2control.png",plutabn_boxplot, width=20, height=13, units = "cm") 

m1<-glm(prop.abn.2C~Treatment, data=abnP2C)
m2<-glm(prop.abn.2C~Treatment+MatePair, data=abnP2C)

AIC(m1, m2) #including mate pair is a better fit
summary(m2) #0.0504 treatment 100ppb sig for abn!!
anova(m2)

#t.tests for rel2C data#####

#blastula-- no sig
x <- abnB2C %>% filter(Treatment=="1000ppb")
t.test(x$prop.abn.2C)

#gastrula
y <- abnG2C %>% filter(Treatment=="1000ppb")
t.test(y$prop.abn.2C)

#pluteus
z <- abnP2C %>% filter(Treatment=="1000ppb")
t.test(z$prop.abn.2C)

#ALL ABN####
#look at all abnormalities across stages at once
abn <- read.csv("Dose_exposure/morphAbn_Summary_data.csv")

abn$prop.abn.bl <- abn$abn_blastula/abn$total_blastula
abn$prop.abn.ga <- abn$abn_gastrula/abn$total_gastrula
abn$prop.abn.pl <- abn$abn_pluteus/abn$total_pluteus

mod <- lm(prop.abn.pl~Treatment+MatePair,data=abn)
anova(mod)

abn.norm <- abn %>%
  group_by(MatePair) %>%
  mutate(
    control_value = prop.abn.pl[Treatment == "Control" & row_number() == which(Treatment == "Control")],
    normalized_value = prop.abn.pl-control_value
  ) %>%
  ungroup()

abn.norm$Treatment <- factor(abn.norm$Treatment,levels=c("Control","100ppb","500ppb","1000ppb"))
abn.norm$MatePair<- factor(x=abn.norm$MatePair, labels=c(
  "1" ="Pair 1",
  "2" = "Pair 2",
  "3" = "Pair 3",
  "4" = "Pair 4"
))

anova(lm(prop.abn.bl~Treatment+as.factor(MatePair),data=abn.norm))
anova(lm(normalized_value~Treatment+MatePair,data=abn.norm))
anova(lm(prop.abn.ga~Treatment+as.factor(MatePair),data=abn.norm))
anova(lm(prop.abn.pl~Treatment+MatePair,data=abn.norm))

# figures in paper####
bl <- ggplot(abn.norm, aes(x=Treatment, y=prop.abn.bl, color=MatePair)) +
  geom_point(size=4)+
  geom_hline(yintercept = 0, color="white")+
  #geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Blastula', x = 'Nonylphenol Concentration',
       y = 'Proportion of Abnormal Blastula') + theme_box()+
  scale_color_brewer(palette="Dark2")

ga <- ggplot(abn.norm, aes(x=Treatment, y=prop.abn.ga, color=MatePair)) +
  geom_point(size=4)+
  #geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Gastrula', x = 'Nonylphenol Concentration',
       y = 'Proportion of Abnormal Gastrula') + theme_box()+
  scale_color_brewer(palette="Dark2")


pl <- ggplot(abn.norm, aes(x=Treatment, y=prop.abn.pl, color=MatePair)) +
  geom_point(size=4)+
  #geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Pluteus', x = 'Nonylphenol Concentration',
       y = 'Proportion of Abnormal Pluteus') + theme_box()+
  scale_color_brewer(palette="Dark2")


full.abn<-bl + ga + pl + plot_layout(guides="collect")
ggsave("all.abn-nolines.png",full.abn, width=40, height=20, units = "cm") 
