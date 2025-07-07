##DOSE EXPOSURE STUDY SCRIPT ###
#code for analyzing size data for blastula, gastrula and pluteus stages
#Madison Armstrong 
#Last Modified 7/7/2025

#with the defaults
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)

#figure theme to keep things consistent
theme_box <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 12), 
          panel.border = element_rect(color = "black", fill = NA, size = 1)) 
}

###BLASTULA####
#Import the blastula measurements, select and rename useful columns, and
#add a column for experiments 1 and 2
blast_measure<-read.csv("Dose_exposure/Blastula_measurements.csv") 
#view(blast_measure)
blast_measure %>% dplyr::select(BeakerID, MatePair, Treatment, Horizontal.mm.) %>%
  rename(size=Horizontal.mm., Beaker=BeakerID) %>%
  mutate(Experiment=ifelse(MatePair==1 | MatePair == 2, 1, 2))-> blast_measure
  
#Factor Treatment and MatePair so the anova tests work
blast_measure$Treatment<-factor(blast_measure$Treatment, levels = c("Control","100ppb","500ppb","1000ppb"))
blast_measure$MatePair<-factor(blast_measure$MatePair, levels=c("1", "2", "3", "4"))

#take averages for blastula size
blast_measure_avg<-blast_measure %>% dplyr::group_by(Treatment, MatePair) %>% 
  dplyr::summarise(mean = mean(size), sd = sd(size))

write.csv(blast_measure_avg, "Dose_exposure/blast_avg.csv")
#Box plots for treatment, mate pair, and experiment
size<-
  ggplot(blast_measure_avg, aes(x=Treatment, y=mean, color=MatePair)) +
  geom_point(size=2)+
  geom_line(aes(color= MatePair, group=MatePair))+
  geom_point(aes(group = MatePair,color = MatePair),size = 2.3) +
  geom_errorbar(aes(x = Treatment, ymin=mean-sd, ymax=mean+sd, color = MatePair), 
                width = .1, linewidth=0.1) + 
  labs(title = '', x = 'Nonylphenol Concentration', y = 'Blastula Size (nm)')+
  theme_box() +
  scale_color_brewer(palette="Dark2")

ggsave("Dose_exposure/figs/Blastula Size/blast_size.png",size, width=18, height=16, units = "cm") 

blast<-glm(size~Treatment*MatePair, data=blast_measure) #not on means
#relationship of size and treatment with mate pair as a random effect
summary(blast) 
anova(blast)

subblast<-glm(size~Treatment*MatePair, data=sub_blast_measure) 
summary(subblast) 

anova(subblast)

#####GASTRULA####
gast_measure<-read.csv("Dose_exposure/Gastrula_measurements.csv") 
gast_measure %>% dplyr::select(BeakerID, MatePair, Treatment, Height.mm., Stomach_Length.mm.) %>%
  mutate(Experiment=ifelse(MatePair==1 | MatePair==2, 1, 2))-> gast_measure

#Factor Treatment and MatePair so the anova tests work
gast_measure$Treatment<-factor(gast_measure$Treatment, levels = c("Control","100ppb","500ppb","1000ppb"))
gast_measure$MatePair<-factor(gast_measure$MatePair, levels=c("1", "2", "3", "4"))

#subset diff measurements into diff datasheets to take averages
#height
Hgast_measure<-gast_measure %>% 
  select(BeakerID, MatePair, Treatment, Height.mm.) %>% 
  drop_na()
Hgastsum<-Hgast_measure %>% dplyr::group_by(Treatment, MatePair) %>% 
  dplyr::summarise(Heightmean = mean(Height.mm.), sd = sd(Height.mm.)) %>% 
  drop_na(sd) #this removes any datapoints that only had one individual aka no sd

write.csv(Hgastsum, "Dose_exposure/gastheight_avg.csv")
#stomach length
SLgast_measure<-gast_measure %>% 
  select(BeakerID, MatePair, Treatment, Stomach_Length.mm.) %>% 
  drop_na()
SLgastsum<-SLgast_measure %>% dplyr::group_by(Treatment, MatePair) %>% 
  dplyr::summarise(SLmean = mean(Stomach_Length.mm.), sd = sd(Stomach_Length.mm.))%>% 
  drop_na(sd) #this removes any datapoints that only had one individual aka no sd
write.csv(SLgastsum, "Dose_exposure/gastSL_avg.csv")
#ggplots!!

gast.height<-
ggplot(Hgastsum, aes(x=Treatment, y=Heightmean, color=MatePair)) +
  geom_point(size=2)+
  geom_line(aes(color= MatePair, group=MatePair))+
  geom_point(aes(group = MatePair,color = MatePair),size = 2.3) +
  geom_errorbar(aes(x = Treatment, ymin=Heightmean-sd, ymax=Heightmean+sd, color = MatePair), 
                width = .1, linewidth=0.1) + 
  labs(title = '', x = 'Nonylphenol Concentration', y = 'Gastrula Height (nm)')+
  theme_box()+
  scale_color_brewer(palette="Dark2")

gast.stomach<-
ggplot(SLgastsum, aes(x=Treatment, y=SLmean, color=MatePair)) +
  geom_point(size=2)+
  geom_line(aes(color= MatePair, group=MatePair))+
  geom_point(aes(group = MatePair,color = MatePair),size = 2.3) +
  geom_errorbar(aes(x = Treatment, ymin=SLmean-sd, ymax=SLmean+sd, color = MatePair), 
                width = .1, linewidth=0.1) + 
  labs(title = '', x = 'Nonylphenol Concentration', y = 'Gastrula Stomach Length (nm)')+
  theme_box()+
  scale_color_brewer(palette="Dark2")

fullgast.figs <- grid.arrange(gast.height, gast.stomach, nrow = 2)

ggsave("fullgast.figs.png",fullgast.figs, width=16, height=18, units = "cm") 

gast<-glm(Height.mm.~Treatment*MatePair, data=gast_measure) 
anova(gast)
gast2<-glm(Stomach_Length.mm.~Treatment*MatePair, data=gast_measure) 
anova(gast2)


####PLUTEUS####
plut_measure<-read.csv("Dose_exposure/pluteus7dpf_measurements.csv") 

plut_measure <- mutate(plut_measure, AL_mm = ifelse(is.na(AL1_mm), ifelse(is.na(AL2_mm), NA, AL2_mm), 
                                                  ifelse(is.na(AL2_mm), AL1_mm, (AL1_mm + AL2_mm)/2)))

plut_measure %>% dplyr::select(BeakerID, MatePair, Treatment, AL_mm, BL_mm, stomach_mm) %>%
  mutate(Experiment=ifelse(MatePair==1 | MatePair==2, 1, 2))-> plut_measure

#Factor Treatment and MatePair so the anova tests work
plut_measure$Treatment<-factor(plut_measure$Treatment, levels = c("Control","100ppb","500ppb","1000ppb"))
plut_measure$MatePair<-factor(plut_measure$MatePair, levels=c("1", "2", "3", "4"))

#subset diff measurements into diff datasheets to take averages
BLplut_measure<-plut_measure %>% 
  select(BeakerID, MatePair, Treatment, BL_mm) %>% 
  drop_na()
BLplutsum<-BLplut_measure %>% dplyr::group_by(Treatment, MatePair) %>% 
  dplyr::summarise(BLmean = mean(BL_mm), sd = sd(BL_mm)) %>% 
  drop_na(sd) #this removes any datapoints that only had one individual aka no sd

write.csv(BLplutsum, "Dose_exposure/plutBL_avg.csv")

ALplut_measure<-plut_measure %>% 
  select(BeakerID, MatePair, Treatment, AL_mm) %>% 
  drop_na()
ALplutsum<-ALplut_measure %>% dplyr::group_by(Treatment, MatePair) %>% 
  dplyr::summarise(ALmean = mean(AL_mm), sd = sd(AL_mm))%>% 
  drop_na(sd) #this removes any datapoints that only had one individual aka no sd

write.csv(ALplutsum, "Dose_exposure/plutAL_avg.csv")

SAplut_measure<-plut_measure %>% 
  select(BeakerID, MatePair, Treatment, stomach_mm) %>% 
  drop_na()
SAplutsum<-SAplut_measure %>% dplyr::group_by(Treatment, MatePair) %>% 
  dplyr::summarise(SAmean = mean(stomach_mm), sd = sd(stomach_mm))%>% 
  drop_na(sd) #this removes any datapoints that only had one individual aka no sd

write.csv(SAplutsum, "Dose_exposure/plutSA_avg.csv")

#now let's make some graphs!!
plutBL<-
  ggplot(BLplutsum, aes(x=Treatment, y=BLmean, color=MatePair)) +
  geom_point(size=2)+
  geom_line(aes(color= MatePair, group=MatePair))+
  geom_point(aes(group = MatePair,color = MatePair),size = 2.3) +
  geom_errorbar(aes(x = Treatment, ymin=BLmean-sd, ymax=BLmean+sd, color = MatePair), 
                width = .1, linewidth=0.1) + 
  labs(title = '', x = 'Nonylphenol Concentration', y = 'Pluteus Body Length (nm)')+
  theme_box()+
    scale_color_brewer(palette="Dark2")

plutAL<-ggplot(ALplutsum, aes(x=Treatment, y=ALmean, color=MatePair)) +
  geom_point(size=2)+
  geom_line(aes(color= MatePair, group=MatePair))+
  geom_point(aes(group = MatePair,color = MatePair),size = 2.3) +
  geom_errorbar(aes(x = Treatment, ymin=ALmean-sd, ymax=ALmean+sd, color = MatePair), 
                width = .1, linewidth=0.1) + 
  labs(title = '', x = 'Nonylphenol Concentration', y = 'Pluteus Arm Length (nm)')+
  theme_box()+
  scale_color_brewer(palette="Dark2")

plutSA<-
  ggplot(SAplutsum, aes(x=Treatment, y=SAmean, color=MatePair)) +
  geom_point(size=2)+
  geom_line(aes(color= MatePair, group=MatePair))+
  geom_point(aes(group = MatePair,color = MatePair),size = 2.3) +
  geom_errorbar(aes(x = Treatment, ymin=SAmean-sd, ymax=SAmean+sd, color = MatePair), 
                width = .1, linewidth=0.1) + 
  labs(title = '', x = 'Nonylphenol Concentration', y = 'Pluteus Stomach Area (nm^2)')+
  theme_box()+
  scale_color_brewer(palette="Dark2")

fullplut.figs <- grid.arrange(plutBL, plutAL, plutSA, nrow = 3)

ggsave("fullplut.figs.png",fullplut.figs, width=16, height=18, units = "cm") 

##stats
#Plut BL stats w mp2
plut1<-glm(BL_mm~Treatment*MatePair, data=plut_measure) 
anova(plut1)
#plut AL stats
plut2<-glm(AL_mm~Treatment*MatePair, data=plut_measure) 
anova(plut2)
#plut SA stats
plut3<-glm(stomach_mm~Treatment*MatePair, data=plut_measure)
summary(plut3)
anova(plut3)

#some beakers only have one individual that was measured so I need to remove those for stats
beaker_counts <- table(sub_plut_measure$BeakerID)
single_beakers <- names(beaker_counts[beaker_counts == 1]) #B_B10
# Identify beakers where only one type of measurement is recorded
measurement_counts <- aggregate(cbind(BL_mm, AL_mm, stomach_mm) ~ BeakerID, data = sub_plut_measure, 
                                function(x) sum(!is.na(x)))  # Count non-NA values
# Find beakers where only one type of measurement is recorded
single_measurement_beakers <- measurement_counts$BeakerID[rowSums(measurement_counts[, -1] > 0) == 1]
# Combine both filters
beakers_to_remove <- unique(c(single_beakers, single_measurement_beakers)) #still just beaker 10
# Subset the data to remove those beakers
sub_plut_measure <- sub_plut_measure[sub_plut_measure$BeakerID !="B_B10",]

#Plut BL stats
subplut<-glm(BL_mm~Treatment*MatePair, data=sub_plut_measure) 
summary(subplut) 
anova(subplut)
#plut AL stats
subplut2<-glm(AL_mm~Treatment*MatePair, data=sub_plut_measure) 
summary(subplut2) 
anova(subplut2)
#plut SA stats
subplut3<-glm(stomach_mm~Treatment*MatePair, data=sub_plut_measure) 
summary(subplut3) 
anova(subplut3)
