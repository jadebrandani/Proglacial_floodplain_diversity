library(car)
library(rstatix)
library(tidyr)

##### data import ----
biomass <- read.csv("data/Biomass_data.csv")
biomass$ID<-as.factor(biomass$ID)
biomass$Glacier<-as.factor(biomass$Glacier)
biomass$Stream<-as.factor(biomass$Stream)

### ANOVA ------
#In order to test how Season - Glacier- and Stream type have an influence on biomass indicators we can use an ANOVA
# Bacterial Abundance
res.aov <- aov(biomass$Bacterial_Abundance  ~ Season + Glacier +Stream, data = biomass)
summary(res.aov)
TukeyHSD(res.aov)
# BCP
biomass$BCP<-biomass$BCP*1000
res.aov <- aov(biomass$BCP  ~ Season +Glacier +Stream, data = biomass)
summary(res.aov)
TukeyHSD(res.aov)
# CHLA
res.aov <- aov(biomass$Chla  ~ Season + Glacier +Stream, data = biomass)
summary(res.aov)
TukeyHSD(res.aov)
# EPS
res.aov <- aov(biomass$EPS ~ Season + Glacier +Stream, data = biomass)
summary(res.aov)
TukeyHSD(res.aov)

#EPS_ per cell
biomass$BCP_cell<-(biomass$BCP/biomass$Bacterial_Abundance)*1000000
res.aov <- aov(biomass$BCP_cell  ~ Season + Glacier +Stream, data = biomass)
summary(res.aov)
t.test(biomass$BCP_cell ~ biomass$Stream, alt="two.sided", conf=0.95, var.eq=F, paired=F)
