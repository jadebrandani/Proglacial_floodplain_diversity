#libraries 
library(ggpubr)
library(ggrepel)
library(rstatix)
library(gridExtra)
library(Hmisc)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(stringr)
library(readr)
library(data.table)
library(permute)
library(tabula)
library(tidyr)
######
# Here we are working with a file containing all water physio chemical information "water". 

#### PCA   ----
Stream <-c("#73DFFF","#FFA300")
#Val Roseg 
water$Is_kryal<-as.factor(water$Is_kryal)
VAR_water<-subset(water,water$Glacier=="ValRoseg")
VAR.param<-VAR_water[,c(5:13,28)]
var.pca <- PCA(VAR.param, graph = FALSE)
library("factoextra")
fviz_pca_biplot(var.pca, repel = TRUE,
                habillage=VAR_water$Is_kryal,
                palette=Stream,
                col.var="black",
                geom.ind = "point",
                label = "all",
                pointsize = 8,
                labelsize = 6,
                arrowsize = 0.5,axes = c(1,2),
                axes.linetype = "dashed", 
                title= "")+
  scale_shape_manual(values=c(15,15))+
  theme(text = element_text(size = 18),legend.position = "none")+
  xlim(-4, 4.5)+
  ylim (-4, 4.5)
# Valsorey 
SOY_water<-subset(water,water$Glacier=="Valsorey")
SOY.param<-SOY_water[,c(5:13,28)]
soy.pca <- PCA(SOY.param, graph = FALSE)
library("factoextra")
fviz_pca_biplot(soy.pca, repel = TRUE,
                habillage=SOY_water$Is_kryal,
                palette=Stream,
                col.var="black",
                geom.ind = "point",
                label = "all",
                pointsize = 8,
                labelsize = 6,
                arrowsize = 0.5,axes = c(1,2),
                axes.linetype = "dashed", 
                title= "")+
  scale_shape_manual(values=c(15,15))+
  theme(text = element_text(size = 18),legend.position = "none")+
  xlim(-4, 4.5)+
  ylim (-4, 4.5)
# Otemma 
OTE_water<-subset(water,water$Glacier=="Otemma")
OTE.param<-OTE_water[,c(5:13,28)]
ote.pca <- PCA(OTE.param, graph = FALSE)
fviz_pca_biplot(ote.pca, repel = TRUE,
                habillage=OTE_water$Is_kryal,
                palette=Stream,
                col.var="black",
                geom.ind = "point",
                label = "all",
                pointsize = 8,
                labelsize = 6,
                arrowsize = 0.5,axes = c(1,2),
                axes.linetype = "dashed", 
                title= "")+
  scale_shape_manual(values=c(15,15))+
  theme(text = element_text(size = 18),legend.position = "none")+
  xlim(-4, 4.5)+
  ylim (-4, 4.5)

#### ANOVA -------
#Testing for each parameter () Season, Glacier or Stream type) explain the variation in the following water physiochemical parameters
temp.aov <- aov(water$Temp~ water$Season + water$Glacier +water$Is_kryal, data = water)
summary(temp.aov)
cond.aov <- aov(water$Conductivity~ water$Season + water$Glacier +water$Is_kryal, data = water)
summary(cond.aov)
DOC.aov <- aov(water$DOC~ water$Season + water$Glacier +water$Is_kryal, data = water)
summary(DOC.aov)
DIN.aov <- aov(water$DIN~ water$Season + water$Glacier +water$Is_kryal, data = water)
summary(DIN.aov)
SRP.aov <- aov(water$SRP~ water$Season + water$Glacier +water$Is_kryal, data = water)
summary(SRP.aov)


#### Longitudinal Gradient analysis only on Glacier-fed samples -------
water_GFS <-water [which(water$Is_kryal=="Glacier-fed"),]
water_GFS <-water_GFS [which(water_GFS$DIstance_category!="NEUTRAL"),]

library(tidyr)

# By Glacier 
#Valsorey
water_GFS_SOY<- water_GFS [which(water_GFS$Glacier=="Valsorey"),]
mydata.long <- water_GFS_SOY[c(5:11,24:27, 30)] %>%
  pivot_longer(-DIstance_category, names_to = "variables", values_to = "value")
mydata.long %>% sample_n(6)

stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ DIstance_category) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

water_GFS_SOY_summary<- water_GFS_SOY[,c(5:11,24:27)] %>%
  group_by( water_GFS_SOY$DIstance_category) %>%
  get_summary_stats( type = "mean_sd")

#Otemma
water_GFS_OTE<- water_GFS [which(water_GFS$Glacier=="Otemma"),]
mydata.long <- water_GFS_OTE[c(5:11,24:27, 30)] %>%
  pivot_longer(-DIstance_category, names_to = "variables", values_to = "value")
mydata.long %>% sample_n(6)

stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ DIstance_category) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

water_GFS_OTE_summary<- water_GFS_OTE [,c(5:11,24:27)] %>%
  group_by(water_GFS_OTE$DIstance_category) %>%
  get_summary_stats( type = "mean_sd")

# VAL ROSEG 
water_GFS_VAR<- water_GFS [which(water_GFS$Glacier=="ValRoseg"),]
mydata.long <- water_GFS_VAR[c(5:11,24:27, 30)] %>%
  pivot_longer(-DIstance_category, names_to = "variables", values_to = "value")
mydata.long 

stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ DIstance_category) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

water_GFS_VAR_summary<- water_GFS_VAR[,c(5:11,24:27)] %>%
  group_by( water_GFS_VAR$DIstance_category) %>%
  get_summary_stats( type = "mean_sd")

