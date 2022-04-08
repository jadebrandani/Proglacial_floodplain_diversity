#Libraries ----
library(ggpubr)
install.packages("gridExtra")
library(ggrepel)
library(gridExtra)
library(FactoMineR)
library(factoextra)
library(ggplot2); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("cluster"); packageVersion("cluster")
library("ape"); packageVersion("ape")
library("vegan"); packageVersion("vegan")
library(reshape)
library("reshape2")
install.packages("tabula")
library(gridExtra)
library(RColorBrewer)
library(stringr)
library(phyloseq); packageVersion('phyloseq')
library("readr")
library(data.table)
library(reshape2)
library(betapart)
library(ade4)
library(lattice)
library(permute)
library(tabula)
library("geosphere") 

Stream <-c("#73DFFF","#FFA300")
### Data Import ----
### 16S Data Import 
#16S data is already pre-filtered. Singleton have been removed and we also removed Chloroplast and mitochondria. 
#The data is in **absolute abundance** (read count * Bacterial abundance)
sDNA_16S_ASV <- read.csv("data/16S_table.csv", row.names = 1, header = T)
sDNA_16S_ASV<-as.matrix(sDNA_16S_ASV)
row.names.remove1<- c("VAR_61", "VAR_62") 
sDNA_16S_ASV <-sDNA_16S_ASV [,!(colnames(sDNA_16S_ASV) %in% row.names.remove1)]

Taxonomy <- read.csv("data/16S_Taxonomy.csv", row.names=1, header=T) #Taxonomy object 
Taxonomy <- as.matrix(Taxonomy)
names(Taxonomy)<-make.names(names(Taxonomy))

Metadata <- read.csv("data/16S_metadata.csv", row.names = 1) #Metadata
Metadata$Season_pairs<-as.character(Metadata$Season_pairs)
row.names.remove<- c("OTE_27", "OTE_48", "VAR_61", "VAR_62") #removing samples that were not sequenced
Metadata <-Metadata [!(row.names(Metadata) %in% row.names.remove),]
rm(row.names.remove)
Metadata$ID<-row.names(Metadata)

OTU = phyloseq::otu_table(sDNA_16S_ASV, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(Taxonomy)
META<- phyloseq::sample_data(Metadata)
data_16S <- phyloseq(OTU, TAX, META) # Phyloseq object 

# 18S data import ####### This data has NOT been pre-filtered before hand 

#OTU object 
OTU_18S_RAW <- read.delim("data/18S_otu_table.txt", row.names = 1, header = T)
OTU_18S_RAW <-as.matrix(OTU_18S_RAW)
row.names.remove1<- c("VAR_61", "VAR_62") 
OTU_18S_RAW  <-OTU_18S_RAW  [,!(colnames(OTU_18S_RAW) %in% row.names.remove1)]

#Taxonomy object 
Tax_18S_RAW <- read.delim("data/18S_taxonomy_NEW.txt", row.names=1, header=T)
Tax_18S_RAW <- as.matrix(Tax_18S_RAW)
names(Tax_18S_RAW)<-make.names(names(Tax_18S_RAW))

Metadata_18S <- read.csv("data/18S_Metadata.csv", row.names = 1)
Metadata_18S$Season_pairs<-as.character(Metadata_18S$Season_pairs)
row.names.keep<-colnames(OTU_18S_RAW)
Metadata_18S<-Metadata_18S [(row.names(Metadata_18S) %in% row.names.keep),]
row.names.remove<- c("OTE_27", "OTE_48", "VAR_61", "VAR_62") #removing samples that were not sequenced
Metadata_18S<-Metadata_18S [!(row.names(Metadata_18S) %in% row.names.remove),]

rows_18S<-row.names(Metadata_18S)
Subset_16S_match18S<- Metadata[row.names(Metadata)%in% rows_18S,c(10,11)]
Metadata_18S<-merge(Metadata_18S, Subset_16S_match18S, by=0)
row.names(Metadata_18S)<-Metadata_18S$Row.names

OTU = phyloseq::otu_table(OTU_18S_RAW, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(Tax_18S_RAW)
META<- phyloseq::sample_data(Metadata_18S)
data_18S_RAW <- phyloseq(OTU, TAX, META) 

#filtering 18S ------
  #Subset by taxonomy (Selecting Photo-autotrophs exclusively)
  data_18S_filtered = subset_taxa(data_18S_RAW,
                                  D_2=="Chloroplastida"| 
                                    D_3=="Ochrophyta"| 
                                    D_2=="Cryptomonadales")

#Removing OTUs with only read count = 1 
Low_read_18S<- prune_taxa(taxa_sums(data_18S_filtered)>1, data_18S_filtered)
#Removing the OTU with only one observation
Low_read_18S<-phyloseq::filter_taxa (Low_read_18S, function(x){sum(x > 0) > 1}, prune = TRUE)
# We remove sample with less than 1 OTU
data_18S<-prune_samples(sample_sums(Low_read_18S)>1, Low_read_18S)

############################################## Alpha diversity #####################################################-------

############################# LATERAL DIMENSION (TRIBUTARY vs. GLACIER-FED) ###############################----
######  16S 
Diversity_16S<- as.data.frame(diversity(t(sDNA_16S_ASV), index="shannon")) 
colnames(Diversity_16S)<-c("shannon")
Diversity_16S$Richness <- specnumber(t(sDNA_16S_ASV))
Diversity_16S$Evenness <- Diversity_16S$shannon/log(Diversity_16S$Richness)
Diversity_16S$Glacier<-Metadata$Glacier
Diversity_16S$Stream<-Metadata$Is_kryal
Diversity_16S$Campaign<-Metadata$Campaign
Diversity_16S$Distance_snout<-Metadata$DIstance_snout

library(tidyverse)
library(rstatix)   
library(ggpubr)


######  RICHNESS
richness.test <- Diversity_16S %>%
  group_by(Glacier) %>%
  t_test(Richness ~ Stream) %>%
  adjust_pvalue(method = "bonferroni") %>%
  mutate(y.position = 4000)
richness.test


ggplot(data=Diversity_16S)+
  geom_boxplot(aes(x=Stream, y=Richness, fill=Stream), color="black")+
  scale_fill_manual(values = Stream)+
  xlab("") + ylab("Richness") +
  facet_wrap(~Glacier)+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  ylim(500, 4000)+
  stat_pvalue_manual(richness.test, label = "p.adj")

######  EVENNESS
evenness.test <- Diversity_16S %>%
  group_by(Glacier) %>%
  t_test(Evenness ~ Stream) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(y.position = 0.95)
evenness.test

ggplot(data=Diversity_16S)+
  geom_boxplot(aes(x=Stream, y=Evenness, fill=Stream), color="black")+
  scale_fill_manual(values = Stream)+
  xlab("") + ylab("Richness") +
  facet_wrap(~Glacier)+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  ylim(0.7,0.95)+
  stat_pvalue_manual(evenness.test, label = "p.adj")


##### LONGITUDINAL DIMENSION ########
ggplot(data=Diversity_16S, aes(x=Distance_snout, y=Richness, shape = Glacier))+
  geom_point(color="#73DFFF", size = 5)+
  geom_smooth(method = "lm", se=FALSE, color="black", linetype = "dashed")+
  scale_shape_manual(values=c(19, 17, 15))+
  facet_grid(~Glacier, scales="free")+
  stat_regline_equation(label.x = c(10, 700), 
                        aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")))+
  xlab("") + ylab("Richness") +
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  ylim(500, 4000)+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3, label.y=800
  )+
  xlab("Distance from Snout Glacier (m)")


#### 18S 
sDNA_18S = as(otu_table(data_18S), "matrix")
Meta_18S = as(sam_data(data_18S), "matrix")
Meta_18S<-as.data.frame(Meta_18S)

Diversity_18S<- as.data.frame(diversity(t(sDNA_18S), index="shannon")) 
colnames(Diversity_18S)<-c("shannon")
Diversity_18S$Richness <- specnumber(t(sDNA_18S))
Diversity_18S$Evenness <- Diversity_18S$shannon/log(Diversity_18S$Richness)
Diversity_18S$Glacier<-Meta_18S$Glacier
Diversity_18S$Stream<-as.factor(Meta_18S$Stream)
Diversity_18S$Campaign<-Meta_18S$Campaign

###### RICHNESS
richness.test <- Diversity_18S %>%
  group_by(Glacier) %>%
  t_test(Richness ~ Stream) %>%
  adjust_pvalue(method = "bonferroni") %>%
  mutate(y.position = 100)
richness.test

ggplot(data=Diversity_18S)+
  geom_boxplot(aes(x=Stream, y=Richness, fill=Stream), color="black")+
  scale_fill_manual(values = Stream)+
  xlab("") + ylab("Richness") +
  facet_wrap(~Glacier)+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  ylim(0, 100)+
  stat_pvalue_manual(richness.test, label = "p.adj")

######  EVENNESS
evenness.test <- Diversity_18S %>%
  group_by(Glacier) %>%
  t_test(Evenness ~ Stream) %>%
  adjust_pvalue(method = "bonferroni") %>%
  mutate(y.position =1)
evenness.test

ggplot(data=Diversity_18S)+
  geom_boxplot(aes(x=Stream, y=Evenness, fill=Stream), color="black")+
  scale_fill_manual(values = Stream)+
  xlab("") + ylab("Evenness") +
  facet_wrap(~Glacier)+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  ylim(0,1)+
  stat_pvalue_manual(evenness.test, label = "p.adj")

####  VENN BAR PLOT  -------
#Split by glacier floodplain for both 16S and 18S
OTE_16S <-subset_samples(data_16S, Glacier == "Otemma")
VAR_16S  <-subset_samples(data_16S, Glacier == "Val Roseg")
SOY_16S<-subset_samples(data_16S, Glacier == "Valsorey")
OTE_18S <-subset_samples(data_18S, Glacier == "Otemma")
VAR_18S  <-subset_samples(data_18S, Glacier == "Val Roseg")
SOY_18S<-subset_samples(data_18S, Glacier == "Valsorey")

#Unsing the MicEco package to extract the values for both unweighted and abundance-weighted
library(MicEco)
OTE_Unweighted<-ps_venn(OTE_16S,
        group= "Is_kryal",
        fraction = 0,
        weight = F,
        type = "percent",
        relative = TRUE,
        plot = T)

OTE_Weighted<-ps_venn(OTE_16S,
               group= "Is_kryal",
               fraction = 0,
               weight = T,
               type = "percent",
               relative = TRUE,
               plot = T)
#repeated for other floodplain and for 18S dataset


# All results obtained from ps_venn is entered manually in a csv file and upload later as "Venn_16S" and "Venn_18S" 
## Importing venn results in a table and plotting as bar chart
### 16S
Venn_16S <- read.csv("data/Venn_16S.csv", header = T)
Venn_16S$Stream <- factor(Venn_16S$Stream , levels = c("Glacier-fed", "Tributary", "Shared"))

ggplot(Venn_16S,
       aes(x=Glacier, y=Proportion, fill=Stream))+
  geom_bar(stat="identity", position = "fill", color="black")+
  facet_wrap(~Type)+
  scale_y_continuous(labels=scales::percent, breaks = c(0,0.2,0.4, 0.6, 0.8,1))+
  scale_fill_manual(values = c("#73DFFF", "#FFA300", "lavender"))+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(strip.background =element_rect(fill="white"))+
  ggtitle("16S")

### 18S
Venn_18S <- read.csv("data/Venn_18S.csv", header = T)
str(Venn_18S$Stream)
Venn_18S$Stream <- factor(Venn_18S$Stream , levels = c("Glacier-fed", "Tributary", "Shared"))

ggplot(Venn_18S,
       aes(x=Glacier, y=Proportion, fill=Stream))+
  geom_bar(stat="identity", position = "fill", color="black")+
  facet_wrap(~Type)+
  scale_y_continuous(labels=scales::percent, breaks = c(0,0.2,0.4, 0.6, 0.8,1))+
  scale_fill_manual(values = c("#73DFFF", "#FFA300", "lavender"))+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(strip.background =element_rect(fill="white"))+
  ggtitle("18S")



#### TAXONOMIC COMPOSITION BAR PLOT -----
### 16S
glom_Order <- tax_glom(data_16S, taxrank = 'Order')
glom_Order # should list # taxa as # order
data_glom_Order<- psmelt(glom_Order) 
data_glom_Order$Order <- as.character(data_glom_Order$Order)

Test<-data_glom_Order%>%group_by(Sample, Order)%>%summarise(Sum=sum(Abundance), .groups='drop') 
Test$Order <-as.character(Test$Order)
Test$Order [Test$Order=="uncultured"]<-"Other"

Test_Fam14 = Test%>% group_by(Order)%>% summarise(total=sum(Sum))%>%top_n(14, total)
Test2 = Test
Test2$Sum[Test2$Sum==0] <- NA
Test2<- Test2[complete.cases(Test2),]
Test2$Order[!(Test2$Order %in% Test_Fam14$Order)] <- "Other"

Test3<- Test2 %>% group_by(Sample)%>% mutate(PercentFam=Sum/sum(Sum)*100)

Metadata$Sample<-rownames(Metadata)
Test4<-merge (Test3, Metadata, by="Sample")

unique(Test4$Order)

Test4$Order <- factor(Test4$Order, levels = c("Blastocatellales","Burkholderiales","Chitinophagales", "Cytophagales",
                                              "Flavobacteriales","Gemmatales", "Gemmatimonadales", "Nitrospirales",
                                              "Pirellulales","Rhizobiales" ,"Rhodobacterales","Sphingomonadales",  "Verrucomicrobiales",
                                              "Vicinamibacterales", "Other"))

colors_15= c("#3D518C","#b8def4","#06BEE1", "#b2df8a","#33a02c",
               "#F8E16C","#ff7f00","#fbc888","#e31a1c", "#fb9a99",   
               "#6a3d9a","#bab5e3","#ee2f65","#68473C",  "#c8ced0")

ggplot(data=Test4, aes(x=Sample, y=PercentFam, fill=Order)) + 
  facet_grid(Glacier~Is_kryal, scales = "free")+
  geom_bar(aes(), stat="identity", position="stack")+
  theme(legend.position="bottom", panel.background = element_rect(fill = 'white', color="grey60"),
        axis.title.x=element_blank(), axis.text.x=element_blank())+
  guides(fill=guide_legend(nrow=5))+scale_fill_manual(values =colors_15)


### 18S 
D_4_Gloom <- tax_glom(data_18S, taxrank = 'D_4')
D_4_Gloom # should list # taxa as # order
D_4_Gloom_melt<- psmelt(D_4_Gloom) 
D_4_Gloom_melt$D_4 <- as.character(D_4_Gloom_melt$D_4)
Test<-D_4_Gloom_melt%>%group_by(Sample, D_4)%>%summarise(Sum=sum(Abundance), .groups='drop') 

D_4_top9 = Test%>% group_by(D_4)%>% summarise(total=sum(Sum))%>%top_n(6, total)
Test2 = Test
Test2$Sum[Test2$Sum==0] <- NA
Test2<- Test2[complete.cases(Test2),]
Test2$D_4[!(Test2$D_4 %in% D_4_top9$D_4)] <- "Other"

Test3<- Test2 %>% group_by(Sample)%>% mutate(PercentFam=Sum/sum(Sum)*100)

Metadata_18S$Sample<-row.names(Metadata_18S)
Test4<-merge (Test3, Metadata_18S, by="Sample")
unique(Test4$D_4)

Test4$D_4 <- factor(Test4$D_4, levels = c("Chlorophyceae","Chrysophyceae", "Diatomea",
                                          "Phragmoplastophyta", "Trebouxiophyceae" , "Ulvophyceae","Other"))

library("RColorBrewer")
getPalette = colorRampPalette(brewer.pal(8, "Paired"))

ggplot(data=Test4, aes(x=Sample, y=PercentFam, fill=D_4)) + 
  facet_grid(Glacier~Stream, scales = "free")+
  geom_bar(aes(), stat="identity", position="stack")+
  theme(legend.position="bottom", panel.background = element_rect(fill = 'white', color="grey60"),
        axis.title.x=element_blank(),axis.text.x=element_text(angle = 90))+
  guides(fill=guide_legend(nrow=5))+
  scale_fill_manual(values=getPalette(7))







