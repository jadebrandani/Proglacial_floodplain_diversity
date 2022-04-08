#Libraries ----
library(ggpubr)
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
library(pairwiseAdonis)

Stream <-c("#73DFFF","#FFA300")

### 16S Data Import -----
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
row.names.remove<- c("OTE_27", "OTE_48", "VAR_61", "VAR_62") #removing samples that were not sequenced and we suspect VAR_61 and VAR_62 to have been exchanged
Metadata <-Metadata [!(row.names(Metadata) %in% row.names.remove),]
rm(row.names.remove)
Metadata$ID<-row.names(Metadata)

OTU = phyloseq::otu_table(sDNA_16S_ASV, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(Taxonomy)
META<- phyloseq::sample_data(Metadata)
data_16S <- phyloseq(OTU, TAX, META)


#NMDS of 16S ALONE using bray-curtis distances 
O1<-ordinate(data_16S, "NMDS", "bray")
plot_ordination(data_16S,O1, type="samples",color= "Is_kryal",shape = "Glacier")+
  geom_point(size=5)+
  scale_color_manual(values = Stream)+
  scale_shape_manual(values=c(19, 17, 15))+
  xlab("NMDS1") + ylab("NMDS2") +ggtitle("NMDS on 16S - 257 samples")+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  guides(color = guide_legend(order=1, override.aes=list(shape=19, size=4)), shape = guide_legend(order=2, override.aes=list(size=4, fill="black")))


######## 18S data import -  This data has NOT been pre-filtered before hand -----
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

#filtering 18S 
#Fungi_18S <-subset_taxa (data_18S_RAW,D_3=="Fungi")
#Fungi_18S_tax_table<- as.data.frame(tax_table(Fungi_18S))
#Fungi_18S_tax_table$D_6 <-as.factor(Fungi_18S_tax_table$D_6)
#levels(Fungi_18S_tax_table$D_6)

#Subset by taxonomy (Photo-autotrophs exclusively)
data_18S_filtered = subset_taxa(data_18S_RAW,
                                D_2=="Chloroplastida"| 
                                  D_3=="Ochrophyta"| 
                                  D_2=="Cryptomonadales")

#Removing OTUs with only read count = 1 
Low_read_18S<- prune_taxa(taxa_sums(data_18S_filtered)>1, data_18S_filtered)
#Removing the OTU with only one observation
Low_read_18S<-phyloseq::filter_taxa (Low_read_18S, function(x){sum(x > 0) > 1}, prune = TRUE)
# We remove sample with less than 1 OTU
data_18S<-prune_samples(sample_sums(Low_read_18S)>1, Low_read_18S) # we are losing 4 samples (SOY_22 - VAR_21 - VAR_41 - VAR_54)


###### BETA diversity ######
sDNA_18S_OTUs_filtered = as(otu_table(data_18S), "matrix")

# Filtering low 18SOTU samples
sDNA_18S_OTUs_filtered2<-as.data.frame(sDNA_18S_OTUs_filtered[,colSums(sDNA_18S_OTUs_filtered>0)>1]) # OTU threshold. We only keep samples that have at least 1 OTU 
setdiff(colnames(sDNA_18S_OTUs_filtered),colnames(sDNA_18S_OTUs_filtered2)) 
Tax_18S_filtered<-Tax_18S_RAW[rownames(sDNA_18S_OTUs_filtered2),]

samples_removed<-setdiff(colnames(sDNA_18S_OTUs_filtered), colnames(sDNA_18S_OTUs_filtered2)) #we have removed 3 samples ( OTE_31, SOY_54, VAR_18)

#18S ordination
MDS_18S<-metaMDS(t(sDNA_18S_OTUs_filtered2))
Metadata_18S2<-Metadata_18S[colnames(sDNA_18S_OTUs_filtered2),]
plot(MDS_18S, disp="sites", type="n")
points(MDS_18S, display="sites", cex=0.8, pch=21, bg="red", select=which(Metadata_18S2$Stream=="Tributary"))
points(MDS_18S, display="sites", cex=0.8, pch=21, bg="blue", select=which(Metadata_18S2$Stream=="Glacier-fed"))

MDS_18S.scores<-scores(MDS_18S)[,1]
Metadata_16S2<-Metadata[colnames(sDNA_18S_OTUs_filtered2),]
# we filter 16S table to keep the same sample as in 18S dataset 
sDNA_16S_ASV2<-as.data.frame(sDNA_16S_ASV[,colnames(sDNA_18S_OTUs_filtered2)])
MDS_16S2<-metaMDS(t(sDNA_16S_ASV2))

ordisurf(MDS_16S2, MDS_18S.scores, knots=4, nlevels=8, col="grey", main="")
ordisurf(MDS_16S2, Metadata_16S2$Chla, knots=4, nlevels=8, col="grey", main="")
points(MDS_16S2, display="sites", cex=1.25, pch=21, bg="red", select=which(Metadata_16S2$Is_kryal=="Tributary"))
points(MDS_16S2, display="sites", cex=1.25, pch=21, bg="blue", select=which(Metadata_16S2$Is_kryal=="Glacier-fed"))


# implementing ordisurf in ggplot
ordisurf<-ordisurf(MDS_16S2, MDS_18S.scores, knots=4, nlevels=8, col="grey", main="")
head(ordisurf)

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}

contour.ordisurf <- extract.xyz(obj = ordisurf)
head(contour.ordisurf)
Metadata_16S2$ID<-row.names(Metadata_16S2)
Metadata_16S3<-data.frame(Metadata_16S2$Glacier, Metadata_16S2$ID, Metadata_16S2$Is_kryal)
colnames(Metadata_16S3)<-c("Glacier", "ID", "Stream")
sample.scores <- as.data.frame(scores(MDS_16S2))
sample.scores$ID<- row.names(sample.scores)
sample.scores.merge<-merge(sample.scores,Metadata_16S3, by="ID")
names(sample.scores.merge)[c(2, 3)] <- c("x", "y")
sample.scores.merge$z<-NA
sample.scores.merge$Glacier<-as.factor(sample.scores.merge$Glacier)
sample.scores.merge$Stream<-as.factor(sample.scores.merge$Stream)

str(sample.scores.merge)

samples_keep<-sample.scores.merge$ID

#FIGURE 3 : ORDI surf NMDS plots ********************
library(ggnewscale)
ggplot() +
  stat_contour(data = contour.ordisurf , aes(x, y, z = z),color = "dimgrey", bins= 9) + 
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  new_scale_color()+
  geom_point(data = sample.scores.merge, aes(x = x, y = y, fill=Stream, shape=Glacier, size= 4),color="black")+
  scale_fill_manual(values=Stream)+
  scale_shape_manual(values=c(21, 24, 22))+
  guides(fill = guide_legend(order=1, override.aes=list(shape=21, size=4)),size = "none", colour="none", shape = guide_legend(order=2, override.aes=list(size=4, fill="black")))+
  xlab("NMDS1") + ylab("NMDS2") +ggtitle("Ordisurf on 16S and 18S - 235 samples")


#NMDS of 18S alone *************************
O2<-ordinate(data_18S, "NMDS", "bray")

plot_ordination(data_18S,O2, type="samples",color= "Stream",shape = "Glacier")+
  geom_point(size=4)+
  scale_color_manual(values = Stream)+
  scale_shape_manual(values=c(19, 17, 15))+
  xlab("NMDS1") + ylab("NMDS2") +ggtitle("NMDS on 18S - 238 samples")+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  guides(color = guide_legend(order=1, override.aes=list(shape=19, size=4)), shape = guide_legend(order=2, override.aes=list(size=4, fill="black")))
#geom_text_repel(aes(NMDS1, NMDS2, label=ID))

###### PERMANOVA on 16S  ---------
bray_full <- phyloseq::distance(data_16S, method="bray")
env<-data.frame(sample_data(data_16S))
#Important test here 
adonis_16S<-adonis(bray_full~ Glacier+Is_kryal + Campaign + Chla, data = env)
adonis_16S$aov.tab$padj<-adjust_pvalue(adonis_16S$aov.tab$`Pr(>F)`,method = "BH")
adonis_16S
disper_full <- betadisper(bray_full, env$Campaign)
permutest(disper_full, pairwise=T)



 pairwise.adonis<-pairwise.adonis(bray_full, env$Glacier)
 pairwise.adonis
 
### PERMANOVA on 18S *************
bray_full_18S <- phyloseq::distance(data_18S, method="bray")
env_18S<-data.frame(sample_data(data_18S))
#Important test here 
adonis_18S<-adonis(bray_full_18S~ Glacier+Stream + Campaign +Chla, data = env_18S)
adonis_18S
capture.output(adonis_18S,file="output/adonis_18S.doc")
adonis_18S<-as.data.frame(adonis_18S$aov.tab)
write.csv(adonis_18S, "output/adonis_18S.csv")

disper_full_18S <- betadisper(bray_full_18S, env_18S$Glacier)
permutest(disper_full, pairwise=T)



#### DISEPRSION ANALYSIS 16S 
#Otemma
OTE <-subset_samples(data_16S, Glacier == "Otemma")
ote_bray<-phyloseq::distance(OTE, method="bray")
env_ote <- data.frame(sample_data(OTE))
adonis2(ote_bray~ Is_kryal, data = env_ote)
disper_ote <- betadisper(ote_bray, env_ote$Is_kryal)
permutest(disper_ote, pairwise=T)
boxplot(disper_ote, main ="Otemma", col=Stream)

#Valsorey
SOY <-subset_samples(data_16S, Glacier == "Valsorey")
soy_bray<-phyloseq::distance(SOY, method="bray")
env_soy <- data.frame(sample_data(SOY))
adonis2(soy_bray~ Is_kryal, data = env_soy)
disper_soy <- betadisper(soy_bray, env_soy$Is_kryal)
permutest(disper_soy, pairwise=T)

plot(disper_soy, main ="Valsorey", col=Stream)
boxplot(disper_soy, main ="Valsorey", col=Stream)

#Val Roseg 
VAR <-subset_samples(data_16S, Glacier == "Val Roseg")
VAR_prune <- prune_samples(!(sample_names(VAR) %in% to_remove), VAR)
var_bray<-phyloseq::distance(VAR_prune, method="bray")
env_var <- data.frame(sample_data(VAR_prune))
disper_var <- betadisper(var_bray, env_var$Is_kryal)
permutest(disper_var, pairwise=T)
plot(disper_var, main ="Val Roseg", col=Stream)
boxplot(disper_var, main ="Val Roseg", col=Stream)

adonis2(var_bray~ Is_kryal, data = env_var)

#Combining all to have one figure with boxplots
df<-env[,c(7,8,19,20)]
df$ID<- rownames(df)
disper_ote_melt<-reshape2::melt(as.matrix(disper_ote$distances), na.rm=TRUE)
colnames(disper_ote_melt)<-c("ID","NA","Distance")
df2<-merge(df,disper_ote_melt,by=c("ID"),all=TRUE)
disper_soy_melt<-reshape2::melt(as.matrix(disper_soy$distances),na.rm=TRUE)
colnames(disper_soy_melt)<-c("ID","NA","Distance")
df3<-merge(df2,disper_soy_melt,by=c("ID"),all=TRUE)
disper_var_melt<-reshape2::melt(as.matrix(disper_var$distances),na.rm=TRUE)
colnames(disper_var_melt)<-c("ID","NA","Distance")
df4<-merge(df3,disper_var_melt,by=c("ID"),all=TRUE)
library(tidyr)
df5<-unite(df4, centroid, c(Distance.x, Distance.y, Distance), remove=FALSE)
df5$centroid<-str_replace_all(df5$centroid, "NA_","")
df5$centroid<-str_replace_all(df5$centroid, "_NA","")

df_final<-df5[,c(1:5,7)]
df_final$centroid<-as.numeric(df_final$centroid)
df_final$centroid<-round(df_final$centroid,digits=3)
rm(df,df2,df3,df4,df5)

#DISPERSION 16S 
library(ggpubr)
p<-ggplot(df_final, aes(x=Is_kryal,y=centroid, fill=Is_kryal))+
  geom_boxplot(outlier.colour="black", outlier.shape=16)+
  facet_wrap(~Glacier,1)+
  #scale_color_manual(values=c( "coral1", "cyan3"), labels=c("Glacier-fed","Tributary"))+
  guides(fill=guide_legend(""))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="bottom", 
        legend.box = "horizontal",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks.x = element_blank())+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size =12, face="bold" ))+
  ggtitle("Dispersion Analysis - 16S")+ 
  labs(y="Distance to centroids")+ 
  scale_fill_manual(values=Stream)+
  stat_compare_means(method="t.test")
stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x=1.5, label.y=0.83)

p


### Comparing average BC value between glacier 
bray_16S<- phyloseq::distance(data_16S, method="bray")
library(reshape2)
bray_melt<-reshape2::melt(as.matrix(bray_16S))
bray_melt = bray_melt %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  distinct(value, .keep_all = TRUE)%>%
  mutate_if(is.factor, as.character)
colnames(bray_melt)<-c("Var1", "Var2", "BC")
bray_melt$Glacier_Var1<-as.factor(Metadata$Glacier[match(bray_melt$Var1, Metadata$ID)])
bray_melt$Glacier_Var2<-as.factor(Metadata$Glacier[match(bray_melt$Var2, Metadata$ID)])

bray_melt = bray_melt %>%
  filter(as.character(Glacier_Var1) != as.character(Glacier_Var2))
bray_melt = bray_melt %>% mutate(GL =
                     case_when(Glacier_Var1 == "Otemma" & Glacier_Var2 == "Val Roseg" ~ "OTE_VAR",
                     Glacier_Var1 == "Val Roseg" & Glacier_Var2 == "Otemma" ~ "OTE_VAR",
                     Glacier_Var1 == "Otemma" & Glacier_Var2 == "Valsorey" ~ "OTE_SOY",
                     Glacier_Var1 == "Valsorey" & Glacier_Var2 == "Otemma" ~ "OTE_SOY",
                     Glacier_Var1 == "Val Roseg" & Glacier_Var2 == "Valsorey" ~ "SOY_VAR",
                     Glacier_Var1 == "Valsorey" & Glacier_Var2 == "Val Roseg" ~ "SOY_VAR"))



bray_melt %>%
  group_by(GL) %>%
  get_summary_stats(BC, type = "mean_sd")

res.aov <- aov(bray_melt$BC~ bray_melt$GL, data = bray_melt)
summary(res.aov)

TukeyHSD(res.aov, "bray_melt$GL")


############################ LONGITUDINAL DIMENSION #####################################################-----
data_16S_GFS<-subset_samples(data_16S, Is_kryal=="Glacier-fed")
meta_GFS<-data.frame(phyloseq::sample_data(data_16S_GFS))

### Creating the dissimilarity matrix 
bray_GFS<- phyloseq::distance(data_16S_GFS, method="bray")
bray_GFS_melt<-reshape2::melt(as.matrix(bray_GFS))
bray_GFS_melt = bray_GFS_melt %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
## Creating the distance matrix 
dist_mat <- distm(GFS_coordinate, fun = distGeo) 
row.names(dist_mat)<-row.names(GFS_coordinate)
colnames(dist_mat)<-row.names(GFS_coordinate)

dist_mat_melt<-reshape2::melt(as.matrix(dist_mat))
dist_mat_melt = dist_mat_melt %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)

# Matching and re-arranging mbray-curtis matrix and distance matrix 
combine<-cbind(bray_GFS_melt, dist_mat_melt$value)
colnames(combine)<- c("Site1", "Site2", "BC", "Distance")

combine$GL_Site1 <- as.factor(meta_GFS$Glacier[match(combine$Site1, meta_GFS$ID)])
combine$GL_Site2 <- as.factor(meta_GFS$Glacier[match(combine$Site2, meta_GFS$ID)])
combine$GL_Season1 <- as.factor(meta_GFS$Campaign[match(combine$Site1, meta_GFS$ID)])
combine$GL_Season2 <- as.factor(meta_GFS$Campaign[match(combine$Site2, meta_GFS$ID)])

combine_filter = combine %>%
  filter(as.character(GL_Site1) == as.character(GL_Site2)) %>%
  mutate_if(is.factor, as.character)
combine_filter2 = combine_filter %>%
  filter(as.character(GL_Season1) == as.character(GL_Season2)) %>%
  mutate_if(is.factor, as.character)
combine_filter2$Similarity<-1-combine_filter2$BC #Transformation in similarity

ggplot(combine_filter2,aes(x=Distance,y=Similarity, shape= GL_Site1))+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
  geom_point(color="#73DFFF", size = 3.5)+
  scale_shape_manual(values=c(19, 17, 15))+
  ylim(0,1)+
  facet_grid(GL_Site1~GL_Season1)+
  ggtitle("Distance_decay")+
  geom_smooth(linetype ="dashed", color="black", method = "lm", se=F, formula = y ~ x)+
  stat_regline_equation(label.x = c(10, 700), 
                        label.y = c(0.95),
                        aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")))+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3, label.y=0.75
  )+
  xlab("Pairwise distance (m) ")