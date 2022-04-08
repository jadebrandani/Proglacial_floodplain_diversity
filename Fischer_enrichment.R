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
library(btools)
library(geosphere)
library(ade4)
library(lattice)
library(permute)
library(tabula)

#data import
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
data_16S <- phyloseq(OTU, TAX, META)
#prevalence filtering ASV must be present at least 5 times to stay in the dataset 
data_16S_prevalence_filtering <- phyloseq::filter_taxa(data_16S, function(x){sum(x > 0) >= 5}, prune = TRUE)

Pres_16S_test<-as.data.frame(otu_table(data_16S_prevalence_filtering))

# Wee need to transform the matrix in a presence absence 
Pres_16S_test[Pres_16S_test > 0] <- 1   
Pres_16S_test<-as.data.frame(t(Pres_16S_test))
Pres_16S_test$Stream<-Metadata$Is_kryal[match(rownames(Pres_16S_test),Metadata$ID)]

Subset<- as.data.frame(Pres_16S_test[,c(5,20)])

Subset$Presence<- (Subset[,1] == 1)
Subset$Stream_boolean <- (Subset$Stream == "Glacier-fed")

Subset_table<-as.data.frame(table(Subset$Presence, Subset$Stream_boolean))
colnames(Subset_table)<-c("Presence", "GFS", "freq")

df <- xtabs(freq~Presence+GFS, data = Subset_table)
print(df)
fisher.result <- fisher.test(df)
print(c(fisher.result$estimate,fisher.result$p.value))

# Fischer's exact test net 
Fischer_result <- data.frame()
# loop through the scales and each variable
for(i in 1:ncol(Pres_16S_test[,c(1:24323)])) {
  
  #Contigency table 
  Subset<- as.data.frame(Pres_16S_test[,c(i,24324)])
  Subset$Presence<-(Subset[,1] == 1)
  Subset$Stream_boolean <-(Subset$Stream == "Glacier-fed")
  
  Subset_table<-as.data.frame(table(Subset$Presence, Subset$Stream_boolean))
  colnames(Subset_table)<-c("Presence", "GFS", "freq")
  xtabs <- xtabs(freq~Presence+GFS, data = Subset_table)
  fisher.result <- fisher.test(xtabs)
  
  ## capture summary stats
  OR <- fisher.result$estimate
  p.value <- fisher.result$p.value
  
  # create temporary data frame
  df <-data.frame(asv = i, OddRatio = OR[1],p.value = p.value [1],stringsAsFactors = F)
  
  # bind rows of temporary data frame to the results data frame
  Fischer_result<- rbind(Fischer_result, df) 
  }

##
ASV_names<-as.data.frame(colnames(Pres_16S_test))
colnames(ASV_names)<-"ID"
Fischer_name<-Fischer_result
Fischer_name$ASV<-ASV_names$ID[match(Fischer_name$asv,rownames(ASV_names))]

Fischer_name$p.adj <-p.adjust(Fischer_name$p.value, method="fdr", n = length(Fischer_name$p.value))

hist(Fischer_name$p.value)

ASV.to.keep <- Fischer_name$ASV

Taxonomy_prevalence_filtered<-as.data.frame(tax_table(data_16S_prevalence_filtering))
Taxonomy.subset <- Taxonomy_prevalence_filtered[(row.names(Taxonomy_prevalence_filtered) %in% ASV.to.keep),]

Fischer.taxonomy<-merge(Fischer_name, Taxonomy.subset, by.x = "ASV", by.y = 0)
Fischer.taxonomy.signif<- Fischer.taxonomy[Fischer.taxonomy$p.adj<= 0.05,]
Fischer.taxonomy.signif$Group<- (Fischer.taxonomy.signif$OddRatio > 1)
Fischer.taxonomy.signif<-Fischer_name
Fischer.taxonomy.signif$Group<- (Fischer_name$OddRatio > 1)

rm(Fischer.taxonmy)

# Fischer's exact test net 
pos_enrichment = data.frame()
for (Genus in unique(Fischer.taxonomy.signif$Genus)){
  if (Genus != ''){
    # Pos enrichment
    pos_table = table(Fischer.taxonomy.signif$Genus == Genus, Fischer.taxonomy.signif$Group == 'TRUE')
    ftest = fisher.test(pos_table, alternative = 'greater')
    p = ftest$p.value
    low_or = ftest$conf.int[1]
    high_or = ftest$conf.int[2]
    or = ftest$estimate
    pos_enrichment = rbind(pos_enrichment, data.frame(Genus=Genus,p=p,low_or=low_or,high_or=high_or,or=or))
  }
}


pos_enrichment$p.adj <-p.adjust(pos_enrichment$p, method="fdr", n = length(pos_enrichment$p))

pos_enrichement.signif<-pos_enrichment[pos_enrichment$p.adj<=0.05,]

# plot 
ggplot()+
  geom_pointrange(data=pos_enrichement.signif, mapping=aes(x=or, y=Genus, xmin=low_or, xmax= high_or))
###################

# We want to know the relative abundance of the identified taxa 
#create a vector with the name of isolated ASV 
Genus_25 <- pos_enrichement.signif$Genus
data_16S_pf_RelativeAbundance <- transform_sample_counts(data_16S_prevalence_filtering, function(x) x/sum(x))

genus_RA<-tax_glom(data_16S_pf_RelativeAbundance, taxrank="Genus", NArm=TRUE)

genus_RA_subset<- subset_taxa(genus_RA, Genus %in% Genus_25)

genus_RA_subset@sam_data
genus_RA_subset.df<-psmelt(genus_RA_subset)

#Overall mean - NOT SPECIFIC TO GFS ( but possible to do when using the sumamrize function )
Genus_summary<-genus_RA_subset.df %>%
  group_by(Genus)%>%
  summarize(mean_abund = mean(Abundance, na.rm=TRUE)) 
head(Genus_summary)

Genus_summary$mean_abund_precent<-Genus_summary$mean_abund *100
ggplot(data=Genus_summary, aes(x=mean_abund_precent, y=Genus))+
  geom_point()

row.names(pos_enrichement.signif)<-pos_enrichement.signif$Genus
row.names(Genus_summary)<-Genus_summary$Genus

Merge_enrichment<-merge(pos_enrichement.signif, Genus_summary, by=0)
Merge_enrichment<-Merge_enrichment[,c(1,3:7,9,10)]
colnames(Merge_enrichment)<-c("Genus", "p", "low_or", "high_or", "OR", "p_adj", "mean_abun", "mean_abund_percent")
Merge_enrichment$Genus<-as.factor(Merge_enrichment$Genus)
Merge_enrichment$abund_category<- if(Merge_enrichment$mean_abund_percent< 1){
  (" < 1%")
}else if (Merge_enrichment$mean_abund_percent< 2.5){
  ("> 2%")
}else ("> 5%")


ggplot(data=Merge_enrichment, aes(x=OR, y=Genus))+
  #geom_pointrange(mapping=aes(x=OR, y=Genus, xmin=low_or, xmax= high_or,size=mean_abund_percent ))
  geom_point(aes(size=mean_abund_percent))+
  lims(c("<1%",))
xlim(0, 40)+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        panel.grid.major = element_line(colour = "grey50"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))


ggplot()+
  geom_pointrange(data=Merge_enrichment, mapping=aes(x=OR, y=Genus, xmin=low_or, xmax= high_or))





