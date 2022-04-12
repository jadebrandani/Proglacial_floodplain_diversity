library(igraph)
library(tidyverse)
library(ggplot2)
setwd("~/Documents/ENSEMBLE /0-CHAPTER1/3-Data_Analysis/Clean_Analysis/floodplain_diversity")

########################################## DATA PREPROCESING ##############################################

###################################### GLACIER - FED STREAM  ##############################################

### Data import
GFS_links <- read.csv("data/Network/GFS_Fungi_links.csv", row.names = 1) 
GFS_nodes <- read.csv("data/Network/GFS_Fungi_nodes.csv", row.names = 1) 

nrow(GFS_links); nrow(unique(GFS_links[,c("from", "to")])) #Checking if data frames have to be - simplified -
# if the numbers differ then it means that we have cases in the data where there ar multiple links between the same two nodes
GFS_links<-aggregate(GFS_links[,3], GFS_links[,-3], sum)
GFS_links<-GFS_links[order(GFS_links$from, GFS_links$to),]
colnames(GFS_links)[6] <- "weight"

GFS_net<-graph_from_data_frame(d=GFS_links, vertices=GFS_nodes, directed=F)  #Creating the graph object 

### Data Filtering
#STEP 1: Keep only the positiveweight   
GFS_net<- delete.edges(GFS_net,which(E(GFS_net)$weight<0))
GFS_net<-delete.vertices(GFS_net, V(GFS_net)[igraph::degree(GFS_net)==0]) #removal of nodes with 0 degrees (not connected to any other nodes)
length(E(GFS_net)) ; length(V(GFS_net))
#STEP 2: Filter on interaction type <<<<<<< Keep = YES includes Bacteria - Photographs AND Bacteria - Fungi  
GFS_net.phot.prok<- delete.edges(GFS_net,which(E(GFS_net)$Keep !="YES"))
GFS_net.phot.prok<-delete.vertices(GFS_net.phot.prok, V(GFS_net.phot.prok)[degree(GFS_net.phot.prok)==0])
length(E(GFS_net.phot.prok)); length(V(GFS_net.phot.prok))

#STEP 3: TOP 10 % of interaction strenght 
weight.list<-as.data.frame(E(GFS_net.phot.prok)$weight)
colnames(weight.list)<- "weight"
weight.list<-sort(weight.list$weight)
cut.off<-0.9
ref <- round(length(weight.list)*cut.off)
w = weight.list[ref]
GFS_net.phot.prok.sp <- delete_edges(GFS_net.phot.prok, E(GFS_net.phot.prok)[weight<w])
GFS_net.phot.prok.sp<-delete.vertices(GFS_net.phot.prok.sp, V(GFS_net.phot.prok.sp)[degree(GFS_net.phot.prok.sp)<1])
length(E(GFS_net.phot.prok.sp)) ; length(V(GFS_net.phot.prok.sp))

#STEP 4: Adding color for type 
V(GFS_net.phot.prok.sp)[V(GFS_net.phot.prok.sp)$Type2=="Bacteria"]$color <- "blue"
V(GFS_net.phot.prok.sp)[V(GFS_net.phot.prok.sp)$Type2=="Phototrophs"]$color <- "green"
V(GFS_net.phot.prok.sp)[V(GFS_net.phot.prok.sp)$Type2=="Fungi"]$color <- "red"

V(GFS_net.phot.prok.sp)[V(GFS_net.phot.prok.sp)$Type2=="Bacteria"]$label <- NA
V(GFS_net.phot.prok.sp)[V(GFS_net.phot.prok.sp)$Type2=="Phototrophs"]$label  <- V(GFS_net.phot.prok.sp)[V(GFS_net.phot.prok.sp)$Type2=="Phototrophs"]$Order
V(GFS_net.phot.prok.sp)[V(GFS_net.phot.prok.sp)$Type2=="Fungi"]$label  <- NA

### COMMUNITY DETECTION -----
GRDY_GFS<-cluster_fast_greedy(GFS_net.phot.prok.sp)
length(GRDY_GFS); sizes(GRDY_GFS); modularity(GRDY_GFS)
plot(GFS_net.phot.prok.sp, mark.groups=communities(GRDY_GFS), vertex.size=3, vertex.label=NA, main = "GRDY_GFS") #fast plot
# Cleaning step to remove clusters < 5 nodes
Small <- which(table(GRDY_GFS$membership) < 5)
Keep_GFS = V(GFS_net.phot.prok.sp)[!(GRDY_GFS$membership %in% Small)] ## Which nodes should be kept
GFS_net.phot.prok.sp2  = induced_subgraph(GFS_net.phot.prok.sp, Keep_GFS)## Get subgraph & plot

GRDY_GFS2 = cluster_louvain(GFS_net.phot.prok.sp2)
length(GRDY_GFS2); sizes(GRDY_GFS2); modularity(GRDY_GFS2)

plot(GFS_net.phot.prok.sp2, 
     mark.groups=communities(GRDY_GFS2),
     vertex.size=5, main = "GRDY_GFS_CLEANED_10%", vertex.label.dist=1, 
     layout=layout_with_fr, layout=layout_with_fr,niter = 500)

# Topology 
length(E(GFS_net.phot.prok.sp2)) ; length(V(GFS_net.phot.prok.sp2))
diameter(GFS_net.phot.prok.sp2) # Diameter
mean_distance(GFS_net.phot.prok.sp2)# Mean distance
edge_density(GFS_net.phot.prok.sp2, loops=FALSE)



############################################  TRIBUTARY Network   ###################################################

### Data import
TRIB_links <- read.csv("data/Network/Trib_Fungi_links.csv", row.names = 1) 
TRIB_nodes <- read.csv("data/Network/Trib_Fungi_nodes.csv", row.names = 1) 
nrow(TRIB_links); nrow(unique(TRIB_links[,c("from", "to")])) #Checking if data frames have to be - simplified -
# if the numbers differ then it means that we have cases in the data where there ar multiple links between the same two nodes
TRIB_links<-aggregate(TRIB_links[,3], TRIB_links[,-3], sum)
TRIB_links<-TRIB_links[order(TRIB_links$from, TRIB_links$to),]
colnames(TRIB_links)[8] <- "weight"
 
TRIB_net<-graph_from_data_frame(d=TRIB_links, vertices=TRIB_nodes, directed=F)  # Creating the graph object
length(E(TRIB_net)); length(V(TRIB_net))

### DATA FILTERING
#STEP 1: Keep only the positive
TRIB_net<- delete.edges(TRIB_net,which(E(TRIB_net)$weight<0))
TRIB_net<-delete.vertices(TRIB_net, V(TRIB_net)[igraph::degree(TRIB_net)==0]) #removal of nodes with 0 degrees (not connected to any other nodes)

#STEP 2: Filter on interaction type <<<<<<< need to create a few files with the different option  
TRIB_net.phot.prok<- delete.edges(TRIB_net,which(E(TRIB_net)$Keep !="YES"))
TRIB_net.phot.prok<-delete.vertices(TRIB_net.phot.prok, V(TRIB_net.phot.prok)[degree(TRIB_net.phot.prok)==0])
length(E(TRIB_net.phot.prok)); length(V(TRIB_net.phot.prok))

#STEP 3: TOP 10 % of interaction strength 
weight.list<-as.data.frame(E(TRIB_net.phot.prok)$weight)
colnames(weight.list)<- "weight"
weight.list<-sort(weight.list$weight)
cut.off<-0.9
ref <- round(length(weight.list)*cut.off)
w = weight.list[ref]
TRIB_net.phot.prok.sp <- delete_edges(TRIB_net.phot.prok, E(TRIB_net.phot.prok)[weight<w])
TRIB_net.phot.prok.sp<-delete.vertices(TRIB_net.phot.prok.sp, V(TRIB_net.phot.prok.sp)[degree(TRIB_net.phot.prok.sp)<1])
length(E(TRIB_net.phot.prok.sp)) ; length(V(TRIB_net.phot.prok.sp))

V(TRIB_net.phot.prok.sp)[V(TRIB_net.phot.prok.sp)$Type2=="Phototrophs"]$color <- "green"
V(TRIB_net.phot.prok.sp)[V(TRIB_net.phot.prok.sp)$Type2=="Bacteria"]$color <- "blue"
V(TRIB_net.phot.prok.sp)[V(TRIB_net.phot.prok.sp)$Type2=="Fungi"]$color <- "red"

V(TRIB_net.phot.prok.sp)[V(TRIB_net.phot.prok.sp)$Type2=="Bacteria"]$label <- V(TRIB_net.phot.prok.sp)[V(TRIB_net.phot.prok.sp)$Type2=="Bacteria"]$Family
V(TRIB_net.phot.prok.sp)[V(TRIB_net.phot.prok.sp)$Type2=="Phototrophs"]$label  <- V(TRIB_net.phot.prok.sp)[V(TRIB_net.phot.prok.sp)$Type2=="Phototrophs"]$Genus
V(TRIB_net.phot.prok.sp)[V(TRIB_net.phot.prok.sp)$Type2=="Fungi"]$label  <- V(TRIB_net.phot.prok.sp)[V(TRIB_net.phot.prok.sp)$Type2=="Fungi"]$Genus


### COMMUNITY DETECTION -----
GRDY_TRIB<-cluster_fast_greedy(TRIB_net.phot.prok.sp)
length(GRDY_TRIB); sizes(GRDY_TRIB); modularity(GRDY_TRIB)
plot(TRIB_net.phot.prok.sp, mark.groups=communities(GRDY_TRIB), vertex.size=3, vertex.label=NA, main = "GRDY_TRIB")
#Cleaning by removing clusters < 5 nodes
Small <- which(table(GRDY_TRIB$membership) < 5)
Keep_TRIB = V(TRIB_net.phot.prok.sp)[!(GRDY_TRIB$membership %in% Small)] ## Which nodes should be kept?
TRIB_net.phot.prok.sp2  = induced_subgraph(TRIB_net.phot.prok.sp, Keep_TRIB)## Get subgraph & plot

GRDY_TRIB2 = cluster_louvain(TRIB_net.phot.prok.sp2)
length(GRDY_TRIB2); sizes(GRDY_TRIB2); modularity(GRDY_TRIB2)

#### PLOT ###
plot(TRIB_net.phot.prok.sp2, 
     mark.groups=communities(GRDY_TRIB2),
     vertex.size=5, main = "Tributary Network", vertex.label.dist=0.5, 
     layout=layout_with_fr,niter = 200)

#Topology 
diameter(TRIB_net.phot.prok.sp2) # Diameter
mean_distance(TRIB_net.phot.prok.sp2)# Mean distance
edge_density(TRIB_net.phot.prok.sp2, loops=FALSE)


############################################### DATA ANALYSIS  ########################################################################
##### Importing already pre-processed file - GFS ---------
GFS_node_edit <- read.csv("data/Network/GFS_node_edit.csv", row.names = 1) 
GFS_links_edit<- read.csv("data/Network/GFS_links_edit.csv", row.names = 1) 
GFS_node_edit%>% dplyr::group_by(Type2) %>% count()

# Creating the graph object 
GFS_FINAL<-graph_from_data_frame(d=GFS_links_edit, vertices=GFS_node_edit, directed=F)

lay_GFS<-layout_with_fr(GFS_FINAL) # fixing the layout

plot(GFS_FINAL, 
     mark.groups=communities(GRDY_GFS2),
     vertex.size=8, main = "Bacteria:Family / Photo: Genus / Fungi : Genus", 
     layout=lay_GFS,niter = 400, vertex.label=V(GFS_FINAL)$label, 
     vertex.label.color= "white", vertex.label.font=2)

# TOPOLOGY
diameter(GFS_FINAL) # Diameter
mean_distance(GFS_FINAL)# Mean distance
edge_density(GFS_FINAL, loops=FALSE)
length(E(GFS_FINAL)) ; length(V(GFS_FINAL))
length(V(GFS_FINAL)[V(GGFS_FINAL)$Type2=="Bacteria"])
length(V(GFS_FINAL)[V(GFS_FINAL)$Type2=="Phototrophs"])
length(V(GFS_FINAL)[V(GFS_FINAL)$Type2=="Fungi"])


##### Importing already pre- processed file - TRIBUTARY --------
TRIB_node_edit <- read.csv("data/Network/Trib_node_edit.csv", row.names = 1) 
TRIB_links_edit<- read.csv("data/Network/Trib_links_edit.csv", row.names = 1) 
TRIB_node_edit%>% dplyr::group_by(Type2) %>% count()


# Creating the graph object 
TRIB_FINAL<-graph_from_data_frame(d=TRIB_links_edit, vertices=TRIB_node_edit, directed=F)

lay_trib<-layout_with_fr(TRIB_FINAL) # fixing the layout

plot(TRIB_FINAL, 
     mark.groups=communities(GRDY_TRIB2),
     vertex.size=7, main = "TRIBUTARY",
     layout=lay_trib,niter = 200, vertex.label=V(TRIB_FINAL)$label, 
     vertex.label.color= "white", vertex.label.font=2)

diameter(TRIB_FINAL) # Diameter
mean_distance(TRIB_FINAL)# Mean distance
edge_density(TRIB_FINAL, loops=FALSE)
length(E(TRIB_FINAL)) ; length(V(TRIB_FINAL))
length(V(TRIB_FINAL)[V(TRIB_FINAL)$Type2=="Bacteria"])
length(V(TRIB_FINAL)[V(TRIB_FINAL)$Type2=="Phototrophs"])
length(V(TRIB_FINAL)[V(TRIB_FINAL)$Type2=="Fungi"])


##### DEGREE BETWEENNESS ANALYSIS #######

## TRIBUTARY
degree.TRIB_df <- sort(igraph::degree(TRIB_FINAL)) %>% as.data.frame() %>% dplyr::rename("degree" = ".") %>% rownames_to_column(var = "Node")
betweenness.TRIB_df <- sort(betweenness(TRIB_FINAL)) %>% as.data.frame() %>% dplyr::rename("betweenness" = ".") %>% rownames_to_column(var = "Node")
hub_score.TRIB_df<-sort(hub_score(TRIB_FINAL)$vector) %>% as.data.frame() %>% dplyr::rename("hub_score" = ".") %>% rownames_to_column(var = "Node")
merge1.TRIB_df <- full_join(degree.TRIB_df, betweenness.TRIB_df, by = "Node")
merge.TRIB_df <- full_join(merge1.TRIB_df, hub_score.TRIB_df, by = "Node")
merge.TRIB_df$Type2 <- V(TRIB_FINAL)$Type2[match(degree.TRIB_df$Node, V(TRIB_FINAL)$name)]
merge.TRIB_df$Order <- V(TRIB_FINAL)$Order[match(degree.TRIB_df$Node, V(TRIB_FINAL)$name)]

#trying to count how many euk and prok per community 
Count<-merge.TRIB_df %>% group_by(Community) %>% count(Type)
Count_recast<-dcast(Count, Community~Type)
library(reshape2)

ggplot(data=merge.TRIB_df, (aes(x=betweenness, y=degree, color= Type2)))+
  geom_point(size=5)+
  geom_text(label=merge.TRIB_df$Order, nudge_x = 0.4, nudge_y = 0.4, 
            check_overlap = T)+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  scale_color_manual(values=c( "blue", "red","green3"))+
  ggtitle("TRIB - Degree /Betweenness")

## GLACIER - FED STREAM 
degree.GFS_df <- sort(igraph::degree(GFS_FINAL)) %>% as.data.frame() %>% dplyr::rename("degree" = ".") %>% rownames_to_column(var = "Node")
betweenness.GFS_df <- sort(betweenness(GFS_FINAL)) %>% as.data.frame() %>% dplyr::rename("betweenness" = ".") %>% rownames_to_column(var = "Node")
hub_score.GFS_df<-sort(hub_score(GFS_FINAL)$vector) %>% as.data.frame() %>% dplyr::rename("hub_score" = ".") %>% rownames_to_column(var = "Node")
merge1.GFS_df <- full_join(degree.GFS_df, betweenness.GFS_df, by = "Node")
merge.GFS_df <- full_join(merge1.GFS_df, hub_score.GFS_df, by = "Node")
merge.GFS_df$Type2 <- V(GFS_FINAL)$Type2[match(degree.GFS_df$Node, V(GFS_FINAL)$name)]
merge.GFS_df$Order <- V(GFS_FINAL)$Order[match(degree.GFS_df$Node, V(GFS_FINAL)$name)]


library(reshape2)
library(ggrepel)

ggplot(data=merge.GFS_df, (aes(x=betweenness, y=degree, color= Type2)))+
  geom_point(size=5)+
  geom_text(label=merge.GFS_df$Order, nudge_x = 0.4, nudge_y = 0.6, 
            check_overlap = T)+
  theme(legend.position="right", legend.box = "vertical",
        panel.background = element_rect(fill = 'white', color="grey60"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  scale_color_manual(values=c( "blue", "red","green3"))+
  ggtitle("GFS - Degree /Betweenness")

####### TAXONOMIC ANALYSIS ------ 
GFS_node_edit$Network<- "GFS"
TRIB_node_edit$Network<- "TRIB"

merge_node<-rbind(GFS_node_edit, TRIB_node_edit)  #merging all nodes together 
prok_colors<-c("#06b1b1","#e57010","#76be37","#644e7e","#e0dd00","#a53b6b","#c5811b","#6e7faf","#B5b8ba", "black", "grey")
# Subdividing into Type2 
merge_node_B<-merge_node[which(merge_node$Type2=="Bacteria"),]
merge_node_P<-merge_node[which(merge_node$Type2=="Phototrophs"),]
merge_node_F<-merge_node[which(merge_node$Type2=="Fungi"),]

##### PHOTOTROPHS ------ ------
merge_node_P_order<- merge_node_P%>%dplyr::group_by(Network, Order)%>%count()
merge_node_P_order$Network<-as.factor(merge_node_P_order$Network)
merge_node_P_order$Order<-as.factor(merge_node_P_order$Order)
merge_node_P_order$n<-as.numeric(merge_node_P_order$n)
colnames(merge_node_P_order)<-c("Network", "Order","Count")
merge_node_P_order<-merge_node_P_order%>%dplyr::group_by(Network) %>% dplyr::mutate(RA=(Count/sum(Count))*100) 
unique(merge_node_P_order$Order)
merge_node_P_order$Order<-factor(merge_node_P_order$Order, levels=c( "Chlorophyceae","Chrysophyceae","Diatomea" ,"Klebsormidiophyceae",
                                                                     "Phragmoplastophyta","Trebouxiophyceae","Ulvophyceae","Xanthophyceae", "Unknown"))
# plotting 
Phototrophs_colors<-c("#335f70","#4ED0C1","#359794","#BFD55D","#DFAC2A","#EDAD5A","#F0852D","#EF956B", "#E45C3A", "#E24D28")

ggplot(merge_node_P_order, aes(x=Network, y=Count, fill=Order))+
  geom_bar(stat="identity",position="stack")+
  theme(legend.position="right", panel.background = element_rect(fill = 'white', color="grey60"),
        axis.title.x=element_blank(),axis.text.x=element_text(size=14))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values =Phototrophs_colors)+
  ggtitle("Phototrophs - Order")+
  ylab("Count")

merge_node_P_order<-merge_node_P_order%>%dplyr::group_by(Network) %>% dplyr::mutate(RA=(Count/sum(Count))*100) 



##### FUNGI ------ ------
## Order
merge_node_F_Order<- merge_node_F%>%dplyr::group_by(Network, Order)%>%count()
merge_node_F_Order<-as_tibble(merge_node_F_Order, stringsAsFactors = FALSE)
merge_node_F_Order$Network<-as.factor(merge_node_F_Order$Network)
merge_node_F_Order$Order<-as.factor(merge_node_F_Order$Order)
merge_node_F_Order$n<-as.numeric(merge_node_F_Order$n)
colnames(merge_node_F_Order)<-c("Network", "Order","Count")

top_bacteria<-merge_node_F_Order%>%dplyr::group_by(Order)%>%dplyr::summarise(Total=sum(Count)) %>% top_n(3,Total) %>% droplevels()
levels(top_bacteria$Order) = c(levels(top_bacteria$Order),"Other")
levels(merge_node_F_Order$Order) = c(levels(merge_node_F_Order$Order),"Other")

merge_node_F_Order$Order[!(merge_node_F_Order$Order %in% top_bacteria$Order)]<-c("Other")

#merge_node_F_Order<-merge_node_F_Order%>%dplyr::group_by(Network) %>% dplyr::mutate(RA=(Count/sum(Count))*100) 

unique(merge_node_F_Order$Order)
merge_node_F_Order$Order<-factor(merge_node_F_Order$Order, levels=c( "Chytridiomycota","Cryptomycota","Zoopagomycota", "Other"))
# plotting 
Fungi_colors<-c("#2cacc9","#7dce82","#e8e288","#ff8360")

ggplot(merge_node_F_Order, aes(x=Network, y=Count, fill=Order))+
  geom_bar(stat="identity",position="stack")+
  theme(legend.position="right", panel.background = element_rect(fill = 'white', color="grey60"),
        axis.title.x=element_blank(),axis.text.x=element_text(size=14))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values =Fungi_colors)+
  ggtitle("Fungi - Order")+
  ylab("Count")


##### BACTERIA ------
## FAMILY 
merge_node_B_family<- merge_node_B%>%dplyr::group_by(Network, Family)%>%count()
merge_node_B_family<-as_tibble(merge_node_B_family, stringsAsFactors = FALSE)
merge_node_B_family$Network<-as.factor(merge_node_B_family$Network)
merge_node_B_family$Family<-as.factor(merge_node_B_family$Family)
merge_node_B_family$n<-as.numeric(merge_node_B_family$n)
colnames(merge_node_B_family)<-c("Network", "Family","Count")

top_bacteria<-merge_node_B_family%>%dplyr::group_by(Family)%>%dplyr::summarise(Total=sum(Count)) %>% top_n(9,Total) %>% droplevels()
#top_bacteria<-as.data.frame(top_bacteria)
str(top_bacteria)
levels(top_bacteria$Family) = c(levels(top_bacteria$Family),"Other")
levels(merge_node_B_family$Family) = c(levels(merge_node_B_family$Family),"Other")

merge_node_B_family$Family[!(merge_node_B_family$Family %in% top_bacteria$Family)]<-c("Other")


#merge_node_B_family<-merge_node_B_family%>%dplyr::group_by(Network) %>% dplyr::mutate(RA=(Count/sum(Count))*100) 

unique(merge_node_B_family$Family)
merge_node_B_family$Family<-factor(merge_node_B_family$Family, levels=c( "Chitinophagaceae","Comamonadaceae","Gemmatimonadaceae","Methylophilaceae","Nitrosomonadaceae",
                                                                         "Pirellulaceae","Pseudomonadaceae","Rhodobacteraceae","Sphingomonadaceae", "Other"))
# plotting 
Bacteria_colors2<-c("#669900","#99CC33","#CCEE66","#006699","#3399CC","#990066","#CC3399","#FF6600", "#FF9900", "#FFCC00")

ggplot(merge_node_B_family, aes(x=Network, y=Count, fill=Family))+
  geom_bar(stat="identity",position="stack")+
  theme(legend.position="right", panel.background = element_rect(fill = 'white', color="grey60"),
        axis.title.x=element_blank(),axis.text.x=element_text(size=14))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values =Bacteria_colors2)+
  ggtitle("Bacteria - Family")+
  ylab("Count")









