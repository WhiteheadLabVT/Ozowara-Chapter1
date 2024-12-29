#Chapter 1 Analysis 
rm(list=ls()) # clears work space
#Install packages---------------------------------------------------------------
library(ggplot2)
library(MASS)
library(reshape2)
library(dplyr)
library(randomForest)
library(glmmTMB) #mixed models
library(car) #mixed models summary tables
library(vegan) #multivariate stats
library(Boruta) #random forest models
library(tidyverse)
library(pls) #PCAs
library(Hmisc)#correlation matrix calculation 
library(corrplot)#plotting correlation matrix 
library(readr)
library(GGally)
library(emmeans) 
library(ggpubr) #multiplotting 
library(viridis) #color pallete 
library(performance) # generating r squareds 

#Read Data Organization and Restructuring---------------------------------------
#Read in Data sheets
#Orchard Level Data (24 obs)
Orchard <- read_csv("Ozowara_et_al_OrchardLevelData.csv")
Orchard$orchard.num <- as.factor(as.character(Orchard$orchard.num))

#Tree Level Data (120 obs)
Tree <- read_csv("Ozowara_et_al_TreeLevelData.csv")
View(Tree_Level_Data)
Tree$orchard.num <- as.factor(as.character(Tree$orchard.num))
Tree$Tree <- as.factor(as.character(Tree$Tree))


#Phenolics Data (359 obs)
d <-read_csv("Ozowara_et_al_FruitLevelData.csv")
d$orchard.num <- as.factor(as.character(d$orchard.num))

#Organizing Phenolics Data for Q1#

#creating totalphen and phenrich
d <- d %>%
  mutate(TotalPhen=rowSums(across(8:41)),
         PhenRich=rowSums(across(8:41)!=0))

#tranforms total phenolics for model analysis
d <- d %>%
  mutate(TotalPhenTrans=((TotalPhen/1000000)+0.0001))

#organizational data as factor 
d <- d %>% 
  mutate_at(c("orchard.num", "orchard.type", "site.code", "Tree", "SampleID"), as.factor)

#separate by tissue type 
d.sk <- filter(d, Tissue=="SKIN")
d.sk <- dplyr::select(d.sk, -(6+which(colSums(d.sk[7:41], na.rm=TRUE) %in% 0)))
d.pu <- filter(d, Tissue=="PULP")
d.pu <- dplyr::select(d.pu, -(6+which(colSums(d.pu[7:41], na.rm=TRUE) %in% 0)))
d.se <- filter(d, Tissue=="SEED")
d.se <- dplyr::select(d.se, -(6+which(colSums(d.se[7:41], na.rm=TRUE) %in% 0)))


#Assign row names
d <- as.data.frame(d)
row.names(d) <- d$SampleID

#For some analyses we need a table with only composition info
d.comp <- d[,8:41]
#and one with just explanatory variables
d.expl <- d[,1:6]

#Also need those by tissue
#skin 
d.comp.sk <- d.sk[,8:41]
d.expl.sk <- d.sk[,1:6]
#pulp
d.comp.pu <- d.pu[,8:37]
d.expl.pu <- d.pu[,1:6]
#Seed
d.comp.se <- d.se[,8:38]
d.expl.se <- d.se[,1:6]

#Condensing to orchard level data for all other questions#
Tree <- Tree %>%
  mutate(avgwgt=(apple.wgt.3-bag.weight)/3) 

Tree <- Tree %>% 
  dplyr::select(-c(5:6))

Tree.sum <- Tree %>%
  group_by(orchard.num) %>%
  summarise_at(c("avgwgt", "Firmness", "SSC", "maturity.index"), mean, na.rm = TRUE)

#categorical latitude data to be used for RF and pest questiosn 
Orchard <- Orchard %>% 
  mutate(lat_cat=cut(Latitude, breaks=c(-Inf, 42, Inf), labels=c("low", "high")))

c <- left_join(Orchard, Tree.sum, by="orchard.num")

#dividing chem data by tissue type to the orchard level  
SkinD <- filter(d, Tissue=="SKIN")
PulpD <- filter(d, Tissue=="PULP")
SeedD <- filter(d, Tissue=="SEED")

#Skin to orchard 
SkinD <- SkinD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich","TotalPhenTrans")
               , mean, na.rm = TRUE)
SkinD <- left_join(SkinD, Tree.sum, by="orchard.num")
SkinD <- left_join(SkinD, Orchard, by="orchard.num")

#Pulp to orchard 
PulpD <- PulpD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich","TotalPhenTrans")
               , mean, na.rm = TRUE)
PulpD <- left_join(PulpD, Tree.sum, by="orchard.num")
PulpD <- left_join(PulpD, Orchard, by="orchard.num")
shapiro.test(PulpD$TotalPhen)

#Seed to orchard 
SeedD <- SeedD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich", "TotalPhenTrans")
               , mean, na.rm = TRUE)
SeedD <- left_join(SeedD, Tree.sum, by="orchard.num")
SeedD <- left_join(SeedD, Orchard, by="orchard.num")


#Q1: How do management systems interact with broad climatic changes across latitude?----
#Analysis: Using GLMMs with tree level data to explore how measured qualities interact
#across latitude 
#Q1-A: Physical Quality------------------------------------------------------------
#binding latitude to the 120 obs data set#
TreeLat <- left_join(Tree, Orchard %>% select(orchard.num, Latitude, lat_cat), 
                     by = c("orchard.num"))

#SSC
Lat1<- glmmTMB(SSC ~ orchard.type*Latitude + maturity.index + (1|site.code/orchard.num), data=TreeLat)
summary(Lat1)
Anova(Lat1) 

plot(SSC ~ maturity.index, data=TreeLat)

#calculate effect size
emmeans(Lat1, ~ Latitude+maturity.index, at = list(Latitude = c(34:48)), type= "reponse")

#Firmness 
Lat2<- glmmTMB(Firmness ~ orchard.type*Latitude + maturity.index + (1|site.code/orchard.num), data=TreeLat)
summary(Lat2)
Anova(Lat2)

#calculate effect size
emmeans(Lat2, ~ Latitude, at = list(Latitude = c(34:48)), type= "reponse")

#Average Weight 
Lat3<- glmmTMB(avgwgt ~ orchard.type*Latitude + maturity.index + (1|site.code/orchard.num), data=TreeLat)
summary(Lat3)
Anova(Lat3) 
#orchard.type:Latitude 9.6577  1   0.001886 **

#calculate effect size
emmeans(Lat3, ~ Latitude|orchard.type, at = list(Latitude = c(34:48)), type= "reponse")


###Investigating Latitude and Average weight 
#Splitting latitude in high and low 
Tree_low <- filter(TreeLat, Latitude<42)
Tree_high <- filter(TreeLat, Latitude>42)

#Low Latitude Sites 
avg_low<- glmmTMB(avgwgt ~ orchard.type + maturity.index +(1|site.code/orchard.num), data=Tree_low)
summary(avg_low)
Anova(avg_low) 

#High Latitude Sites 
avg_hgh<- glmmTMB(avgwgt ~ orchard.type + maturity.index + (1|site.code/orchard.num), data=Tree_high)
summary(avg_hgh)
Anova(avg_hgh)

#Q1-B: Fruit Chemistry-------------------------------------------------------------
#binding latitude to the 359 obs data set 
ChemLat <- left_join(d, Orchard %>% select(orchard.num, Latitude), 
                     by = c("orchard.num"))
ChemLat <- left_join(ChemLat, TreeLat %>% select(orchard.num, Tree, maturity.index), 
                     by = c("orchard.num", "Tree"))

#Total Phenolics with beta distribution 
tp1 <- glmmTMB(TotalPhenTrans ~ orchard.type*Tissue*Latitude + maturity.index +
                 (1|site.code/orchard.num), data=ChemLat, family=beta_family(link="logit"))
summary(tp1)
Anova(tp1)

#binding latitudinal data to tissue data sets
d.sk <- left_join(d.sk, Orchard[,c(2,4)], by="orchard.num")
d.sk <- left_join(d.sk, TreeLat %>% select(orchard.num, Tree, maturity.index), 
                  by = c("orchard.num", "Tree"))

d.pu <- left_join(d.pu, Orchard[,c(2,4)], by="orchard.num")
d.pu <- left_join(d.pu, TreeLat %>% select(orchard.num, Tree, maturity.index), 
                  by = c("orchard.num", "Tree"))


d.se <- left_join(d.se, Orchard[,c(2,4)], by="orchard.num")
d.se <- left_join(d.se, TreeLat %>% select(orchard.num, Tree, maturity.index), 
                  by = c("orchard.num", "Tree"))

###total phenolics###
#skin
tp.sk <- glmmTMB(TotalPhenTrans ~ orchard.type*Latitude + maturity.index + (1|site.code/orchard.num), 
                 data=d.sk, family=beta_family(link="logit"))
summary(tp.sk)
Anova(tp.sk)
#nothing 

#pulp
tp.pu <- glmmTMB(TotalPhenTrans ~ orchard.type*Latitude + maturity.index + (1|site.code/orchard.num), 
                 data=d.pu, family=beta_family(link="logit"))
summary(tp.pu)
Anova(tp.pu)
#nothing 

#seed
tp.se <- glmmTMB(TotalPhenTrans ~ orchard.type*Latitude + maturity.index +(1|site.code/orchard.num/Tree), 
                 data=d.se, family=beta_family(link="logit"))
summary(tp.se)
Anova(tp.se)

#calculate effect size
emmeans(tp.se, ~maturity.index, type="response")

###phenolics richness###
pr1 <- glmmTMB(PhenRich~ orchard.type*Latitude*Tissue + maturity.index + (1|site.code/orchard.num), data=ChemLat, 
               family=poisson(link="log"))
summary(pr1)
Anova(pr1)


#skin
pr.sk <- glmmTMB(PhenRich ~ orchard.type*Latitude+maturity.index+(1|site.code/orchard.num), data=d.sk, 
                 family=poisson(link="log"))
summary(pr.sk)
Anova(pr.sk)


#pulp
pr.p <- glmmTMB(PhenRich ~ orchard.type*Latitude+maturity.index+(1|site.code/orchard.num), data=d.pu, 
                family=poisson(link="log"))
summary(pr.p)
Anova(pr.p)


#seed
pr.se <- glmmTMB(PhenRich ~ orchard.type*Latitude+maturity.index+(1|site.code/orchard.num), data=d.se, 
                 family=poisson(link="log"))
summary(pr.se)
Anova(pr.se)

plot(PhenRich ~ maturity.index, data=d.se)

#calculate effect size
emmeans(pr.se,~maturity.index, type="response")

#Q1-C: NMDS + PERMANOVA---------------------------------------------------------
#left joining the "expl" by latitude  

d.expl <- left_join(d.expl, TreeLat %>% select(orchard.num, Tree, maturity.index,lat_cat), 
          by = c("orchard.num", "Tree"))

d.expl.sk <- left_join(d.expl.sk, TreeLat %>% select(orchard.num, Tree, maturity.index, lat_cat), 
                    by = c("orchard.num", "Tree"))

d.expl.pu <- left_join(d.expl.pu, TreeLat %>% select(orchard.num, Tree, maturity.index,lat_cat), 
                       by = c("orchard.num", "Tree"))

d.expl.se <- left_join(d.expl.se, TreeLat %>% select(orchard.num, Tree, maturity.index, lat_cat), 
                       by = c("orchard.num", "Tree"))

###NMDS for Skin 
m.NMDS.sk <- metaMDS(d.comp.sk, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.sk

plot(m.NMDS.sk, type="t")

#Extract the site scores (coordinates of points in ordination space)
nmds_scores <- vegan::scores(m.NMDS.sk, display = "sites")

# Convert to data frame
data.scores <- as.data.frame(nmds_scores)
data.scores$SampleID <- rownames(data.scores)

# Merge NMDS scores with metadata
merged_data <- merge(data.scores, d.expl.sk, by = "row.names")

names(merged_data)
# Plot with ggplot2
sk.nmds <- ggplot(merged_data, aes(x = NMDS1, y = NMDS2, color=orchard.type, shape= lat_cat)) +
  geom_point(size = 2) +
  geom_jitter(size = 2, width = 0.5, height = 0.5) +
  theme_classic() +
  scale_color_manual(values = c("#3b0f70", "#de4968"), name = "Management System")+
  scale_shape_manual(values = c(16, 17), name = "Latitudinal Category") +  
  stat_ellipse(mapping = aes(group = orchard.type), level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "Skin",
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme(legend.position = "bottom")

#For some reason ggplot does not want to plot the shape and color at the same time
#this code supresses the warning and prints anyways 
suppressWarnings(print(sk.nmds))

sk.nmds  <- sk.nmds  + guides(color = "none")


#PERMANOVA can test whether the visualized differences are significant
m.perm1 <- adonis2(d.comp.sk~orchard.type*lat_cat+maturity.index, data=d.expl.sk)
m.perm1


###NMDS for Pulp 
m.NMDS.pu <- metaMDS(d.comp.pu, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.pu


#Extract the site scores (coordinates of points in ordination space)
nmds_scores <- vegan::scores(m.NMDS.pu, display = "sites")

# Convert to data frame
data.scores <- as.data.frame(nmds_scores)
data.scores$SampleID <- rownames(data.scores)

# Merge NMDS scores with metadata
merged_data <- merge(data.scores, d.expl.pu, by = "row.names")

# Plot with ggplot2
pu.nmds <- ggplot(merged_data, aes(x = NMDS1, y = NMDS2, color=orchard.type, shape= lat_cat)) +
  geom_point(size = 2) +
  geom_jitter(size = 2, width = 0.5, height = 0.5) +
  theme_classic() +
  scale_color_manual(values = c("#3b0f70", "#de4968"), name = "Management System")+
  scale_shape_manual(values = c(16, 17), name = "Latitudinal Category") +  
  stat_ellipse(mapping = aes(group = orchard.type), level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "Pulp",
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme(legend.position = "bottom")


suppressWarnings(print(pu.nmds))

#PERMANOVA
m.perm2 <- adonis2(d.comp.pu~orchard.type*lat_cat+maturity.index, data=d.expl.pu)
m.perm2


###NMDS for Seed 
m.NMDS.se <- metaMDS(d.comp.se, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.se

nmds_scores <- vegan::scores(m.NMDS.se, display = "sites")

# Convert to data frame
data.scores <- as.data.frame(nmds_scores)
data.scores$SampleID <- rownames(data.scores)

# Merge NMDS scores with metadata
merged_data <- merge(data.scores, d.expl.se, by = "row.names")

# Plot with ggplot2
se.nmds <- ggplot(merged_data, aes(x = NMDS1, y = NMDS2, color=orchard.type, shape= lat_cat)) +
  geom_point(size = 2) +
  geom_jitter(size = 2, width = 0.5, height = 0.5) +
  theme_classic() +
  scale_color_manual(values = c("#3b0f70", "#de4968"), name = "Management System")+
  scale_shape_manual(values = c(16, 17), name = "Latitudinal Category") +  
  stat_ellipse(mapping = aes(group = orchard.type), level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "Seed",
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme(legend.position = "bottom")


suppressWarnings(print(se.nmds))

se.nmds  <- se.nmds  + guides(color = "none")


#PERMANOVA
m.perm3 <- adonis2(d.comp.se~orchard.type*lat_cat+maturity.index, data=d.expl.se)
m.perm3



##joining all for multiplot 
ggarrange(sk.nmds, pu.nmds, se.nmds, nrow = 1, ncol = 3, labels = c("A", "B", "C"))
ggsave("fig5.png", width=20, height=24, units="cm", dpi=600)






#Q1-C: Random Forest------------------------------------------------------------
#Management Sorting 
#Skin 
m1.rf.sk <- randomForest(d.comp.sk,d.expl.sk$orchard.type, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.sk$importance
varImpPlot(m1.rf.sk)
MDSplot(m1.rf.sk, d.expl.sk$orchard.type)

m1.rf.b.sk <- Boruta(d.comp.sk,d.expl.sk$orchard.type)
m1.rf.b.sk
plot(m1.rf.b.sk,las = 2, cex.axis = 0.7)  

getSelectedAttributes(m1.rf.b.sk) 
attStats(m1.rf.b.sk)

##Running MANOVAS (needs to be fixed)
d.comp.sk.sel <- data.matrix(d.comp.sk[,getSelectedAttributes(m1.rf.b.sk)])
m1.man.sk <-manova(d.comp.sk.sel ~ d.expl.sk$orchard.type)
summary(m1.man.sk)  

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.sk)  

#some quick plots of all of them
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.sk.sel))){
  d.temp=d.comp.sk.sel[,i]
  plot(d.temp ~ d.expl.sk$orchard.type, ylab=colnames(d.comp.sk.sel)[i])
}
dev.off()


###PULP###
m1.rf.pu <- randomForest(d.comp.pu,d.expl.pu$orchard.type, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.pu$importance
varImpPlot(m1.rf.pu)
MDSplot(m1.rf.pu, d.expl.pu$orchard.type)

m1.rf.b.pu <- Boruta(d.comp.pu,d.expl.pu$orchard.type)
m1.rf.b.pu
plot(m1.rf.b.pu,las = 2, cex.axis = 0.7)  

getSelectedAttributes(m1.rf.b.pu) 
attStats(m1.rf.b.pu)

##Running MANOVAS
d.comp.pu.sel <- data.matrix(d.comp.pu[,getSelectedAttributes(m1.rf.b.pu)])
m1.man.pu <-manova(d.comp.pu.sel ~ d.expl.pu$orchard.type)
summary(m1.man.pu)  

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.pu)  

#some quick plots of all of them
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.pu.sel))){
  d.temp=d.comp.pu.sel[,i]
  plot(d.temp ~ d.expl.pu$orchard.type, ylab=colnames(d.comp.pu.sel)[i])
}
dev.off()

###SEEDS###
m1.rf.se <- randomForest(d.comp.se,d.expl.se$orchard.type, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.se$importance
varImpPlot(m1.rf.se)
MDSplot(m1.rf.se, d.expl.se$orchard.type)


m1.rf.b.se <- Boruta(d.comp.se,d.expl.se$orchard.type)
m1.rf.b.se
plot(m1.rf.b.se,las = 2, cex.axis = 0.7)  

getSelectedAttributes(m1.rf.b.se)
attStats(m1.rf.b.se)

#Running MANOVAS 
d.comp.se.sel <- data.matrix(d.comp.se[,getSelectedAttributes(m1.rf.b.se)])
m1.man.se <- manova(d.comp.se.sel ~ d.expl.se$orchard.type)
summary(m1.man.se)  

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.se)  

par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.se.sel))){
  d.temp=d.comp.se.sel[,i]
  plot(d.temp ~ d.expl.se$orchard.type, ylab=colnames(d.comp.se.sel)[i])
}
dev.off()

###Random Forest for Latitude ###
###skin
m1.rf.sk <- randomForest(d.comp.sk,d.expl.sk$lat_cat, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.sk$importance
varImpPlot(m1.rf.sk)
MDSplot(m1.rf.sk, d.expl.sk$lat_cat)

#using boruta 
m1.rf.b.sk <- Boruta(d.comp.sk,d.expl.sk$lat_cat)
m1.rf.b.sk

#plot 
plot(m1.rf.b.sk,las = 2, cex.axis = 0.7)  
getSelectedAttributes(m1.rf.b.sk) 
attStats(m1.rf.b.sk)


#performing MANOVA 
d.comp.sk.sel <- data.matrix(d.comp.sk[,getSelectedAttributes(m1.rf.b.sk)])
m1.man.sk <- manova(d.comp.sk.sel ~ d.expl.sk$lat_cat)
summary(m1.man.sk) 

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.sk)  

#some quick plots of all of them
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.sk.sel))){
  d.temp=d.comp.sk.sel[,i]
  plot(d.temp ~ d.expl.sk$lat_cat, ylab=colnames(d.comp.sk.sel)[i])
}
dev.off()

###PULP###
m1.rf.pu <- randomForest(d.comp.pu,d.expl.pu$lat_cat, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.pu$importance
varImpPlot(m1.rf.pu)
MDSplot(m1.rf.pu, d.expl.pu$lat_cat)

m1.rf.b.pu <- Boruta(d.comp.pu,d.expl.pu$lat_cat)
m1.rf.b.pu
plot(m1.rf.b.pu, las = 2, cex.axis = 0.7)  

getSelectedAttributes(m1.rf.b.pu) 
attStats(m1.rf.b.pu)

##Running MANOVAS
d.comp.pu.sel <- data.matrix(d.comp.pu[,getSelectedAttributes(m1.rf.b.pu)])
m1.man.pu <-manova(d.comp.pu.sel ~ d.expl.pu$lat_cat)
summary(m1.man.pu)

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.pu)  

#some quick plots of all of them
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.pu.sel))){
  d.temp=d.comp.pu.sel[,i]
  plot(d.temp ~ d.expl.pu$lat_cat, ylab=colnames(d.comp.pu.sel)[i])
}
dev.off()

###SEEDS###
m1.rf.se <- randomForest(d.comp.se,d.expl.se$lat_cat, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.se$importance
varImpPlot(m1.rf.se)
MDSplot(m1.rf.se, d.expl.se$lat_cat)


m1.rf.b.se <- Boruta(d.comp.se,d.expl.se$lat_cat)
m1.rf.b.se
plot(m1.rf.b.se, las = 2, cex.axis = 0.7)  
getSelectedAttributes(m1.rf.b.se) 
attStats(m1.rf.b.se)

#Running MANOVAS 
d.comp.se.sel <- data.matrix(d.comp.se[,getSelectedAttributes(m1.rf.b.se)])
m1.man.se <- manova(d.comp.se.sel ~ d.expl.se$lat_cat)
summary(m1.man.se)  

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.se)  

par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.se.sel))){
  d.temp=d.comp.se.sel[,i]
  plot(d.temp ~ d.expl.se$lat_cat, ylab=colnames(d.comp.se.sel)[i])
}
dev.off()


#Q2: Which abiotic factors are the most important drivers of fruit quality?---------
#Analysis: principle components analysis followed by PC regression with each variable
#create a vector for all of the PCs

p_clim <- dplyr::select(c,c("Prox.Water","elevation", 
                            "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI"))


#calculate principal components
results <- prcomp(p_clim, scale = TRUE)

#this gives you the % variance explained
summary(results)  

#display principal components
results$x

results$rotation
#display the first six scores
head(results$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results,choices=1:2, cex = 1,
       col = c(NA, "BLACK"),
       scale = FALSE, xlabs = rep("*", 24), 
       xlab = "Principal Component 1", 
       ylab = "Principal Component 2", 
       main = "PCA Biplot",
       cex.lab = 1.2, 
       cex.axis = 1.1, 
       cex.main = 1.5,
       xlim = c(-1.5, 1.5), # Adjust these limits to zoom in on the x-axis
       ylim = c(-1.75, 1.75))

###structuring figure for export###

#Extract loadings
loadings <- results$rotation[, 1:2]  # Replace 'rotation' if using a different object structure

#adding jitter to arrows 
jittered_loadings <- loadings * 2
jittered_loadings[, 1] <- jitter(jittered_loadings[, 1], amount = 0.2)
jittered_loadings[, 2] <- jitter(jittered_loadings[, 2], amount = 0.15)

#Plot
pca <-plot(results$x[, 1:2], type = "n", 
     xlab = "Principal Component 1", 
     ylab = "Principal Component 2", 
     xlim = c(-1.7, 1.7), 
     ylim = c(-1.7, 1.7))
pca <- pca + arrows(0, 0, jittered_loadings[, 1], jittered_loadings[, 2], length = 0.1, col = "black")
pca




#calculate total variance explained by each principal component
summary(results)$importance
summary(results)$importance[2,]

var_explained = results$sdev^2 / sum(results$sdev^2)

df <- data.frame(PC=1:7, var_explained=var_explained)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

#bind the results to the data frame 
pc_clim <- as.data.frame(results$x)

#correlation matrix  
p_clim1 <- dplyr::select(c,c("Prox.Water","elevation", 
                             "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip",
                             "Szn.UVI"))

#The first matrix shows the correlation coefficients between the variables  
#the second matrix shows the corresponding p-values.
rcorr(as.matrix(p_clim1))

#now let's visualize this 
corrplot(cor(p_clim1))

testRes = cor.mtest(p_clim1, conf.level = 0.95)

## add all p-values
corrplot(cor(p_clim1), p.mat = testRes$p, insig = 'p-value', sig.level = -1)

## add significant level stars
corrplot(cor(p_clim1), p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')

#Q2-A: Physical Quality---------------------------------------------------------
##Linear models for PC x quality 
pc_clim_phys <- cbind(pc_clim, c)

###Firmness###
p1 <- glmmTMB(Firmness ~ orchard.type+ maturity.index + PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p1)

plot(Firmness ~ PC1, data=pc_clim_phys)
plot(Firmness ~ PC2, data=pc_clim_phys)
plot(Firmness ~ PC3, data=pc_clim_phys)

results$rotation

###SSC###
P2 <- glmmTMB(SSC ~ orchard.type + maturity.index + PC1 + PC2 + PC3 + PC4 + (1|site.code), data= pc_clim_phys)
summary(P2)

plot(SSC ~ PC1, data=pc_clim_phys)
plot(SSC ~ PC2, data=pc_clim_phys)

results$rotation

###AVGWGT
p3 <- glmmTMB(avgwgt ~ orchard.type + maturity.index + PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p3)

plot(avgwgt ~ PC3, data=pc_clim_phys)

results$rotation


#Q2-B: Fruit Chemistry----------------------------------------------------------
#Linear models for PC X total phenolics and phenolic richness 
#SKIN#
pc.sk_clim <- cbind(pc_clim, SkinD)

#TotalPhen
pc.sk.tp <- glmmTMB(TotalPhenTrans ~ orchard.type+ maturity.index + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.sk_clim, family=beta_family (link="logit"))
summary(pc.sk.tp)

#PhenRich
pc.sk.pr <- glmmTMB(PhenRich ~ orchard.type + maturity.index + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.sk_clim)
summary(pc.sk.pr)

#PULP#
pc.pu_clim <- cbind(pc_clim, PulpD)

##TotalPhen###
pc.pu.tp <- glmmTMB(TotalPhenTrans ~ orchard.type + maturity.index + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim, family=beta_family (link="logit"))
summary(pc.pu.tp)

#PhenRich 
pc.pu.pr <- glmmTMB(PhenRich ~ orchard.type + maturity.index + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim)
summary(pc.pu.pr)

#SEED#
pc.se_clim <- cbind(pc_clim, SeedD)

#TotalPhen
pc.se.tp <- glmmTMB(TotalPhenTrans ~ orchard.type + maturity.index + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.se_clim, family=beta_family (link="logit"))
summary(pc.se.tp)

#PhenRich 
pc.se.pr <- glmmTMB(PhenRich ~ orchard.type + maturity.index + PC1 + PC2 + PC3 + PC4 +
                      (1|site.code), data=pc.se_clim)
summary(pc.se.pr)


#Dummy Coded Corroplot----------------------------------------------------------

cheese <- read.csv("Cheese.csv")
names(cheese)


p_mgmt <- dplyr::select(cheese,c("orchard.type", "Com_Mul", "Cover_Crops", "Weed_Mats", "Herbicides", "Acres",
                            "Mowing", "Cultivation","Prox.Water","elevation", 
                            "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI"))



p_mgmt <- p_mgmt[-c(17), ]

#The first matrix shows the correlation coefficients between the variables  
#the second matrix shows the corresponding p-values.
rcorr(as.matrix(p_mgmt))

#now let's visualize this 
corrplot(cor(p_mgmt))

testRes = cor.mtest(p_mgmt, conf.level = 0.95)

## add all p-values
corrplot(cor(p_mgmt), p.mat = testRes$p, insig = 'p-value', sig.level = -1)

## add significant level stars
corrplot(cor(p_mgmt), p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')


#Q3: Which specific management practices are the most important drivers of fruit quality?----
#For this analysis, we'll be running GLMMs against every mgmt variable. 
#Q3-A: Physical Quality---------------------------------------------------------
#removing orchard number 17 from analysis 
c <- c[-c(17), ]
view(c)

#SSC
mgmt1 <- glmmTMB(SSC ~ orchard.type+ Cultivation + Herbicides + Com_Mul + Mowing +
                   Cover_Crops  + (1|site.code), 
                 data=c)
summary(mgmt1)
Anova(mgmt1)

## Calculate the effect size
emmeans(mgmt1,~ Herbicides   , type="response") 
emmeans(mgmt1,~ Mowing   , type="response") 
emmeans(mgmt1,~ Cover_Crops   , type="response") 


#Average Weight
mgmt2 <- glmmTMB(avgwgt ~ orchard.type + Cultivation + Herbicides + Com_Mul + 
                   Mowing +Cover_Crops + (1|site.code), data=c)
summary(mgmt2)
Anova(mgmt2)

## Calculate the effect size
emmeans(mgmt2,~Mowing, type="response")

#Firmness 
mgmt3 <- glmmTMB(Firmness ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Cover_Crops + (1|site.code), 
                 data=c)
summary(mgmt3)
Anova(mgmt3)

## Calculate the effect size
emmeans(mgmt3,~Herbicides, type="response")

#Maturity Index 
mgmt4 <- glmmTMB(maturity.index ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Cover_Crops + (1|site.code), 
                 data=c)
summary(mgmt4)
Anova(mgmt4)

#calculate the effect size 
emmeans(mgmt4,~Cultivation, type="response")
emmeans(mgmt4,~Herbicides, type="response")
emmeans(mgmt4,~Cover_Crops, type="response")

#Q3-B: Fruit Chemistry----------------------------------------------------------
###Removing Row 17 from analysis because orchard did not fill out this portion of survey data
SkinD <- SkinD[-c(9), ]
PulpD <- PulpD[-c(9), ]
SeedD <- SeedD[-c(9), ]

#Dropping outputs from Weed_mats (only conventional), firemgmt (1 orchard), and acres 

###Skin Total Phenolics###
mgmt.sk <- glmmTMB(TotalPhenTrans ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing 
                   + Cover_Crops + (1|site.code), 
                   data=SkinD, family=beta_family(link="logit"))
summary(mgmt.sk)
Anova(mgmt.sk)

## Calculate the effect size
emmeans(mgmt.sk,  ~ Herbicides, type="response")
emmeans(mgmt.sk,  ~ Com_Mul, type="response")

###Pulp Total Phenolics###
mgmt.pu <- glmmTMB(TotalPhenTrans ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                     Cover_Crops + (1|site.code), 
                   data=PulpD, family=beta_family(link="logit"))
summary(mgmt.pu)
Anova(mgmt.pu)

## Calculate the effect size
emmeans(mgmt.pu, ~ Cultivation, type="response")
emmeans(mgmt.pu, ~ Cover_Crops,type="response")

###Seed Total Phenolics###
mgmt.se <- glmmTMB(TotalPhenTrans ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                     Cover_Crops + (1|site.code), 
                   data=SeedD, family=beta_family(link="logit"))
summary(mgmt.se)
Anova(mgmt.se)

###Skin Phenolic Richness### 
mgmt.pr.sk <- glmmTMB(PhenRich ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                        Cover_Crops + (1|site.code), 
                      data=SkinD)
summary(mgmt.pr.sk)
Anova(mgmt.pr.sk)

###Pulp Phenolic Richness### 
mgmt.pr.pu <- glmmTMB(PhenRich ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                        Cover_Crops + (1|site.code), 
                      data=PulpD)
summary(mgmt.pr.pu)
Anova(mgmt.pr.pu)

## Calculate the effect size
emmeans(mgmt.pr.pu,~Cultivation, type= "response")


###Seed Phenolic Richness### 
mgmt.pr.se <- glmmTMB(PhenRich ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                        Cover_Crops + (1|site.code), 
                      data=SeedD)
summary(mgmt.pr.se)
Anova(mgmt.pr.se)

## Calculate the effect size
emmeans(mgmt.pr.se,~Com_Mul, type= "response")


#Figure 2-------------------------------------------- 
#Physical Traits over Latitude
#SSC
ag1= ggplot(TreeLat, aes(x=Latitude, y=SSC, color=orchard.type)) +
  geom_point(size= 1, position=position_jitterdodge(jitter.width=.2))+
  ylab ("Soluble Sugar Content (°Bx)") +
  xlab ("Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3b0f70", "#de4968"),name="Management System")+
  scale_x_continuous(breaks = seq(34, 49, by = 3))+
  guides(color = "none") 

#Firmness
ag2= ggplot(TreeLat, aes(x=Latitude, y=Firmness, color=orchard.type)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Firnmness (N)") +
  xlab ("Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3b0f70", "#de4968"),name="Management System")+ 
  scale_x_continuous(breaks = seq(34, 49, by = 3))+ 
  guides(color = "none") 


#Average Weight 
ag3= ggplot(TreeLat, aes(x=Latitude, y=avgwgt, color=orchard.type)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Average Weight (g)") +
  xlab ("Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3b0f70", "#de4968"),name="Management System")+
  scale_x_continuous(breaks = seq(34, 49, by = 3))+ 
  guides(color = "none") 

ag4= ggplot(Tree_low, aes(x=Latitude, y=avgwgt, color=orchard.type)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Average Weight (g)") +
  xlab ("Low (<42°) Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3b0f70", "#de4968"),name="Management System")+
  scale_x_continuous(breaks = seq(34, 42, by = 2))+ 
  guides(color = "none") 

ag5= ggplot(Tree_high, aes(x=Latitude, y=avgwgt, color=orchard.type)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Average Weight (g)") +
  xlab ("High (>42°) Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3b0f70", "#de4968"),name="Management System")+
  scale_x_continuous(breaks = seq(44, 49, by = 2))+ 
  guides(color = "none") 



ggarrange(ag1, ag2, ag3, ag4, ag5, nrow = 2, ncol = 3, labels = c("A", "B", "C", "D", "E"))
ggsave("Figure2.jpeg", width=30, height=20, units="cm", dpi=600)


#Figure 3------------------------------------
#Phenolic Richness and Total Phenolics Plot
###Q1-B Phenolic Richness by Tissue Type over Latitude 
ctp = ggplot(ChemLat, aes(x=Latitude, y=TotalPhenTrans, color=Tissue)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Total Phenolics (PPW)") +
  xlab ("Latitude")+
  geom_smooth(method=lm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +
  scale_color_manual(values=c("#22a884", "#fe9f6d", "#414487"),name="Tissue")+
  scale_x_continuous(breaks = seq(34, 49, by = 3))+ 
  guides(shape = guide_legend(override.aes = list(shape = c(16, 17))))+
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "bottom")


cpr = ggplot(ChemLat, aes(x=Latitude, y=PhenRich, color=Tissue)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Phenolic Richness") +
  xlab ("Latitude")+
  geom_smooth(method=lm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +
  scale_color_manual(values=c("#22a884", "#fe9f6d", "#414487"),name="Tissue")+
  scale_x_continuous(breaks = seq(34, 49, by = 3))


cpr  <- cpr  + guides(color = "none")



#mgmt and tissue 
opr = ggplot(ChemLat, aes(x = Tissue, y = PhenRich, color = orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(size= 1,position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  ylab("Phenolic Richness") +
  xlab("") +
  theme_classic() +
  scale_color_manual(values = c("#3b0f70", "#de4968"), name = "Tissue")+
  labs(color = "Tissue")+
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "bottom")

opr <- opr + guides(color = "none")


#mgmt and tissue 
otp = ggplot(ChemLat, aes(x = Tissue, y = TotalPhenTrans, color = orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(size= 1, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  ylab("Total Phenolics (PDW)") +
  xlab("") +
  theme_classic() +
  scale_color_manual(values = c("#3b0f70", "#de4968"), name = "Tissue")+
  labs(color = "Tissue")+
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "bottom")




ggarrange(otp, ctp, opr, cpr, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
ggsave("Figure3.jpeg", width=20, height=24, units="cm", dpi=600)

#Figure 4-----------------------------------------------------------------
#Management practices Plots
###SCC###
#creating data frame to plot all together 
brix <- data.frame(c$SSC,c$avgwgt,c$Firmness, c$Herbicides,c$Mowing, c$Cover_Crops, c$Cultivation, c$Com_Mul )
brix <- brix[-c(17), ]
brix <- data.frame(
  Treatment = rep(c("c.Herbicides", "c.Com_Mul","c.Cultivation","c.Cover_Crops",  "c.Mowing"), each = nrow(brix)),
  Presence = c(brix$c.Herbicides, brix$c.Mowing, brix$c.Cover_Crops, brix$c.Cultivation, brix$c.Com_Mul),
  SSC = rep(brix$c.SSC, 5),
  avgwgt = rep(brix$c.avgwgt, 5),
  Firmness = rep(brix$c.Firmness, 5)
)

#SSC
sscplot = ggplot(brix, aes(x = Treatment, y = SSC, color = Presence)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("SSC (°Bx) ") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  scale_color_manual(values = c("#22a884", "#fe9f6d"), name = "Usage")+
  xlab(" ")+
  theme(
    axis.text.y = element_text(size = 8))+
  guides(fill = FALSE)

sscplot  <- sscplot  + guides(color = "none")

#weight
wgtplot = ggplot(brix, aes(x = Treatment, y = avgwgt, color = Presence)) +
  geom_boxplot(outlier.shape=NA) +
  ylab("Average Weight (g) ") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  scale_color_manual(values = c("#22a884", "#fe9f6d"), name = "Usage")+
  xlab("")+
  theme(
    axis.text.y = element_text(size = 8))+
  guides(fill = FALSE)

wgtplot  <- wgtplot  + guides(color = "none")

#Firmness 
firplot = ggplot(brix, aes(x = Treatment, y = Firmness, color = Presence)) +
  geom_boxplot(outlier.shape=NA) +
  ylab("Firmness (N) ") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  scale_color_manual(values = c("#22a884", "#fe9f6d"), name = "Usage")+
  xlab(" ")+
  theme(
    axis.text.y = element_text(size = 8))+
  guides(fill = FALSE)

firplot  <- firplot  + guides(color = "none")


#####Skin Total Phenolics##

#creating data frame to plot all together 
sk <- data.frame(SkinD$Herbicides,SkinD$Com_Mul,
 SkinD$Cultivation, SkinD$Cover_Crops, SkinD$Mowing, SkinD$TotalPhenTrans, SkinD$PhenRich)

sk <- sk[-c(9), ]


sk <- data.frame(
  Treatment = rep(c("sk$SkinD.Herbicides", "sk$SkinD.Com_Mul", "sk$SkinD.Cultivation",
                    "sk$SkinD.Cover_Crops", "sk$SkinD.Mowing"), each = nrow(sk)),
  Presence = c(sk$SkinD.Herbicides, sk$SkinD.Com_Mul,sk$SkinD.Cultivation, sk$SkinD.Cover_Crops, sk$SkinD.Mowing),
  TP = rep(sk$SkinD.TotalPhenTrans, 5),
  PR = rep(sk$SkinD.PhenRich,5)
)

SKTP = ggplot(sk, aes(x = Treatment, y = TP, color = Presence)) +
  geom_boxplot(outlier.shape=NA) +
  ylab("Skin Total Phenolics (PDW)") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  scale_color_manual(values = c("#22a884", "#fe9f6d"), name = "Usage")+
  xlab("")+
  theme(
    axis.text.y = element_text(size = 8))+
  guides(fill = FALSE)

SKTP  <- SKTP  + guides(color = "none")


SKPR = ggplot(sk, aes(x = Treatment, y = PR, color = Presence)) +
  geom_boxplot(outlier.shape=NA) +
  ylab("Skin Phenolic Richness") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  scale_color_manual(values = c("#22a884", "#fe9f6d"), name = "Usage")+
  xlab("Cover Crops       Herbicides        Mowing")+
  theme(
    axis.text.y = element_text(size = 8))+
  guides(fill = FALSE)

SKPR  <- SKPR  + guides(color = "none")


#####Pulp Total Phenolics##
##creating data frame to plot all together 
PU <- data.frame(PulpD$Herbicides,PulpD$Com_Mul,
                 PulpD$Cultivation, PulpD$Cover_Crops, PulpD$Mowing, PulpD$TotalPhenTrans, PulpD$PhenRich)
PU <- PU[-c(9), ]

PU <- data.frame(
  Treatment = rep(c("PU$PulpD.Herbicides", "PU$PulpD.Com_Mul", "PU$PulpD.Cultivation",
                    "PU$PulpD.Cover_Crops", "PU$PulpD.Mowing"), each = nrow(PU)),
  Presence = c(PU$PulpD.Herbicides, PU$PulpD.Com_Mul,PU$PulpD.Cultivation, PU$PulpD.Cover_Crops, PU$PulpD.Mowing),
  TP = rep(PU$PulpD.TotalPhenTrans, 5),
  PR = rep(PU$PulpD.PhenRich,5)
)

PUTP = ggplot(PU, aes(x = Treatment, y = TP, color = Presence)) +
  geom_boxplot(outlier.shape=NA) +
  ylab("Pulp Total Phenolics (PDW)") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  scale_color_manual(values = c("#22a884", "#fe9f6d"), name = "Usage")+
  xlab("")+
  theme(
    axis.text.y = element_text(size = 8))+
  guides(fill = FALSE)
PUTP  <- PUTP  + guides(color = "none")


PUPR = ggplot(PU, aes(x = Treatment, y = PR, color = Presence)) +
  geom_boxplot(outlier.shape=NA) +
  ylab("Pulp Phenolic Richness") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  scale_color_manual(values = c("#22a884", "#fe9f6d"), name = "Usage")+
  xlab("")+
  theme(
    axis.text.y = element_text(size = 8))+
  guides(fill = FALSE)

PUPR  <- PUPR  + guides(color = "none")


#Seed Phenolic Richness
SE <- data.frame(SeedD$Herbicides,SeedD$Com_Mul,
                 SeedD$Cultivation, SeedD$Cover_Crops, SeedD$Mowing, SeedD$TotalPhenTrans, SeedD$PhenRich)
SE <- SE[-c(9), ]

SE <- data.frame(
  Treatment = rep(c("SE$SeedD.Herbicides", "SE$SeedD.Com_Mul", "SE$SeedD.Cultivation",
                    "SE$SeedD.Cover_Crops", "SE$SeedD.Mowing"), each = nrow(SE)),
  Presence = c(SE$SeedD.Herbicides, SE$SeedD.Com_Mul,SE$SeedD.Cultivation, SE$SeedD.Cover_Crops, SE$SeedD.Mowing),
  TP = rep(SE$SeedD.TotalPhenTrans, 5),
  PR = rep(SE$SeedD.PhenRich,5)
)

SETP = ggplot(SE, aes(x = Treatment, y = TP, color = Presence)) +
  geom_boxplot(outlier.shape=NA) +
  ylab("Seed Total Phenolics (PDW)") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  scale_color_manual(values = c("#22a884", "#fe9f6d"), name = "Usage")+
  xlab("")+
  theme(
    axis.text.y = element_text(size = 8))+
  guides(fill = FALSE)

SETP  <- SETP  + guides(color = "none")


SEPR = ggplot(PU, aes(x = Treatment, y = PR, color = Presence)) +
  geom_boxplot(outlier.shape=NA) +
  ylab("Seed Phenolic Richness") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  scale_color_manual(values = c("#22a884", "#fe9f6d"), name = "Usage")+
  xlab("")+
  theme(
    axis.text.y = element_text(size = 8))+
  guides(fill = FALSE)

SEPR  <- SEPR  + guides(color = "none")



#create one large plot
ggarrange(sscplot, wgtplot, firplot, SKTP, SKPR, PUTP, PUPR, SETP, SEPR
          ,nrow = 3, ncol = 3, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "F"))
ggsave("Figure4.jpeg", width=20, height=20, units="cm", dpi=600)



#Aggregating Data for tables---------------------------------------------------
d.pu <- left_join(d.pu, Orchard[,c(2,40)], by="orchard.num")
d.sk <- left_join(d.sk, Orchard[,c(2,40)], by="orchard.num")
d.se <- left_join(d.se, Orchard[,c(2,40)], by="orchard.num")



aggregate(d, A~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, B~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, C~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, E~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, F~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, G~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, H~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, PB1~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, gentistic_acid~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, chlorogenic_acid~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, syringic_acid~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, U2~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, K~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, isoquercitin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, L~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, quercetin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, p_coumaric_acid~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, cyanidin_Galactoside~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, caffeic_acid~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, ePicatechin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, I~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, rutin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, reynoutrin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, quercitrin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, phloretin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, ferulic_acid~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, catechin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, PB2~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, U1~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, J~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, M~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, hyperin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, avicularin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, phloridzin~orchard.type+Tissue, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))





aggregate(d, A~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, B~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, C~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, E~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, F~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, G~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, H~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, PB1~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, gentistic_acid~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, chlorogenic_acid~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, syringic_acid~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, U2~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, K~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, isoquercitin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, L~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, quercetin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, p_coumaric_acid~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, cyanidin_Galactoside~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, caffeic_acid~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, ePicatechin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, I~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, rutin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, reynoutrin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, quercitrin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, phloretin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, ferulic_acid~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, catechin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, PB2~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, U1~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, J~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, M~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, hyperin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, avicularin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
aggregate(d, phloridzin~orchard.type, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))



aggregate(d.pu, caffeic_acid~lat_cat, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))





#prctice 
#Finding r squared values for models--------------------------------------------

#for figure 2

#Subsetting data into Organic and Conventional 
TreeO <- filter(TreeLat, orchard.type == "Organic")
TreeC <- filter(TreeLat, orchard.type == "Conventional")

#Fit separate linear models
ssco <- lm(SSC ~ Latitude, data = TreeO)
sscc <- lm(SSC ~ Latitude, data = TreeC)
firo <- lm(Firmness ~ Latitude, data = TreeO)
firc <- lm(Firmness ~ Latitude, data = TreeC)
wgto <- lm(avgwgt ~ Latitude, data = TreeO)
wgtc <- lm(avgwgt ~ Latitude, data = TreeC)

#extract R-squared 
summary(ssco)$r.squared
summary(sscc)$r.squared
summary(firo)$r.squared
summary(firc)$r.squared
summary(wgto)$r.squared
summary(wgtc)$r.squared

#extract p values
summary(ssco)$coefficients
summary(sscc)$coefficients
summary(firo)$coefficients
summary(firc)$coefficients
summary(wgto)$coefficients
summary(wgtc)$coefficients

#sub again for low 
LowO <- filter(Tree_low, orchard.type == "Organic")
LowC <- filter(Tree_low, orchard.type == "Conventional")

lwo <- lm(avgwgt ~ Latitude, data = LowO)
lwc <- lm(avgwgt ~ Latitude, data = LowC)

summary(lwo)$r.squared
summary(lwc)$r.squared

summary(lwo)$coefficients
summary(lwc)$coefficients

#sub again for high 
HighO <- filter(Tree_high, orchard.type == "Organic")
HighC <- filter(Tree_high, orchard.type == "Conventional")

hwo <- lm(avgwgt ~ Latitude, data = HighO)
hwc <- lm(avgwgt ~ Latitude, data = HighC)

summary(hwo)$r.squared
summary(hwc)$r.squared

summary(hwo)$coefficients
summary(hwc)$coefficients

view(ChemLat)


#for figure 3#
#Subsetting data into tissue 
SkinR <- filter(ChemLat, Tissue == "SKIN")
PulpR <- filter(ChemLat, Tissue == "PULP")
SeedR <- filter(ChemLat, Tissue == "SEED")


s1 <- lm(TotalPhenTrans ~ Latitude, data = SkinR)
s2 <- lm(PhenRich ~ Latitude, data = SkinR)

p1 <- lm(TotalPhenTrans ~ Latitude, data = PulpR)
p2 <- lm(PhenRich ~ Latitude, data = PulpR)

d1 <- lm(TotalPhenTrans ~ Latitude, data = SeedR)
d2 <- lm(PhenRich ~ Latitude, data = SeedR)


summary(s1)$r.squared
summary(p1)$r.squared
summary(d1)$r.squared

summary(s1)$coefficients
summary(p1)$coefficients
summary(d1)$coefficients

summary(s2)$r.squared
summary(p2)$r.squared
summary(d2)$r.squared


summary(s2)$coefficients
summary(p2)$coefficients
summary(d2)$coefficients




#messing around ----------------------------------------------------------------

M <- left_join(Tree, Orchard %>% select(orchard.num, Latitude, lat_cat, AHD, GHD), 
                     by = c("orchard.num"))

M$AHD <- as.Date(M$AHD, format = "%Y%m%d")
view(M)

#SSC
m1<- glmmTMB(SSC ~ AHD + maturity.index + (1|site.code/orchard.num), data=M)
summary(m1)
Anova(m1) 

plot(SSC ~ AHD, data=M)
plot(SSC ~ Latitude, data=M)

m1<- glmmTMB(avgwgt ~ AHD + (1|site.code/orchard.num), data=M)
summary(m1)
Anova(m1) 

plot(avgwgt ~ AHD, data=M)
plot(avgwgt ~ Latitude, data=M)
