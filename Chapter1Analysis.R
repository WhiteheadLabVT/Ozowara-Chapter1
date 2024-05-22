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



#Read Data Organization and Restructuring---------------------------------------
#Read in Data sheets
#Orchard Level Data (24 obs)
Orchard <- read_csv("Ozowara_et_al_OrchardLevelData.csv")
Orchard$orchard.num <- as.factor(as.character(Orchard$orchard.num))

#Tree Level Data (120 obs)
Tree <- read_csv("Ozowara_et_al_TreeLevelData.csv")
View(Tree_Level_Data)
Tree$orchard.num <- as.factor(as.character(Tree$orchard.num))

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

#creating pest index
Orchard <- Orchard %>%
  mutate(Pest.Index=rowSums(across(21:38))/18)

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
TreeLat <- left_join(Tree, Orchard[,c(2,4)], by="orchard.num")
#SSC
Lat1<- glmmTMB(SSC ~ orchard.type*Latitude + (1|site.code/orchard.num), data=TreeLat)
summary(Lat1)
Anova(Lat1) 

#calculate effect size
emmeans(Lat1, ~ Latitude, at = list(Latitude = c(34:48)), type= "reponse")

#Firmness 
Lat2<- glmmTMB(Firmness ~ orchard.type*Latitude + (1|site.code/orchard.num), data=TreeLat)
summary(Lat2)
Anova(Lat2)

#calculate effect size
emmeans(Lat2, ~ Latitude, at = list(Latitude = c(34:48)), type= "reponse")

#Average Weight 
Lat3<- glmmTMB(avgwgt ~ orchard.type*Latitude + (1|site.code/orchard.num), data=TreeLat)
summary(Lat3)
Anova(Lat3) 

#calculate effect size
emmeans(Lat3, ~ Latitude|orchard.type, at = list(Latitude = c(34:48)), type= "reponse")


#Maturity Index 
Lat4<- glmmTMB(maturity.index ~ orchard.type*Latitude + (1|site.code/orchard.num), 
               data=TreeLat)
summary(Lat4)
Anova(Lat4)  

#calculate effect size
emmeans(Lat4, ~ Latitude, at = list(Latitude = c(34:48)), type= "reponse")


###Investigating Latitude and Average weight 
#Splitting latitude in high and low 
Tree_low <- filter(TreeLat, Latitude<42)
Tree_high <- filter(TreeLat, Latitude>42)

#Low Latitude Sites 
avg_low<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/orchard.num), data=Tree_low)
summary(avg_low)
Anova(avg_low) 

#High Latitude Sites 
avg_hgh<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/orchard.num), data=Tree_high)
summary(avg_hgh)
Anova(avg_hgh)

## Calculate the effect size
emmeans(avg_hgh,pairwise~orchard.type, type="response")

#Q1-B: Fruit Chemistry-------------------------------------------------------------
#binding latitude to the 359 obs data set 
ChemLat <- left_join(d, Orchard[,c(2,4,40)], by="orchard.num")

#Total Phenolics with beta distribution 
tp1 <- glmmTMB(TotalPhenTrans ~ orchard.type*Tissue*Latitude + 
                 (1|site.code/orchard.num/Tree), data=ChemLat, family=beta_family(link="logit"))
summary(tp1)
Anova(tp1)

#binding latitudinal data to tissue data sets
d.sk <- left_join(d.sk, Orchard[,c(2,4)], by="orchard.num")
d.pu <- left_join(d.pu, Orchard[,c(2,4)], by="orchard.num")
d.se <- left_join(d.se, Orchard[,c(2,4)], by="orchard.num")

###total phenolics###
#skin
tp.sk <- glmmTMB(TotalPhenTrans ~ orchard.type*Latitude + (1|site.code/orchard.num), 
                 data=d.sk, family=beta_family(link="logit"))
summary(tp.sk)
Anova(tp.sk)
#nothing 

#pulp
tp.pu <- glmmTMB(TotalPhenTrans ~ orchard.type*Latitude + (1|site.code/orchard.num), 
                 data=d.pu, family=beta_family(link="logit"))
summary(tp.pu)
Anova(tp.pu)
#nothing 

#seed
tp.se <- glmmTMB(TotalPhenTrans ~ orchard.type*Latitude + (1|site.code/orchard.num/Tree), 
                 data=d.se, family=beta_family(link="logit"))
summary(tp.se)
Anova(tp.se)

#calculate effect size
emmeans(tp.se, pairwise~orchard.type)


###phenolics richness###
pr1 <- glmmTMB(PhenRich~ orchard.type*Latitude*Tissue + (1|site.code/orchard.num/Tree), data=ChemLat, 
               family=poisson(link="log"))
summary(pr1)
Anova(pr1)

#pulp
pr.p <- glmmTMB(PhenRich ~ orchard.type*Latitude+(1|site.code/orchard.num), data=d.pu,
                family=poisson(link="log"))
summary(pr.p)
Anova(pr.p)

#skin
pr.sk <- glmmTMB(PhenRich ~ orchard.type*Latitude+(1|site.code/orchard.num), data=d.sk,
                 family=poisson(link="log"))
summary(pr.sk)
Anova(pr.sk)

#seed
pr.se <- glmmTMB(PhenRich ~ orchard.type*Latitude+(1|site.code/orchard.num), data=d.se,
                 family=poisson(link="log"))
summary(pr.se)
Anova(pr.se)

#calculate effect size
emmeans(pr.se, ~ Latitude, at = list(Latitude = c(34:48)), type= "reponse")


#Q1-C: Which compounds distinguish fruits based on management systems?----------
#Analysis: NMDS and random forest
#Q1-C: NMDS + PERMANOVA---------------------------------------------------------
#left joining the "expl" by latitude  
d.expl <- left_join(d.expl, Orchard[,c(2,40)], by="orchard.num")
d.expl.sk <- left_join(d.expl.sk, Orchard[,c(2,40)], by="orchard.num")
d.expl.pu <- left_join(d.expl.pu, Orchard[,c(2,40)], by="orchard.num")
d.expl.se <- left_join(d.expl.se, Orchard[,c(2,40)], by="orchard.num")

###NMDS for Skin 
m.NMDS.sk <- metaMDS(d.comp.sk, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.sk

plot(m.NMDS.sk, type="t")


#for plotting, need to add columns that give values for 
#colors and symbols we want in plot
d.expl.sk$Color <- recode_factor(d.expl.sk$lat_cat,
                                 low="red", high="blue")
d.expl.sk$Symbol <- recode_factor(d.expl.sk$orchard.type,
                                  Organic=2, Conventional=1)
d.expl.sk$Symbol <- as.numeric(as.character(d.expl.sk$Symbol))

d.expl.sk$lat_cat = as.factor(d.expl.sk$lat_cat)


plot(m.NMDS.sk, type="n") #plots the ordination axes only
points(m.NMDS.sk, pch=d.expl.sk$Symbol,
       col=as.character(d.expl.sk$Color), cex = 0.8)     
ordiellipse(m.NMDS.sk, d.expl.sk$Symbol, conf = 0.95)

#PERMANOVA can test whether the visualized differences are significant
m.perm1 <- adonis2(d.comp.sk~orchard.type*lat_cat, data=d.expl.sk)
m.perm1


###NMDS for Pulp 
m.NMDS.pu <- metaMDS(d.comp.pu, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.pu

plot(m.NMDS.pu, type="t")

#for plotting, need to add columns that give values for 
#colors and symbols we want in plot
d.expl.pu$Color <- recode_factor(d.expl.pu$lat_cat,
                                 low="red", high="blue")
d.expl.pu$Symbol <- recode_factor(d.expl.pu$orchard.type,
                                  Organic=2, Conventional=1)
d.expl.pu$Symbol <- as.numeric(as.character(d.expl.pu$Symbol))

d.expl.pu$lat_cat = as.factor(d.expl.pu$lat_cat)


NMDSpu=plot(m.NMDS.pu, type="n") #plots the ordination axes only
points(m.NMDS.pu, pch=d.expl.pu$Symbol,
       col=as.character(d.expl.pu$Color), cex = 0.8)     
ordiellipse(m.NMDS.pu, d.expl.pu$Symbol, conf = 0.95)
NMDSpu

#PERMANOVA
m.perm2 <- adonis2(d.comp.pu~orchard.type*lat_cat, data=d.expl.pu)
m.perm2

###NMDS for Seed 
m.NMDS.se <- metaMDS(d.comp.se, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.se

plot(m.NMDS.se, type="t")

#for plotting, need to add columns that give values for 
#colors and symbols we want in plot
d.expl.se$Color <- recode_factor(d.expl.se$orchard.type,
                                 Organic="red", Conventional="blue")
d.expl.se$Symbol <- recode_factor(d.expl.se$lat_cat,
                                  high=2, low=1)
d.expl.se$Symbol <- as.numeric(as.character(d.expl.se$Symbol))

d.expl.se$lat_cat = as.factor(d.expl.se$lat_cat)


NMDSse=plot(m.NMDS.se, type="n") #plots the ordination axes only
points(m.NMDS.se, pch=d.expl.se$Symbol,
       col=as.character(d.expl.se$Color), cex = 0.8)     
ordiellipse(m.NMDS.se, d.expl.se$Symbol, conf = 0.95)
NMDSse

#PERMANOVA
m.perm3 <- adonis2(d.comp.se~orchard.type*lat_cat, data=d.expl.se)
m.perm3

#Q1-C: Random Forest------------------------------------------------------------
###Skin was dropped 

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

##Running MANOVAS (needs to be fixed)
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
m1.rf.sk <- randomForest(d.comp.sk,d.expl.sk$lat_cat, d.expl.sk$orchard.type, importance=TRUE, 
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
biplot(results,choices=1:2, cex = 0.75,
       col = c('darkblue', 'BLACK'),
       scale = FALSE, xlabs = rep("*", 24))



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

#correlation matrix including Latitude and Longitude 
p_clim1 <- dplyr::select(c,c("Prox.Water","elevation", 
"Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip",
"Szn.UVI","Latitude", "Longitude"))

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
p1 <- glmmTMB(Firmness ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p1)

plot(Firmness ~ PC1, data=pc_clim_phys)
plot(Firmness ~ PC2, data=pc_clim_phys)
plot(Firmness ~ PC3, data=pc_clim_phys)

results$rotation

###SSC###
P2 <- glmmTMB(SSC ~ orchard.type + PC1 + PC2 + PC3 + PC4 + (1|site.code), data= pc_clim_phys)
summary(P2)

plot(SSC ~ PC1, data=pc_clim_phys)
plot(SSC ~ PC2, data=pc_clim_phys)

results$rotation

###AVGWGT
p3 <- glmmTMB(avgwgt ~ orchard.type + PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p3)

plot(avgwgt ~ PC2, data=pc_clim_phys)

results$rotation

###Maturity Index 
p4 <- glmmTMB(maturity.index ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p4)

plot(maturity.index ~ PC1, data=pc_clim_phys)

results$rotation

#Q2-B: Fruit Chemistry----------------------------------------------------------
#Linear models for PC X total phenolics and phenolic richness 
#SKIN#
pc.sk_clim <- cbind(pc_clim, SkinD)

#TotalPhen
pc.sk.tp <- glmmTMB(TotalPhenTrans ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.sk_clim, family=beta_family (link="logit"))
summary(pc.sk.tp)

#PhenRich
pc.sk.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.sk_clim)
summary(pc.sk.pr)

#PULP#
pc.pu_clim <- cbind(pc_clim, PulpD)

##TotalPhen###
pc.pu.tp <- glmmTMB(TotalPhenTrans ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim, family=beta_family (link="logit"))
summary(pc.pu.tp)

#PhenRich 
pc.pu.pr <- glmmTMB(PhenRich ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim)
summary(pc.pu.pr)

#SEED#
pc.se_clim <- cbind(pc_clim, SeedD)

#TotalPhen
pc.se.tp <- glmmTMB(TotalPhenTrans ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.se_clim, family=beta_family (link="logit"))
summary(pc.se.tp)

#PhenRich 
pc.se.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 +
                      (1|site.code), data=pc.se_clim)
summary(pc.se.pr)

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
emmeans(mgmt1,pairwise~Herbicides, type="response")
emmeans(mgmt1,pairwise~Mowing, type="response")
emmeans(mgmt1,pairwise~Cover_Crops, type="response")


#Average Weight
mgmt2 <- glmmTMB(avgwgt ~ orchard.type + Cultivation + Herbicides + Com_Mul + 
Mowing +Cover_Crops + (1|site.code), data=c)
summary(mgmt2)
Anova(mgmt2)

## Calculate the effect size
emmeans(mgmt2,pairwise~Mowing, type="response")

#Firmness 
mgmt3 <- glmmTMB(Firmness ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
Cover_Crops + (1|site.code), 
                 data=c)
summary(mgmt3)
Anova(mgmt3)

## Calculate the effect size
emmeans(mgmt3,pairwise~Herbicides, type="response")

#Maturity Index 
mgmt4 <- glmmTMB(maturity.index ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
Cover_Crops + (1|site.code), 
                 data=c)
summary(mgmt4)
Anova(mgmt4)

#calculate the effect size 
emmeans(mgmt4,pairwise~Cultivation, type="response")
emmeans(mgmt4,pairwise~Herbicides, type="response")
emmeans(mgmt4,pairwise~Cover_Crops, type="response")

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
emmeans(mgmt.sk, pairwise ~ Herbicides)
emmeans(mgmt.sk, pairwise ~ Com_Mul)

###Pulp Total Phenolics###
mgmt.pu <- glmmTMB(TotalPhenTrans ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
 Cover_Crops + (1|site.code), 
                   data=PulpD, family=beta_family(link="logit"))
summary(mgmt.pu)
Anova(mgmt.pu)

## Calculate the effect size
emmeans(mgmt.pu, pairwise ~ Cultivation)
emmeans(mgmt.pu, pairwise ~ Cover_Crops)

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
emmeans(mgmt.pr.pu,pairwise~Cultivation)


###Seed Phenolic Richness### 
mgmt.pr.se <- glmmTMB(PhenRich ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
Cover_Crops + (1|site.code), 
                      data=SeedD)
summary(mgmt.pr.se)
Anova(mgmt.pr.se)

## Calculate the effect size
emmeans(mgmt.pr.se,pairwise~Com_Mul)


#Total Phenolics & PhenRich Plot by Mgmt--------------------------------------------------------------------------- 
#Skin total Phenolics 
sk.tp.ag=ggplot(d.sk, aes(x=orchard.type, y=TotalPhen/1000, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Skin Total Phenolics (ug/g)") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")


sk.tp.ag=ggplot(d.sk, aes(x=orchard.type, y=TotalPhen/1000000, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Skin Total Phenolics (pdw)") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")


#remove legend for multiplot
sk.tp.ag <- sk.tp.ag + guides(color = "none")


#skin phenolic richness 
sk.pr.ag=ggplot(d.sk, aes(x=orchard.type, y=PhenRich, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Skin Phenolic Richness") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

#remove legend for multiplot
sk.pr.ag <- sk.pr.ag + guides(color = "none")


#pulp total phenolics 
pu.tp.ag=ggplot(d.pu, aes(x=orchard.type, y=TotalPhen/1000, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Pulp Total Phenolics (ug/g)") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

#remove legend for multiplot
pu.tp.ag <- pu.tp.ag + guides(color = "none")


#pulp phenolic richness 
pu.pr.ag=ggplot(d.pu, aes(x=orchard.type, y=PhenRich, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Pulp Phenolic Richness") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")


#remove legend for multiplot
pu.pr.ag <- pu.pr.ag + guides(color = "none")


#seed total phenolics 
se.tp.ag=ggplot(d.se, aes(x=orchard.type, y=TotalPhen/1000, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Seed Total Phenolics (ug/g)") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

#remove legend for multiplot
se.tp.ag <- se.tp.ag + guides(color = "none")

#seed phenolic richness 
se.pr.ag=ggplot(d.se, aes(x=orchard.type, y=PhenRich, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Seed Phenolic Richness") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

#remove legend for multiplot
se.pr.ag <- se.pr.ag + guides(color = "none")


ggarrange(sk.pr.ag, sk.tp.ag, pu.pr.ag,
          pu.tp.ag, se.pr.ag,se.tp.ag,nrow = 2, ncol = 3, labels = c("A", "B", "C", "D", "E", "F"))
ggsave("chemmgmt.png", width=10, height=16, units="cm", dpi=600)


#Physical Traits over Latitude-------------------------------------------- 
#SSC
ag1= ggplot(TreeLat, aes(x=Latitude, y=SSC, color=orchard.type)) +
  geom_point(size= 1, position=position_jitterdodge(jitter.width=.2))+
  ylab ("Soluble Sugar Content") +
  xlab ("Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_continuous(breaks = seq(34, 49, by = 3))+ 
  guides(shape = guide_legend(override.aes = list(shape = c(16, 17))))+
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "top")


#Firmness
ag2= ggplot(TreeLat, aes(x=Latitude, y=Firmness, color=orchard.type)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Firnmness (N)") +
  xlab ("Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+ 
  scale_x_continuous(breaks = seq(34, 49, by = 3))+ 
  guides(color = "none")  


#Average Weight 
ag3= ggplot(TreeLat, aes(x=Latitude, y=avgwgt, color=orchard.type)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Average Weight (g)") +
  xlab ("Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_continuous(breaks = seq(34, 49, by = 3))+ 
  guides(color = "none") 

#Maturity Index 
ag4= ggplot(TreeLat, aes(x=Latitude, y=maturity.index, color=orchard.type)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("CSI Value") +
  xlab ("Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+  
  scale_x_continuous(breaks = seq(34, 49, by = 3))+ 
  guides(color = "none") 


#avg high 
ag5= ggplot(Tree_low, aes(x=Latitude, y=avgwgt, color=orchard.type)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Average Weight (g)") +
  xlab ("Low Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_continuous(breaks = seq(34, 41, by = 2))+ 
  guides(color = "none") 

#avg high 
ag6= ggplot(Tree_high, aes(x=Latitude, y=avgwgt, color=orchard.type)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Average Weight (g)") +
  xlab ("High Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +  
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_continuous(breaks = seq(42, 49, by = 2))+ 
guides(color = "none") 



ggarrange(ag1, ag3, ag2, ag4, ag5, ag6, nrow = 3, ncol = 2, labels = c("A", "B", "C", "D", "E", "F"))
ggsave("fig1.png", width=20, height=24, units="cm", dpi=600)


#Phenolic Richness and Total Phenolics Plot------------------------------------

###Q1-B Phenolic Richness by Tissue Type over Latitude 
ctp = ggplot(ChemLat, aes(x=Latitude, y=TotalPhenTrans, color=Tissue)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Total Phenolics") +
  xlab ("Latitude")+
  geom_smooth(method=lm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +
  scale_color_manual(values=c("#3EBCD2", "#9A607F", "darkgreen"),name="Tissue")+
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
  scale_color_manual(values=c("#3EBCD2", "#9A607F", "darkgreen"),name="Tissue")+
  scale_x_continuous(breaks = seq(34, 49, by = 3))
  

cpr  <- cpr  + guides(color = "none")



#mgmt and tissue 
opr = ggplot(ChemLat, aes(x = Tissue, y = PhenRich, color = orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(size= 1,position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  ylab("Phenolic Richness") +
  xlab("") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Tissue")+
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
  ylab("Total Phenolics") +
  xlab("") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Tissue")+
  labs(color = "Tissue")+
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "bottom")


#seed total phenolics 
se.tp.ag=ggplot(d.se, aes(x=orchard.type, y=TotalPhenTrans, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab ("Seed Total Phenolics") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
 guides(color = "none")



se.pr.l = ggplot(SeedD, aes(x=Latitude, y=PhenRich, color=orchard.type)) +
  geom_point(size= 1,position=position_jitterdodge(jitter.width=.2))+
  ylab (" Seed Phenolic Richness") +
  xlab ("Latitude")+
  geom_smooth(method=lm ,alpha = .15,aes(fill = NULL))+
  theme_classic() +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="mgmt")+
  scale_x_continuous(breaks = seq(34, 49, by = 3))+
  guides(color = "none")


ggarrange(otp, ctp, opr, cpr, se.tp.ag , se.pr.l, nrow = 3, ncol = 2, labels = c("A", "B", "C", "D", "E", "F"))
ggsave("fig2.png", width=20, height=24, units="cm", dpi=600)

#Management practices Plots-----------------------------------------------------------------

#SSC
sxh = ggplot(c, aes(x = Herbicides, y = SSC, fill = Herbicides)) +
  geom_boxplot(outlier.shape=NA)+
  geom_smooth(method=glm, se=FALSE)+
   ylab("SSC") +
  xlab("Herbicides") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Herbicides")+  
  guides(fill = FALSE)  

sxm = ggplot(c, aes(x = Mowing, y = SSC, fill = Mowing)) +
  geom_boxplot(outlier.shape=NA)+
  geom_smooth(method=glm, se=FALSE)+
  ylab("SSC") +
  xlab("Mowing") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Herbicides")+  
  guides(fill = FALSE)  

sxcc = ggplot(c, aes(x = Cover_Crops, y = SSC, fill = Cover_Crops)) +
  geom_boxplot(outlier.shape=NA)+
  geom_smooth(method=glm, se=FALSE)+
  ylab("SSC") +
  xlab("Cover Crops") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Herbicides")+  
  guides(fill = FALSE)  

#Average Weight
axc = ggplot(c, aes(x = Mowing, y = avgwgt, fill = Mowing)) +
  geom_boxplot(outlier.shape=NA)+
  geom_smooth(method=glm, se=FALSE)+    ylab("Average Weight (g)") +
  xlab("Mowing") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
  guides(fill = FALSE)  

#Firmness
fxh = ggplot(c, aes(x = Herbicides, y = Firmness, fill = Herbicides)) +
  geom_boxplot(outlier.shape=NA)+
  geom_smooth(method=glm, se=FALSE)+ylab("Firmness (N)") +
  xlab("Herbicides") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
  guides(fill = FALSE)  

#Maturity Index 
mxh = ggplot(c, aes(x = Herbicides, y = maturity.index, fill = Herbicides)) +
  geom_boxplot(outlier.shape=NA)+
  geom_smooth(method=glm, se=FALSE)+  ylab("CSI Value") +
  xlab("Herbicides") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
  guides(fill = FALSE)  

mxcul = ggplot(c, aes(x = Cultivation, y = maturity.index, fill = Cultivation)) +
  geom_boxplot(outlier.shape=NA)+
  geom_smooth(method=glm, se=FALSE)+  ylab("CSI Value") +
  xlab("Cultivation") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
  guides(fill = FALSE)  

mxcc = ggplot(c, aes(x = Cover_Crops, y = maturity.index, fill = Cover_Crops)) +
  geom_boxplot(outlier.shape=NA)+
  geom_smooth(method=glm, se=FALSE)+  ylab("CSI Value") +
  xlab("Cover Crops") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
  guides(fill = FALSE)  

ggarrange(mxh, mxcul, mxcc, axc, fxh, sxh, sxm, sxcc
          ,nrow = 4, ncol = 2, labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
ggsave("mgmt_prac.png", width=16, height=20, units="cm", dpi=600)




#Skin Total Phenolics 
stphxh = ggplot(SkinD, aes(x = Herbicides, y = TotalPhenTrans, fill = Herbicides)) +
  geom_boxplot(outlier.shape=NA)+
  geom_smooth(method=glm, se=FALSE)+ 
  ylab("CSI Value") +  ylab("Skin Total Phenolics PDW") +
  xlab("Herbicides") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
  guides(fill = FALSE)  

stpxwm = ggplot(SkinD, aes(x = Com_Mul, y = TotalPhenTrans, fill = Com_Mul)) +
  geom_boxplot(outlier.shape=NA)+
  ylab("Skin Total Phenolics PDW") +
  xlab("Weed Mats") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
  guides(fill= FALSE)


#Pulp Total Phenolics 
pxc = ggplot(PulpD, aes(x = Cultivation, y = TotalPhenTrans, fill = Cultivation)) +
  geom_boxplot(outlier.shape=NA)+
  ylab("Pulp Total Phenolics PDW") +
  xlab("Cultivation") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
guides(fill= FALSE)

pxcc = ggplot(PulpD, aes(x = Cover_Crops, y = TotalPhenTrans, fill = Cover_Crops)) +
  geom_boxplot(outlier.shape=NA)+
  ylab("Pulp Total Phenolics PDW") +
  xlab("Cover Crops") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
  guides(fill= FALSE)



#Pulp Phenolic Richness
prcm = ggplot(PulpD, aes(x = Cultivation, y = TotalPhenTrans, fill = Cultivation)) +
  geom_boxplot(outlier.shape=NA)+
  ylab("Pulp Phenolic Richness") +
  xlab("Cultivation") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
  guides(fill= FALSE)

#Seed Phenolic Richness
secm = ggplot(SeedD, aes(x = Com_Mul, y = TotalPhenTrans, fill = Com_Mul)) +
  geom_boxplot(outlier.shape=NA)+
  ylab("Seed Phenolic Richness") +
  xlab("Compost/Mulching") +
  theme_classic() +
  scale_color_manual(values = c("#3EBCD2", "#9A607F"), name = "Management System")+
  guides(fill= FALSE)






ggarrange(stphxh, stpxwm, pxc, pxcc, prcm, secm, 
          nrow = 3, ncol = 2, labels = c("A", "B", "C", "D", "E", "F"))
ggsave("mgmt_prac_chem.png", width=16, height=20, units="cm", dpi=600)


#Map Plot----------------------------------------------------------------------
install.packages("usmap")
library(usmap)
library(ggplot2)

# Your data transformation code
c1 <- data.frame(
  Orchard = c$orchard.num,
  Latitude = c$Latitude,
  Longitude = c$Longitude,
  Otype = c$orchard.type)

#transform the points 
us_map <- usmap_transform(data = c1, input_names = c("Longitude", "Latitude"))

c2 <- data.frame(
  Geometry = us_map$geometry)

us_map <- plot_usmap(include = c("CA", "OR", "WA"))

# Add your latitude and longitude points
map_with_points <- us_map +
  geom_point(data = c2, color = "red", size = 3)

print(map_with_points)


