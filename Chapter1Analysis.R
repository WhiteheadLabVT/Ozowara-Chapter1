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
#Read Data Organization and Restructuring---------------------------------------
#Read in Data sheets
#Orchard Level Data (24 obs)
Orchard <- read_csv("Orchard Level Data.csv")
Orchard$orchard.num <- as.factor(as.character(Orchard$orchard.num))

#Tree Level Data (120 obs)
Tree <- read_csv("Tree Level Data.csv")
View(Tree_Level_Data)
Tree$orchard.num <- as.factor(as.character(Tree$orchard.num))

#Phenolics Data (359 obs)
d <-read_csv("FruitLevelData_revised - Sheet1.csv")
d$orchard.num <- as.factor(as.character(d$orchard.num))

#Organizing Phenolics Data for Q1#

#creating totalphen and phenrich
d <- d %>%
  mutate(TotalPhen=rowSums(across(8:41)),
         PhenRich=rowSums(across(8:41)!=0))

#testing 
shapiro.test(d$TotalPhen)
#W = 0.83175, p-value < 2.2e-16, not normal
shapiro.test(d$PhenRich)
#W = 0.94445, p-value = 2.358e-10, not normal


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
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)
SkinD <- left_join(SkinD, Tree.sum, by="orchard.num")
SkinD <- left_join(SkinD, Orchard, by="orchard.num")

shapiro.test(SkinD$TotalPhen)
#W = 0.84362, p-value = 0.001668, not normal 
hist(SkinD$TotalPhen)
shapiro.test(SkinD$PhenRich)#normal!

#Pulp to orchard 
PulpD <- PulpD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)
PulpD <- left_join(PulpD, Tree.sum, by="orchard.num")
PulpD <- left_join(PulpD, Orchard, by="orchard.num")
shapiro.test(PulpD$TotalPhen)
#W = 0.75627, p-value = 6.158e-05, not normal 
shapiro.test(PulpD$PhenRich)#normal!


#Seed to orchard 
SeedD <- SeedD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)
SeedD <- left_join(SeedD, Tree.sum, by="orchard.num")
SeedD <- left_join(SeedD, Orchard, by="orchard.num")
shapiro.test(SeedD$TotalPhen)
#W = 0.8686, p-value = 0.004939, so close but not normal 
shapiro.test(SeedD$PhenRich)#normal!

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
#orchard.type           0.5016  1     0.4788    
#Latitude              18.2814  1  1.906e-05 ***
#orchard.type:Latitude  0.2196  1     0.6393  

## Calculate the effect size
ems <- emmeans(Lat1, c("Latitude","orchard.type"), infer = c(T, T))
eff_size(ems, sigma = sigma(Lat1), edf = df.residual(Lat1))
plot(ems)
#high level of uncertainty 


#Firmness 
Lat2<- glmmTMB(Firmness ~ orchard.type*Latitude + (1|site.code/orchard.num), data=TreeLat)
summary(Lat2)
Anova(Lat2)
#orchard.type           0.0828  1     0.7735    
#Latitude              34.5700  1  4.112e-09 ***
#orchard.type:Latitude  0.9233  1     0.3366 


## Calculate the effect size
ems <- emmeans(Lat2, c("Latitude","orchard.type"), infer = c(T, T))
eff_size(ems, sigma = sigma(Lat2), edf = df.residual(Lat2))
plot(ems)

#Average Weight 
Lat3<- glmmTMB(avgwgt ~ orchard.type*Latitude + (1|site.code/orchard.num), data=TreeLat)
summary(Lat3)
#orchard.typeOrganic           142.408     45.320   3.142  0.00168 **
Anova(Lat3) 
#orchard.type          0.1535  1   0.695206   
#Latitude              0.3713  1   0.542288   
#orchard.type:Latitude 9.7210  1   0.001822 **


#Maturity Index 
Lat4<- glmmTMB(maturity.index ~ orchard.type*Latitude + (1|site.code/orchard.num), 
               data=TreeLat)
summary(Lat4)
Anova(Lat4)  
#orchard.type          0.0047  1    0.94544  
#Latitude              3.8790  1    0.04889 *
#orchard.type:Latitude 0.9598  1    0.32724 


###Investigating Latitude and Average weight 
#Splitting latitude in high and low 
Tree_low <- filter(TreeLat, Latitude<42)
Tree_high <- filter(TreeLat, Latitude>42)

#Low
avg_low<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/orchard.num), data=Tree_low)
summary(avg_low)
Anova(avg_low) 
#orchard.type 2.5824  1     0.1081

avg_hgh<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/orchard.num), data=Tree_high)
summary(avg_hgh)
Anova(avg_hgh)
#orchard.type 3.6794  1    0.05509 .


#Q1-B: Fruit Chemistry-------------------------------------------------------------
#binding latitude to the 359 obs data set 
ChemLat <- left_join(d, Orchard[,c(2,4,40)], by="orchard.num")

#Total Phenolics with beta distribution 
tp1 <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Tissue*Latitude + 
                 (1|site.code/orchard.num/Tree), data=ChemLat, family=beta_family(link="logit"))
summary(tp1)
Anova(tp1)
#orchard.type                   4.4968  1    0.03396 *  
#Tissue                       278.1257  2    < 2e-16 ***
#orchard.type:Tissue            8.3481  2    0.01539 *  

#strong effects of tissue and interaction between tissue and orchard type
#splitting by tissue
d.sk <- left_join(d.sk, Orchard[,c(2,4)], by="orchard.num")
d.pu <- left_join(d.pu, Orchard[,c(2,4)], by="orchard.num")
d.se <- left_join(d.se, Orchard[,c(2,4)], by="orchard.num")

tp.sk <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Latitude + (1|site.code/orchard.num), 
                 data=d.sk, family=beta_family(link="logit"))
summary(tp.sk)
Anova(tp.sk)
#nothing 

tp.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Latitude + (1|site.code/orchard.num), 
                 data=d.pu, family=beta_family(link="logit"))
summary(tp.pu)
Anova(tp.pu)
#nothing 


tp.se <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Latitude + (1|site.code/orchard.num/Tree), 
                 data=d.se, family=beta_family(link="logit"))
summary(tp.se)
Anova(tp.se)
#orchard.type          8.0185  1    0.00463 **

#visualize this 
ggplot(d.se, aes(x=Latitude, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher in conventional 


#phenolics richness
pr1 <- glmmTMB(PhenRich~ orchard.type*Latitude*Tissue + (1|site.code/orchard.num/Tree), data=ChemLat, 
               family=poisson(link="log"))
summary(pr1)
Anova(pr1)
#Latitude:Tissue                7.6070  2    0.02229 *  
#Tissue                       798.6850  2    < 2e-16 ***

#per tissue type 
pr.p <- glmmTMB(PhenRich ~ orchard.type*Latitude+(1|site.code/orchard.num), data=d.pu,
                family=poisson(link="log"))
summary(pr.p)
Anova(pr.p)
#nothing 

pr.sk <- glmmTMB(PhenRich ~ orchard.type*Latitude+(1|site.code/orchard.num), data=d.sk,
                 family=poisson(link="log"))
summary(pr.sk)
Anova(pr.sk)
#nothing

pr.se <- glmmTMB(PhenRich ~ orchard.type*Latitude+(1|site.code/orchard.num), data=d.se,
                 family=poisson(link="log"))
summary(pr.se)
Anova(pr.se)
#Latitude              6.0320  1    0.01405 *


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
#Dimensions: 2 
#Stress:     0.06984612 
#Stress type 1, weak ties

m.NMDS.sk <- metaMDS(d.comp.sk, distance = "jaccard", trymax=100, autotransform =FALSE)
m.NMDS.sk
#Dimensions: 2 
#Stress:     0.07062257 
#Stress type 1, weak ties

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
#                      Df SumOfSqs      R2      F Pr(>F)
#orchard.type           1   0.0331 0.00322 0.3847  0.855
#lat_cat                1   0.1297 0.01260 1.5065  0.164
#orchard.type:lat_cat   1   0.1432 0.01391 1.6624  0.134
#Residual             116   9.9901 0.97028              
#Total                119  10.2962 1.00000              


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
#                      Df SumOfSqs      R2      F Pr(>F)  
#orchard.type           1   0.6803 0.02462 3.0582  0.040 *
#lat_cat                1   0.3623 0.01311 1.6286  0.172  
#orchard.type:lat_cat   1   0.7898 0.02858 3.5504  0.015 *
#Residual             116  25.8034 0.93370                
#Total                119  27.6357 1.00000  


###NMDS for Seed 
m.NMDS.se <- metaMDS(d.comp.se, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.se

plot(m.NMDS.se, type="t")

#for plotting, need to add columns that give values for 
#colors and symbols we want in plot
d.expl.se$Color <- recode_factor(d.expl.se$lat_cat,
                                 low="red", high="blue")
d.expl.se$Symbol <- recode_factor(d.expl.se$orchard.type,
                                  Organic=2, Conventional=1)
d.expl.se$Symbol <- as.numeric(as.character(d.expl.se$Symbol))

d.expl.se$lat_cat = as.factor(d.expl.se$lat_cat)


NMDSse=plot(m.NMDS.se, type="n") #plots the ordination axes only
points(m.NMDS.se, pch=d.expl.pu$Symbol,
       col=as.character(d.expl.pu$Color), cex = 0.8)     
ordiellipse(m.NMDS.se, d.expl.pu$Symbol, conf = 0.95)
NMDSse

#PERMANOVA
m.perm3 <- adonis2(d.comp.se~orchard.type*lat_cat, data=d.expl.se)
m.perm3
#                     Df SumOfSqs      R2      F Pr(>F)  
#orchard.type           1   0.6419 0.03179 3.9284  0.017 *
#lat_cat                1   0.5534 0.02741 3.3871  0.030 *
#orchard.type:lat_cat   1   0.2044 0.01012 1.2508  0.234  
#Residual             115  18.7902 0.93067                
#Total                118  20.1899 1.00000 


#Q1-C: Random Forest------------------------------------------------------------
###PULP###

m1.rf.pu <- randomForest(d.comp.pu,d.expl.pu$orchard.type, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.pu$importance
varImpPlot(m1.rf.pu)
MDSplot(m1.rf.pu, d.expl.pu$orchard.type)

m1.rf.b.pu <- Boruta(d.comp.pu,d.expl.pu$orchard.type)
m1.rf.b.pu
plot(m1.rf.b.pu,las = 2, cex.axis = 0.7)  
# 2 attributes confirmed important: G

getSelectedAttributes(m1.rf.b.pu) 

##Running MANOVAS (needs to be fixed)
d.comp.pu.sel <- data.matrix(d.comp.pu[,getSelectedAttributes(m1.rf.b.pu)])
m1.man.pu <-manova(d.comp.pu.sel ~ d.expl.pu$orchard.type)
summary(m1.man.pu)  #overall significance for MANOVA
#d.expl.pu$orchard.type   1 0.038823   2.3629      2    117 0.09863 .

summary.aov(m1.man.pu)  #follow-up ANOVAs for each individual compound
# g= 0.0353 *

#some quick plots of all of them
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.pu.sel))){
  d.temp=d.comp.pu.sel[,i]
  plot(d.temp ~ d.expl.pu$orchard.type, ylab=colnames(d.comp.pu.sel)[i])
}
dev.off()

#G is higher in conventional


###SEEDS###
m1.rf.se <- randomForest(d.comp.se,d.expl.se$orchard.type, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.se$importance
varImpPlot(m1.rf.se)
MDSplot(m1.rf.se, d.expl.se$orchard.type)


m1.rf.b.se <- Boruta(d.comp.se,d.expl.se$orchard.type)
m1.rf.b.se
plot(m1.rf.b.se,las = 2, cex.axis = 0.7)  #important variables (better than shadow) are in green
getSelectedAttributes(m1.rf.b.se) #lists all important ones
#[1] "E"           "PB2"         "ePicatechin" 

#Running MANOVAS 
d.comp.se.sel <- data.matrix(d.comp.se[,getSelectedAttributes(m1.rf.b.se)])
m1.man.se <- manova(d.comp.se.sel ~ d.expl.se$orchard.type)
summary(m1.man.se)  
#d.expl.se$orchard.type   1  0.136    6.034      3    115 0.000747 ***


summary.aov(m1.man.se)  
#E= 0.003836 **
#PB2= 0.001805 **
#epicatechin = 0.0006019 ***


par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.se.sel))){
  d.temp=d.comp.se.sel[,i]
  plot(d.temp ~ d.expl.se$orchard.type, ylab=colnames(d.comp.se.sel)[i])
}
dev.off()
#o:E
#PB2, EPI, C 

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
#12 important attributes 

#plot 
plot(m1.rf.b.sk,las = 2, cex.axis = 0.7)  
getSelectedAttributes(m1.rf.b.sk) #lists all important ones

# [1] "A"                "F"                "G"               
#[4] "catechin"         "chlorogenic_acid" "PB2"             
#[7] "ePicatechin"      "U1"               "U2"              
#[10] "I"                "K"                "phloridzin"  

#performing MANOVA 
d.comp.sk.sel <- data.matrix(d.comp.sk[,getSelectedAttributes(m1.rf.b.sk)])
m1.man.sk <- manova(d.comp.sk.sel ~ d.expl.sk$lat_cat)
summary(m1.man.sk) 
#d.expl.sk$lat_cat   1 0.48521   8.4043     12    107 5.317e-11 ***


#follow-up ANOVAs for each individual compound
summary.aov(m1.man.sk)  
#A= 0.003673 
#F= 0.002651 
#G= 0.004627 
#catechin= 0.0001634 
#pb2= 0.00372 
#ePicatechin= 0.01994 
#u2= 0.001417 
#i= 0.01639 
#k= 0.000237 


#some quick plots of all of them
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.sk.sel))){
  d.temp=d.comp.sk.sel[,i]
  plot(d.temp ~ d.expl.sk$lat_cat, ylab=colnames(d.comp.sk.sel)[i])
}
dev.off()

#high: F, G, PB2, EPI, U2, I, K
#LOW: A, 


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
#[1] "E"            "F"            "G"            "H"           
#[5] "caffeic_acid" "phloridzin" 

##Running MANOVAS
d.comp.pu.sel <- data.matrix(d.comp.pu[,getSelectedAttributes(m1.rf.b.pu)])
m1.man.pu <-manova(d.comp.pu.sel ~ d.expl.pu$lat_cat)
summary(m1.man.pu)  #overall significance for MANOVA
#d.expl.pu$lat_cat   1 0.17619   4.0279      6    113 0.001077 **

summary.aov(m1.man.pu)  #follow-up ANOVAs for each individual compound
#E= 0.08855 
#F= 0.001154 
#G= 1.867e-05
#H= 0.02089 
#phloridzin= 0.0375 

#some quick plots of all of them
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.pu.sel))){
  d.temp=d.comp.pu.sel[,i]
  plot(d.temp ~ d.expl.pu$lat_cat, ylab=colnames(d.comp.pu.sel)[i])
}
dev.off()

#HIGH:all 

###SEEDS###
m1.rf.se <- randomForest(d.comp.se,d.expl.se$lat_cat, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.se$importance
varImpPlot(m1.rf.se)
MDSplot(m1.rf.se, d.expl.se$lat_cat)


m1.rf.b.se <- Boruta(d.comp.se,d.expl.se$lat_cat)
m1.rf.b.se
plot(m1.rf.b.se, las = 2, cex.axis = 0.7)  #important variables (better than shadow) are in green
getSelectedAttributes(m1.rf.b.se) #lists all important ones
#[1] "PB2"         "ePicatechin" "U2"          "I"           "J"          
#[6] "reynoutrin"  "quercetin" 

#Running MANOVAS 
d.comp.se.sel <- data.matrix(d.comp.se[,getSelectedAttributes(m1.rf.b.se)])
m1.man.se <- manova(d.comp.se.sel ~ d.expl.se$lat_cat)
summary(m1.man.se)  
#d.expl.se$lat_cat   1 0.22307   5.3595      6    112 6.753e-05 ***

summary.aov(m1.man.se)  
#PB2= 8.717e-05
#epicatechin = 1.703e-05
#U2 4.526e-06
#I= 3.415e-06
#J= 3.426e-05
#reynoutrin= 1.585e-06


par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.se.sel))){
  d.temp=d.comp.se.sel[,i]
  plot(d.temp ~ d.expl.se$lat_cat, ylab=colnames(d.comp.se.sel)[i])
}
dev.off()
#HIGH: PB2, EPI, U2, I, J, REY
#LOW: NONE 




#Q2: Which abiotic factors are the most important drivers of fruit quality?---------
#Analysis: principle components analysis followed by PC regression with each variable
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
#  PC1     PC2     PC3     PC4     PC5     PC6     PC7 
#0.61114 0.15131 0.12861 0.07630 0.02307 0.00957 0.00000 

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

#PC's 5, 6, 7, 8, 9 are very low/ not worth looking at


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
##Linear models PC x quality 
pc_clim_phys <- cbind(pc_clim, c)


###Firmness###
p1 <- glmmTMB(Firmness ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p1)
#orchard.typeOrganic  0.04264    1.27173   0.034 0.973252    
#PC1                  1.60404    0.32095   4.998  5.8e-07 ***
#PC2                  1.56234    0.64399   2.426 0.015264 *  
#PC3                  2.36819    0.70735   3.348 0.000814 ***
#PC4                 -1.70559    0.89457  -1.907 0.056573 . 


plot(Firmness ~ PC1, data=pc_clim_phys)
plot(Firmness ~ PC2, data=pc_clim_phys)
plot(Firmness ~ PC3, data=pc_clim_phys)

results$rotation

###SSC###
P2 <- glmmTMB(SSC ~ orchard.type + PC1 + PC2 + PC3 + PC4 + (1|site.code), data= pc_clim_phys)
summary(P2)
#orchard.typeOrganic -0.08014    0.35712  -0.224  0.82244    
#PC1                 -1.11783    0.15404  -7.257 3.97e-13 ***
#PC2                  1.05801    0.33345   3.173  0.00151 ** 
#PC3                  0.28125    0.33223   0.847  0.39725    
#PC4                  0.72345    0.43259   1.672  0.09445 . 


plot(SSC ~ PC1, data=pc_clim_phys)
plot(SSC ~ PC2, data=pc_clim_phys)

results$rotation

###AVGWGT
p3 <- glmmTMB(avgwgt ~ orchard.type + PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p3)
#orchard.typeOrganic   0.1284     8.2181   0.016  0.98754    
#PC1                  -3.2341     2.5251  -1.281  0.20027    
#PC2                  15.4784     5.3450   2.896  0.00378 ** 
#PC3                   1.1532     5.4339   0.212  0.83194    
#PC4                 -10.1635     7.1234  -1.427  0.15364   

plot(avgwgt ~ PC2, data=pc_clim_phys)


results$rotation

###Maturity Index 
p4 <- glmmTMB(maturity.index ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p4)
#orchard.typeOrganic  0.018648   0.098212   0.190  0.84941    
#PC1                 -0.302249   0.094930  -3.184  0.00145 ** 
#PC2                  0.276450   0.192868   1.433  0.15175    
#PC3                 -0.007687   0.185102  -0.042  0.96687    
#PC4                 -0.374573   0.263331  -1.422  0.15490  

plot(maturity.index ~ PC1, data=pc_clim_phys)

results$rotation

#Q2-B: Fruit Chemistry----------------------------------------------------------
#SKIN#
pc.sk_clim <- cbind(pc_clim, SkinD)

#TotalPhen
pc.sk.tp <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.sk_clim, family=beta_family (link="logit"))
summary(pc.sk.tp)
#nothing 

#PhenRich
pc.sk.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.sk_clim)
summary(pc.sk.pr)
#nothing 

results$rotation

#PULP#
pc.pu_clim <- cbind(pc_clim, PulpD)

##TotalPhen###
pc.pu.tp <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim, family=beta_family (link="logit"))
summary(pc.pu.tp)

#PhenRich 
pc.pu.pr <- glmmTMB(PhenRich ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim)
summary(pc.pu.pr)


#SEED#
pc.se_clim <- cbind(pc_clim, SeedD)

#TotalPhen
pc.se.tp <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.se_clim, family=beta_family (link="logit"))
summary(pc.se.tp)

#PhenRich 
pc.se.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 +
                      (1|site.code), data=pc.se_clim)
summary(pc.se.pr)


#Q3: Which specific management practices are the most important drivers of fruit quality?----
#For this analysis, we'll be running GLMMs against every mgmt variable. 
#Q3-A: Physical Quality---------------------------------------------------------
#SSC
mgmt1 <- glmmTMB(SSC ~ orchard.type+ Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt1)
Anova(mgmt1)
#Herbicides   49.8749  1  1.639e-12 ***
  
ggplot(c, aes(x=Herbicides, y=SSC, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()

##higher ssc in organic orchards without herbicides 
##higher ssc in conventional orchards with herbicides 



#avgwgt
mgmt2 <- glmmTMB(avgwgt ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt2)
Anova(mgmt2)
#Acres        9.7082  1   0.001835 **
#Cultivation  4.0286  1   0.044734 * 
  
  
ggplot(c, aes(x=Acres, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

#steady increase in weigth as acregae increases 

ggplot(c, aes(x=Cultivation, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher weights in con without cult


#Firmness 
mgmt3 <- glmmTMB(Firmness ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt3)
Anova(mgmt3)
#Herbicides   4.8742  1    0.02726 *
  
ggplot(c, aes(x=Herbicides, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher firmness wihtout herbicides, higher in con 

#Maturity Index 
mgmt4 <- glmmTMB(maturity.index ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt4)
Anova(mgmt4)
#Herbicides   20.4531  1  6.111e-06 ***
#orchard.type 20.6444  1  5.530e-06 ***
#Acres         3.4818  1   0.062046 .  


ggplot(c, aes(x=Herbicides, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#organic orchards that used herbicides had higher maturity values 


ggplot(c, aes(x=Acres, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#lower maturity in smaller acreage for conventional but increases past organic at higher 


#Q3-B: Fruit Chemistry----------------------------------------------------------
#Total Phenolics
#Skin 
mgmt.sk <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                     Weed_Mats + Cover_Crops + Acres + (1|site.code), 
                   data=SkinD, family=beta_family(link="logit"))
summary(mgmt.sk)
Anova(mgmt.sk)
#Herbicides   14.5890  1  0.0001337 ***
#Weed_Mats     3.1533  1  0.0757759 .  
#orchard.type  5.3281  1  0.0209846 *  


ggplot(SkinD, aes(x=Herbicides, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher total phen in herbicides, higher in organic 
ggplot(SkinD, aes(x=Weed_Mats, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#only coventional used weed_mats and was higehr 


#Pulp 
mgmt.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                     Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                   data=PulpD, family=beta_family(link="logit"))
summary(mgmt.pu)
Anova(mgmt.pu)
#orchard.type 4.7614  1   0.029106 * 
#Cultivation  6.7080  1   0.009598 **
#Cover_Crops  7.1191  1   0.007627 **
  
ggplot(PulpD, aes(x=Cultivation, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher in conventional with and higher in organic without without 

ggplot(PulpD, aes(x=Cover_Crops, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher total phen when using cover crops 


#Seed
mgmt.se <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                     Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                   data=SeedD, family=beta_family(link="logit"))
summary(mgmt.se)
Anova(mgmt.se)
#Cover_Crops    8.7445  1   0.003105 ** 
#Acres          8.4539  1   0.003643 ** 
#orchard.type  57.0345  1  4.282e-14 ***
  
ggplot(SeedD, aes(x=Cover_Crops, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher in conventional, with and without 

ggplot(SeedD, aes(x=Acres, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#lowers in conventional, increases in organic 



#Phenolics Richness 
#Skin 
mgmt.pr.sk <- glmmTMB(PhenRich ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                        Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                      data=SkinD)
summary(mgmt.pr.sk)
Anova(mgmt.pr.sk)
#orchard.type 10.9965  1  0.0009128 ***
#Herbicides   20.6498  1  5.514e-06 ***
#Mowing        3.2551  1  0.0712025 .  
#Weed_Mats    12.2388  1  0.0004681 ***
#Cover_Crops   3.9363  1  0.0472553 *  
#Acres         6.5262  1  0.0106297 *  

ggplot(SkinD, aes(x=Herbicides, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher in organic that used, lower when not sued 


ggplot(SkinD, aes(x=Weed_Mats, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#no organic reported using weed mats, higehr in organic when not used 

ggplot(SkinD, aes(x=Cover_Crops, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher in organic not using cover crps

ggplot(SkinD, aes(x=Acres, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decreased as acreage increased 


#Pulp 
mgmt.pr.pu <- glmmTMB(PhenRich ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                        Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                      data=PulpD)
summary(mgmt.pr.pu)
Anova(mgmt.pr.pu)
#Cultivation   13.5405  1  0.0002335 ***
#Com_Mul        6.6978  1  0.0096531 ** 
  
ggplot(PulpD, aes(x=Cultivation, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher in roanic not using it, higher in conventoanl using it 


ggplot(PulpD, aes(x=Com_Mul, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higehr in conventional not using, highe rin organic using

#Seed
mgmt.pr.se <- glmmTMB(PhenRich ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                        Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                      data=SeedD)
summary(mgmt.pr.se)
Anova(mgmt.pr.se)
#orchard.type  6.3516  1   0.011728 *  
#Cultivation   8.5015  1   0.003548 ** 
#Herbicides    8.6565  1   0.003259 ** 
#Mowing        4.6069  1   0.031844 *  
#Weed_Mats     4.0724  1   0.043590 *  
#Cover_Crops   5.3724  1   0.020458 *  
#Acres         5.7112  1   0.016857 * 


ggplot(SeedD, aes(x=Cultivation, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher values when not used, higher in organic

ggplot(SeedD, aes(x=Herbicides, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher in conventional when not used, higher in organic when used

ggplot(SeedD, aes(x=Mowing, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher in conventional when used, not ysed not reported for conventional 

ggplot(SeedD, aes(x=Weed_Mats, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#higher in conventional when not used 

ggplot(SeedD, aes(x=Cover_Crops, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#lower in organic when used, higher in conventional when not sued 

ggplot(SeedD, aes(x=Acres, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#starts higher in coventional and decrease oast organic, both decrease 

#Multiple plot function----------------------------------------------------------------------------------------------
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#plots----------------------------------------------------------------
#Q1-A
#Physical traits by management level  
ag1= ggplot(TreeLat, aes(x=orchard.type, y=SSC, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Soluble Sugar Content (%)") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"))

ag2= ggplot(TreeLat, aes(x=orchard.type, y=Firmness, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Firmness (N)") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic()+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"))

ag3= ggplot(TreeLat, aes(x=orchard.type, y=avgwgt, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Average Weight (g)") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

ag4= ggplot(TreeLat, aes(x=orchard.type, y=maturity.index, color="#3EBCD2", "#9A607F")) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Cornell Starch-Iodine Value") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() 

multiplot(ag1,ag2,ag3,ag3, cols=2)
  
#Q1-B: Total Phenolics and Phenolic Richness by Managment System 
sk.tp.ag=ggplot(d.sk, aes(x=orchard.type, y=TotalPhen/1000, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Skin Total Phenolics (ug/g)") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

sk.pr.ag=ggplot(d.sk, aes(x=orchard.type, y=PhenRich, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Skin Phenolic Richness") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
#pulp
pu.tp.ag=ggplot(d.pu, aes(x=orchard.type, y=TotalPhen/1000, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Pulp Total Phenolics (ug/g)") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

pu.pr.ag=ggplot(d.pu, aes(x=orchard.type, y=PhenRich, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Pulp Phenolic Richness") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
#seed
se.tp.ag=ggplot(d.se, aes(x=orchard.type, y=TotalPhen/1000, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Seed Total Phenolics (ug/g)") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

se.pr.ag=ggplot(d.se, aes(x=orchard.type, y=PhenRich, color=orchard.type)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Seed Phenolic Richness") +
  xlab ("Management System")+
  geom_smooth(method=glm, se=FALSE)+
  theme_classic() +  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

multiplot(sk.tp.ag, sk.pr.ag, pu.tp.ag, pu.pr.ag, se.tp.ag, se.pr.ag, cols=3)

#Q1-A: Physical Traits over Latitude 
la1= ggplot(TreeLat, aes(x=Latitude, y=SSC, color=orchard.type)) +
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Soluble Sugar Content") +
  xlab ("Latitude")+
  geom_smooth(method=lm)+
  theme_bw()+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

la2= ggplot(TreeLat, aes(x=Latitude, y=Firmness, color=orchard.type)) +
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Firmness") +
  xlab ("Latitude")+
  geom_smooth(method=lm)+
  theme_bw()+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

la3= ggplot(TreeLat, aes(x=Latitude, y=avgwgt, color=orchard.type)) +
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Average Weight") +
  xlab ("Latitude")+
  geom_smooth(method=lm)+
  theme_bw()+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

la4= ggplot(TreeLat, aes(x=Latitude, y=maturity.index, color=orchard.type)) +
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Maturity Index") +
  xlab ("Latitude")+
  geom_smooth(method=lm)+
  theme_bw()+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")


multiplot(la1,la2,la3,la4, cols=2)


#Q1-A: Average Weight split by Latitudinal Categories 
avglow <- ggplot(Tree_low, aes(x=orchard.type, y=avgwgt, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  xlab ("Latitude < 42 ") +
  ylab ("Average Weight") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("Conventional", "Organic"))
avglow

avghigh <- ggplot(Tree_high, aes(x=orchard.type, y=avgwgt, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  xlab ("Latitude > 42 ") +
  ylab ("Average Weight") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("Conventional", "Organic"))
avghigh

multiplot(avglow, avghigh, cols=2)


###Q1-B Phenolic Richness by Tissue Type over Latitude 
ggplot(ChemLat, aes(x=Latitude, y=TotalPhen/1000, color=Tissue, shape=orchard.type)) +
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Total Phenolics ug/g") +
  xlab ("Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_bw()+
  scale_color_manual(values=c("#3EBCD2", "#9A607F", "darkgreen"),name="Tissue")+
  scale_shape_manual(values=c(16, 17), name="Management System")+
guides(shape = guide_legend(override.aes = list(shape = c(16, 17))))







##Q1-B Phenolic Richness by Tissue Type over Latitude 
ggplot(ChemLat, aes(x=Latitude, y=PhenRich, color=Tissue), shape=orchard.type) +
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Phenolic Richness") +
  xlab ("Latitude")+
  geom_smooth(method=glm ,alpha = .15,aes(fill = NULL))+
  theme_bw()+
  scale_color_manual(values=c("#3EBCD2", "#9A607F", "darkgreen"),name="Tissue")

#Map Plot----------------------------------------------------------------------
install.packages("usmap")
library(usmap)
library(ggplot2)

# Your data transformation code
c1 <- data.frame(
  Orchard = c$orchard.num,
  Latitude = c$Latitude,
  Longitude = c$Longitude,
  Otype = c$orchard.type
)


c2 <- usmap_transform(
  c1,
  input_names = c("Longitude", "Latitude"),
  output_names = c("x", "y")
)

c2 <- c2 %>%
  mutate(Orchard = ifelse(duplicated(select(c2, x)), 
"Both", as.character(Otype)))


# Create the US map
p1 <- plot_usmap(include = c("CA", "OR", "WA")) +
  geom_point(data = c2, aes(x = x, y = y, color = Otype), 
             size = 3, 
alpha = 0.7, position = position_jitter(width = 10, height = 10)) +
  labs(
    title = "          Participating Orchards",
    caption = "                                 August 2023"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(hjust = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    panel.background = element_rect(fill = "white"),  # Set background color
    panel.grid = element_blank(),                      # Remove grid lines
    panel.grid.major = element_blank()                  # Remove major grid lines
  ) +
  scale_color_manual(values = c("darkgreen", "#9A607F", "#3EBCD2"), name = "Management System")+
  scale_x_continuous(labels = c("-123°W", "-120°W", "-117°W", "-114°W", "-111°W"), expand = c(0, 0)) +  # Set breaks and labels for x-axis
  scale_y_continuous(labels = c("34°N", "38°N", "42°N", "46°N"), expand = c(0, 0))    # Set breaks and labels for y-axis
p1

ggsave("output_plot.png", plot = p1, device = "png", bg = "transparent")


