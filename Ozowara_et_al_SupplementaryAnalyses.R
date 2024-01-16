#Chapter 1 Analysis 
setwd("C:\\Users\\xozow\\OneDrive\\Documents\\Dissertation_Data\\Ozowara-Chapter1")
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
#Read in Data sheets#

#Orchard Level Data (24 obs)
Orchard <- read_csv("Orchard Level Data.csv")
View(Orchard)
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

#Q1-C: Which compounds distinguish fruits based on management systems?----------
#Analysis: NMDS and random forest

#left joining the "expl" by latitude  
d.expl <- left_join(d.expl, Orchard[,c(2,40)], by="orchard.num")
d.expl.sk <- left_join(d.expl.sk, Orchard[,c(2,40)], by="orchard.num")
d.expl.pu <- left_join(d.expl.pu, Orchard[,c(2,40)], by="orchard.num")
d.expl.se <- left_join(d.expl.se, Orchard[,c(2,40)], by="orchard.num")

###NMD for Skin 
m.NMDS.sk <- metaMDS(d.comp.sk, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.sk

plot(m.NMDS.sk, type="t")


# Extract NMDS coordinates
skin_NMDS <- vegan::scores(m.NMDS.sk,display="sites")
skin_NMDS1 <- vegan::scores(m.NMDS.sk,display="species")

# Create a data frame for ggplot2
skin_NMDS.df1 <- as.data.frame(skin_NMDS)
skin_NMDS.df2 <- as.data.frame(skin_NMDS1)

# Create NMDS plot with ggplot2
ggplot(skin_NMDS.df1, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = sites)) +
  theme_minimal() +
  ggtitle("NMDS Plot")






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
 

###NMD for Pulp 
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

#permanova
m.perm2 <- adonis2(d.comp.pu~orchard.type*lat_cat, data=d.expl.pu)
m.perm2
#                      Df SumOfSqs      R2      F Pr(>F)  
#orchard.type           1   0.6803 0.02462 3.0582  0.040 *
#lat_cat                1   0.3623 0.01311 1.6286  0.172  
#orchard.type:lat_cat   1   0.7898 0.02858 3.5504  0.015 *
#Residual             116  25.8034 0.93370                
#Total                119  27.6357 1.00000  


###NMD for Seed 
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


m.perm3 <- adonis2(d.comp.se~orchard.type*lat_cat, data=d.expl.se)
m.perm3
#                     Df SumOfSqs      R2      F Pr(>F)  
#orchard.type           1   0.6419 0.03179 3.9284  0.017 *
#lat_cat                1   0.5534 0.02741 3.3871  0.030 *
#orchard.type:lat_cat   1   0.2044 0.01012 1.2508  0.234  
#Residual             115  18.7902 0.93067                
#Total                118  20.1899 1.00000 


###Random Forest for Orchard Type###
###skin
m1.rf.sk <- randomForest(d.comp.sk,d.expl.sk$orchard.type, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.sk$importance
varImpPlot(m1.rf.sk)
MDSplot(m1.rf.sk, d.expl.sk$orchard.type)

#using boruta 
m1.rf.b.sk <- Boruta(d.comp.sk,d.expl.sk$orchard.type)
m1.rf.b.sk
# 9 attributes confirmed important: A, cyanidin_Galactoside, ePicatechin, F, G and 4 more;

#plot 
plot(m1.rf.b.sk,las = 2, cex.axis = 0.7)  
getSelectedAttributes(m1.rf.b.sk) #lists all important ones

#[1] "A"                    "F"                    "G"                    "PB1"                 
#[5] "cyanidin_Galactoside" "ePicatechin"          "U1"                   "I"                   
#[9] "K" 

#performing MANOVA 
d.comp.sk.sel <- data.matrix(d.comp.sk[,getSelectedAttributes(m1.rf.b.sk)])
m1.man.sk <- manova(d.comp.sk.sel ~ d.expl.sk$orchard.type)
summary(m1.man.sk) 
#d.expl.sk$orchard.type   1 0.23711   3.7986      9    110 0.0003348 ***

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.sk)  
#A= 0.09424 .
#PB1 =0.05021 .
#gala = 0.08316 .

#some quick plots of all of them
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.sk.sel))){
  d.temp=d.comp.sk.sel[,i]
  plot(d.temp ~ d.expl.sk$orchard.type, ylab=colnames(d.comp.sk.sel)[i])
}
dev.off()

#all three marginally higher in organic 

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



#Q4: Which pest or diseases presence has the most significant affect on fruit quality-----
#visualize pest data in indices 

p <- pivot_longer(data=c, cols=21:38, names_to="pest", values_to="index")

ggplot(p, aes(x=pest, y=index))+
  geom_boxplot()+
  geom_jitter(width=0.2, height=0.1)

#removing everything that doesn't pass 3 on index 
#what were left with: 
#pests: aphids, apple maggots, coddling moth, 
#disease: fireblight, powdery mildew, 

#correlation matrix
p_pest <- dplyr::select(c,c(Aphids, Apple.Maggots, Codling.Moth, 
                            Powdery.mildew, Fire.Blight))

rcorr(as.matrix(p_pest))
corrplot(cor(p_pest))

testRes = cor.mtest(p_pest, conf.level = 0.95)

#add all p-values
corrplot(cor(p_pest), p.mat = testRes$p, insig = 'p-value', sig.level = -1)

#add significant level stars
corrplot(cor(p_pest), p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')


#strong relationships between: 
#aphinds and apple maggots
#apple maggots and coddling moth 

#Q4-A: Physical Quality---------------------------------------------------------
#Examining Relationships of Pest.Index, orchard.type, and latitude 
#pest index 
#ssc
ssc.index <- glmmTMB(SSC ~ orchard.type*Pest.Index+ (1|site.code), 
                     data=c)
summary(ssc.index)
Anova(ssc.index)

#firmness 
firm.index <- glmmTMB(Firmness ~ orchard.type*Pest.Index+ (1|site.code), 
                      data=c)
summary(firm.index)
Anova(firm.index)
#nothing 

#avg wgt 
wgt.index <- glmmTMB(avgwgt ~ orchard.type*Pest.Index+ (1|site.code), 
                     data=c)
summary(wgt.index)
Anova(wgt.index)

ggplot(c, aes(x=Pest.Index, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

#maturity  
mat.index <- glmmTMB(maturity.index ~ orchard.type*Pest.Index+ (1|site.code), 
                     data=c)
summary(mat.index)
Anova(mat.index)
#Pest.Index              11.4270  1  0.0007238 ***

ggplot(c, aes(x=Pest.Index, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


###responses to species###
###SSC###
ssc_pest <- glmmTMB(SSC ~ orchard.type + Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Fire.Blight+ 
                      (1|site.code), 
                    data=c)
summary(ssc_pest)
Anova(ssc_pest)
#Fire.Blight    5.3607  1    0.02060 *
#orchard.type   5.6523  1    0.01743 * (organic)


ggplot(c, aes(x=Fire.Blight, y=SSC, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#lower SSC with lower reported fireblight

###avgwgt###
wgt_pest <- glmmTMB(avgwgt ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Fire.Blight+ (1|site.code), 
                    data=c)
summary(wgt_pest)
Anova(wgt_pest)
#Apple.Maggots  4.8605  1    0.02748 *
#Powdery.mildew 3.8439  1    0.04993 *
#Fire.Blight    5.6336  1    0.01762 *


ggplot(c, aes(x=Apple.Maggots, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#weight dcreased as apple maggots increase 

ggplot(c, aes(x=Powdery.mildew, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avgwgt decrease as pow mil increases 

ggplot(c, aes(x=Fire.Blight, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avgwgt decrease as fireblight increases 

ggplot(c, aes(x=Codling.Moth, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avgwgt higher in orchards without codling moth, organic weight higher with high pressure 

###firmness###
frm_pest <- glmmTMB(Firmness ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Fire.Blight+  (1|site.code), 
                    data=c)
summary(frm_pest)
Anova(frm_pest)
#Aphids         4.0144  1    0.04511 *
#Apple.Maggots  3.5262  1    0.06040 .
#Powdery.mildew 3.3710  1    0.06636 .

ggplot(c, aes(x=Aphids, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#frimness inreased as aphids increased in conventional orchards 

ggplot(c, aes(x=Apple.Maggots, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increase in conventional, as apple maggots increased
#decrease in organic 

ggplot(c, aes(x=Powdery.mildew, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decrease as pow mil increase

###maturity###
mat_pest <- glmmTMB(maturity.index ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Fire.Blight+ (1|site.code), 
                    data=c)
summary(mat_pest)
Anova(mat_pest)
#Powdery.mildew 5.1418  1    0.02336 *

ggplot(c, aes(x=Powdery.mildew, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decrease in maturity as pow mil increases 


#Q4-B: Fruit Chemistry ---------------------------------------------------------
#Total Phenolics
#Skin 
pest.index.tp.sk <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Pest.Index +(1|site.code), 
                            data=SkinD, family=beta_family(link="logit"))
summary(pest.index.tp.sk)
Anova(pest.index.tp.sk)
#nothing 

pest.tp.sk <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+ (1|site.code), 
                      data=SkinD, family=beta_family(link="logit"))
summary(pest.tp.sk)
Anova(pest.tp.sk)
#Powdery.mildew 5.0525  1    0.02459 *
#Fire.Blight    3.2125  1    0.07308 .

ggplot(SkinD, aes(x=Powdery.mildew, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher total phen in orchards with higher pow mol 

ggplot(SkinD, aes(x=Fire.Blight, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher total phen in orchards with higher fire.blight


#Pulp
pest.index.tp.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Pest.Index +(1|site.code), 
                            data=PulpD, family=beta_family(link="logit"))
summary(pest.index.tp.pu)
Anova(pest.index.tp.pu)
#nothing 

pest.tp.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+ (1|site.code), 
                      data=PulpD, family=beta_family(link="logit"))
summary(pest.tp.pu)
Anova(pest.tp.pu)
#orchard.type   5.1942  1    0.02266 * 
#Aphids         4.6777  1    0.03056 * 
#Apple.Maggots  4.0006  1    0.04548 * 
#Powdery.mildew 7.5087  1    0.00614 **
#Fire.Blight    6.5499  1    0.01049 * 

ggplot(PulpD, aes(x=Aphids, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher total phen in conventional , steady decrease 

ggplot(PulpD, aes(x=Apple.Maggots, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#total phen decrease as applemag increase, steadily higher in conventional 

ggplot(PulpD, aes(x=Powdery.mildew, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increases

ggplot(PulpD, aes(x=Fire.Blight, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increases, higher in conventional 

#Seed 
pest.index.tp.se <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Pest.Index +(1|site.code), 
                            data=SeedD, family=beta_family(link="logit"))
summary(pest.index.tp.se)
Anova(pest.index.tp.se)


pest.tp.se <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+ (1|site.code), 
                      data=SeedD, family=beta_family(link="logit"))
summary(pest.tp.se)
Anova(pest.tp.se)
#orchard.type   13.5904  1  0.0002273 ***


##PhenRich##
#skin 
pest.index.pr.sk <- glmmTMB(PhenRich~ orchard.type*Pest.Index +(1|site.code), 
                            data=SkinD)
summary(pest.index.pr.sk)
Anova(pest.index.pr.sk)


pest.pr.sk <- glmmTMB(PhenRich ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+ (1|site.code), 
                      data=SkinD)
summary(pest.pr.sk)
Anova(pest.pr.sk)
#Powdery.mildew 14.2281  1  0.0001619 ***

ggplot(SkinD, aes(x=Powdery.mildew, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increase as powmil increases 

#pulp
pest.index.pr.pu <- glmmTMB(PhenRich~ orchard.type*Pest.Index +(1|site.code), 
                            data=PulpD)
summary(pest.index.pr.pu)
Anova(pest.index.pr.pu)


pest.pr.pu <- glmmTMB(PhenRich ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+(1|site.code), 
                      data=PulpD)
summary(pest.pr.pu)
Anova(pest.pr.pu)
#Aphids          7.0735  1   0.007823 ** 
#Apple.Maggots   8.1773  1   0.004242 ** 
#Powdery.mildew 19.5452  1  9.825e-06 ***

ggplot(PulpD, aes(x=Apple.Maggots, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increase in conventional, decrease in organic 

ggplot(PulpD, aes(x=Aphids, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#o starts low eend shigh, conventional starts high and ends low 

ggplot(PulpD, aes(x=Powdery.mildew, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#sharp increase in o, starts low ends high, both increae 


#Seed
pest.index.pr.se <- glmmTMB(PhenRich~ orchard.type*Pest.Index +(1|site.code), 
                            data=SeedD)
summary(pest.index.pr.se)
Anova(pest.index.pr.se)


pest.pr.se <- glmmTMB(PhenRich ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+(1|site.code), 
                      data=SeedD)
summary(pest.pr.se)
Anova(pest.pr.se)
#orchard.type   13.7879  1  0.0002046 ***
#Aphids          4.4591  1  0.0347154 *  
#Codling.Moth   19.2622  1  1.139e-05 ***


ggplot(SeedD, aes(x=Aphids, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#c starts high ends low, o reverse 

ggplot(SeedD, aes(x=Codling.Moth, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decrease c,increase o



#Q5: How does fruit quality compare to total phenolics and phenolic richness--------
#skin 
sk.tp.qual <- glmmTMB((TotalPhen/1000000)+0.0001~ orchard.type+SSC+Firmness+avgwgt+maturity.index + (1|site.code)
                      , data=SkinD,family=beta_family(link="logit"))
summary(sk.tp.qual)
Anova(sk.tp.qual)
#nothing 

sk.pr.qual <- glmmTMB(PhenRich~ orchard.type+SSC+Firmness+avgwgt+maturity.index + (1|site.code)
                      ,data=SkinD)
summary(sk.pr.qual)
Anova(sk.pr.qual)
#nothing 


#pulp 
pu.tp.qual <- glmmTMB((TotalPhen/1000000)+0.0001~ orchard.type+SSC+Firmness+avgwgt+maturity.index + (1|site.code)
                      , data=PulpD,family=beta_family(link="logit"))
summary(pu.tp.qual)
Anova(pu.tp.qual)
#nothing 

pu.pr.qual <- glmmTMB(PhenRich~ orchard.type+SSC+Firmness+avgwgt+maturity.index + (1|site.code)
                      ,data=PulpD)
summary(pu.pr.qual)
Anova(pu.pr.qual)
#SSC            3.7418  1    0.05307 .
#avgwgt         4.7708  1    0.02895 *
#maturity.index 2.9161  1    0.08770 .

ggplot(PulpD, aes(x=SSC, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#SSC increase as phen rich decreases 

ggplot(PulpD, aes(x=avgwgt, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avg wgt increase as phen rich increasses 


ggplot(PulpD, aes(x=maturity.index, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#maturity increase as phen rich decreases 


#pulp avgwgt (figure is not correct)
ggplot(PulpD, aes(x=avgwgt, y=Latitude, color=orchard.type)) + 
  geom_smooth(method = "lm") +
  geom_point() + 
  scale_y_continuous(
    "Latitude", 
    sec.axis = sec_axis(~ . * .25, name = "Phenolic Richness")
  )+
  theme_classic()

#seed
se.tp.qual <- glmmTMB((TotalPhen/1000000)+0.0001~ orchard.type+SSC+Firmness+avgwgt+maturity.index + (1|site.code)
                      , data=SeedD,family=beta_family(link="logit"))
summary(se.tp.qual)
Anova(se.tp.qual)
#nothing 

se.pr.qual <- glmmTMB(PhenRich~ orchard.type+SSC+Firmness+avgwgt+maturity.index + (1|site.code)
                      ,data=SeedD)
summary(se.pr.qual)
Anova(se.pr.qual)
#SSC            10.1002  1   0.001483 **
#avgwgt          5.0848  1   0.024136 * 


ggplot(SeedD, aes(x=SSC, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#SSC increase as phen rich decrease 

ggplot(SeedD, aes(x=avgwgt, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#AVG wgt increase as phen rich decreases 



#Q2: Which abiotic factors are the most important drivers of fruit quality?---------
###original analysis with latitude and longitude as PCs

#Analysis: principle components analysis followed by PC regression with each variable
p_clim <- dplyr::select(c,c("Prox.Water", "Longitude","Latitude","elevation", 
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
#    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9 
#0.57948 0.19358 0.11433 0.07249 0.02452 0.00894 0.00572 0.00094 0.00000 

df <- data.frame(PC=1:9, var_explained=var_explained)

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
#The first matrix shows the correlation coefficients between the variables  
#the second matrix shows the corresponding p-values.
rcorr(as.matrix(p_clim))

#now let's visualize this 
corrplot(cor(p_clim))

testRes = cor.mtest(p_clim, conf.level = 0.95)

## add all p-values
corrplot(cor(p_clim), p.mat = testRes$p, insig = 'p-value', sig.level = -1)

## add significant level stars
corrplot(cor(p_clim), p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')


#PC's 5, 6, 7, 8, 9 are very low/ not worth looking at



#Q2-A: Physical Quality---------------------------------------------------------
##Linear models PC x quality 
pc_clim_phys <- cbind(pc_clim, c)

hist(c$Firmness)
hist(c$SSC)
hist(c$avgwgt)
hist(c$maturity.index)

###Firmness###

p1 <- glmmTMB(Firmness ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p1)
#PC1                  1.57530    0.27341   5.762 8.33e-09 ***
#PC3                  2.70855    0.61936   4.373 1.22e-05 ***


plot(Firmness ~ PC1, data=pc_clim_phys)
plot(Firmness ~ PC3, data=pc_clim_phys)

results$rotation

###SSC###
P2 <- glmmTMB(SSC ~ orchard.type + PC1 + PC2 + PC3 + PC4 + (1|site.code), data= pc_clim_phys)
summary(P2)
#PC1                 -1.01359    0.15615  -6.491 8.51e-11 ***
#PC2                 -0.58994    0.28553  -2.066   0.0388 *  
#PC3                  0.87971    0.37121   2.370   0.0178 *


plot(SSC ~ PC1, data=pc_clim_phys)
plot(SSC ~ PC2, data=pc_clim_phys)
plot(SSC ~ PC3, data=pc_clim_phys)

results$rotation

###AVGWGT
p3 <- glmmTMB(avgwgt ~ orchard.type + PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p3)
#PC2                  -8.2380     4.2160  -1.954   0.0507 .  
#PC3                  10.2852     5.2833   1.947   0.0516 .  
#PC4                 -12.4581     6.0755  -2.051   0.0403 *   

plot(avgwgt ~ PC2, data=pc_clim_phys)
plot(avgwgt ~ PC3, data=pc_clim_phys)
plot(avgwgt ~ PC4, data=pc_clim_phys)

results$rotation

###Maturity Index 
p4 <- glmmTMB(maturity.index ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p4)
#PC1                 -0.26665    0.08381  -3.181  0.00147 ** 
#PC4                 -0.44814    0.20812  -2.153  0.03129 *  


plot(maturity.index ~ PC1, data=pc_clim_phys)
plot(maturity.index ~ PC4, data=pc_clim_phys)

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
#PC1                 -0.08577    0.05167  -1.660 0.096964 .  
#PC2                  0.27957    0.08336   3.354 0.000797 ***  

plot(TotalPhen ~ PC2, data=pc.pu_clim)


#PhenRich 
pc.pu.pr <- glmmTMB(PhenRich ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim)
summary(pc.pu.pr)
#PC1           0.3597     0.1899    1.89  0.05823 .  
#PC2           0.9112     0.3412    2.67  0.00758 ** 


plot(PhenRich ~ PC1, data=pc.pu_clim)
plot(PhenRich ~ PC2, data=pc.pu_clim)
results$rotation


#SEED#
pc.se_clim <- cbind(pc_clim, SeedD)

#TotalPhen
pc.se.tp <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.se_clim, family=beta_family (link="logit"))
summary(pc.se.tp)
#PC2          0.128315   0.055966    2.29   0.0219 *  

plot(TotalPhen ~ PC2, data=pc.se_clim)

results$rotation


#PhenRich 
pc.se.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 +
                      (1|site.code), data=pc.se_clim)
summary(pc.se.pr)
#PC1           0.3597     0.1899    1.89  0.05823 .  
#PC2           0.9112     0.3412    2.67  0.00758 **


plot(PhenRich ~ PC1, data=pc.se_clim)
plot(PhenRich ~ PC2, data=pc.se_clim)

results$rotation

#Q1-C: Which compounds distinguish fruits based on management systems?----------
#Analysis: NMDS and random forest
###NMDS###
m.NMDS <- metaMDS(d.comp, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS

plot(m.NMDS, type="t")

#Data:     wisconsin(sqrt(d.comp)) 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.135198  
#Stress type 1, weak ties
#Best solution was not repeated after 100 tries

#for plotting, need to add columns that give values for 
#colors and symbols we want in plot
d.expl$Color <- recode_factor(d.expl$orchard.type,
                              Organic="red", Conventional="blue")
d.expl$Symbol <- recode_factor(d.expl$Tissue,
                               SKIN=2, PULP=3, SEED=1)
d.expl$Symbol <- as.numeric(as.character(d.expl$Symbol))

d.expl$orchard.type = as.factor(d.expl$orchard.type)


plot(m.NMDS, type="n") #plots the ordination axes only
points(m.NMDS, pch=d.expl$Symbol,
       col=as.character(d.expl$Color), cex = 0.8)     
ordiellipse(m.NMDS, d.expl$Symbol, conf = 0.95)




#PERMANOVA can test whether the visualized differences are significant
m.perm <- adonis2(d.comp~orchard.type*Tissue, data=d.expl)
m.perm
#orchard.type          1    0.450 0.00386   2.5040  0.024 *  
#Tissue                2   51.790 0.44450 144.1093  0.001 ***
#orchard.type:Tissue   2    0.843 0.00724   2.3465  0.008 ** 
#Residual            353   63.431 0.54440                    
#Total               358  116.515 1.00000  



###Random Forest###
###skin
m1.rf.sk <- randomForest(d.comp.sk,d.expl.sk$orchard.type, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.sk$importance
varImpPlot(m1.rf.sk)
MDSplot(m1.rf.sk, d.expl.sk$orchard.type)

#using boruta 
m1.rf.b.sk <- Boruta(d.comp.sk,d.expl.sk$orchard.type)
m1.rf.b.sk
# 9 attributes confirmed important: A, cyanidin_Galactoside, ePicatechin, F, G and 4 more;

#plot 
plot(m1.rf.b.sk)  
getSelectedAttributes(m1.rf.b.sk) #lists all important ones

#[1] "A"                    "F"                    "G"                    "PB1"                 
#[5] "cyanidin_Galactoside" "ePicatechin"          "U1"                   "I"                   
#[9] "K" 

#performing MANOVA 
d.comp.sk.sel <- data.matrix(d.comp.sk[,getSelectedAttributes(m1.rf.b.sk)])
m1.man.sk <- manova(d.comp.sk.sel ~ d.expl.sk$orchard.type)
summary(m1.man.sk) 
#d.expl.sk$orchard.type   1 0.23711   3.7986      9    110 0.0003348 ***

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.sk)  
#A= 0.09424 .
#PB1 =0.05021 .
#gala = 0.08316 .

#some quick plots of all of them
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.sk.sel))){
  d.temp=d.comp.sk.sel[,i]
  plot(d.temp ~ d.expl.sk$orchard.type, ylab=colnames(d.comp.sk.sel)[i])
}
dev.off()

#all three marginally higher in organic 

###PULP###

m1.rf.pu <- randomForest(d.comp.pu,d.expl.pu$orchard.type, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.pu$importance
varImpPlot(m1.rf.pu)
MDSplot(m1.rf.pu, d.expl.pu$orchard.type)

m1.rf.b.pu <- Boruta(d.comp.pu,d.expl.pu$orchard.type)
m1.rf.b.pu
plot(m1.rf.b.pu)  
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
plot(m1.rf.b.se)  #important variables (better than shadow) are in green
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

#Q1-D: Which compounds distinguish fruits based on latitude?----------
#Analysis: NMDS and random forest

#left joining the "expl" by latitude for RF 
d.expl <- left_join(d.expl, Orchard[,c(2,40)], by="orchard.num")
d.expl.sk <- left_join(d.expl.sk, Orchard[,c(2,40)], by="orchard.num")
d.expl.pu <- left_join(d.expl.pu, Orchard[,c(2,40)], by="orchard.num")
d.expl.se <- left_join(d.expl.se, Orchard[,c(2,40)], by="orchard.num")


###NMDS###
m.NMDS_2 <- metaMDS(d.comp, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS_2
plot(m.NMDS_2, type="t")

d.expl$Color <- recode_factor(d.expl$lat_cat,
                              low="red", high="blue")
d.expl$Symbol <- recode_factor(d.expl$Tissue,
                               SKIN=2, PULP=3, SEED=1)
d.expl$Symbol <- as.numeric(as.character(d.expl$Symbol))

d.expl$lat_cat = as.factor(d.expl$lat_cat)


plot(m.NMDS_2, type="n") #plots the ordination axes only
points(m.NMDS_2, pch=d.expl$Symbol,
       col=as.character(d.expl$Color), cex = 0.8)     
ordiellipse(m.NMDS_2, d.expl$Symbol, conf = 0.95)




#PERMANOVA can test whether the visualized differences are significant
m.perm <- adonis2(d.comp~lat_cat*Tissue, data=d.expl)
m.perm
#lat_cat          1    0.390 0.00334   2.1747  0.041 *  
#Tissue           2   51.788 0.44448 144.5504  0.001 ***
#lat_cat:Tissue   2    1.102 0.00946   3.0766  0.002 ** 
#Residual       353   63.235 0.54272                    
#Total          358  116.515 1.00000  



###Random Forest###
###skin
m1.rf.sk <- randomForest(d.comp.sk,d.expl.sk$lat_cat, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.sk$importance
varImpPlot(m1.rf.sk)
MDSplot(m1.rf.sk, d.expl.sk$lat_cat)

#using boruta 
m1.rf.b.sk <- Boruta(d.comp.sk,d.expl.sk$lat_cat)
m1.rf.b.sk
#12 important attributes 

#plot 
plot(m1.rf.b.sk)  
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
plot(m1.rf.b.pu)  

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
plot(m1.rf.b.se)  #important variables (better than shadow) are in green
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

