#phenolic analysis
rm(list=ls()) # clears work space
#Install packages-------------------------------------------------------------
library(ggplot2)
library(MASS)
library(reshape2)
library(dplyr)
library(randomForest)
library(glmmTMB) #mixed models
library(car) #mixed models summary tables
library(vegan) #multivariate stats
library(Boruta) #random forest models
library(readr)
library(pls) #PCAs
library(Hmisc)#correlation matrix calculation 
library(corrplot)#plotting correlation matrix 
library(car)

######Data Organization and Restructuring###--------------------------------------------


#Fruit Level Data--------------------------------------------------------------
#This data set is going to be used to answer one portion of Q1. This data set, 
#contains the data for each known compound as well as the mgmt system, and site
#code

d <- read_csv("chem_all_dat.csv")
View(d)
d <- d %>%
  mutate(TotalPhen=rowSums(across(8:29)),
         PhenRich=rowSums(across(8:29)!=0))

d <- d %>% 
  mutate_at(c("Orchard.num", "orchard.type", "site.code", "Tree", "SampleID"), as.factor)
View(d)


#logit transformation worked the best so were going to use that 
d <- d %>%
  mutate(TotalPhentrans= log10(d$TotalPhen))
d$TotalPhentrans[d$TotalPhentrans==-Inf] <- NA
na.omit(d$TotalPhentrans)


shapiro.test(d$TotalPhentrans)
#W = 0.86448, p-value < 2.2e-16
hist(d$TotalPhen)
hist(d$TotalPhentrans)
#not good but looks better so we'll mess around with this 
View(d)


#separate by tissue type 
d.sk <- filter(d, Tissue=="SKIN")
d.sk <- dplyr::select(d.sk, -(6+which(colSums(d.sk[7:32], na.rm=TRUE) %in% 0)))
d.pu <- filter(d, Tissue=="PULP")
d.pu <- dplyr::select(d.pu, -(6+which(colSums(d.pu[7:32], na.rm=TRUE) %in% 0)))
d.se <- filter(d, Tissue=="SEED")
d.se <- dplyr::select(d.se, -(6+which(colSums(d.se[7:32], na.rm=TRUE) %in% 0)))

#Assign row names
d <- as.data.frame(d)
row.names(d) <- d$SampleID

#For some analyses we need a table with only composition info
d.comp <- d[,8:29]
#and one with just explanatory variables
d.expl <- d[,1:6]


#Also need those by tissue
#skin 
d.comp.sk <- d.sk[,8:29]
d.expl.sk <- d.sk[,1:6]
#pulp
d.comp.pu <- d.pu[,8:29]
d.expl.pu <- d.pu[,1:6]
#Seed
d.comp.se <- d.se[,8:29]
d.expl.se <- d.se[,1:6]

#Orchard & Tree Level Data----------------------------------------
#These data sets will be condensed to the orchard level (24 observations)
#We're going to create four different data sets to analyze parts of the fruit: 
#Whole Fruit, Skin, Seed, and Pulp

###Orchard and Tree level data will be the same for all###
Orchard <- read_csv("Orchard Level Data.csv")
View(Orchard)

#convert weights based on averages and bag weights
Tree <- read_csv("Tree Level Data.csv")
View(Tree)
Tree<- Tree %>%
  mutate(avgwgt=(apple.wgt.3-bag.weight)/3)%>%
dplyr::select(-c(3:4))

Tree <- Tree %>%
  group_by(orchard.num) %>%
  summarise_at(c("Firmness", "SSC", "maturity.index", "avgwgt")
               , mean, na.rm = TRUE)
View(Tree)

#Whole Fruit--------------------------------------------------------------------
Fruit0 <- read_csv("Fruit Level Data.csv")
Fruit0 <- Fruit0 %>%
  mutate(TotalPhen=rowSums(across(8:29)),
         PhenRich=rowSums(across(8:29)!=0))

Fruit0 <- Fruit0 %>%
  mutate(TotalPhentrans= log10(Fruit0$TotalPhen))
Fruit0$TotalPhentrans[Fruit0$TotalPhentrans==-Inf] <- NA
na.omit(Fruit0$TotalPhentrans)
view(Fruit0$TotalPhentrans)

Fruit0 <- Fruit0 %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich","TotalPhentrans")
               , mean, na.rm = TRUE)

CCD <- left_join(Tree, Orchard, by="orchard.num")
CCD <- left_join(CCD, Fruit0, by="orchard.num")
View(CCD)

#Skin---------------------------------------------------------------------------
#condense chem data by tissue and orchard 
Fruit1 <- read_csv("Fruit Level Data.csv")
Fruit1 <- Fruit1 %>%
  mutate(TotalPhen=rowSums(across(8:29)),
         PhenRich=rowSums(across(8:29)!=0))

Fruit1 <- Fruit1 %>%
  mutate(TotalPhentrans= log10(Fruit1$TotalPhen))
Fruit1$TotalPhentrans[Fruit1$TotalPhentrans==-Inf] <- NA
na.omit(Fruit1$TotalPhentrans)



Fruit1 <- Fruit1 %>%
dplyr::filter(Tissue != "PULP") %>%
dplyr::filter(Tissue != "SEED") 


Fruit1 <- Fruit1 %>%
group_by(orchard.num) %>%
summarise_at(c("TotalPhen", "PhenRich", "TotalPhentrans")
               , mean, na.rm = TRUE)

View(Fruit1)

#combine all data sets to 120 values 
SkinD <- left_join(Tree, Orchard, by="orchard.num")
SkinD <- left_join(SkinD, Fruit1, by="orchard.num")
View(SkinD)

#Pulp---------------------------------------------------------------------------
Fruit2 <- read_csv("Fruit Level Data.csv")
Fruit2 <- Fruit2 %>%
  mutate(TotalPhen=rowSums(across(8:29)),
         PhenRich=rowSums(across(8:29)!=0))

Fruit2 <- Fruit2 %>%
  mutate(TotalPhentrans= log10(Fruit2$TotalPhen))
Fruit2$TotalPhentrans[Fruit2$TotalPhentrans==-Inf] <- NA
na.omit(Fruit2$TotalPhentrans)


Fruit2 <- Fruit2 %>%
  dplyr::filter(Tissue != "SKIN") %>%
  dplyr::filter(Tissue != "SEED") 


Fruit2 <- Fruit2 %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich", "TotalPhentrans")
               , mean, na.rm = TRUE)

View(Fruit2)

PulpD <- left_join(Tree, Orchard, by="orchard.num")
PulpD <- left_join(PulpD, Fruit2, by="orchard.num")
View(PulpD)
#Seed#--------------------------------------------------------------------------
Fruit3 <- read_csv("Fruit Level Data.csv")
Fruit3 <- Fruit3 %>%
  mutate(TotalPhen=rowSums(across(8:29)),
         PhenRich=rowSums(across(8:29)!=0))

Fruit3 <- Fruit3 %>%
  mutate(TotalPhentrans= log10(Fruit3$TotalPhen))
Fruit3$TotalPhentrans[Fruit3$TotalPhentrans==-Inf] <- NA
na.omit(Fruit3$TotalPhentrans)


Fruit3 <- Fruit3 %>%
  dplyr::filter(Tissue != "SKIN") %>%
  dplyr::filter(Tissue != "PULP") 

Fruit3 <- Fruit3 %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich", "TotalPhentrans")
               , mean, na.rm = TRUE)

View(Fruit3)

#combine all data sets to 120 values 
SeedD <- left_join(Tree, Orchard, by="orchard.num")
SeedD <- left_join(SeedD, Fruit3, by="orchard.num")
View(SeedD)



#####Data Analysis-------------------------------------------------------------------



#How do management systems (organic vs conventional) shape chemical comp and richness?-------
#Analyses: one linear mixed models for total phenolics and phenolic richness for the 
#whole fruit and for each tissue type 

###total phenolics###
tp1 <- glmmTMB(TotalPhentrans~ orchard.type*Tissue + (1|site.code/Orchard.num/Tree), data=d)
summary(tp1)
Anova(tp1)
#orchard.type          1.6935  1    0.19314    
#Tissue              321.4999  2    < 2e-16 ***
#orchard.type:Tissue   5.1914  2    0.07459 . 

#total p per tissue type 
tp.p <- glmmTMB(TotalPhentrans ~ orchard.type + (1|site.code/Orchard.num/Tree), data=d.pu)
summary(tp.p)
Anova(tp.p)

tp.sk <- glmmTMB(TotalPhentrans ~ orchard.type + (1|site.code/Orchard.num/Tree), data=d.sk)
summary(tp.sk)
Anova(tp.sk)

tp.se <- glmmTMB(TotalPhentrans ~ orchard.type + (1|site.code/Orchard.num/Tree), data=d.se)
summary(tp.se)
Anova(tp.se) #p=0.006587 **


###Q1B:phenolics richness###
pr1 <- glmmTMB(PhenRich~ orchard.type*Tissue + (1|site.code/Orchard.num/Tree), data=d, family= poisson)
summary(pr1)
Anova(pr1)

#orchard.type           2.2590  1    0.16471    
#Tissue              1279.1586  2    < 2e-16 ***
#orchard.type:Tissue    6.4470  2    0.03982 *  

diagnose(pr1)
shapiro.test(resid(pr1))
hist(resid(pr1))
#= 0.98958, p-value = 0.01156


#phen rich per tissue type 
pr.p <- glmmTMB(PhenRich ~ orchard.type+(1|site.code/Orchard.num/Tree), data=d.pu, family= poisson)
summary(pr.p)
Anova(pr.p)
#not sig 

pr.sk <- glmmTMB(PhenRich ~ orchard.type+(1|site.code/Orchard.num/Tree), data=d.sk, family= poisson)
summary(pr.sk)
Anova(pr.sk)
#not sig 

pr.se <- glmmTMB(PhenRich ~ orchard.type+(1|site.code/Orchard.num/Tree), data=d.se, family= poisson)
summary(pr.se)
Anova(pr.se)
#orchard.type  0.01726 *
#Which compounds distinguish fruits raised in their respective management systems?----------
#Analysis: NMDS and random forest
###NMDS###
#someone online said to change the distance due to 0 values 
m.NMDS <- metaMDS(d.comp, distance = "euclidean", trymax=100, autotransform =TRUE)
m.NMDS
plot(m.NMDS, type="t")

#Dimensions: 2 
#Stress:     0.1048574 (fair?)
#Stress type 1, weak ties
#Best solution was repeated 1 time in 56 tries

#for plotting, need to add columns that give values for 
#colors and symbols we want in plot
d.expl$Color <- recode_factor(d.expl$orchard.type,
                              Organic="red", Conventional="blue")
d.expl$Symbol <- recode_factor(d.expl$Tissue,
                               SKIN=8, PULP=1, SEED=2)
d.expl$Symbol <- as.numeric(as.character(d.expl$Symbol))

d.expl$orchard.type = as.factor(d.expl$orchard.type)


plot(m.NMDS, type="n") #plots the ordination axes only
points(m.NMDS, pch=d.expl$Symbol,
       col=as.character(d.expl$Color), cex = 0.8)     

#PERMANOVA can test whether the visualized differences are significant
##this portion doesn't work because of 0s 
m.perm <- adonis2(d.comp~orchard.type*Tissue, data=d.expl)
m.perm

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
#epicatechin, pb1, syr.acid significant 

#plot 
plot(m1.rf.b.sk)  
getSelectedAttributes(m1.rf.b.sk) #lists all important ones

#performing MANOVA 
d.comp.sk.sel <- data.matrix(d.comp.sk[,getSelectedAttributes(m1.rf.b.sk)])
m1.man.sk <- manova(d.comp.sk.sel ~ d.expl.sk$orchard.type)
summary(m1.man.sk) 
#overall significance for MANOVA p= 0.3281

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.sk)  
#pb1 = 0.07191 .
#syr.acid= 0.5659
#epi= 0.4211

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
plot(m1.rf.b.pu)  
# No attributes deemed important.

getSelectedAttributes(m1.rf.b.pu) 


##Running MANOVAS (needs to be fixed)
d.comp.pu.sel <- data.matrix(d.comp.pu[,getSelectedAttributes(m1.rf.b.pu)])
m1.man.pu <-manova(d.comp.pu.sel ~ d.expl.pu$orchard.type)
summary(m1.man.pu)  #overall significance for MANOVA
summary.aov(m1.man.pu)  #follow-up ANOVAs for each individual compound
#0.9745

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
plot(m1.rf.b.se)  #important variables (better than shadow) are in green
getSelectedAttributes(m1.rf.b.se) #lists all important ones
#"pb2"        "syr.acid"   "reynoutrin" "U2" 

#Running MANOVAS 
d.comp.se.sel <- data.matrix(d.comp.se[,getSelectedAttributes(m1.rf.b.se)])
m1.man.se <- manova(d.comp.se.sel ~ d.expl.se$orchard.type)
summary(m1.man.se)  
#0.008044 **
summary.aov(m1.man.se)  
#ppb2=0.004867 **
#syr.acid = 0.0003586 ***
#reynoutrin = 0.8561
#U2 =  0.006278 **

par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.se.sel))){
  d.temp=d.comp.se.sel[,i]
  plot(d.temp ~ d.expl.se$orchard.type, ylab=colnames(d.comp.se.sel)[i])
}
dev.off()
#values higher in conventional orchards 

#How do management systems interact with broad climatic changes across latitude?-----

###whole fruit###
#Phen Rich 
pr.la <- glmmTMB(PhenRich~ orchard.type*Latitude + (1|site.code), data=CCD)
summary(pr.la)
Anova(pr.la)
#orchard.type:Latitude 3.8925  1     0.0485 *

pr.lo <- glmmTMB(PhenRich~ orchard.type*Longitude + (1|site.code), data=CCD)
summary(pr.lo)
Anova(pr.lo)
#orchard.type           3.2739  1    0.07039 .
#Longitude              4.4290  1    0.03533 *

pr.el <- glmmTMB(PhenRich~ orchard.type*elevation + (1|site.code), data=CCD)
summary(pr.el)
Anova(pr.el)
#elevation              10.4298  1    0.00124 **

#Total Phen 
tp.la <- glmmTMB(TotalPhentrans~ orchard.type*Latitude + (1|site.code), data=CCD)
summary(tp.la)
Anova(tp.la)

tp.lo <- glmmTMB(TotalPhentrans~ orchard.type*Longitude, data=CCD)
summary(tp.lo)
Anova(tp.lo)
#Longitude              11.3098  1   0.000771 ***

tp.el <- glmmTMB(TotalPhentrans~ orchard.type*elevation, data=CCD)
summary(tp.el)
Anova(tp.el)
#elevation              6.9692  1   0.008293 **


#Which abiotic factors are the most important drivers of fruit chem?------------

#Analysis: principle components analysis followed by PC regression with each variable
#followed by a correlation matrix for visual interpretation 
#1 for whole fruit, skin, pulp, seed

###WHOLE FRUIT###
p_clim <- dplyr::select(CCD,c("Prox.Water", "Longitude","Latitude","elevation", 
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
biplot(results,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


#calculate total variance explained by each principal component
summary(results)$importance
summary(results)$importance[2,]

var_explained = results$sdev^2 / sum(results$sdev^2)

df <- data.frame(PC=1:9, var_explained=var_explained)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


#PC's  7, 8, 9 are very low/ not worth looking at

##Linear models PC x quality 

pc_clim <- as.data.frame(results$x)
pc_clim <- cbind(pc_clim, CCD)


##TotalPhen###

p1 <- glmmTMB(TotalPhentrans ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + (1|site.code), data=pc_clim)
summary(p1)
#PC2                 -0.127781   0.029778   -4.29 1.78e-05 ***     

plot(TotalPhen ~ PC2, data=pc_clim)

results$rotation

p2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + (1|site.code), data=pc_clim)
summary(p2)
#PC2  0.000175 ***    

plot(PhenRich ~ PC2, data=pc_clim)

results$rotation

###correlation matrix 
p_phenclim <- dplyr::select(CCD,c("Prox.Water", "Longitude","Latitude","elevation", 
 "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI", 
 "TotalPhentrans", "PhenRich"))

rcorr(as.matrix(p_phenclim))

corrplot(cor(p_phenclim))

###SKIN###
psk_clim <- dplyr::select(SkinD,c("Prox.Water", "Longitude","Latitude","elevation", 
"Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI"))

#calculate principal components
results.sk <- prcomp(psk_clim, scale = TRUE)

#this gives you the % variance explained
summary(results.sk)  

#display principal components
results.sk$x

results$rotation
#display the first six scores
head(results.sk$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results.sk,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


#calculate total variance explained by each principal component
summary(results.sk)$importance
summary(results.sk)$importance[2,]

var_explained.sk = results.sk$sdev^2 / sum(results.sk$sdev^2)

df <- data.frame(PC=1:9, var_explained.sk=var_explained.sk)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained.sk)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


#PC's  6, 7, 8, 9 are very low/ not worth looking at

##Linear models PC x quality 

pc.sk_clim <- as.data.frame(results.sk$x)
pc.sk_clim <- cbind(pc.sk_clim, SkinD)


##TotalPhen###
psk1 <- glmmTMB(TotalPhen ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc.sk_clim)
summary(psk1)
#PC2  0.00102 **     

plot(TotalPhen ~ PC2, data=pc_clim)

results.sk$rotation

psk2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc.sk_clim)
summary(psk2)
#nothing   


###correlation matrix 
pc.sk_phenclim <- dplyr::select(SkinD,c("Prox.Water", "Longitude","Latitude","elevation", 
"Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI", 
"TotalPhen", "PhenRich"))

rcorr(as.matrix(pc.sk_phenclim))
corrplot(cor(pc.sk_phenclim))

###PULP###
p.pu_clim <- dplyr::select(PulpD,c("Prox.Water", "Longitude","Latitude","elevation", 
"Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI"))

#calculate principal components
results.pu <- prcomp(p.pu_clim, scale = TRUE)

#this gives you the % variance explained
summary(results.pu)  

#display principal components
results.pu$x

results.pu$rotation
#display the first six scores
head(results.pu$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results.pu,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


#calculate total variance explained by each principal component
summary(results.pu)$importance
summary(results.pu)$importance[2,]

var_explained.pu = results.pu$sdev^2 / sum(results.pu$sdev^2)

df <- data.frame(PC=1:9, var_explained.pu=var_explained.pu)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained.pu)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


#PC's  6, 7, 8, 9 are very low/ not worth looking at

##Linear models PC x quality 

pc.pu_clim <- as.data.frame(results.pu$x)
pc.pu_clim <- cbind(pc.pu_clim, PulpD)


##TotalPhen###
ppu1 <- glmmTMB(TotalPhen ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc.pu_clim)
summary(ppu1)
#PC1                  -1238.2      622.5  -1.989  0.04669 *  
#PC2                  -4267.0     1082.2  -3.943 8.05e-05 ***
#PC3                  -4064.2     1410.4  -2.882  0.00396 **

plot(TotalPhen ~ PC1, data=pc.pu_clim)
plot(TotalPhen ~ PC2, data=pc.pu_clim)
plot(TotalPhen ~ PC3, data=pc.pu_clim)

results.pu$rotation

ppu2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc.pu_clim)
summary(ppu2)
#PC2                 -0.86158    0.24963  -3.451 0.000558 ***
#PC3                 -0.63019    0.32536  -1.937 0.052757 .     

plot(PhenRich ~ PC2, data=pc.pu_clim)
plot(PhenRich ~ PC3, data=pc.pu_clim)

results.pu$rotation

###correlation matrix 
p.pu_phenclim <- dplyr::select(PulpD,c("Prox.Water", "Longitude","Latitude","elevation", 
"Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI", 
"TotalPhen", "PhenRich"))

rcorr(as.matrix(p.pu_phenclim))

corrplot(cor(p.pu_phenclim))



###SEED###
p.se_clim <- dplyr::select(PulpD,c("Prox.Water", "Longitude","Latitude","elevation", 
"Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI"))

#calculate principal components
results.se <- prcomp(p.se_clim, scale = TRUE)

#this gives you the % variance explained
summary(results.se)  

#display principal components
results.se$x

results.se$rotation
#display the first six scores
head(results.se$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results.se,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


#calculate total variance explained by each principal component
summary(results.se)$importance
summary(results.se)$importance[2,]

var_explained.se = results.se$sdev^2 / sum(results.se$sdev^2)

df <- data.frame(PC=1:9, var_explained.se=var_explained.se)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained.se)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


#PC's  5, 6, 7, 8, 9 are very low/ not worth looking at

##Linear models PC x quality 
pc.se_clim <- as.data.frame(results.se$x)
pc.se_clim <- cbind(pc.se_clim, SeedD)


##TotalPhen###
p.se1 <- glmmTMB(TotalPhen ~orchard.type+ PC1 + PC2 + PC3 + PC4, data=pc.se_clim)
summary(p.se1)
#PC2                   -17200       6114  -2.813   0.0049 **

plot(TotalPhen ~ PC2, data=pc.se_clim)
results.se$rotation

p.se2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc.se_clim)
summary(p.se2)
#PC1                   0.3052     0.1202   2.539  0.01111 *  
#PC2                  -1.2048     0.2112  -5.704 1.17e-08 ***
#PC3                  -0.5668     0.2741  -2.068  0.03863 *    

plot(PhenRich ~ PC1, data=pc.se_clim)
plot(PhenRich ~ PC2, data=pc.se_clim)
plot(PhenRich ~ PC3, data=pc.se_clim)

results.se$rotation

###correlation matrix 
p.se_phenclim <- dplyr::select(PulpD,c("Prox.Water", "Longitude","Latitude","elevation", 
 "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI", 
                                       "TotalPhen", "PhenRich"))

rcorr(as.matrix(p.se_phenclim))

corrplot(cor(p.se_phenclim))

#Which specific management practices are the most important drivers of fruit chem?----

#Analysis: principle components analysis followed by PC regression with each variable
#followed by a correlation matrix for visual interpretation 
#1 for whole fruit, skin, pulp, seed

###WHOLE FRUIT###
p_mgmt <- dplyr::select(CCD,c("Cultivation","Herbicides","Com_Mul", "Mowing",
                            "Weed_Mats", "Cover_Crops","Fire_Mgmt","Acres"))

#calculate principal components
results1 <- prcomp(p_mgmt, scale = TRUE)

#this gives you the % variance explained
summary(results1)  

#display principal components
results1$x

results1$rotation
#display the first six scores
head(results1$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results1,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


#calculate total variance explained by each principal component
summary(results1)$importance
summary(results1)$importance[2,]

var_explained1 = results1$sdev^2 / sum(results1$sdev^2)

df1 <- data.frame(PC=1:8, var_explained1=var_explained1)

#create scree plot
ggplot(df1, aes(x=PC, y=var_explained1)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


##Linear models PC x quality 
#all sig

pc_mgmt <- as.data.frame(results1$x)
pc_mgmt <- cbind(pc_mgmt, CCD)

###TotalPhen###
p3 <- glmmTMB(TotalPhen ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 
+ PC8 +(1|site.code), data=pc_mgmt)
summary(p3)
#PC1                  -8693.5     2698.2  -3.222  0.00127 **
#PC6                 -20437.0     2704.7  -7.556 4.15e-14 ***

plot(TotalPhen ~ PC1, data=pc_mgmt)
plot(TotalPhen ~ PC6, data=pc_mgmt)

###PhenRich###
p4 <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 
              + PC8 +(1|site.code), data=pc_mgmt)
summary(p4)
#PC1                 -0.34937    0.14959  -2.335 0.019518 *  
#PC5                 -0.87040    0.36264  -2.400 0.016388 *  
#PC6                 -0.38873    0.15781  -2.463 0.013764 *  
#PC8                 -1.02755    0.26770  -3.838 0.000124 ***

plot(TotalPhen ~ PC1, data=pc_mgmt)
plot(TotalPhen ~ PC5, data=pc_mgmt)
plot(TotalPhen ~ PC6, data=pc_mgmt)
plot(TotalPhen ~ PC8, data=pc_mgmt)

###correlation matrix 
p_phenmgmt <- dplyr::select(CCD,c("Cultivation","Herbicides","Com_Mul", "Mowing",
"Weed_Mats", "Cover_Crops","Fire_Mgmt","Acres", "TotalPhen", "PhenRich"))


rcorr(as.matrix(p_phenmgmt))

corrplot(cor(p_phenmgmt))

###SKIN###
p.sk_mgmt <- dplyr::select(SkinD,c("Prox.Water", "Longitude","Latitude","elevation", 
  "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI"))

#calculate principal components
results.sk1 <- prcomp(p.sk_mgmt, scale = TRUE)

#this gives you the % variance explained
summary(results.sk1)  

#display principal components
results.sk1$x

results$rotation
#display the first six scores
head(results.sk1$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results.sk1,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


#calculate total variance explained by each principal component
summary(results.sk1)$importance
summary(results.sk1)$importance[2,]

var_explained.sk1 = results.sk1$sdev^2 / sum(results.sk1$sdev^2)

df <- data.frame(PC=1:9, var_explained.sk1=var_explained.sk1)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained.sk1)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


#PC's  6, 7, 8, 9 are very low/ not worth looking at

##Linear models PC x quality 

pc.sk_mgmt <- as.data.frame(results.sk$x)
pc.sk_mgmt <- cbind(pc.sk_mgmt, SkinD)


##TotalPhen###
psk1 <- glmmTMB(TotalPhen ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc.sk_mgmt)
summary(psk1)
#PC2  0.00102 **     

plot(TotalPhen ~ PC2, data=pc.sk_mgmt)

results.sk$rotation

psk2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc.sk_mgmt)
summary(psk2)
#nothing   


###correlation matrix 
pc.sk_phenclim <- dplyr::select(SkinD,c("Prox.Water", "Longitude","Latitude","elevation", 
                                        "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI", 
                                        "TotalPhen", "PhenRich"))

rcorr(as.matrix(pc.sk_phenclim))
corrplot(cor(pc.sk_phenclim))

###PULP###
p.pu_clim <- dplyr::select(PulpD,c("Prox.Water", "Longitude","Latitude","elevation", 
                                   "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI"))

#calculate principal components
results.pu <- prcomp(p.pu_clim, scale = TRUE)

#this gives you the % variance explained
summary(results.pu)  

#display principal components
results.pu$x

results.pu$rotation
#display the first six scores
head(results.pu$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results.pu,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


#calculate total variance explained by each principal component
summary(results.pu)$importance
summary(results.pu)$importance[2,]

var_explained.pu = results.pu$sdev^2 / sum(results.pu$sdev^2)

df <- data.frame(PC=1:9, var_explained.pu=var_explained.pu)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained.pu)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


#PC's  6, 7, 8, 9 are very low/ not worth looking at

##Linear models PC x quality 

pc.pu_clim <- as.data.frame(results.pu$x)
pc.pu_clim <- cbind(pc.pu_clim, PulpD)


##TotalPhen###
ppu1 <- glmmTMB(TotalPhen ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc.pu_clim)
summary(ppu1)
#PC1                  -1238.2      622.5  -1.989  0.04669 *  
#PC2                  -4267.0     1082.2  -3.943 8.05e-05 ***
#PC3                  -4064.2     1410.4  -2.882  0.00396 **

plot(TotalPhen ~ PC1, data=pc.pu_clim)
plot(TotalPhen ~ PC2, data=pc.pu_clim)
plot(TotalPhen ~ PC3, data=pc.pu_clim)

results.pu$rotation

ppu2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc.pu_clim)
summary(ppu2)
#PC2                 -0.86158    0.24963  -3.451 0.000558 ***
#PC3                 -0.63019    0.32536  -1.937 0.052757 .     

plot(PhenRich ~ PC2, data=pc.pu_clim)
plot(PhenRich ~ PC3, data=pc.pu_clim)

results.pu$rotation

###correlation matrix 
p.pu_phenclim <- dplyr::select(PulpD,c("Prox.Water", "Longitude","Latitude","elevation", 
                                       "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI", 
                                       "TotalPhen", "PhenRich"))

rcorr(as.matrix(p.pu_phenclim))

corrplot(cor(p.pu_phenclim))



###SEED###
p.se_clim <- dplyr::select(PulpD,c("Prox.Water", "Longitude","Latitude","elevation", 
                                   "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI"))

#calculate principal components
results.se <- prcomp(p.se_clim, scale = TRUE)

#this gives you the % variance explained
summary(results.se)  

#display principal components
results.se$x

results.se$rotation
#display the first six scores
head(results.se$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results.se,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


#calculate total variance explained by each principal component
summary(results.se)$importance
summary(results.se)$importance[2,]

var_explained.se = results.se$sdev^2 / sum(results.se$sdev^2)

df <- data.frame(PC=1:9, var_explained.se=var_explained.se)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained.se)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


#PC's  5, 6, 7, 8, 9 are very low/ not worth looking at

##Linear models PC x quality 
pc.se_clim <- as.data.frame(results.se$x)
pc.se_clim <- cbind(pc.se_clim, SeedD)


##TotalPhen###
p.se1 <- glmmTMB(TotalPhen ~orchard.type+ PC1 + PC2 + PC3 + PC4, data=pc.se_clim)
summary(p.se1)
#PC2                   -17200       6114  -2.813   0.0049 **

plot(TotalPhen ~ PC2, data=pc.se_clim)
results.se$rotation

p.se2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc.se_clim)
summary(p.se2)
#PC1                   0.3052     0.1202   2.539  0.01111 *  
#PC2                  -1.2048     0.2112  -5.704 1.17e-08 ***
#PC3                  -0.5668     0.2741  -2.068  0.03863 *    

plot(PhenRich ~ PC1, data=pc.se_clim)
plot(PhenRich ~ PC2, data=pc.se_clim)
plot(PhenRich ~ PC3, data=pc.se_clim)

results.se$rotation

###correlation matrix 
p.se_phenclim <- dplyr::select(PulpD,c("Prox.Water", "Longitude","Latitude","elevation", 
                                       "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI", 
                                       "TotalPhen", "PhenRich"))

rcorr(as.matrix(p.se_phenclim))

corrplot(cor(p.se_phenclim))


#Which pest or diseases presence has the most significant affect on fruit chem-----

names(c)

p_pest <- dplyr::select(c,c("Anthracnose","European Canker","Bullseye Rot", "Powdery mildew",     
                            "Apple scab","Root Rot","Fire Blight","Apple Maggots","Codling Moth","Aphids","Tree Borer",
                            "Cedar Apple Rust","Bitter Rot","Leaf Roller","Horned Caterpillars","Trhips","Stemble",
                            "Scale", "Pest_Index"))


#calculate principal components
results2 <- prcomp(p_pest, scale = TRUE)

#this gives you the % variance explained
summary(results1)  

#display principal components
results1$x

results1$rotation
#display the first six scores
head(results1$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results1,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


#calculate total variance explained by each principal component
summary(results1)$importance
summary(results1)$importance[2,]

var_explained1 = results1$sdev^2 / sum(results1$sdev^2)

df1 <- data.frame(PC=1:8, var_explained1=var_explained1)

#create scree plot
ggplot(df1, aes(x=PC, y=var_explained1)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


##Linear models PC x quality 

pc_mgmt <- as.data.frame(results1$x)
pc_mgmt <- cbind(pc_mgmt, c
                 
                 
                 








#How does fruit quality compare to total phenolics and phenolic richness--------
#Figures------------------------------------------------------------------------ 

#How do management systems (organic vs conventional) shape chemical comp and richness?
#whole fruit 
p1 = ggplot(d, aes(x = orchard.type, y = TotalPhentrans)) +
  geom_point(aes(color = orchard.type)) +
  geom_smooth(method=glm, se=FALSE)+
  facet_wrap(~orchard.type)+
  scale_color_viridis_d()
p1

ggplot(d, aes(x = Tissue, y = PhenRich)) +
  geom_boxplot(aes(color = orchard.type)) +
  geom_smooth(method=glm, se=FALSE)+
  facet_wrap(~orchard.type)+
  scale_color_viridis_d()


ggplot(d, aes(x = Tissue, y = pb2)) +
  geom_boxplot(aes(color = orchard.type)) +
  geom_smooth(method=glm, se=FALSE)+
  facet_wrap(~orchard.type)+
  scale_color_viridis_d()

#Tissue
sk.p1 <- ggplot(d.sk, aes(x=orchard.type, y=TotalPhentrans, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Total Phenolics") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("Conventional", "Organic"))
sk.p1

pu.p1 <- ggplot(d.pu, aes(x=orchard.type, y=TotalPhentrans, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Total Phenolics") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("Conventional", "Organic"))
pu.p1

se.p1 <- ggplot(d.se, aes(x=orchard.type, y=TotalPhentrans, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Total Phenolics") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("Conventional", "Organic"))
se.p1



multiplot(se.p1,sk.p1, pu.p1)

#Which compounds distinguish fruits raised in their respective management systems?























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





