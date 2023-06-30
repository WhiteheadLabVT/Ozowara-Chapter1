#phenol analysis
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
###Data Organization###
###Data Organization 
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


#separate by tissue type 
d.sk <- filter(d, Tissue=="SKIN")
d.sk <- dplyr::select(d.sk, -(6+which(colSums(d.sk[7:29], na.rm=TRUE) %in% 0)))
d.pu <- filter(d, Tissue=="PULP")
d.pu <- dplyr::select(d.pu, -(6+which(colSums(d.pu[7:29], na.rm=TRUE) %in% 0)))
d.se <- filter(d, Tissue=="SEED")
d.se <- dplyr::select(d.se, -(6+which(colSums(d.se[7:29], na.rm=TRUE) %in% 0)))

#Assign row names
d <- as.data.frame(d)
row.names(d) <- d$SampleID

#For some analyses we need a table with only composition info
d.comp <- d[,8:29]
#and one with just explanatory variables
d.expl <- d[,1:5]


#Also need those by tissue
#skin 
d.comp.sk <- d.sk[,8:29]
d.expl.sk <- d.sk[,1:5]
#pulp
d.comp.pu <- d.pu[,8:28]
d.expl.pu <- d.pu[,1:5]
#Seed
d.comp.se <- d.se[,8:29]
d.expl.se <- d.se[,1:5]


#Orchard Level Data----------------------------------------
#These data sets will be condensed to the orchard level (24 observations)
#We're going to create four different data sets to analyze parts of the fruit: 
#Whole Fruit, Skin, Seed, and Pulp

###Orchard and Tree level data will be the same for all###
Orchard <- read_csv("Orchard Level Data.csv")
View(Orchard)

#convert weights based on averages and bag weights
Tree <- read_csv("Tree Level Data.csv")
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
  mutate(TotalPhen=rowSums(across(5:26)),
         PhenRich=rowSums(across(5:26)!=0))

Fruit0 <- Fruit0 %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)

CCD <- left_join(Tree, Orchard, by="orchard.num")
CCD <- left_join(CCD, Fruit0, by="orchard.num")
View(CCD)

#Skin---------------------------------------------------------------------------
#condense chem data by tissue and orchard 
Fruit1 <- read_csv("Fruit Level Data.csv")
Fruit1 <- Fruit1 %>%
  mutate(TotalPhen=rowSums(across(5:26)),
         PhenRich=rowSums(across(5:26)!=0))

Fruit1 <- Fruit1 %>%
dplyr::filter(Tissue != "PULP") %>%
dplyr::filter(Tissue != "SEED") 


Fruit1 <- Fruit1 %>%
group_by(orchard.num) %>%
summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)

View(Fruit1)

#combine all data sets to 120 values 
SkinD <- left_join(Tree, Orchard, by="orchard.num")
SkinD <- left_join(SkinD, Fruit1, by="orchard.num")
View(SkinD)

#Pulp---------------------------------------------------------------------------
Fruit2 <- read_csv("Fruit Level Data.csv")
Fruit2 <- Fruit2 %>%
  mutate(TotalPhen=rowSums(across(5:26)),
         PhenRich=rowSums(across(5:26)!=0))

Fruit2 <- Fruit2 %>%
  dplyr::filter(Tissue != "SKIN") %>%
  dplyr::filter(Tissue != "SEED") 


Fruit2 <- Fruit2 %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)

View(Fruit2)

PulpD <- left_join(Tree, Orchard, by="orchard.num")
PulpD <- left_join(PulpD, Fruit2, by="orchard.num")
PulpD <- c[-25,]
View(PulpD)
#Seed#--------------------------------------------------------------------------
Fruit3 <- read_csv("Fruit Level Data.csv")
Fruit3 <- Fruit3 %>%
  mutate(TotalPhen=rowSums(across(5:26)),
         PhenRich=rowSums(across(5:26)!=0))

Fruit3 <- Fruit3 %>%
  dplyr::filter(Tissue != "SKIN") %>%
  dplyr::filter(Tissue != "PULP") 

Fruit3 <- Fruit3 %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)

View(Fruit3)

#combine all data sets to 120 values 
SeedD <- left_join(Tree, Orchard, by="orchard.num")
SeedD <- left_join(SeedD, Fruit3, by="orchard.num")
SeedD <- c[-25,]
View(SeedD)




#Data Transformation------------------------------------------------------------
#Total Phenol data is skewed weird and yields a high Shapiro score 
#Need to transform this data so that we can work with it 

#first we'll look at the data 
shapiro.test(d$TotalPhen)
#W = 0.59801, p-value < 2.2e-16
hist(d$TotalPhen)
#skewed to 0 on the positive side 

#logit transformation 
d <- d %>%
  mutate(TotalPhentrans= log10(d$TotalPhen))
d$TotalPhentrans[d$TotalPhentrans==-Inf] <- NA
na.omit(d$TotalPhentrans)
view(d$TotalPhentrans)

shapiro.test(d$TotalPhentrans)
#W = 0.86448, p-value < 2.2e-16
hist(d$TotalPhentrans)
###didnt seem to do anything but were going to use this 



#Asin sqrt transformation 
d <- d %>%
  mutate(TotalPhenArc= sqrt(d$TotalPhen))
d$TotalPhenArc[d$TotalPhenArc==0] <- NA
na.omit(d$TotalPhenArc)
view(d$TotalPhenArc)

#shapiro test this 
shapiro.test(d$TotalPhenArc)
#W = 0.9281, p-value = 9.022e-12

#Histogram for the transformation 
hist(d$TotalPhenArc)
#did not make data normal 

#Inverse transformation
d <- d %>%
  mutate(TotalPhenInv= 1/(d$TotalPhen))
d$TotalPhenInv[d$TotalPhenInv==Inf] <- NA
na.omit(d$TotalPhenInv)
view(d$TotalPhenInv)

#shapiro test this 
shapiro.test(d$TotalPhenInv)
#W = 0.091764, p-value < 2.2e-16

#Histogram for the transformation 
hist(d$TotalPhenInv)
#did not make data normal 

#PCA----------------------------------------------------------------------------
install.packages("pls")
library(pls)
#second we're just going to look at predictor variables and their affect on chem
#create practice data set 

p1 <- d
view(p1)
#looking at just repsonse variables 
p_chem <- dplyr::select(p1, c(8:29))

#calculate principal components
results <- prcomp(p_chem, scale = TRUE)

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
       scale = TRUE, xlabs = rep("*", 358))

#calculate total variance explained by each principal component
summary(results)$importance
summary(results)$importance[2,]

var_explained = results$sdev^2 / sum(results$sdev^2)

df <- data.frame(PC=1:22, var_explained=var_explained)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


#correlation matrix-----------------------------------------------------------
#we're going to round out looking at the influence of our variables by creating a,
#correlation matrix. This will help to pin point which interactions to look at, 
#during analysis 

#were going to use the same practice set we used for PCA 
library(Hmisc)
#The first matrix shows the correlation coefficients between the variables  
#the second matrix shows the corresponding p-values.
rcorr(as.matrix(p_chem))
#any value below 0.05 is not a statistically significant relationship 
# need to sort through and organzie these values 

#now let's visualize this 
library(corrplot)
corrplot(cor(p_chem))


###again with data set C with all values 
view(c)
chem_clim <- dplyr::select(c, c(2:4, 7:26))

rcorr(as.matrix(chem_clim))
corrplot(cor(chem_clim))





###NMDS-------------------------------------------------------------------------
#someone online said to change the distance
m.NMDS <- metaMDS(d.comp,distance = "euclidean", k=2, trymax=100, autotransform =TRUE)
m.NMDS

#for plotting, need to add columns that give values for 
#colors and symbols we want in plot
d.expl$Color <- recode_factor(d.expl$orchard.type,
                              Organic="red", Conventional="blue")
d.expl$Symbol <- recode_factor(d.expl$Tissue,
                               SK=8, PU=1, SE=2)
d.expl$Symbol <- as.numeric(as.character(d.expl$Symbol))

d.expl$orchard.type = as.factor(d.expl$orchard.type)


plot(m.NMDS, type="n") #plots the ordination axes only
#cant get to plot 
points(m.NMDS, pch=d.expl$Symbol,
       col=as.character(d.expl$Color), cex = 0.8)     



#PERMANOVA can test whether the visualized differences are significant

m.perm <- adonis2(d.comp~orchard.type*Tissue, data=d.expl)
m.perm



#PERMANOVA can test whether the visualized differences are significant

m.perm <- adonis2(d.comp~Dom.Status*Tissue, data=d.expl)
m.perm

#Random Forest --------------------------------------------------
###skin
m1.rf.sk <- randomForest(d.comp.sk,d.expl.sk$orchard.type, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.sk$importance
varImpPlot(m1.rf.sk)
MDSplot(m1.rf.sk, d.expl.sk$orchard.type)

#using boruta 
m1.rf.b.sk <- Boruta(d.comp.sk,d.expl.sk$orchard.type)
m1.rf.b.sk
#plot 
plot(m1.rf.b.sk)  
getSelectedAttributes(m1.rf.b.sk) #lists all important ones
#"pb1"         "syr.acid"    "epicatechin"

#performing MANOVA 
d.comp.sk.sel <- data.matrix(d.comp.sk[,getSelectedAttributes(m1.rf.b.sk)])
m1.man.sk <- manova(d.comp.sk.sel ~ d.expl.sk$orchard.type)
summary(m1.man.sk) 
#overall significance for MANOVA p= 0.2883

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.sk)  
#pb1 = 0.3935

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
#important variables: chl.acid, PhenRich, TotalPhen, U1
#lists all important ones
getSelectedAttributes(m1.rf.b.pu) 
#"chl.acid"  "U1"

##Running MANOVAS 
d.comp.pu.sel <- data.matrix(d.comp.pu[,getSelectedAttributes(m1.rf.b.pu)])
m1.man.pu <- manova(d.comp.pu.sel ~ d.expl.pu$orchard.type) #doesnt work for some reason 
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
#pb2"        "syr.acid"   "reynoutrin" "U2"         "PhenRich" 
#Running MANOVAS 
d.comp.se.sel <- data.matrix(d.comp.se[,getSelectedAttributes(m1.rf.b.se)])
m1.man.se <- manova(d.comp.se.sel ~ d.expl.se$orchard.type)
summary(m1.man.se)  #overall significance for MANOVA
summary.aov(m1.man.se)  
#p=0.02016 *

par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.se.sel))){
  d.temp=d.comp.se.sel[,i]
  plot(d.temp ~ d.expl.se$orchard.type, ylab=colnames(d.comp.se.sel)[i])
}
dev.off()



#Question 1---------------------------------------------------------------------
#How do management systems interact with broad variation in abiotic conditions across latitude to shape:
#a) fruit chemical composition, b) fruit chemical diversity

###Q1A: total phenolics###
tp1 <- glmmTMB(TotalPhentrans~ orchard.type*Tissue + (1|site.code/Orchard.num/Tree), data=d)
summary(tp1)
Anova(tp1)
#orchard.type          1.6935  1    0.19314    
#Tissue              321.4999  2    < 2e-16 ***
#orchard.type:Tissue   5.1914  2    0.07459 . 

diagnose(tp1)
shapiro.test(resid(tp1))
hist(resid(tp1))

#residuals not good
#data needs to be transformed using either a logit or arcsin method 



#total p per tissue type 
tp.p <- glmmTMB(TotalPhentrans ~ orchard.type + (1|site.code/Orchard.num/Tree), data=d.pu)
summary(tp.p)
Anova(tp.p)

tp.sk <- glmmTMB(TotalPhentrans ~ orchard.type + (1|site.code/Orchard.num/Tree), data=d.sk)
summary(tp.sk)
Anova(tp.sk)

tp.se <- glmmTMB(TotalPhentrans ~ orchard.type + (1|site.code/Orchard.num/Tree), data=d.se)
summary(tp.se)
Anova(tp.se)
#orchard.type 7.3822  1   0.006587 **


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


#How do management systems interact with broad abiotic conditions
#Q1A
##latitude##
tpl1 <- glmmTMB(TotalPhen~ orchard.type*Latitude + (1|site.code/orchard.num), data=c)
summary(tpl1)
Anova(tpl1)

##elevation## 
tpe1 <- glmmTMB(TotalPhen~ orchard.type*elevation + (1|site.code/orchard.num), data=c)
summary(tpe1)
Anova(tpe1)

##longitude## 
tplo1 <- glmmTMB(TotalPhen~ orchard.type*Longitude + (1|site.code/orchard.num), data=c)
summary(tplo1)
Anova(tplo1)


#Q1B 
##longitude##
prlo1 <- glmmTMB(PhenRich~ orchard.type*Longitutde + (1|site.code/Orchard.num), data=d)
summary(prlo1)
Anova(prlo1)

#per tissue 
pr.sk <- glmmTMB(PhenRich~ orchard.type*Longitutde + (1|site.code/Orchard.num), data=d.sk)
summary(pr.sk)
Anova(pr.sk)
#no sig 
pr.se <- glmmTMB(PhenRich~ orchard.type*Longitutde + (1|site.code/Orchard.num), data=d.se)
summary(pr.se)
Anova(pr.se)

pr.pu <- glmmTMB(PhenRich~ orchard.type*Longitutde + (1|site.code/Orchard.num), data=d.pu)
summary(pr.pu)
Anova(pr.pu)
#Longitutde              4.3082  1    0.03793 *



##latitude##
prl1 <- glmmTMB(PhenRich~ orchard.type*Latitude + (1|site.code/Orchard.num), data=d)
summary(prl1)
Anova(prl1)
#orchard.type:Latitude 2.8000  1    0.09426 .

#per tissue 
pr.sk1 <- glmmTMB(PhenRich~ orchard.type*Latitude + (1|site.code/Orchard.num), data=d.sk)
summary(pr.sk1)
Anova(pr.sk1)
#no sig 

pr.se1 <- glmmTMB(PhenRich~ orchard.type*Latitude + (1|site.code/Orchard.num), data=d.se)
summary(pr.se1)
Anova(pr.se1)
#orchard.type          6.8775  1   0.008729 **
#Latitude              5.0659  1   0.024401 * 

pr.pu1 <- glmmTMB(PhenRich~ orchard.type*Latitude + (1|site.code/Orchard.num), data=d.pu)
summary(pr.pu1)
Anova(pr.pu1)
#orchard.type:Latitude 2.9541  1    0.08566 .



pre1 <- glmmTMB(PhenRich~ orchard.type*elevation + (1|site.code/Orchard.num), data=d)
summary(pre1)
Anova(pre1)

#elevation              14.6850  1  0.0001271 ***

#Question 1 Figures------------------------------------------------------------

tphen <- ggplot(d, aes(x=Tissue, y=TotalPhen, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Total Phenolics (ug per g") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("pulp", "seeds", "skin"))
tphen



phenrich <- ggplot(d, aes(x=Tissue, y=PhenRich, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Phenolic Richness") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("pulp", "seeds", "skin"))
phenrich


#Question 2---------------------------------------------------------------------
#Q1B Which abiotic factors are the most important drivers of fruit chemistry
###total phenolics###
p1<- glmmTMB(TotalPhen ~ Szn.Total.Precip*orchard.type + (1|site.code/onum), data=c)
summary(pl)
anova(p1)

p2<- glmmTMB(TotalPhen ~ Prox.Water*orchard.type + (1|site.code/onum), data=c)
summary(p2)
anova(p2)


p3<- glmmTMB(TotalPhen ~ Szn.Temp.Avg*orchard.type + (1|site.code/onum), data=c)
summary(p3)
anova(p3)

#Question 2 Figures-------------------------------------------------------------
#Question 3---------------------------------------------------------------------
#Question 3 Figures-------------------------------------------------------------