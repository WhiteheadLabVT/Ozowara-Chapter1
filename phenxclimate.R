#phenolic analysis--------------------------------------------------------------
rm(list=ls()) # clears work space
#install packages 
library(ggplot2)
library(MASS)
library(reshape2)
library(dplyr)
library(randomForest)
library(glmmTMB) #mixed models
library(car) #mixed models summary tables
library(vegan) #multivariate stats
library(Boruta) #random forest models

#Data organization--------------------------------------------------------------
library(readr)
d <- read_csv("phenol_clim - Sheet1.csv")
View(d)

##restructure data 
d <- d %>%
  mutate(TotalPhen=rowSums(across(17:38)),
         PhenRich=rowSums(across(17:38)!=0))

d <- d %>% 
  mutate_at(10:15, as.factor)
View(d)

#seperate by tissue type 
d.sk <- filter(d, Tissue=="SKIN")
d.sk <- select(d.sk, -(15+which(colSums(d.sk[17:38], na.rm=TRUE) %in% 0)))
d.pu <- filter(d, Tissue=="PULP")
d.pu <- select(d.pu, -(15+which(colSums(d.pu[17:38], na.rm=TRUE) %in% 0)))
d.se <- filter(d, Tissue=="SEED")
d.se <- select(d.se, -(15+which(colSums(d.se[17:38], na.rm=TRUE) %in% 0)))


d <- as.data.frame(d)
row.names(d) <- d$SampleID

#For some analyses we need a table with only composition info
#and one with just explanatory variables
d.comp <- d[,17:38]
d.expl <- d[,c(1:15, 39:40)]


#Also need those by tissue
###fix this 
d.comp.sk <- d.sk[,17:38]
d.expl.sk <- d.sk[,c(1:15, 37:38)]
d.comp.pu <- d.pu[,17:37]
d.expl.pu <- d.pu[,c(1:15, 36:37)]
d.comp.se <- d.se[,17:37]
d.expl.se <- d.se[,c(1:15, 37:38)]



#How do phenol differ across tissue types and mgmt systems----------------------
###total phenols###
tp1 <- glmmTMB(TotalPhen~ orchard.type*Tissue + (1|site.code/Orchard.num/Tree), data=d)
summary(tp1)
Anova(tp1)
#orchard.type        8.4285  1   0.003694 **
#Tissue              4.9692  2   0.083359 . 
#orchard.type:Tissue 0.3021  2   0.859807 

diagnose(tp1)
shapiro.test(resid(tp1))
hist(resid(tp1))
#residuals not good :(


#total p per tissue type 
tp.p <- glmmTMB(TotalPhen ~ orchard.type + (1|site.code/Orchard.num/Tree), data=d.pu, family=beta_family)
summary(tp.p)
Anova(tp.p)
#not sig 

tp.sk <- glmmTMB(TotalPhen ~ orchard.type + (1|site.code/Orchard.num/Tree), data=d.sk)
summary(tp.sk)
Anova(tp.sk)

#not sig 
tp.se <- glmmTMB(TotalPhen ~ orchard.type + (1|site.code/Orchard.num/Tree), data=d.se)
summary(tp.se)
Anova(tp.se)

p2 <- ggplot(d, aes(x=Tissue, y=TotalPhen, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Total Phenolics (ug per g") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("pulp", "seeds", "skin"))
p2


###phen rich ###
pr1 <- glmmTMB(PhenRich~ orchard.type*Tissue + (1|site.code/Orchard.num/Tree), data=d)
summary(pr1)
Anova(pr1)

#orchard.type           2.2590  1    0.16471    
#Tissue              1279.1586  2    < 2e-16 ***
#orchard.type:Tissue    5.4814  2    0.08499 .   

diagnose(pr1)
shapiro.test(resid(pr1))
hist(resid(pr1))
#= 0.98958, p-value = 0.01156


#phen rich per tissue type 
pr.p <- lm(PhenRich ~ orchard.type, data=d.pu)
summary(pr.p)
Anova(pr.p)
#not sig 

pr.sk <- lm(PhenRich ~ orchard.type, data=d.sk)
summary(pr.sk)
Anova(pr.sk)
#not sig 
pr.sd <- lm(PhenRich ~ orchard.type, data=d.se)
summary(pr.sd)
Anova(pr.sd)
#orchard.type 0.02436 *


p1 <- ggplot(d, aes(x=Tissue, y=PhenRich, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Phenolic Richness") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("pulp", "seeds", "skin"))
p1

#Q1A How do management systems interact with broad abiotic conditions to shape fruit quality------
###total phenolics x climate###

##latitude##
tpl1 <- glmmTMB(TotalPhen~ orchard.type*Latitude + (1|site.code/Orchard.num), data=d)
summary(tpl1)
Anova(tpl1)

#orchard.type          0.0060  1     0.9382
#Latitude              2.2062  1     0.1375
#orchard.type:Latitude 0.0001  1     0.9928

#per tissue 
tp.sk1 <- glmmTMB(TotalPhen~ orchard.type*Latitude + (1|site.code/Orchard.num), data=d.sk)
summary(tp.sk1)
Anova(tp.sk1)

tp.se1 <- glmmTMB(TotalPhen~ orchard.type*Latitude + (1|site.code/Orchard.num), data=d.se)
summary(tp.se1)
Anova(tp.se1)
#orchard.type          8.8497  1   0.002931 **

tp.pu1 <- glmmTMB(TotalPhen~ orchard.type*Latitude + (1|site.code/Orchard.num), data=d.pu)
summary(tp.pu1)
Anova(tp.pu1)


##elevation## 
tpe1 <- glmmTMB(TotalPhen~ orchard.type*elevation + (1|site.code/Orchard.num), data=d)
summary(tpe1)
Anova(tpe1)
#orchard.type           2.9601  1    0.08534 .
#elevation              5.7922  1    0.01610 *
#orchard.type:elevation 1.3073  1    0.25289  

#per tissue 
tp.sk2 <- glmmTMB(TotalPhen~ orchard.type*elevation + (1|site.code/Orchard.num), data=d.sk)
summary(tp.sk2)
Anova(tp.sk2)
#elevation              9.5094  1   0.002044 **

tp.se2 <- glmmTMB(TotalPhen~ orchard.type*elevation + (1|site.code/Orchard.num), data=d.se)
summary(tp.se2)
Anova(tp.se2)
#elevation              31.2623  1  2.254e-08 ***
  

tp.pu2 <- glmmTMB(TotalPhen~ orchard.type*elevation + (1|site.code/Orchard.num), data=d.pu)
summary(tp.pu2)
Anova(tp.pu2)

#elevation              17.7032  1  2.582e-05 ***
#orchard.type:elevation  3.0118  1    0.08266 .  



##longitude## 
tplo1 <- glmmTMB(TotalPhen~ orchard.type*Longitutde + (1|site.code/Orchard.num), data=d)
summary(tplo1)
Anova(tplo1)
#orchard.type            0.0047  1     0.9453
#Longitutde              0.0158  1     0.8998
#orchard.type:Longitutde 0.0000  1     0.9992


#per tissue 
tp.sk3 <- glmmTMB(TotalPhen~ orchard.type*Longitutde + (1|site.code/Orchard.num), data=d.sk)
summary(tp.sk3)
Anova(tp.sk3)

tp.se3 <- glmmTMB(TotalPhen~ orchard.type*Longitutde + (1|site.code/Orchard.num), data=d.se)
summary(tp.se3)
Anova(tp.se3)
#orchard.type            9.2674  1   0.002333 **


tp.pu3 <- glmmTMB(TotalPhen~ orchard.type*Longitutde + (1|site.code/Orchard.num), data=d.pu)
summary(tp.pu3)
Anova(tp.pu3)


###phen rich x climate### 

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

###corrolation matrix-----------------------------------------------------------

#edit to do a true corrolation matrix with true numeric values 

library(corrplot)
d.cor = cor(d)
corrplot(d.cor)
library("Hmisc")
#creates the corrolation matrix 
d.rcorr = rcorr(as.matrix(d))
d.rcorr
#generates p-values between varibales 
d.coeff = d.rcorr$r
d.p = d.rcorr$P #use this corroplot

#corrplot will show the significance of relationships between varibales.
#look for colored tiles on the dependent variables 
corrplot(d.p)
###theres a lot going on :(





















###NMDS-------------------------------------------------------------------------

m.NMDS <- metaMDS(d.comp,distance = "bray", trymax=100, autotransform =TRUE)
m.NMDS



#for plotting, need to add columns that give values for 
#colors and symbols we want in plot
d.expl$Color <- recode_factor(d.expl$Dom.Status,
                              domesticated="#3EBCD2", wild="#9A607F")
d.expl$Symbol <- recode_factor(d.expl$Tissue,
                               SK=8, PU=1, SE=2)
d.expl$Symbol <- as.numeric(as.character(d.expl$Symbol))


plot(m.NMDS, type="n") #plots the ordination axes only
points(m.NMDS, pch=d.expl$Symbol,
       col=as.character(d.expl$Color), cex = 0.8)     
#text(m.NMDS, pos = 4, cex = 0.3, display = "sites")  #if you want to see sample names


#PERMANOVA can test whether the visualized differences are significant

m.perm <- adonis2(d.comp~Dom.Status*Tissue, data=d.expl)
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



#Q1B Which abiotic factors are the most important drivers of fruit chemistry----
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
