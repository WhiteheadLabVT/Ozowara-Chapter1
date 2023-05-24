#phenolic analysis--------------------------------------------------------------
rm(list=ls()) # clears work space
#install packages 
library(ggplot2)
library(MASS)
library(reshape2)
library(tidyverse)
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
  mutate(TotalPhen=rowSums(across(where(is.numeric) & !trans.cinnamic.acid)),
         PhenRich=rowSums(across(where(is.numeric) & !trans.cinnamic.acid & !TotalPhen)!=0))
d <- d %>% 
  mutate_at(c("Plant", "Tissue", "orchard.type"), as.factor)

d$site.code <- d$site.code %>% as.factor() %>% as.numeric()
d$orchard.type <- d$orchard.type %>% as.factor() %>% as.numeric()
d$Tissue <- d$Tissue %>% as.factor() %>% as.numeric()
d$Plant <- d$Plant %>% as.factor() %>% as.numeric()

ga=  d$`gallic acid`
#REMOVE COLLUMNS WITH NO VALUES 
d %>% select(-ga)

d %>% select(-glucoside)

view(d)

#Divide by tissue type, for some analyses we will want to look at them separate
d.sk <- filter(d, Tissue=="SKIN")
d.sk <- select(d.sk, -(13+which(colSums(d.sk[14:36], na.rm=TRUE) %in% 0)))
d.pu <- filter(d, Tissue=="PULP")
d.pu <- select(d.pu, -(13+which(colSums(d.pu[14:36], na.rm=TRUE) %in% 0)))
d.se <- filter(d, Tissue=="SEED")
d.se <- select(d.se, -(13+which(colSums(d.se[14:36], na.rm=TRUE) %in% 0)))

#For some analyses we need a table with only composition info
#and one with just explanatory variables

d.comp <- d[,15:36]
row.names(d.comp) <- d$Orchard_Num


d.expl <- d[,1:14]
row.names(d.expl) <- d$Orchard_Num

#Also need those by tissue

d.comp.sk <- d.sk[,15:33]
row.names(d.comp.sk) <- d.sk$Orchard_Num

d.expl.sk <- d.sk[,1:14]
row.names(d.expl.sk) <- d.sk$Orchard_Num

d.comp.pu <- d.pu[,15:33]
row.names(d.comp.pu) <- d.pu$Orchard_Num

d.expl.pu <- d.pu[,1:14]
row.names(d.expl.pu) <- d.pu$Orchard_Num

d.comp.se <- d.se[,15:33]
row.names(d.comp.se) <- d.se$Orchard_Num

d.expl.se <- d.se[,1:14]
row.names(d.expl.se) <- d.se$Orchard_Num

###corrolation matrix-----------------------------------------------------------
library(corrplot)
d.cor = cor(d)
corrplot(d.cor)
install.packages("Hmisc")
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




















###-----------------------------------------------------------------------------

#Q2: Does overall chemical composition differ across tissues and
# depending on domestication status
#----------------------

#First, use an NMDS to visualize differences across sample types



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
#important variables: quercetin 
getSelectedAttributes(m1.rf.b.sk) #lists all important ones
#procyanidin B1" "syringic acid"  "epicatechin"   

#performing MANOVA 
d.comp.sk.sel <- data.matrix(d.comp.sk[,getSelectedAttributes(m1.rf.b.sk)])
m1.man.sk <- manova(d.comp.sk.sel ~ d.expl.sk$orchard.type)
summary(m1.man.sk) 
#overall significance for MANOVA p= 0.2883

#follow-up ANOVAs for each individual compound
summary.aov(m1.man.sk)  
#pb1 = 0.06155 .

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
#important variables: coumaric and reyonutrin 
#lists all important ones
getSelectedAttributes(m1.rf.b.pu) 
#phloridzin

##Running MANOVAS 
d.comp.pu.sel <- data.matrix(d.comp.pu[,getSelectedAttributes(m1.rf.b.pu)])
m1.man.pu <- manova(d.comp.pu.sel ~ d.expl.pu$orchard.type) #doesnt work for some reason 
summary(m1.man.pu)  #overall significance for MANOVA
summary.aov(m1.man.pu)  #follow-up ANOVAs for each individual compound


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
#"gentistic acid" "procyanidin B2" "syringic acid"  "epicatechin"    "reynoutrin"  

#Running MANOVAS 
d.comp.se.sel <- data.matrix(d.comp.se[,getSelectedAttributes(m1.rf.b.se)])
m1.man.se <- manova(d.comp.se.sel ~ d.expl.se$orchard.type)
summary(m1.man.se)  #overall significance for MANOVA
#P= 0.0008877 ***
summary.aov(m1.man.se)  


par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.se.sel))){
  d.temp=d.comp.se.sel[,i]
  plot(d.temp ~ d.expl.se$orchard.type, ylab=colnames(d.comp.se.sel)[i])
}
dev.off()


