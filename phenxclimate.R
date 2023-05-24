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
d <- read_csv("~/Dissertation Data/phenol_clim - Sheet1.csv")
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

#Random Forest Analysis---------------------------------------------------------
#change collumn names for analysis 
names(d)
names(d)[names(d) == 'site.code'] <- 'site.code'
names(d)[names(d) == 'Tissue'] <- 'Tissue'
names(d)[names(d) == 'gentistic acid'] <- 'gent.acid'
names(d)[names(d) == 'glucoside'] <- 'glucoside'
names(d)[names(d) == 'epicatechin'] <- 'epicatechin'
names(d)[names(d) == 'hyperin'] <- 'hyperin'
names(d)[names(d) == 'quercitrin'] <- 'quercitrin'
names(d)[names(d) == 'lat'] <- 'lat'
names(d)[names(d) == 'PhenRich'] <- 'PhenRich'
names(d)[names(d) == 'orchard.type'] <- 'orchard.type'
names(d)[names(d) == 'trans.cinnamic.acid'] <- 'trans.cinnamic.acid'
names(d)[names(d) == 'catechin'] <- 'catechin'
names(d)[names(d) == 'caffeic acid??'] <- 'caf.acid'
names(d)[names(d) == 'p-coumaric'] <- 'p.coumaric'
names(d)[names(d) == 'isoquercitin'] <- 'isoquercitin'
names(d)[names(d) == 'phloridzin'] <- 'phloridzin'
names(d)[names(d) == 'elev'] <- 'elev'
names(d)[names(d) == 'Orchard_Num'] <- 'Orchard_Num'
names(d)[names(d) == 'gallic acid'] <- 'gal.acid'
names(d)[names(d) == 'cyanidin galactoside'] <- 'galactoside'
names(d)[names(d) == 'procyanidin B2'] <- 'pb2'
names(d)[names(d) == 'ferulic'] <- 'ferulic'
names(d)[names(d) == 'reynoutrin'] <- 'reynoutrin'
names(d)[names(d) == 'quercetin'] <- 'quercetin'
names(d)[names(d) == 'procyanidin B1'] <- 'pb1'
names(d)[names(d) == 'chlorogenic acid'] <- 'chl.acid'
names(d)[names(d) == 'syringic acid'] <- 'syr.acid'
names(d)[names(d) == 'rutin'] <- 'rutin'
names(d)[names(d) == 'avicularin'] <- 'avicularin'
names(d)[names(d) == 'phloretin'] <- 'phloretin'
names(d)[names(d) == 'TotalPhen'] <- 'TotalPhen'


#make this example reproducible
set.seed(1)

#fit the random forest model
model <- randomForest(
  formula = TotalPhen ~ .,
  data = d
)

#display fitted model
model

#find number of trees that produce lowest test MSE
which.min(model$mse)

#find RMSE of best model
sqrt(model$mse[which.min(model$mse)]) 

#plot the test MSE by number of trees
plot(model)


#produce variable importance plot
varImpPlot(model) 


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



m.NMDS <- metaMDS(d.comp, distance = "bray", trymax=100, autotransform =TRUE)
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

#Random Forest Susan's Method --------------------------------------------------
#random forest can only handle one set of predictors, so we will do
#this separately for each tissue type

#skin--------
m1.rf.sk <- randomForest(d.comp.sk,d.expl.sk$Dom.Status, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.sk$importance
varImpPlot(m1.rf.sk)
MDSplot(m1.rf.sk, d.expl.sk$Dom.Status)

#Boruta is a wrapper algorithm for random forest for feature selection
#It determines which variables are "important" for distinguishing
#among groups by comparing each variable to its "shadow attribute"
#which is basically just taking each variable and shuffling it across
#wild and domesticated, so it could no longer be important

m1.rf.b.sk <- Boruta(d.comp.sk,d.expl.sk$Dom.Status)
m1.rf.b.sk
plot(m1.rf.b.sk)  #important variables (better than shadow) are in green
getSelectedAttributes(m1.rf.b.sk) #lists all important ones


#To follow up feature selection, we typically perform a MANOVA
#on all selected attributes to determine if they are significantly
#different between wild and domesticated

d.comp.sk.sel <- data.matrix(d.comp.sk[,getSelectedAttributes(m1.rf.b.sk)])
m1.man.sk <- manova(d.comp.sk.sel ~ d.expl.sk$Dom.Status)
summary(m1.man.sk)  #overall significance for MANOVA
summary.aov(m1.man.sk)  #follow-up ANOVAs for each individual compound


#some quick plots of all of them
#note most compounds are higher in wild but a few are higher in domesticated
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.sk.sel))){
  d.temp=d.comp.sk.sel[,i]
  plot(d.temp ~ d.expl.sk$Dom.Status, ylab=colnames(d.comp.sk.sel)[i])
}
dev.off()


#pulp-------

m1.rf.pu <- randomForest(d.comp.pu,d.expl.pu$Dom.Status, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.pu$importance
varImpPlot(m1.rf.pu)
MDSplot(m1.rf.pu, d.expl.pu$Dom.Status)

#Boruta is a wrapper algorithm for random forest for feature selection
#It determines which variables are "important" for distinguishing
#among groups by comparing each variable to its "shadow attribute"
#which is basically just taking each variable and shuffling it across
#wild and domesticated, so it could no longer be important

m1.rf.b.pu <- Boruta(d.comp.pu,d.expl.pu$Dom.Status)
m1.rf.b.pu
plot(m1.rf.b.pu)  #important variables (better than shadow) are in green
getSelectedAttributes(m1.rf.b.pu) #lists all important ones


#To follow up feature selection, we typically perform a MANOVA
#on all selected attributes to determine if they are significantly
#different between wild and domesticated

d.comp.pu.sel <- data.matrix(d.comp.pu[,getSelectedAttributes(m1.rf.b.pu)])
m1.man.pu <- manova(d.comp.pu.sel ~ d.expl.pu$Dom.Status)
summary(m1.man.pu)  #overall significance for MANOVA
summary.aov(m1.man.pu)  #follow-up ANOVAs for each individual compound


#some quick plots of all of them
#note most compounds are higher in wild but a few are higher in domesticated
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.pu.sel))){
  d.temp=d.comp.pu.sel[,i]
  plot(d.temp ~ d.expl.pu$Dom.Status, ylab=colnames(d.comp.pu.sel)[i])
}
dev.off()


#seeds--------------

m1.rf.se <- randomForest(d.comp.se,d.expl.se$Dom.Status, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.se$importance
varImpPlot(m1.rf.se)
MDSplot(m1.rf.se, d.expl.se$Dom.Status)

#Boruta is a wrapper algorithm for random forest for feature selection
#It determines which variables are "important" for distinguishing
#among groups by comparing each variable to its "shadow attribute"
#which is basically just taking each variable and shuffling it across
#wild and domesticated, so it could no longer be important

m1.rf.b.se <- Boruta(d.comp.se,d.expl.se$Dom.Status)
m1.rf.b.se
plot(m1.rf.b.se)  #important variables (better than shadow) are in green
getSelectedAttributes(m1.rf.b.se) #lists all important ones


#To follow up feature selection, we typically perform a MANOVA
#on all selected attributes to determine if they are significantly
#different between wild and domesticated

d.comp.se.sel <- data.matrix(d.comp.se[,getSelectedAttributes(m1.rf.b.se)])
m1.man.se <- manova(d.comp.se.sel ~ d.expl.se$Dom.Status)
summary(m1.man.se)  #overall significance for MANOVA
summary.aov(m1.man.se)  #follow-up ANOVAs for each individual compound


#some quick plots of all of them
#note most compounds are higher in wild but a few are higher in domesticated
par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.se.sel))){
  d.temp=d.comp.se.sel[,i]
  plot(d.temp ~ d.expl.se$Dom.Status, ylab=colnames(d.comp.se.sel)[i])
}
dev.off()


