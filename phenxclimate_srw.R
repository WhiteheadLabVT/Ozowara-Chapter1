#SRW: Great start to this!!
#Next steps:
#1: review the changes I made to the data cleaning/organization
#2: review xls sheet--need to double check data, there are a lot of samples
#with zero values all the way across (i.e. TotalPhen=0) which should not
#be the case. 
#3: update models for Total Phen to try with the beta-distribution (see example
#and notes below). Right now the first one I ran doesn't look great but 
#hopefully fixing the zeros in the dataset might help that. If not I think we 
#should try a transformation of the total phenolics variable
#either a logit transform of the proportion, or an arcsin sqrt transform of the %



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
d <- read_csv("phenol_updated - Sheet1.csv")
View(d)



##restructure data 

#SRW: I changed this--it was summing across all numeric columns to get totalphen
#but you have a lot of other numeric variables you don't want to do that with
#e.g. elevation, fruit size
d <- d %>%
  mutate(TotalPhen=rowSums(across(15:36)),
         PhenRich=rowSums(across(15:36)!=0))


#SRW: you want all the tree ID, orchard ID etc, variables to be factors
d <- d %>% 
  mutate_at(9:13, as.factor)
view(d)

#SRW: suggest you add a "SampleID" column, which will a unique code for each sample
d <- d%>%
  mutate(SampleID=paste(Orchard.num, Tree, Tissue, sep="_"))

#SRW: also suggest a unique code for each tree ID, to avoid the analyses treating
#tree 1 from one orchard the same as tree 1 from another
d <- d %>%
  mutate(Tree=paste(Orchard.num, Tree, sep="_"))


#SRW: don't do this--these should be factors
#convert into numeric frames 
#d$site.code <- d$site.code %>% as.factor() %>% as.numeric()
#d$orchard.type <- d$orchard.type %>% as.factor() %>% as.numeric()

#SRW: changed from 13+ and 14:36
#Divide by tissue type, for some analyses we will want to look at them separate
d.sk <- filter(d, Tissue=="SKIN")
d.sk <- select(d.sk, -(14+which(colSums(d.sk[15:36], na.rm=TRUE) %in% 0)))
d.pu <- filter(d, Tissue=="PULP")
d.pu <- select(d.pu, -(14+which(colSums(d.pu[15:36], na.rm=TRUE) %in% 0)))
d.se <- filter(d, Tissue=="SEED")
d.se <- select(d.se, -(14+which(colSums(d.se[15:36], na.rm=TRUE) %in% 0)))


#SRW: don't do this, tissue and tree should be factors
#convert ro numeric after creating filters
#d$Tissue <- d$Tissue %>% as.factor() %>% as.numeric()
#d$Tree <- d$Tree %>% as.factor() %>% as.numeric()



#SRW: assign rownames once here to avoid having to do it over and over below
#use the sampleID (not orchard ID) so you have a unique code for each sample 
#you can link across composition and explanatory data

d <- as.data.frame(d)
row.names(d) <- d$SampleID


#For some analyses we need a table with only composition info
#and one with just explanatory variables
#SRW: added columns 37:39 to all your explanatory variables for full datasets, 
#SRW: updated column IDs below for each of the skin, pulp, seed datasets
#they all have different numbers of columns because we deleted the zero compounds
#above, so you need to make sure it is using the right column IDs
#SRW: got rid of separate rowname assignments since that was done above

d.comp <- d[,15:36]
d.expl <- d[,c(1:14, 37:39)]


#Also need those by tissue
d.comp.sk <- d.sk[,15:34]
d.expl.sk <- d.sk[,c(1:14, 35:37)]
d.comp.pu <- d.pu[,15:32]
d.expl.pu <- d.pu[,c(1:14, 33:35)]
d.comp.se <- d.se[,15:31]
d.expl.se <- d.se[,c(1:14, 32:34)]


#SRW: your values for total phen and phen rich were summing across all
#numeric columns so they were off, so all your numbers
#are going to be different below
#SRW: also a lot of the factor variables were being treated as numeric
#so that will also have a big impact on your results
#SRW: for your models where you have the different tissue types together
#you all need to include treeID as part of your random effects
#structure of your model (because you have multiple samples per tree)

#How do phenol differ across tissue types and mgmt systems----------------------
###total phenols###
tp1 <- glmmTMB(TotalPhen~ orchard.type*Tissue + (1|site.code/Orchard.num/Tree), data=d)
summary(tp1)
Anova(tp1)

#orchard.type         4.0456  1    0.04429 *  
#Tissue              91.2389  2    < 2e-16 ***
#orchard.type:Tissue  4.0761  2    0.13028   

diagnose(tp1)
shapiro.test(resid(tp1))
hist(resid(tp1))

#SRW: these model are not good, the residuals are way way deviating from 
#normal
hist(d$TotalPhen)
#SRW: the response variable is very far from normal also, we will need to use
#a different distribution (beta distribution) or transform


###beta-distribution model
#beta-distribution only works for proportions, so, before running, you need to 
#your concentration data to a proportion dry weight
#if right now it is ug/g, you would divide by 1,000,000 to get g/g
#Also you have a lot of zero values?  Need to double check this, I would not
#expect zeros in this dataset. For now I am adding a small constant because
#the model won't run with zeros
tp2 <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Tissue + (1|site.code/Orchard.num/Tree), 
               data=d, family=beta_family(link="logit"))
summary(tp2)
Anova(tp2)

diagnose(tp2)
shapiro.test(resid(tp2))
hist(resid(tp2))  
#still really not great! We should consider a transformation. But first lets 
#double check all the data, update zeros, etc. I think there should not be
#so many zeros and that will affect the distribution a lot



#SRW: you would also want site code and orchard number in these models
#you can take out tree ID once you split by tissue since now you only have
#one sample per tree

#SRW: these models below also all have non-normal residuals, I would suggest 
#first checking the data, then re-running the total model above, then updating
#all these models below to include the random effects (like I did for tp.p as 
#an example), and the beta distribution (like I did for the totalPhen as an example)

#total p per tissue type 
tp.p <- glmmTMB(TotalPhen ~ orchard.type + (1|site.code/Orchard.num), data=d.pu)
summary(tp.p)
Anova(tp.p)
#not sig 

tp.sk <- lm(TotalPhen ~ orchard.type, data=d.sk)
summary(tp.sk)
Anova(tp.sk)
#not sig 
tp.sd <- lm(TotalPhen ~ orchard.type, data=d.se)
summary(tp.sd)
Anova(tp.sd)
#orchard.type 0.01425 *

p2 <- ggplot(d, aes(x=Tissue, y=TotalPhen, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Total Phenolics (ug per g") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("pulp", "seeds", "skin"))
p2


###phen rich ###
pr1 <- glmmTMB(PhenRich~ orchard.type*Tissue + (1|site.code/Orchard.num), data=d)
summary(pr1)
Anova(pr1)

#orchard.type           2.2590  1    0.13284    
#Tissue              1279.1586  2    < 2e-16 ***
#orchard.type:Tissue    5.4814  2    0.06453 . 

diagnose(pr1)
shapiro.test(resid(pr1))
hist(resid(pr1))


#SRW: update these models to include the random effects of site and orchard
#also check some diagnostics as above, with diagnose and plotting the 
#histogram of the residuals
#if they look good you are fine, but sometimes for richness data a poisson
#distribution might work better because it is discrete count data

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
#orchard.type 0.02016 *


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


