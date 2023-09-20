#AMAZING, it is really starting to come together!!! Here's my biggest comments

#1) For Q3 and Q4 (and maybe Q2 as well), I think we should include orchard type 
#as another
#explanatory factor (since this is such a key part of the study design) in the 
#main model. This would just be an additive effect, unless you have a strong
#hypothesis for how orchard.type should interact with a particular management
#factor (e.g. herbicides), in which case you include just that two-way interaction.
#Then instead of following those big models with the interactive models you usually 
#did that tested interactions between orchard.type and other significant factors, 
#I would suggest the big overall model is followed by some model simplification 
#(see comments below) that narrow down which variables are important, then you 
#report the p-values and parameter estimates from the simplified version

#2) I really like the addition of Q4, but it did make me wonder whether logically
#then it would also make sense to include some big picture measure of pest pressure
#(like the index you calculated) as part of Q1, just to make things tidy, so we first 
#ask broadly how management system, abiotic conditions, and pest pressure
#interact to shape fruit quality. Then we drill into each of those factors in 
#Q2, Q3, Q4. I did play around with that a little bit just to see what happens
#and that does have some effects. Could be worth including, or at least 
#playing around with a bit to see what it looks like and then we can decide whether 
#that makes sense to include in the paper or just leave with the two-way interaction
#between orchard.type*Latitude

#3) As we discussed when we met there was an error in the code so that it was
#not using the right data for total phenolics for skin, pulp, and seed. See notes
#below on that, it does have a big impact on the results but I think it makes
#more sense now 



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
library("purrr")
library(Hmisc) #correlation matrix
library(readr)
library(GGally)
#Read Data Organization and Restructuring---------------------------------------
#Read in Data sheets#

#Orchard Level Data (24 obs)
Orchard <- read_csv("Orchard Level Data.csv")
Orchard$orchard.num <- as.factor(as.character(Orchard$orchard.num))

#Tree Level Data (120 obs)
Tree <- read_csv("Tree Level Data.csv")
Tree$orchard.num <- as.factor(as.character(Tree$orchard.num))

#Phenolics Data (359 obs)
d <- read_csv("FruitLevelData_revised - Sheet1.csv")
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
#Whole Fruit Orchard Level Data (Will also be used for physical quality analysis)
Tree <- Tree %>%
  mutate(avgwgt=(apple.wgt.3-bag.weight)/3) 

Tree <- Tree %>% 
  dplyr::select(-c(5:6))

Tree.sum <- Tree %>%
  group_by(orchard.num) %>%
  summarise_at(c("avgwgt", "Firmness", "SSC", "maturity.index"), mean, na.rm = TRUE)

c <- d
c <- c %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)

c <- left_join(c, Tree.sum, by="orchard.num")
c <- left_join(c, Orchard, by="orchard.num")

##SRW: here you are averaging across all three tissue types? I am not sure that
#makes sense since they are very different. I would suggest keeping them 
#separate and creating one variable for each. 

#SRW: I wrote the below section before realizing that you had summarized them
#separately below for each tissue dataset, which accomplishes the same thing. 
#So this part is not necessary, just leaving it as an example of another 
#approach using pivot_wider to spread the data so the tissues are separate columns

c2 <- pivot_wider(d[,c(1:5, 42:43)], names_from="Tissue", values_from=c("TotalPhen", "PhenRich"))

c3 <- c2 %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen_SKIN", "PhenRich_SKIN", "TotalPhen_PULP", 
                 "PhenRich_PULP", "TotalPhen_SEED", "PhenRich_SEED")
               , mean, na.rm = TRUE)

c <- left_join(c, c3, by="orchard.num")

#SRW: though then they look a lot less normal...

hist(c$TotalPhen_SKIN)
hist(c$TotalPhen_PULP)
hist(c$TotalPhen_SEED)

shapiro.test(c$TotalPhen)
#normal, W = 0.94176, p-value = 0.1785
shapiro.test(c$PhenRich)


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
###Physical Quality###

#binding latitude to the 120 obs data set#
#SRW: also added pest_index to this so we can play around with that
TreeLat <- left_join(Tree, Orchard[,c(2,4, 39)], by="orchard.num")

#SSC
Lat1<- glmmTMB(SSC ~ orchard.type*Latitude*Pest_Index + (1|site.code/orchard.num), data=TreeLat)
summary(Lat1)
Anova(Lat1) 
#orchard.type           0.5016  1     0.4788    
#Latitude              18.2814  1  1.906e-05 ***
#orchard.type:Latitude  0.2196  1     0.6393  

hist(resid(Lat1))  #looks great
diagnose(Lat1)

ggplot(c, aes(x=Latitude, y=SSC, color=Pest_Index))+
  geom_smooth(method = "lm") +
  geom_point()

#Firmness 
Lat2<- glmmTMB(Firmness ~ orchard.type*Latitude + (1|site.code/orchard.num), data=TreeLat)
summary(Lat2)
Anova(Lat2)
#orchard.type           0.0828  1     0.7735    
#Latitude              34.5700  1  4.112e-09 ***
#orchard.type:Latitude  0.9233  1     0.3366 

hist(resid(Lat2))  #looks great
diagnose(Lat2)

ggplot(c, aes(x=Latitude, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

#Average Weight 
Lat3<- glmmTMB(avgwgt ~ orchard.type*Latitude + (1|site.code/orchard.num), data=TreeLat)
summary(Lat3)
#orchard.typeOrganic           142.408     45.320   3.142  0.00168 **

Anova(Lat3) 
#orchard.type          0.1535  1   0.695206   
#Latitude              0.3713  1   0.542288   
#orchard.type:Latitude 9.7210  1   0.001822 **

#SRW: weird the orchard type looks so significant in the summary but does not
#come out in the Anova...don't really get that...

hist(resid(Lat3))  #looks great
diagnose(Lat3)
##strongest interaction 

#Maturity Index 
#using poisson distribution 
Lat4<- glmmTMB(maturity.index ~ orchard.type*Latitude + (1|site.code/orchard.num), 
               data=TreeLat)
summary(Lat4)
Anova(Lat4)  
#orchard.type          0.0047  1    0.94544  
#Latitude              3.8790  1    0.04889 *
#orchard.type:Latitude 0.9598  1    0.32724 

hist(resid(Lat4)) 
diagnose(Lat4)

ggplot(c, aes(x=Latitude, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

#Figure: Average Weight x Latitude 
plot1 = ggplot(TreeLat, aes(x=Latitude, y=avgwgt, color=orchard.type)) +
  geom_point() +
  ylab ("Average Weight (g)") +
  xlab ("Latitude")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
plot1
#organic weight decreases as latitude increases 
#conventional weight increases as latitude increases 



###Investigating Latitude and Average weight 
#Splitting latitude in high and low 
Tree_low <- filter(TreeLat, Latitude<42)
Tree_high <- filter(TreeLat, Latitude>42)

#Low
avg_low<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/orchard.num), data=Tree_low)
summary(avg_low)
Anova(avg_low) 

avg_hgh<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/orchard.num), data=Tree_high)
summary(avg_hgh)
Anova(avg_hgh)
#orchard.type 3.6794  1    0.05509 .


#Figure: Avg weight boxplots at high and low latitudes 
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

#SRW: what package is this? not working for me with what is listed at the top
multiplot(avglow, avghigh, cols=2)





#Q1-B: Fruit Chemistry-------------------------------------------------------------
#binding latitude to the 359 obs data set 
ChemLat <- left_join(d, Orchard[,c(2,4)], by="orchard.num")

#Total Phenolics with beta distribution 
tp1 <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Tissue*Latitude +
                 (1|site.code/orchard.num/Tree), data=ChemLat, family=beta_family(link="logit"))
summary(tp1)
Anova(tp1)
#orchard.type                   4.4968  1    0.03396 *  
#Tissue                       278.1257  2    < 2e-16 ***
#orchard.type:Tissue            8.3481  2    0.01539 *  

#visualize this 
ggplot(ChemLat, aes(x=Tissue, y=TotalPhen/1000, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#seeds and skin have higher values 
#seeds highest 

ggplot(ChemLat, aes(x=Latitude, y=TotalPhen/1000, color=Tissue))+
  geom_smooth(method = "lm") +
  geom_point()
#over latitude we see skin and seeds increase 
#pulp appears to decrease 

ggplot(ChemLat, aes(x=Latitude, y=TotalPhen/1000, color=Tissue, shape= orchard.type))+
  geom_smooth(method = "lm") +
  geom_point(size=1)


#SRW: for the final plots I would suggest using the proportion dry weight as units
#to match the analyses. 

#SRW: Also, to simplify for the presentation of results in the paper, you could
#just state than an initial set of models revealed strong effects of tissue
#and interactions between tissue and other variables, so you just report the
#results for each tissue individually


#strong effects of tissue and interaction between tissue and orchard type
#splitting by tissue
d.sk <- left_join(d.sk, Orchard[,c(2,4, 39)], by="orchard.num")
d.pu <- left_join(d.pu, Orchard[,c(2,4, 39)], by="orchard.num")
d.se <- left_join(d.se, Orchard[,c(2,4, 39)], by="orchard.num")

tp.sk <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Latitude*Pest_Index + (1|site.code/orchard.num), 
                 data=d.sk, family=beta_family(link="logit"))
summary(tp.sk)
Anova(tp.sk)
#nothing 

tp.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Latitude*Pest_Index + (1|site.code/orchard.num), 
                 data=d.pu, family=beta_family(link="logit"))
summary(tp.pu)
Anova(tp.pu)
#nothing 


tp.se <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Latitude*Pest_Index + (1|site.code/orchard.num/Tree), 
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

#Visualize this 
#latitude by tissue
ggplot(ChemLat, aes(x=Latitude, y=PhenRich, color=Tissue))+
  geom_smooth(method = "lm") +
  geom_point()
#highest richness in skin 

#otype
ggplot(ChemLat, aes(x=Tissue, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()


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

#visualize this 
ggplot(d.se, aes(x=Latitude, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#organic less rich at lower latitudes

#SRW: although it looks in the plot like there might be an interaction between
#orchard type and latitude for seeds it is not significant so you might not need this

###Investigating Latitude Interactions  
#Splitting latitude in high and low 
Seed_low <- filter(d.se, Latitude<42)
Seed_high <- filter(d.se, Latitude>42)

#Low
sd_low<- glmmTMB(PhenRich ~ orchard.type + (1|site.code), data=Seed_low)
summary(sd_low)
Anova(sd_low) 
#orchard.type 3.6992  1    0.05444 .

sd_hgh<- glmmTMB(PhenRich ~ orchard.type + (1|site.code), data=Seed_high)
summary(sd_hgh)
Anova(sd_hgh)
#orchard.type 0.0737  1      0.786

#low lats make up for higher richness 

ggplot(Seed_low, aes(x=Latitude, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher richness in conventional seeds at lower latitudes 


#Q1-C: Which compounds distinguish fruits raised in their respective management systems?----------
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
                               SKIN=8, PULP=1, SEED=2)
d.expl$Symbol <- as.numeric(as.character(d.expl$Symbol))

d.expl$orchard.type = as.factor(d.expl$orchard.type)


plot(m.NMDS, type="n") #plots the ordination axes only
points(m.NMDS, pch=d.expl$Symbol,
       col=as.character(d.expl$Color), cex = 0.8)     


#PERMANOVA can test whether the visualized differences are significant
m.perm <- adonis2(d.comp~orchard.type*Tissue, data=d.expl)
m.perm
#orchard.type          1    0.450 0.00386   2.5040  0.024 *  
#Tissue                2   51.790 0.44450 144.1093  0.001 ***
#orchard.type:Tissue   2    0.843 0.00724   2.3465  0.008 ** 
#Residual            353   63.431 0.54440                    
#Total               358  116.515 1.00000  

#SRW: cool


###Random Forest###
###skin
m1.rf.sk <- randomForest(d.comp.sk,d.expl.sk$orchard.type, importance=TRUE, 
                         proximity=TRUE, oob.prox=TRUE, ntree=2000)
m1.rf.sk$importance
varImpPlot(m1.rf.sk)
MDSplot(m1.rf.sk, d.expl.sk$orchard.type)  #SRW: getting an error here

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

#SRW: differences are small and marginal, but it seems for all three that are
#marginal the organic is higher



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
#par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.pu.sel))){
  d.temp=d.comp.pu.sel[,i]
  plot(d.temp ~ d.expl.pu$orchard.type, ylab=colnames(d.comp.pu.sel)[i])
}
dev.off()

#SRW: G is higher in conventional


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


par(mfrow=c(1,3))
for (i in 1:length(colnames(d.comp.se.sel))){
  d.temp=d.comp.se.sel[,i]
  plot(d.temp ~ d.expl.se$orchard.type, ylab=colnames(d.comp.se.sel)[i])
}
dev.off()

#2 of 3 higher in conventional...so overall a mix


#run GLMMs of the significant compounds against Latitude 

#SRW:  Hmmm, here I don't see the justification to look at the effects of 
#latitude just for those compounds selected by the random forest, since the
#random forest was just for orchard type
#If you really wanted to you could also do a random forest with latitude
#as a predictor to see which compounds change with latidue. This method
#was originally developed for categorical variables, but it does work 
#with continuous ones


#Skin 
A.sk <- glmmTMB(A~ orchard.type*Latitude + 
                   (1|site.code), data=d.sk)
summary(A.sk)
Anova(A.sk)
#orchard.type          3.2881  1    0.06978 .

#visualize this 
ggplot(d.sk, aes(x=Latitude, y=A/1000, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


PB1.sk <- glmmTMB( PB1~ orchard.type*Latitude + 
                     (1|site.code), data=d.sk)
summary(PB1.sk)
Anova(PB1.sk)
#orchard.type:Latitude 10.8982  1  0.0009626 ***

#visualize this 
ggplot(d.sk, aes(x=Latitude, y=PB1/1000, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


gala.sk <- glmmTMB( cyanidin_Galactoside~ orchard.type*Latitude + 
                      (1|site.code), data=d.sk)
summary(gala.sk)
Anova(gala.sk)
#orchard.type          3.2137  1    0.07302 .

#visualize this 
ggplot(d.sk, aes(x=Latitude, y=cyanidin_Galactoside, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


#Pulp
G.pu <- glmmTMB( G~ orchard.type*Latitude + 
                   (1|site.code), data=d.pu)
summary(G.pu)
Anova(G.pu)
#orchard.type           4.8776  1    0.02721 *  
#Latitude              19.7803  1  8.687e-06 ***
#orchard.type:Latitude  4.2743  1    0.03869 *  

#visualize this 
ggplot(d.pu, aes(x=Latitude, y=G/1000, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


#Seed
E.se <- glmmTMB( E~ orchard.type*Latitude + 
                   (1|site.code), data=d.se)
summary(E.se)
Anova(E.se)
#orchard.type          9.9955  1   0.001569 **

ggplot(d.se, aes(x=Latitude, y=E/1000, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


PB2.se <- glmmTMB( PB2~ orchard.type*Latitude + 
                     (1|site.code), data=d.se)
summary(PB2.se)
Anova(PB2.se)
#orchard.type          11.8700  1  0.0005705 ***
#Latitude               4.6970  1  0.0302156 *

ggplot(d.se, aes(x=Latitude, y=PB2/1000, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

epicatechin.se <- glmmTMB( ePicatechin~ orchard.type*Latitude + 
                             (1|site.code), data=d.se)
summary(epicatechin.se)
Anova(epicatechin.se)
#orchard.type          14.7469  1  0.0001229 ***
#Latitude               5.9931  1  0.0143623 *  

ggplot(d.se, aes(x=Latitude, y=ePicatechin/1000, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


U1.se <- glmmTMB( U1~ orchard.type*Latitude + 
                    (1|site.code), data=d.se)
summary(U1.se)
Anova(U1.se)
#orchard.type          8.2778  1   0.004013 **

ggplot(d.se, aes(x=Latitude, y=U1/1000, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


#Q2: Which abiotic factors are the most important drivers of fruit quality?---------
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
biplot(results,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


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

#PC's 5, 6, 7, 8, 9 are very low/ not worth looking at
#Q2-A: Physical Quality---------------------------------------------------------
##Linear models PC x quality 
pc_clim <- as.data.frame(results$x)
pc_clim <- cbind(pc_clim, c)

hist(c$Firmness)
hist(c$SSC)
hist(c$avgwgt)
hist(c$maturity.index)

###Firmness###

p1 <- glmmTMB(Firmness ~ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim)
summary(p1)
#PC1  7.51e-09 ***
#PC3  4.43 9.52e-06 ***  

#SRW: not clear the justification for these interactive models? I might just
#include orchard.type in the main model (no interactions, I don't think 
#you have the power to test all the interactions here and it's not clear to
#me there would be some a priori hypothesis that orchard type should interact
#with PC1 but not PC2, etc.) 

f_pc1 <- glmmTMB(Firmness ~ orchard.type*PC1 + (1|site.code), data=pc_clim)
summary(f_pc1)
#PC1                       1.8549     0.5581   3.324 0.000889 ***

f_pc3 <- glmmTMB(Firmness ~ orchard.type*PC3 + (1|site.code), data=pc_clim)
summary(f_pc3)
#PC3                        2.851      1.561   1.826   0.0679 .  


plot(Firmness ~ PC1, data=pc_clim)
plot(Firmness ~ PC3, data=pc_clim)

#So, to interpret this you go back to look at the loading on each PC axis.
results$rotation

###SSC###
P2 <- glmmTMB(SSC ~ + PC1 + PC2 + PC3 + PC4 + (1|site.code), data= pc_clim)
summary(P2)
#PC1         -1.00949    0.15340   -6.58 4.68e-11 ***
#PC2         -0.58972    0.27900   -2.11   0.0345 *  
#PC3          0.87079    0.36240    2.40   0.0163 *  


s_pc1 <- glmmTMB(SSC ~ orchard.type*PC1 + (1|site.code), data=pc_clim)
summary(s_pc1)
#PC1                     -1.10112    0.21770  -5.058 4.24e-07 ***

s_pc2 <- glmmTMB(SSC ~ orchard.type*PC2 + (1|site.code), data=pc_clim)
summary(s_pc2)
#nothing

s_pc3 <- glmmTMB(SSC ~ orchard.type*PC3 + (1|site.code), data=pc_clim)
summary(s_pc3)
#nothing

plot(SSC ~ PC1, data=pc_clim)  
plot(SSC ~ PC2, data=pc_clim)
plot(SSC ~ PC3, data=pc_clim)

#SRW: WOW, seems to be strongly related to temps, high temps=high sugar

results$rotation

###AVGWGT
p3 <- glmmTMB(avgwgt ~ + PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim)
summary(p3)
#PC2           -8.186      4.121  -1.986   0.0470 *  
#PC3           10.331      5.216   1.981   0.0476 *  
#PC4          -12.449      6.064  -2.053   0.0401 * 


w_pc2 <- glmmTMB(avgwgt ~ orchard.type*PC2 + (1|site.code), data=pc_clim)
summary(w_pc2)
#nothing 

w_pc3 <- glmmTMB(avgwgt ~ orchard.type*PC3 + (1|site.code), data=pc_clim)
summary(w_pc3)
#PC3                       19.002      8.778   2.165   0.0304 *  

w_pc4 <- glmmTMB(avgwgt ~ orchard.type*PC4 + (1|site.code), data=pc_clim)
summary(w_pc4)
#nothing 

plot(avgwgt ~ PC2, data=pc_clim)
plot(avgwgt ~ PC3, data=pc_clim)
plot(avgwgt ~ PC4, data=pc_clim)

results$rotation

###Maturity Index 
p4 <- glmmTMB(maturity.index ~ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim)
summary(p4)
#PC1         -0.26478    0.08360  -3.167  0.00154 ** 
#PC4         -0.43710    0.20314  -2.152  0.03142 *  


m_pc1 <- glmmTMB(maturity.index ~ orchard.type*PC1 + (1|site.code), data=pc_clim)
summary(m_pc1)
#PC1                     -0.27376    0.09559  -2.864  0.00419 ** 

m_pc4 <- glmmTMB(maturity.index ~ orchard.type*PC4 + (1|site.code), data=pc_clim)
summary(m_pc4)
#nothing

plot(maturity.index ~ PC1, data=pc_clim)
plot(maturity.index ~ PC4, data=pc_clim)

results$rotation

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



#Q2-B: Fruit Chemistry----------------------------------------------------------
#Whole Fruit#
#using pc_clim for this portion 

#TotalPhen
pc.c.tp <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                             (1|site.code), data=pc_clim, family=beta_family (link="logit"))

summary(pc.c.tp)
#PC2          0.128315   0.055966    2.29   0.0219 *   


pc.c.tp.pc2 <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*PC2 + 
                         (1|site.code), data=pc_clim, family=beta_family (link="logit"))
summary(pc.c.tp.pc2)
#orchard.typeOrganic     -0.26424    0.10203  -2.590   0.0096 ** 

plot(TotalPhen ~ PC2, data=pc_clim)
plot(TotalPhen ~ as.factor(orchard.type), data=pc_clim)

results$rotation

#SRW: this is pretty interesting! Suggests that once we account for the climate
#variables there is a strong affect of orchard type with conventional 
#having higher phenolics. BUT, its a little strange to interpret since the
#values are averaged across tissues

#SRW: Also, as above, I would have suggested including orchard type in the full
#model rather than doing it separately after with just the significant PC axes
#Still it turns out significant, see above

#Also, more generally, when you have a model where the effect of a predictor
#only shows up once you add in other fixed or random effects, it can be helpful
#to plot the least squared means 
#(i.e. the model predicted effects of orchard type after accounting for
#the effects of all the climate factors)
#you can do that with the emmeans package

library(emmeans)
glmm_emmeans <-emmeans(pc.c.tp.pc2, ~ orchard.type, type="response")
#SRW: need to look into what this error means, I am not sure...

glmm_emmeans
glmm_emmeans <- as.data.frame(glmm_emmeans)

#SRW: so the "response" and SE here can be interpreted as the model predicted
#mean and SE for each orchard type, assuming you held all the climate 
#variables constant. Then you could plot just those means and SE, probably 
#overlaid on top of the raw data points

#SRW: But, again, I am not sure we should include the total fruit average




#PhenRich
pc.c.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                     (1|site.code), data=pc_clim)
summary(pc.c.pr)
#PC1           0.3597     0.1899    1.89  0.05823 .  
#PC2           0.9112     0.3412    2.67  0.00758 ** 


pc.c.pr.pc1 <- glmmTMB(PhenRich ~ orchard.type*PC1 + 
                         (1|site.code), data=pc_clim)
summary(pc.c.pr.pc1)
#nothing 

pc.c.pr.pc2 <- glmmTMB(PhenRich ~ orchard.type*PC2 + 
                         (1|site.code), data=pc_clim)
summary(pc.c.pr.pc2)
#nothing 

plot(PhenRich ~ PC1, data=pc_clim)
plot(PhenRich ~ PC2, data=pc_clim)

results$rotation


#SRW: I realized part way through the section below that it was using the total 
#fruit average for all the analyses that were supposed to be on skin, pulp, and 
#seed (that is why all the model results were the same)



#SKIN#

#SRW: Here was the problem, when you bound these two datsets like this it ended
#up creating a dataframe that had two TotalPhen and PhenRich columns, and I guess
#it was just using the first one for the analysis. I'm surprised it even
#lets you do that in a dataframe object...

#pc.sk_clim <- cbind(pc_clim, SkinD)

#SRW: also, always better to use a join to combine two dataframes...if your rows
#were in different orders with a cbind you could get into trouble

pc.sk_clim <- left_join(pc_clim[,1:10], SkinD, by="orchard.num")


#TotalPhen
pc.sk.tp <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.sk_clim, family=beta_family (link="logit"))
summary(pc.sk.tp)
#PC2          0.128315   0.055966    2.29   0.0219 *  

plot(TotalPhen ~ as.factor(orchard.type), data=pc.sk_clim)

#SRW: No longer any effects for skin phenolics



pc.sk.tp.pc2 <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*PC2 + 
                          (1|site.code), data=pc.sk_clim, family=beta_family (link="logit"))
summary(pc.sk.tp.pc2)
#orchard.typeOrganic     -0.26424    0.10203  -2.590   0.0096 ** 
 
plot(TotalPhen ~ PC2, data=pc.sk_clim)

#PhenRich
pc.sk.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.sk_clim)
summary(pc.sk.pr)
#PC1           0.3597     0.1899    1.89  0.05823 .  
#PC2           0.9112     0.3412    2.67  0.00758 ** 


pc.sk.pr.pc1 <- glmmTMB(PhenRich ~orchard.type* PC1 + 
                          (1|site.code), data=pc.sk_clim)
summary(pc.sk.pr.pc1)
#nothing 

pc.sk.pr.pc2 <- glmmTMB(PhenRich ~orchard.type* PC2 + 
                          (1|site.code), data=pc.sk_clim)
summary(pc.sk.pr.pc2)
#nothing  

plot(PhenRich ~ PC1, data=pc.sk_clim)
plot(PhenRich ~ PC2, data=pc.sk_clim)

results$rotation



#PULP#

#SRW: updating as above
#pc.pu_clim <- cbind(pc_clim, PulpD)

pc.pu_clim <- left_join(pc_clim[,1:10], PulpD, by="orchard.num")




##TotalPhen###
pc.pu.tp <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim, family=beta_family (link="logit"))
summary(pc.pu.tp)
#PC2          0.128315   0.055966    2.29   0.0219 *  

#SRW: strong effects of PC2

pc.pu.tp.pc2 <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type* PC2 + 
                          (1|site.code), data=pc.pu_clim, family=beta_family (link="logit"))
summary(pc.pu.tp.pc2)
#orchard.typeOrganic     -0.26424    0.10203  -2.590   0.0096 ** 


plot(TotalPhen ~ PC2, data=pc.pu_clim)


#PhenRich 
pc.pu.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim)
summary(pc.pu.pr)
#PC1           0.3597     0.1899    1.89  0.05823 .  
#PC2           0.9112     0.3412    2.67  0.00758 ** 


pc.pu.pr.pc2 <- glmmTMB(PhenRich ~orchard.type* PC2 + 
                          (1|site.code), data=pc.pu_clim)
summary(pc.pu.pr.pc2)
#nothing 

pc.pu.pr.pc1 <- glmmTMB(PhenRich ~orchard.type* PC1 + 
                          (1|site.code), data=pc.pu_clim)
summary(pc.pu.pr.pc1)
#nothing 

plot(PhenRich ~ PC1, data=pc.pu_clim)
plot(PhenRich ~ PC2, data=pc.pu_clim)
results$rotation


#SEED#

#SRW: updating as above
#pc.se_clim <- cbind(pc_clim, SeedD)

pc.se_clim <- left_join(pc_clim[,1:10], SeedD, by="orchard.num")


#TotalPhen
pc.se.tp <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.se_clim, family=beta_family (link="logit"))
summary(pc.se.tp)
#PC2          0.128315   0.055966    2.29   0.0219 *  


#SRW: This makes a lot more sense now and is in agreement with what you have 
#above...an effect of orchard type only for seeds, conventional have higher
#phenolic content. My guess is that result you had for "total fruit" was largely
#driven by the seeds. 

pc.se.tp.pc2 <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type* PC2 + (1|site.code)
                        , data=pc.se_clim, family=beta_family (link="logit"))
summary(pc.se.tp.pc2)
#orchard.typeOrganic     -0.26424    0.10203  -2.590   0.0096 ** 

plot(TotalPhen ~ PC2, data=pc.se_clim)

results$rotation


#SRW: Oh crap after the third round of very similar results I just realized that 
#actually all your data subsets that were supposed
#to be separated by tissue actually have the same data in them for TotalPhen
#and PhenRich, probably the average across all three tissues



#PhenRich 
pc.se.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 +
                      (1|site.code), data=pc.se_clim)
summary(pc.se.pr)
#PC1           0.3597     0.1899    1.89  0.05823 .  
#PC2           0.9112     0.3412    2.67  0.00758 **

pc.se.pr.pc1 <- glmmTMB(PhenRich ~ orchard.type*PC1 +
                          (1|site.code), data=pc.se_clim)
summary(pc.se.pr.pc1)
#nothing 

pc.se.pr.pc2 <- glmmTMB(PhenRich ~ orchard.type*PC2 +
                          (1|site.code), data=pc.se_clim)
summary(pc.se.pr.pc2)
#nothing 


plot(PhenRich ~ PC1, data=pc.se_clim)
plot(PhenRich ~ PC2, data=pc.se_clim)

results$rotation






#Q3: Which specific management practices are the most important drivers of fruit quality?----
#For this analysis, we'll be running GLMMs against every mgmt variable. When there's-
# a significant interaction we'll examine it on its analyzed along side orchard.type
#Q3-A: Physical Quality---------------------------------------------------------

#SRW: Again for these models I would suggest including orchard type in the large
#overall model but not following up with the interactive models. If you have a strong
#hypothesis for why there should be an interaction between orchard type and one
#of these other management variables, you could include that interaction in the 
#large model to start with, e.g.

#SSC
mgmt1 <- glmmTMB(SSC ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + 
                   orchard.type*Herbicides + (1|site.code), 
                 data=c, na.action = "na.fail")
summary(mgmt1)
Anova(mgmt1)
#Herbicides  47.9815  1  4.303e-12 ***


#SRW: with such a complex model it is usually recommended to do some sort of 
#model simplification/selection process, since model parameter estimates
#can be really off if you have a bunch of non-significant predictor variables in 
#the model. One approach is just to sequentially drop the least significant 
#variables until you have a "minimum adequate model" in which all predictor 
#variables are significant. This is the approach taken by Crawley in The R Book
#Side note: you should read The R Book Chp. 19 on Mixed Models, it is great

#Alternatively, you can do what they call a "model averaging" where you 
#fit all possible subsets of the full model, then average the parameter estimates
#across multiple models to account for uncertainty. This is often done on just
#a subset of the best-fitting models among all possible models 

#The classic text on this is:
#Burnham,K.P. and Anderson,D.R.2002 Model selection and multimodel inference:a practical information-theoretic approach.2nded.NewYork,Springer-Verlag.

#Also see:
#Dormann, Carsten F., et al. "Model averaging in ecology: A review of Bayesian, information‐theoretic, and tactical approaches for predictive inference." Ecological Monographs 88.4 (2018): 485-504.

#you can implement model averaging in the MuMIn package:

library(MuMIn)
dd <- dredge(mgmt1)  #note you have to specify na.action=na.fail in model above 
ma <- model.avg(dd, subset = delta < 4)
summary(ma)

#This gives you parameter estimates and p-values for all your predictors that are 
#based on a weighted average across all the "component models", which you define
#with the "subset= delta <4" argument, in this case keeping all the models that
#have dAIC <4


ssc_herb <- glmmTMB(SSC ~  orchard.type*Herbicides+ (1|site.code), data=c)
summary(ssc_herb)
Anova(ssc_herb)
#orchard.type            4.0006  1   0.045485 * 
#Herbicides              7.2867  1   0.006947 **


ggplot(c, aes(x=Herbicides, y=SSC, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot() +
  geom_jitter(width=0.1)

#SRW: This plot is interesting but the orchard.type*herbicides interaction does
#not come out as significant in the overall model or in the MuMIn approach.
#I would instead just suggest plotting all the significant main effects as single
#variable, you could still color points by orchard.type but don't separate
#into different box plots if there is no significant interaction, e.g.

ggplot(c, aes(x=Herbicides, y=SSC))+
  geom_boxplot() +
  geom_jitter(width=0.1, aes(x=Herbicides, y=SSC, color=orchard.type))

#Or plot the least squared means as I described above.
library(emmeans)
glmm_emmeans <-emmeans(ma, ~ orchard.type, type="response")
#SRW: Some errors here, I think  this is possible with a model averaged 
#object but I did not figure it out with a quick search

#you could also do it just on the top model
tm <- get.models(dd, 1)[[1]]
summary(tm)

glmm_emmeans <-emmeans(tm, ~ Herbicides, type="response")

glmm_emmeans
glmm_emmeans <- as.data.frame(glmm_emmeans)

ggplot(glmm_emmeans, aes(x=Herbicides, y = emmean))+
  geom_jitter(width=0.1, data=c, aes(x=Herbicides, y=SSC, color=orchard.type)) +
  geom_point () +
  geom_errorbar(aes(x=Herbicides, ymin = emmean-SE, ymax = emmean+SE), width = 0.1)
  
#Still there is not a large effect based on the plot, honestly I am surprised
#how significant the effect of herbicides is given the plots

#So I think we should be cautious with the interpretation but it is super 
#interesting that herbicides come out as such a strong factor
#I was just reading about organic herbicides and I guess there are quite a few
#most are vinegar based or contain a lot of terpenes like clove oil, limonene, etc


#SRW: I did not update the models below except to add orchard.type (just to see)
#but all the same as above would apply


#avgwgt
mgmt2 <- glmmTMB(avgwgt ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt2)
Anova(mgmt2)
#Acres       11.9567  1  0.0005445 ***
#Cultivation  3.9132  1  0.0479077 *  


wgt_acre <- glmmTMB(avgwgt ~orchard.type*Acres+ (1|site.code), data=c)
summary(wgt_acre)
Anova(wgt_acre)
#Acres              9.9260  1    0.00163 **

ggplot(c, aes(x=Acres, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

#Firmness 
mgmt3 <- glmmTMB(Firmness ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt3)
Anova(mgmt3)
#Herbicides  4.5342  1    0.03322 *


firm_herb <- glmmTMB(Firmness ~  orchard.type*Herbicides+ (1|site.code), data=c)
summary(firm_herb)
Anova(firm_herb)
#nothing 

ggplot(c, aes(x=Herbicides, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()

#Maturity Index 
mgmt4 <- glmmTMB(maturity.index ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt4)
Anova(mgmt4)
#Herbicides  4.7331  1    0.02959 *
#Mowing      3.9845  1    0.04592 *
#Weed_Mats   4.9615  1    0.02592 *
#Acres       4.4712  1    0.03447 *

#SRW: orchard type also has a strong effect on maturity

mat_herb <- glmmTMB(maturity.index ~  orchard.type*Herbicides+ (1|site.code), data=c)
summary(mat_herb)
Anova(mat_herb)
#orchard.type:Herbicides  6.6521  1   0.009904 ** 
#organic orchards that used herbicides had higher maturity values 


mat_herb <- glmmTMB(maturity.index ~  orchard.type*Acres+ (1|site.code), data=c)
summary(mat_herb)
Anova(mat_herb)
#Acres              5.5277  1    0.01872 *


ggplot(c, aes(x=Herbicides, y=maturity.index))+
  geom_smooth(method = "lm") +
  geom_boxplot()

ggplot(c, aes(x=Acres, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


#SRW: So both the acreage and the maturity results seem to be driven largely by 
#that one big 400 acre orchard, will have to interpret cautiously

#Also the fact that so many other management variable are related to maturity, 
#plus the fact that maturity is no doubt a key driver of other fruit quality
#variables like sugar content and firmness, I think that it would be worth
#considering the inclusion of fruit maturity as a key co-variate (i.e. predictor
#variable) in the other models, e.g.

#SSC
mgmt1 <- glmmTMB(SSC ~ orchard.type + maturity.index + Cultivation + Herbicides + 
                   Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + 
                   (1|site.code), 
                 data=c, na.action = "na.fail")
summary(mgmt1)
Anova(mgmt1)  

#Oh wow, now orchard type is way significant
emmeans(mgmt1, ~ orchard.type, type="response")
#SRW: organic has higher sugar content once you account for maturity and these
#other factors

#SRW: I also tried putting this into some of the phenolics models just to see 
#and it doesn't seem to change the interpretation there, at least regarding
#the effects of orchard.type which was my biggest question. So I don't think
#it is necessarily worth adding this into all the different models as another
#predictor. Just something to be aware of when it comes to interpretation so we
#don't over-interpret the results



#Q3-B: Fruit Chemistry ---------------------------------------------------------
#Total Phenolics

#SRW: Again I would skip the whole fruit analysis

#whole fruit 
mgmt.c <- glmmTMB(TotalPhen ~ Cultivation + Herbicides + Com_Mul + Mowing +
                    Weed_Mats + Cover_Crops + Acres + (1|site.code), data=c,)
summary(mgmt.c)
Anova(mgmt.c)
#Mowing      3.1708  1    0.07496 .
#Cover_Crops 2.9479  1    0.08599 .
#Acres       4.2656  1    0.03889 * 

#per tissue#

#Skin 
mgmt.sk <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                     Weed_Mats + Cover_Crops + Acres + (1|site.code), 
                   data=SkinD, family=beta_family(link="logit"))
summary(mgmt.sk)
Anova(mgmt.sk)
#Herbicides  3.3487  1    0.06726 .

#orchard type actually comes out significant here
emmeans(mgmt.sk, ~ orchard.type, type="response")

#SRW: I am getting some warnings with these models, not sure the source, probably
#it is just a lot of variables and some model simplification (as described
#above) will help to make sure the ones that are coming out significant really are
diagnose(mgmt.sk)  #this warning is very common


#Pulp 
mgmt.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                     Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                   data=PulpD, family=beta_family(link="logit"))
summary(mgmt.pu)
Anova(mgmt.pu)
#Cover_Crops 3.0087  1    0.08282 .
#Cultivation 4.4598  1    0.03470 *

#SRW: orchard type also comes out here but in the opposite direction?? 
emmeans(mgmt.pu, ~ orchard.type, type="response")


#Seed
mgmt.se <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                     Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                   data=SeedD, family=beta_family(link="logit"))
summary(mgmt.se)
Anova(mgmt.se)
#Nothing 

emmeans(mgmt.se, ~ orchard.type, type="response")



#Phenolics Richness 
#whole fruit 
mgmt.pr.c <- glmmTMB(PhenRich ~ Cultivation + Herbicides + Com_Mul + Mowing +
                       Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                     data=c)
summary(mgmt.pr.c)
Anova(mgmt.pr.c)
#Cultivation   7.6732  1   0.005605 ** 
#Mowing        4.2425  1   0.039425 *  
#Weed_Mats     3.9059  1   0.048116 *   

#per tissue #
#Skin 
mgmt.pr.sk <- glmmTMB(PhenRich ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                        Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                      data=SkinD)
summary(mgmt.pr.sk)
Anova(mgmt.pr.sk)
#Herbicides  3.8902  1   0.048569 * 

#richness higher in organic
emmeans(mgmt.pr.sk, ~ orchard.type, type="response")



mgmt.pr.sk.h <- glmmTMB(PhenRich ~ orchard.type*Herbicides+ (1|site.code), 
                        data=SkinD)
summary(mgmt.pr.sk.h)
Anova(mgmt.pr.sk.h)
#Herbicides              4.7020  1    0.03013 *


#Pulp 
mgmt.pr.pu <- glmmTMB(PhenRich ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                        Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                      data=PulpD)
summary(mgmt.pr.pu)
Anova(mgmt.pr.pu)
#Cultivation  11.8454  1  0.0005781 ***
#Com_Mul       5.7323  1  0.0166559 * 

mgmt.pr.pu.c <- glmmTMB(PhenRich ~ orchard.type*Cultivation+ (1|site.code), 
                        data=PulpD)
summary(mgmt.pr.pu.c)
Anova(mgmt.pr.pu.c)
#Cultivation              3.2555  1    0.07119 .

mgmt.pr.pu.cm <- glmmTMB(PhenRich ~ orchard.type*Com_Mul+ (1|site.code), 
                         data=PulpD)
summary(mgmt.pr.pu.cm)
Anova(mgmt.pr.pu.cm)
#orchard.type:Com_Mul 2.8062  1     0.0939 .

#Visualize this 
ggplot(c, aes(x=Com_Mul, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()

#Seed
mgmt.pr.se <- glmmTMB(PhenRich ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                        Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                      data=SeedD)
summary(mgmt.pr.se)
Anova(mgmt.pr.se)
#Cultivation  8.4221  1   0.003707 ** 
#Herbicides   4.7221  1   0.029778 *  
#Mowing       3.9252  1   0.047567 *  
#Weed_Mats    6.7875  1   0.009180 ** 
#Cover_Crops  6.7150  1   0.009560 ** 
#Acres        9.0503  1   0.002626 ** 

#richness lower in organic
emmeans(mgmt.pr.se, ~ orchard.type, type="response")


#Cultivation 
mgmt.pr.se.cu <- glmmTMB(PhenRich ~ orchard.type*Cultivation+ (1|site.code), 
                         data=SeedD)
summary(mgmt.pr.se.cu)
Anova(mgmt.pr.se.cu)
#orchard.type:Cultivation 2.9180  1    0.08760 .


#Herbicides 
mgmt.pr.se.he <- glmmTMB(PhenRich ~ orchard.type*Herbicides+ (1|site.code), 
                         data=SeedD)
summary(mgmt.pr.se.he)
Anova(mgmt.pr.se.he)
#nothing 

#Mowing 
mgmt.pr.se.m <- glmmTMB(PhenRich ~ orchard.type*Mowing+ (1|site.code), 
                        data=SeedD)
summary(mgmt.pr.se.m)
Anova(mgmt.pr.se.m)
#doesnt work 

#Weed_Mats 
mgmt.pr.se.w <- glmmTMB(PhenRich ~ orchard.type*Weed_Mats+ (1|site.code), 
                        data=SeedD)
summary(mgmt.pr.se.w)
Anova(mgmt.pr.se.w)
#doesnt work 

#Cover_Crops 
mgmt.pr.se.cc <- glmmTMB(PhenRich ~ orchard.type*Cover_Crops+ (1|site.code), 
                         data=SeedD)
summary(mgmt.pr.se.cc)
Anova(mgmt.pr.se.cc)
#nothing 

#Acres 
mgmt.pr.se.a <- glmmTMB(PhenRich ~ orchard.type*Acres+ (1|site.code), 
                        data=SeedD)
summary(mgmt.pr.se.a)
Anova(mgmt.pr.se.a)
#nothing 


#Q4: Which pest or diseases presence has the most significant affect on fruit quality-----

#SRW: I am not sure why but a lot of the variable names in this set of code 
#were not the same as in the datasets I have at this point, probably something to
#do with what I changed earlier. Time is now short for me to get this to you today
#as promised, so I am not going to try to change all that to make these run for me
#BUT, generally I really like having this as a separate question. 
#I would just suggest once you update the sections above, to update this to 
#a parallel approach (including orchard.type, model simplification, etc.)


#visualize pest data in indices 
ggpairs(c, columns=25:42) 

d <- pivot_longer(data=c, cols=27:44, names_to="pest", values_to="index")

ggplot(d, aes(x=pest, y=index, color=orchard.type))+
  geom_boxplot()+
  geom_jitter(width=0.2, height=0.1)

#SRW: added orchard type to this plot to visualize, it's interesting the pest
#pressure does not seem too different between org/conv

#removing everything that doesn't pass 3 on index 
#what were left with: 
#pests: aphids, apple maggots, coddling moth, 
#disease: apple scab, bitter rot, fireblight, powdery mildew, root rot 

#correlation matrix
p_pest <- dplyr::select(c,c("Powdery.mildew",     
                            "Apple.scab","Root.Rot","Fire.Blight","Apple.Maggots","Codling.Moth","Aphids","Bitter.Rot"))
#SRW: updated names in select to make this run

rcorr(as.matrix(p_pest))
corrplot(cor(p_pest))

testRes = cor.mtest(p_pest, conf.level = 0.95)

#add all p-values
corrplot(cor(p_pest), p.mat = testRes$p, insig = 'p-value', sig.level = -1)

#add significant level stars
corrplot(cor(p_pest), p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')

#SRW: some errors in these plot codes I didn't have time to troubleshoot



#strong relationships between: 
#apple scab and root rot 
#apple scab and coddling moth
#aphids and root rot 
#aphids and maggots
#maggots and coddling moth 
#Q4-A: Physical Quality---------------------------------------------------------
###SSC###
ssc_pest <- glmmTMB(SSC ~ Aphids+applemaggots+codmoth+
                      powmil+bitterrot+applescab+rootrot+fireblight+ (1|site.code/orchard.num), 
                    data=c)
summary(ssc_pest)
Anova(ssc_pest)
#fireblight   15.1670  1  9.841e-05 ***


ssc_fire <- glmmTMB(SSC ~ orchard.type*fireblight+ (1|site.code), 
                    data=c)
summary(ssc_fire)
Anova(ssc_fire)
#orchard.type            3.7194  1   0.053782 . 
#fireblight              6.7217  1   0.009525 **


ggplot(c, aes(x=fireblight, y=SSC, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#lower SSC with lower reported fireblight


ssc_pressure <- glmmTMB(SSC ~ orchard.type*Pest_Index+ (1|site.code), 
                        data=c)
summary(ssc_pressure)
Anova(ssc_pressure)
#orchard.type:Pest_Index 3.0044  1    0.08304 .

ggplot(c, aes(x=Pest_Index, y=SSC, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher SSC in orangic orchards with high pest indexes 

###avgwgt###
wgt_pest <- glmmTMB(avgwgt ~ Aphids+applemaggots+codmoth+
                      powmil+bitterrot+applescab+rootrot+fireblight+ (1|site.code), 
                    data=c)
summary(wgt_pest)
Anova(wgt_pest)
#rootrot      8.0041  1   0.004667 **
#fireblight   4.6472  1   0.031103 *

wgt_rootrot <- glmmTMB(avgwgt ~ orchard.type*rootrot+ (1|site.code), 
                       data=c)
summary(wgt_rootrot)
Anova(wgt_rootrot)
#rootrot              5.1811  1    0.02283 *


wgt_fire <- glmmTMB(avgwgt ~ orchard.type*fireblight+ (1|site.code), 
                    data=c)
summary(wgt_fire)
Anova(wgt_fire)
#orchard.type:fireblight 75.8075  1     <2e-16 ***


ggplot(c, aes(x=fireblight, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avgwgt decrease as fireblight increases 

ggplot(c, aes(x=rootrot, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increased rootrot increases avg wgt? 

wgt_pressure <- glmmTMB(avgwgt ~ orchard.type*Pest_Index+ (1|site.code), 
                        data=c)
summary(wgt_pressure)
Anova(wgt_pressure)
#nothing 

###firmness###
frm_pest <- glmmTMB(Firmness ~ Aphids+applemaggots+codmoth+
                      powmil+bitterrot+applescab+rootrot+fireblight+ (1|site.code), 
                    data=c)
summary(frm_pest)
Anova(frm_pest)
#Aphids       5.6389  1    0.01757 *
#applemaggots 3.1469  1    0.07607 .
#codmoth      3.3438  1    0.06746 .


frm_aphid <- glmmTMB(Firmness ~ Aphids*Pest_Index+ (1|site.code), 
                     data=c)
summary(frm_aphid)
Anova(frm_aphid)
#nothing 

frm_amaggots <- glmmTMB(Firmness ~ applemaggots*Pest_Index+ (1|site.code), 
                        data=c)
summary(frm_amaggots)
Anova(frm_amaggots)
#nothing 

frm_codmoth <- glmmTMB(Firmness ~ codmoth*Pest_Index+ (1|site.code), 
                       data=c)
summary(frm_codmoth)
Anova(frm_codmoth)
#nothing 

#pest index 
frm_pressure <- glmmTMB(Firmness ~ orchard.type*Pest_Index+ (1|site.code), 
                        data=c)
summary(frm_pressure)
Anova(frm_pressure)
#nothing 


ggplot(c, aes(x=Aphids, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

ggplot(c, aes(x=applemaggots, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

ggplot(c, aes(x=codmoth, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

###maturity###
mat_pest <- glmmTMB(maturity.index ~ Aphids+applemaggots+codmoth+
                      powmil+bitterrot+applescab+rootrot+fireblight+ (1|site.code), 
                    data=c)
summary(mat_pest)
Anova(mat_pest)
#nothing 

mat_pressure <- glmmTMB(maturity.index ~ orchard.type*Pest_Index+ (1|site.code), 
                        data=c)
summary(mat_pressure)
Anova(mat_pressure)
#Pest_Index              10.9447  1  0.0009387 ***

ggplot(c, aes(x=Pest_Index, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

###pest pressure by mgmt system 
pest_pressure <- glmmTMB(Pest_Index ~ orchard.type + (1|site.code), 
                         data=c)
summary(pest_pressure)
Anova(pest_pressure)
#nothing 
#Q4-B: Fruit Chemistry ---------------------------------------------------------
##TotalPhen##
pest.tp.c <- glmmTMB(TotalPhen ~ Aphids+Apple.Maggots+Codling.Moth+
                       Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                     data=c)
summary(pest.tp.c)
Anova(pest.tp.c)
#Powdery.mildew 4.9524  1    0.02605 *


#per tissue# 
#Skin 
pest.tp.sk <- glmmTMB(TotalPhen ~ Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                      data=SkinD)
summary(pest.tp.sk)
Anova(pest.tp.sk)
#Bitter.Rot     3.0038  1    0.08307 .
#Fire.Blight    4.2174  1    0.04001 *

#Pulp
pest.tp.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                      data=PulpD, family=beta_family(link="logit"))
summary(pest.tp.pu)
Anova(pest.tp.pu)
#Aphids         4.6917  1    0.03031 *
#Apple.Maggots  5.6223  1    0.01773 *
#Powdery.mildew 2.7381  1    0.09798 .
#Root.Rot       2.8451  1    0.09165 .


#Seed 
pest.tp.se <- glmmTMB(TotalPhen ~ Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                      data=SkinD)
summary(pest.tp.se)
Anova(pest.tp.se)
#Bitter.Rot     3.0038  1    0.08307 .
#Fire.Blight    4.2174  1    0.04001 *


##PhenRich##
#Whole Fruit 
pest.pr.c <- glmmTMB(PhenRich ~ Aphids+Apple.Maggots+Codling.Moth+
                       Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                     data=c)
summary(pest.pr.c)
Anova(pest.pr.c)
#Apple.Maggots  4.8184  1    0.02816 *
#Codling.Moth   4.3041  1    0.03802 *

#per tissue# 
#skin 
pest.pr.sk <- glmmTMB(PhenRich ~ Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                      data=SkinD)
summary(pest.pr.sk)
Anova(pest.pr.sk)
#Powdery.mildew 2.9521  1    0.08577 .
#Apple.Maggots  3.0997  1    0.07831 .


#pulp
pest.pr.pu <- glmmTMB(PhenRich ~ Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+(1|site.code), 
                      data=PulpD)
summary(pest.pr.pu)
Anova(pest.pr.pu)
#Aphids          9.7261  1  0.0018167 ** 
#Apple.Maggots  13.5379  1  0.0002338 ***
#Powdery.mildew 12.5320  1  0.0004000 ***

#Seed
pest.pr.se <- glmmTMB(PhenRich ~ Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+(1|site.code), 
                      data=SeedD)
summary(pest.pr.se)
Anova(pest.pr.se)
#Aphids          2.8796  1  0.0897092 .  
#Bitter.Rot      2.7320  1  0.0983590 .  
#Root.Rot       12.4586  1  0.0004161 **
#Q5: How does fruit quality compare to total phenolics and phenolic richness--------
#PhenRich 
tp.qual <- glmmTMB(TotalPhen~ SSC+Firmness+avgwgt+maturity.index + (1|site.code)
                   , data=c)
summary(tp.qual)
Anova(tp.qual)
#Firmness       2.7617  1    0.09654 .

#Firmness 
tp.firm <- glmmTMB(TotalPhen~ orchard.type*Firmness + 
                     (1|site.code), data=c)
summary(tp.firm)
Anova(tp.firm)
#Firmness              5.2726  1    0.02166 *


ggplot(c, aes(x=Firmness, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher total phen in convential as firmness increases 

#PhenRich 
pr.qual <- glmmTMB(PhenRich~ SSC+Firmness+avgwgt+maturity.index + (1|site.code)
                   ,data=c)
summary(pr.qual)
Anova(pr.qual)
#SSC            5.8072  1    0.01596 *
#avgwgt         3.7377  1    0.05320 .


#SSC 
pr.ssc <- glmmTMB(PhenRich~ orchard.type*SSC + 
                    (1|site.code), data=c)
summary(pr.ssc)
Anova(pr.ssc)
#SSC              4.6719  1    0.03066 *

ggplot(c, aes(x=SSC, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#phen rich decrease as SSC decreases 

#avgwgt 
pr.wgt <- glmmTMB(PhenRich~ orchard.type*avgwgt + 
                    (1|site.code), data=c)
summary(pr.wgt)
Anova(pr.wgt)
#nothing 

ggplot(c, aes(x=avgwgt, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#phen rich decrease as avgwgt increases? 



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

#Figures------------------------------------------------------------------------