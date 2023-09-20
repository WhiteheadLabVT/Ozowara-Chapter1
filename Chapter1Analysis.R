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

c <- d
c <- c %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)

c <- left_join(c, Tree.sum, by="orchard.num")
c <- left_join(c, Orchard, by="orchard.num")

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
TreeLat <- left_join(Tree, Orchard[,c(2,4)], by="orchard.num")
#SSC
Lat1<- glmmTMB(SSC ~ orchard.type*Latitude + (1|site.code/orchard.num), data=TreeLat)
summary(Lat1)
Anova(Lat1) 
#orchard.type           0.5016  1     0.4788    
#Latitude              18.2814  1  1.906e-05 ***
#orchard.type:Latitude  0.2196  1     0.6393  

hist(resid(Lat1))  #looks great
diagnose(Lat1)

ggplot(c, aes(x=Latitude, y=SSC, color=orchard.type))+
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

hist(resid(Lat3))  #looks great
diagnose(Lat3)
##strongest interaction 

#Maturity Index 
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

#Q1-D: Which compounds distinguish fruits based on latitude?----------
#Analysis: NMDS and random forest
#converting latitude into categorical variable
Orchard <- Orchard %>% 
  mutate(lat_cat=cut(Latitude, breaks=c(-Inf, 42, Inf), labels=c("low", "high")))

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
                               SKIN=8, PULP=1, SEED=2)
d.expl$Symbol <- as.numeric(as.character(d.expl$Symbol))

d.expl$lat_cat = as.factor(d.expl$lat_cat)


plot(m.NMDS_2, type="n") #plots the ordination axes only
points(m.NMDS_2, pch=d.expl$Symbol,
       col=as.character(d.expl$Color), cex = 0)



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

#Q3: Which specific management practices are the most important drivers of fruit quality?----
#For this analysis, we'll be running GLMMs against every mgmt variable. When there's-
# a significant interaction we'll examine it on its analyzed along side orchard.type
#Q3-A: Physical Quality---------------------------------------------------------
#SSC
mgmt1 <- glmmTMB(SSC ~ orchard.type+ Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt1)
Anova(mgmt1)
#Herbicides   49.0400  1  2.508e-12 ***

ggplot(c, aes(x=Herbicides, y=SSC, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()

#avgwgt
mgmt2 <- glmmTMB(avgwgt ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt2)
Anova(mgmt2)
#Acres        8.6584  1   0.003256 **
#Cultivation  3.9567  1   0.046687 * 

ggplot(c, aes(x=Acres, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

ggplot(c, aes(x=Cultivation, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()

#Firmness 
mgmt3 <- glmmTMB(Firmness ~ orchard.type + Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt3)
Anova(mgmt3)
#Herbicides  4.5342  1    0.03322 *

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

ggplot(c, aes(x=Herbicides, y=maturity.index))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#organic orchards that used herbicides had higher maturity values 


ggplot(c, aes(x=Acres, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


#Q3-B: Fruit Chemistry ---------------------------------------------------------

#Q4: Which pest or diseases presence has the most significant affect on fruit quality-----
#visualize pest data in indices 
ggpairs(c, columns=25:42) 

d <- pivot_longer(data=c, cols=25:42, names_to="pest", values_to="index")

ggplot(d, aes(x=pest, y=index))+
  geom_boxplot()+
  geom_jitter(width=0.2, height=0.1)

#removing everything that doesn't pass 3 on index 
#what were left with: 
#pests: aphids, apple maggots, coddling moth, 
#disease: apple scab, bitter rot, fireblight, powdery mildew, root rot 

#correlation matrix
p_pest <- dplyr::select(c,c(Aphids, Apple.Maggots, Codling.Moth, 
Powdery.mildew, Bitter.Rot, Apple.scab, Root.Rot, Fire.Blight))

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
#apple scab and root rot 
#apple scab and coddling moth
#aphids and root rot 
#aphids and maggots
#maggots and coddling moth 
#Q4-A: Physical Quality---------------------------------------------------------
###SSC###
ssc_pest <- glmmTMB(SSC ~ orchard.type + Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ 
                      (1|site.code), 
                    data=c)
summary(ssc_pest)
Anova(ssc_pest)
#Fire.Blight    17.4889  1   2.89e-05 ***

ggplot(c, aes(x=Fire.Blight, y=SSC, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#lower SSC with lower reported fireblight


#pest pressure index 
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
wgt_pest <- glmmTMB(avgwgt ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                    data=c)
summary(wgt_pest)
Anova(wgt_pest)
#orchard.type    2.9807  1   0.084263 . 
#Codling.Moth    2.7342  1   0.098221 . 
#Root.Rot       10.0361  1   0.001535 **
#Fire.Blight     6.9323  1   0.008465 **

ggplot(c, aes(x=Root.Rot, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increased rootrot increases avg wgt? 

ggplot(c, aes(x=Fire.Blight, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avgwgt decrease as fireblight increases 

ggplot(c, aes(x=Codling.Moth, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avgwgt higher in orchards without codling moth, organic weight higher with high pressure 

#pest pressure index 
wgt_pressure <- glmmTMB(avgwgt ~ orchard.type*Pest_Index+ (1|site.code), 
                        data=c)
summary(wgt_pressure)
Anova(wgt_pressure)
#nothing 

###firmness###
frm_pest <- glmmTMB(Firmness ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+  (1|site.code), 
                    data=c)
summary(frm_pest)
Anova(frm_pest)
#orchard.type   12.8200  1  0.0003429 ***
#Aphids          5.1133  1  0.0237430 *  
#Apple.Maggots   3.6107  1  0.0574072 .  


ggplot(c, aes(x=Aphids, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

ggplot(c, aes(x=Apple.Maggots, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

#pest index 
frm_pressure <- glmmTMB(Firmness ~ orchard.type*Pest_Index+ (1|site.code), 
                        data=c)
summary(frm_pressure)
Anova(frm_pressure)
#nothing 

###maturity###
mat_pest <- glmmTMB(maturity.index ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                    data=c)
summary(mat_pest)
Anova(mat_pest)
#nothing 

#pest pressure 
mat_pressure <- glmmTMB(maturity.index ~ orchard.type*Pest_Index+ (1|site.code), 
                        data=c)
summary(mat_pressure)
Anova(mat_pressure)
#Pest_Index              10.9447  1  0.0009387 ***


ggplot(c, aes(x=Pest_Index, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


#Q4-B: Fruit Chemistry ---------------------------------------------------------
#Total Phenolics
#Skin 
pest.tp.sk <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                      data=SkinD, family=beta_family(link="logit"))
summary(pest.tp.sk)
Anova(pest.tp.sk)
#Powdery.mildew 2.8202  1    0.09309 .
#Bitter.Rot     4.1312  1    0.04210 *

ggplot(SkinD, aes(x=Bitter.Rot, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher total phen in orchards with higher bitter rot 

ggplot(SkinD, aes(x=Powdery.mildew, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher total phen in orchards with higher pow mil  


#Pulp
pest.tp.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                      data=PulpD, family=beta_family(link="logit"))
summary(pest.tp.pu)
Anova(pest.tp.pu)
#orchard.type   5.8612  1   0.015479 * 
#Aphids         9.3802  1   0.002193 **
#Apple.Maggots  7.6171  1   0.005782 **
#Root.Rot       3.5093  1   0.061025 . 
#Fire.Blight    4.1346  1   0.042014 * 

ggplot(PulpD, aes(x=Aphids, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher total phen in conventional 

ggplot(PulpD, aes(x=Apple.Maggots, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#total phen decrease as applemag increase, steadily higher in conventional 

ggplot(PulpD, aes(x=Root.Rot, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decreases, higher in conventional 

ggplot(PulpD, aes(x=Fire.Blight, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increases, higher in conventional 


#Seed 
pest.tp.se <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                      data=SeedD, family=beta_family(link="logit"))
summary(pest.tp.se)
Anova(pest.tp.se)
#orchard.type   6.9068  1   0.008587 **
#Root.Rot       9.2816  1   0.002315 **

ggplot(SeedD, aes(x=Root.Rot, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decreases, higher in conventional 


##PhenRich##
#skin 
pest.pr.sk <- glmmTMB(PhenRich ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+ (1|site.code), 
                      data=SkinD)
summary(pest.pr.sk)
Anova(pest.pr.sk)
#Apple.Maggots  3.9307  1    0.04741 *
#Codling.Moth   5.2264  1    0.02225 *
#Powdery.mildew 3.2477  1    0.07152 .
#Bitter.Rot     3.3729  1    0.06628 .

ggplot(SkinD, aes(x=Apple.Maggots, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increase o, decrease c 

ggplot(SkinD, aes(x=Codling.Moth, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increase o, decrease c 

#pulp
pest.pr.pu <- glmmTMB(PhenRich ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+(1|site.code), 
                      data=PulpD)
summary(pest.pr.pu)
Anova(pest.pr.pu)
#Aphids          9.8163  1  0.0017297 ** 
#Apple.Maggots  13.6250  1  0.0002232 ***
#Powdery.mildew 11.0199  1  0.0009014 ***

ggplot(PulpD, aes(x=Apple.Maggots, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decrease, higher in c, sharoer in o 

ggplot(PulpD, aes(x=Aphids, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#o starts low ends higher than c 

ggplot(PulpD, aes(x=Powdery.mildew, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#sharp increase in o, starts low ends high 


#Seed
pest.pr.se <- glmmTMB(PhenRich ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Bitter.Rot+Apple.scab+Root.Rot+Fire.Blight+(1|site.code), 
                      data=SeedD)
summary(pest.pr.se)
Anova(pest.pr.se)
#Aphids          3.2302  1  0.0722926 .  
#Bitter.Rot      2.9140  1  0.0878156 .  
#Root.Rot       12.9458  1  0.0003206 ***


ggplot(SeedD, aes(x=Aphids, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#c starts high ends low, o reverse 

ggplot(SeedD, aes(x=Bitter.Rot, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decrease c,increase o

ggplot(SeedD, aes(x=Root.Rot, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decreases, higher in c 

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

ggplot(PulpD, aes(x=avgwgt, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

ggplot(PulpD, aes(x=maturity.index, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


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

ggplot(SeedD, aes(x=avgwgt, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


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