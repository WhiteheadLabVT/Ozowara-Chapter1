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

#Data Organization and Restructuring--------------------------------------------
#Fruit level data 
d <- read_csv("Fruit Level Data.csv")

d <- d %>%
  mutate(TotalPhen=rowSums(across(8:29)),
         PhenRich=rowSums(across(8:29)!=0))

d <- d %>% 
  mutate_at(c("orchard.num", "orchard.type", "site.code", "Tree", "SampleID"), as.factor)
View(d)

#separate by tissue type 
d.sk <- filter(d, Tissue=="SKIN")
d.sk <- dplyr::select(d.sk, -(6+which(colSums(d.sk[7:31], na.rm=TRUE) %in% 0)))
d.pu <- filter(d, Tissue=="PULP")
d.pu <- dplyr::select(d.pu, -(6+which(colSums(d.pu[7:31], na.rm=TRUE) %in% 0)))
d.se <- filter(d, Tissue=="SEED")
d.se <- dplyr::select(d.se, -(6+which(colSums(d.se[7:31], na.rm=TRUE) %in% 0)))

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
d.comp.pu <- d.pu[,8:28]
d.expl.pu <- d.pu[,1:6]
#Seed
d.comp.se <- d.se[,8:29]
d.expl.se <- d.se[,1:6]

#Condensing to orchard level data 
Orchard <- read_csv("Orchard Level Data.csv")
Orchard$orchard.num <- as.factor(as.character(Orchard$orchard.num))
applemaggots= c$`Apple Maggots`
codmoth = c$`Codling Moth`
powmil= c$`Powdery mildew`
bitterrot= c$`Bitter Rot`
applescab=c$`Apple scab`
rootrot= c$`Root Rot`
fireblight= c$`Fire Blight`

#convert weights based on averages and bag weights
Tree <- read_csv("Tree Level Data.csv")
Tree$orchard.num <- as.factor(as.character(Tree$orchard.num))

Tree<- Tree %>%
  mutate(avgwgt=(apple.wgt.3-bag.weight)/3)%>%
dplyr::select(-c(3:4))

Tree <- Tree %>%
  group_by(orchard.num) %>%
  summarise_at(c("Firmness", "SSC", "maturity.index", "avgwgt")
               , mean, na.rm = TRUE)
View(Tree)

#whole fruit 
Fruit <- d 
Fruit$orchard.num <- as.factor(as.character(Fruit$orchard.num))

Fruit <- Fruit %>%
  group_by(orchard.num, Tissue) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)

#Join all data 
c <- left_join(Tree, Orchard, by="orchard.num")
c <- left_join(c, Fruit, by="orchard.num")
View(c)

#Split by tissue 
SkinD <- filter(c, Tissue=="SKIN")
PulpD <- filter(c, Tissue=="PULP")
SeedD <- filter(c, Tissue=="SEED")


#How do management systems interact with broad climatic changes across latitude?-----
#total phenolics
tp1 <- glmmTMB(TotalPhen ~ orchard.type*Latitude + (1|site.code/orchard.num), data=c)
summary(tp1)
Anova(tp1)

#total p per tissue type 
tp.p <- glmmTMB(TotalPhen ~ orchard.type*Latitude + (1|site.code/orchard.num), data=SkinD)
summary(tp.p)
Anova(tp.p)

tp.sk <- glmmTMB(TotalPhen ~ orchard.type*Latitude + (1|site.code/orchard.num), data=PulpD)
summary(tp.sk)
Anova(tp.sk)

tp.se <- glmmTMB(TotalPhen ~ orchard.type*Latitude + (1|site.code/orchard.num), data=SeedD)
summary(tp.se)
Anova(tp.se) 

#phenolics richness
pr1 <- glmmTMB(PhenRich~ orchard.type*Latitude + (1|site.code/orchard.num), data=c)
summary(pr1)
Anova(pr1)

#phen rich per tissue type 
pr.p <- glmmTMB(PhenRich ~ orchard.type*Latitude+(1|site.code/orchard.num), data=SkinD)
summary(pr.p)
Anova(pr.p)
#orchard.type:Latitude 3.6179  1    0.05716 .


pr.sk <- glmmTMB(PhenRich ~ orchard.type*Latitude+(1|site.code/orchard.num), data=PulpD)
summary(pr.sk)
Anova(pr.sk)
#not sig 

pr.se <- glmmTMB(PhenRich ~ orchard.type*Latitude+(1|site.code/orchard.num), data=SeedD)
summary(pr.se)
Anova(pr.se)
#orchard.type          6.1841  1    0.01289 *
#Latitude              4.1285  1    0.04217 *

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

#Which abiotic factors are the most important drivers of fruit chem?------------

#Analysis: principle components analysis followed by PC regression with each variable
#1 for whole fruit, skin, pulp, seed

#can technically use the same first PCA for all 4 different analyses because
#the values of the selected variables are the same

###Principle Components Analysis###
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
dev.off()
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
#PC's 6, 7, 8, 9 are very low/ not worth looking at


###Principle Components Regression###

#Whole Fruit#
pc_clim <- as.data.frame(results$x)
pc_clim <- cbind(pc_clim, CCD)

#TotalPhen
p1 <- glmmTMB(TotalPhentrans ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + 
                (1|site.code), data=pc_clim)
summary(p1)
#PC2 p= 2.25e-05 ***     

plot(TotalPhen ~ PC2, data=pc_clim)
results$rotation
#Total Phen in the whole fruit increases as temps, UVI and precip increase. 
#Total Phen in the whole fruit increases as lat, elevation, and proximity to 
#water decrease
#not specific to orchard type


#PhenRich
p2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + 
                (1|site.code), data=pc_clim)
summary(p2)
#PC2 p= 0.00019 ***    

plot(PhenRich ~ PC2, data=pc_clim)
results$rotation
#Phen Rich in the whole fruit increases as temps, UVI and precip increase. 
#Phen Rich in the whole fruit increases as lat, elevation, and proximity to 
#water decrease
#not specific to ochrad type

#SKIN#
pc.sk_clim <- as.data.frame(results$x)
pc.sk_clim <- cbind(pc.sk_clim, SkinD)

#TotalPhen
p.sk1 <- glmmTMB(TotalPhentrans ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5
                 +(1|site.code), data=pc.sk_clim)
summary(p.sk1)#nothing

#PhenRich
psk2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + 
                  (1|site.code), data=pc.sk_clim)
summary(psk2)#nothing


#PULP#
pc.pu_clim <- as.data.frame(results$x)
pc.pu_clim <- cbind(pc.pu_clim, PulpD)

##TotalPhen###
ppu1 <- glmmTMB(TotalPhentrans ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + 
                  (1|site.code), data=pc.pu_clim)
summary(ppu1)
#orchard.typeOrganic p= 0.0706 . 
#PC2 p= 3.86e-05 ***

plot(TotalPhentrans ~ PC2, data=pc.pu_clim)
results$rotation
#Total Phen in the pulp increases as temps, UVI and precip increase. 
#Total Phen in the pulp increases as lat, elevation, and proximity to 
#water decrease
#slight significance in organic 

ppu2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 + 
                  (1|site.code), data=pc.pu_clim)
summary(ppu2)
#PC2 p= 0.000558 ***
#PC3 p= 0.052757 .     

plot(PhenRich ~ PC2, data=pc.pu_clim)
plot(PhenRich ~ PC3, data=pc.pu_clim)
results$rotation

#Phen Rich in the pulp increases as temps, UVI and precip increase. 
#Phen Rich in the pulp increases as lat, elevation, and proximity to 
#water decrease
#not specific to orchard type


#SEED#
pc.se_clim <- as.data.frame(results$x)
pc.se_clim <- cbind(pc.se_clim, SeedD)

#TotalPhen
p.se1 <- glmmTMB(TotalPhentrans ~orchard.type+ PC1 + PC2 + PC3 + PC4 +PC5
                 , data=pc.se_clim)
summary(p.se1)
#orchard.typeOrganic p= 0.000994 ***
#PC2 p= 7.16e-08 ***

plot(TotalPhentrans ~ PC2, data=pc.se_clim)
results$rotation
#total phen increase specific to organic 

p.se2 <- glmmTMB(PhenRich ~orchard.type+ PC1 + PC2 + PC3 + PC4 + PC5 +
                   (1|site.code), data=pc.se_clim)
summary(p.se2)
#orchard.typeOrganic p=0.0133 *
#PC1 0.0050 **   
#PC2 4.49e-10 ***
#PC3 0.0197 *
#PC5 0.0596 . 

plot(PhenRich ~ PC1, data=pc.se_clim)
plot(PhenRich ~ PC2, data=pc.se_clim)
plot(PhenRich ~ PC3, data=pc.se_clim)
plot(PhenRich ~ PC5, data=pc.se_clim)
results$rotation

#phen rich increase in organic systems 


#Which specific management practices are the most important drivers of fruit chem?----
#Total Phenolics
#whole fruit 
mgmt1 <- glmmTMB(TotalPhen ~ Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt1)
Anova(mgmt1)
#nothing 

#per tissue 
#Skin 
mgmt1.sk <- glmmTMB(TotalPhen ~ Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=SkinD)
summary(mgmt1.sk)
Anova(mgmt1.sk)

#pulp 
mgmt1.pu <- glmmTMB(TotalPhen ~ Cultivation + Herbicides + Com_Mul + Mowing +
                      Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                    data=PulpD)
summary(mgmt1.pu)
Anova(mgmt1.pu)

#seed
mgmt1.se <- glmmTMB(TotalPhen ~ Cultivation + Herbicides + Com_Mul + Mowing +
                      Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                    data=SeedD)
summary(mgmt1.se)
Anova(mgmt1.se)



#Phenolics Richness 
#whole fruit 
mgmt2 <- glmmTMB(PhenRich ~ Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt2)
Anova(mgmt2)
#nothing 

#per tissue 
#Skin 
mgmt2.sk <- glmmTMB(PhenRich ~ Cultivation + Herbicides + Com_Mul + Mowing +
                      Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                    data=SkinD)
summary(mgmt2.sk)
Anova(mgmt2.sk)
#Herbicides  5.7083  1    0.01688 *

sk.herb <- glmmTMB(PhenRich ~ orchard.type*Herbicides+ (1|site.code), 
                    data=SkinD)
summary(sk.herb)
Anova(sk.herb)
#nothing 


#pulp 
mgmt2.pu <- glmmTMB(PhenRich ~ Cultivation + Herbicides + Com_Mul + Mowing +
                      Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                    data=PulpD)
summary(mgmt2.pu)
Anova(mgmt2.pu)
#Cultivation 10.2296  1   0.001382 **
#Mowing       3.0575  1   0.080363 . 

pu.cult <- glmmTMB(PhenRich ~ orchard.type*Cultivation+ (1|site.code), 
                   data=PulpD)
summary(pu.cult)
Anova(pu.cult)
#Cultivation              7.7980  1    0.00523 **

ggplot(PulpD, aes(x=Cultivation, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#lower richness when cultivation used
#need to explore further 

#NANS 
pu.mow <- glmmTMB(PhenRich ~ orchard.type*Mowing+ (1|site.code), 
                   data=PulpD)
summary(pu.mow)
Anova(pu.mow)


#seed
mgmt2.se <- glmmTMB(PhenRich ~ Cultivation + Herbicides + Com_Mul + Mowing +
                      Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                    data=SeedD)
summary(mgmt2.se)
Anova(mgmt2.se)
#Cultivation 7.8545  1   0.005069 **
#Herbicides  8.0635  1   0.004517 **
#Weed_Mats   9.4595  1   0.002101 **
#Fire_Mgmt   6.7139  1   0.009567 **

se.cult <- glmmTMB(PhenRich ~ orchard.type*Cultivation+ (1|site.code), 
                   data=SeedD)
summary(se.cult)
Anova(se.cult)
#orchard.type             5.7312  1    0.01667 *
#Cultivation              2.8729  1    0.09008 .


se.herb <- glmmTMB(PhenRich ~ orchard.type*Herbicides+ (1|site.code), 
                   data=SeedD)
summary(se.herb)
Anova(se.herb)
#orchard.type            6.2786  1    0.01222 *

#NANS 
se.weedm <- glmmTMB(PhenRich ~ orchard.type*Weed_Mats+ (1|site.code), 
                   data=SeedD)
summary(se.weedm)
Anova(se.weedm)


ggplot(SeedD, aes(x=Cultivation, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#lower richness when cultivation used 

ggplot(SeedD, aes(x=Herbicides, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#weird

#Which pest or diseases presence has the most significant affect on fruit chem?-----
#first part explained in the PhysicalxClimate Code

#GLMMs models

#Whole Fruit 
pr.pest1 <- glmmTMB(PhenRich ~ Aphids+applemaggots+codmoth+
powmil+bitterrot+applescab+rootrot+fireblight+ (1|site.code), 
                    data=c)
summary(pr.pest1)
Anova(pr.pest1)
#nothing 

#per tissue 
#skin 
pr.pest1.sk <- glmmTMB(PhenRich ~ Aphids+applemaggots+codmoth+
                      powmil+bitterrot+applescab+rootrot+fireblight+ (1|site.code), 
                    data=SkinD)
summary(pr.pest1.sk)
Anova(pr.pest1.sk)


#pulp
pr.pest1.pu <- glmmTMB(PhenRich ~ Aphids+applemaggots+codmoth+
                         powmil+bitterrot+applescab+rootrot+fireblight+ (1|site.code), 
                       data=PulpD)
summary(pr.pest1.pu)
Anova(pr.pest1.pu)

#How does fruit quality compare to total phenolics and phenolic richness--------

#SSC x Chemsitry 
ssc.tp <- glmmTMB(TotalPhentrans~ orchard.type*SSC + (1|site.code/orchard.num), data=CCD)
summary(ssc.tp)
Anova(ssc.tp)

ssc.pr <- glmmTMB(PhenRich~ orchard.type*SSC + (1|site.code/orchard.num), data=CCD)
summary(ssc.pr)
Anova(ssc.pr)
#orchard.type:SSC 3.5753  1    0.05864 .


#avgwgt x Chemsitry 
avg.tp <- glmmTMB(TotalPhentrans~ orchard.type*avgwgt + (1|site.code/orchard.num), data=CCD)
summary(avg.tp)
Anova(avg.tp)

avg.pr <- glmmTMB(PhenRich~ orchard.type*avgwgt + (1|site.code/orchard.num), data=CCD)
summary(avg.pr)
Anova(avg.pr)
#avgwgt              5.2503  1    0.02194 *


#Firmness x Chemsitry 
frm.tp <- glmmTMB(TotalPhentrans~ orchard.type*Firmness + (1|site.code/orchard.num), data=CCD)
summary(frm.tp)
Anova(frm.tp)

frm.pr <- glmmTMB(PhenRich~ orchard.type*Firmness + (1|site.code/orchard.num), data=CCD)
summary(frm.pr)
Anova(frm.pr)

#maturity.index x Chemsitry 
mi.tp <- glmmTMB(TotalPhentrans~ orchard.type*maturity.index + (1|site.code/orchard.num), data=CCD)
summary(mi.tp)
Anova(mi.tp)

mi.pr <- glmmTMB(PhenRich~ orchard.type*maturity.index + (1|site.code/orchard.num), data=CCD)
summary(mi.pr)
Anova(mi.pr)


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


ggplot(d, aes(x = Tissue, y = TotalPhentrans)) +
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





#How do management systems interact with broad climatic changes across latitude?

#phen rich 
ggplot(CCD, aes(x=Latitude, y=PhenRich, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("PhenRich") +
  xlab ("Latitude")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

ggplot(CCD, aes(x=elevation, y=PhenRich, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("PhenRich") +
  xlab ("Elevation (m)")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")


ggplot(CCD, aes(x=Longitude, y=PhenRich, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("PhenRich") +
  xlab ("Elevation (m)")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

ggplot(CCD, aes(x=Prox.Water, y=PhenRich, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("PhenRich") +
  xlab ("Elevation (m)")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

#total phenolics 
ggplot(CCD, aes(x=Latitude, y=TotalPhentrans, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("Total Phenolics") +
  xlab ("Latitude")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

ggplot(CCD, aes(x=elevation, y=TotalPhentrans, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("PhenRich") +
  xlab ("Elevation (m)")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

ggplot(CCD, aes(x=Longitude, y=TotalPhentrans, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("PhenRich") +
  xlab ("Elevation (m)")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

ggplot(CCD, aes(x=Prox.Water, y=TotalPhentrans, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("PhenRich") +
  xlab ("Elevation (m)")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")


ggplot(pc_clim, aes(x=PC2, y=TotalPhentrans, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("PhenRich") +
  xlab ("Elevation (m)")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")


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





