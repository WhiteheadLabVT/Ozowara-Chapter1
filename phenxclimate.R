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

#Read in Data sheets 
Orchard <- read_csv("Orchard Level Data.csv")
Orchard$orchard.num <- as.factor(as.character(Orchard$orchard.num))

Tree <- read_csv("Tree Level Data.csv")
Tree$orchard.num <- as.factor(as.character(Tree$orchard.num))

d <- read_csv("FruitLevelData_revised - Sheet1.csv")
d$orchard.num <- as.factor(as.character(d$orchard.num))

#creating totalphen and phenrich
d <- d %>%
  mutate(TotalPhen=rowSums(across(8:41)),
         PhenRich=rowSums(across(8:41)!=0))

#testing 
shapiro.test(d$TotalPhen)#not normal
shapiro.test(d$PhenRich)#not normal


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


#Condensing to orchard level data 
#convert weights based on averages and bag weights
Tree<- Tree %>%
  mutate(avgwgt=(apple.wgt.3-bag.weight)/3)%>%
dplyr::select(-c(3:4)) #removing tree & otype collumns  

#summarizing measured traits 
Tree <- Tree %>%
  group_by(orchard.num) %>%
  summarise_at(c("Firmness", "SSC", "maturity.index", "avgwgt")
               , mean, na.rm = TRUE)

#Whole Fruit Orchard Level Data 
c <- d
c <- c %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)
c <- left_join(c, Tree, by="orchard.num")
c <- left_join(c, Orchard, by="orchard.num")
shapiro.test(c$TotalPhen)#normal!
shapiro.test(c$PhenRich)#normal!


#dividing chem data by tissue type to the orchard level  
SkinD <- filter(d, Tissue=="SKIN")
PulpD <- filter(d, Tissue=="PULP")
SeedD <- filter(d, Tissue=="SEED")

#Skin to orchard 
SkinD <- SkinD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)
SkinD <- left_join(SkinD, Tree, by="orchard.num")
SkinD <- left_join(SkinD, Orchard, by="orchard.num")
shapiro.test(SkinD$TotalPhen)#normal!
shapiro.test(SkinD$PhenRich)#normal!

#Pulp to orchard 
PulpD <- PulpD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)
PulpD <- left_join(PulpD, Tree, by="orchard.num")
PulpD <- left_join(PulpD, Orchard, by="orchard.num")
shapiro.test(PulpD$TotalPhen)#not normal?
shapiro.test(PulpD$PhenRich)#normal!


#Seed to orchard 
SeedD <- SeedD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich")
               , mean, na.rm = TRUE)
SeedD <- left_join(SeedD, Tree, by="orchard.num")
SeedD <- left_join(SeedD, Orchard, by="orchard.num")
shapiro.test(SeedD$TotalPhen)#normal!
shapiro.test(SeedD$PhenRich)#normal!


#How do management systems interact with broad climatic changes across latitude?-----
d <- left_join(d, Orchard[,c(2,4)], by="orchard.num")

#Total Phenolics with beta distribution 
tp1 <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Tissue*Latitude + 
(1|site.code/orchard.num/Tree), data=d, family=beta_family(link="logit"))
summary(tp1)
Anova(tp1)
#orchard.type                   4.4968  1    0.03396 *  
#Tissue                       278.1257  2    < 2e-16 ***
#orchard.type:Tissue            8.3481  2    0.01539 *  

#visualize this 
ggplot(d, aes(x=Tissue, y=TotalPhen/1000, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()
#seeds and skin have higher values 
#seeds highest 

ggplot(d, aes(x=Latitude, y=TotalPhen/1000, color=Tissue))+
  geom_smooth(method = "lm") +
  geom_point()
#over latitude we see skin and seeds increase 
#pulp appears to decrease 

ggplot(d, aes(x=Latitude, y=TotalPhen/1000, color=Tissue, shape= orchard.type))+
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
pr1 <- glmmTMB(PhenRich~ orchard.type*Latitude*Tissue + (1|site.code/orchard.num/Tree), data=d, 
               family=poisson(link="log"))
summary(pr1)
Anova(pr1)
#Latitude:Tissue                7.6070  2    0.02229 *  
#Tissue                       798.6850  2    < 2e-16 ***

#Visualize this 
#latitude by tissue
ggplot(d, aes(x=Latitude, y=PhenRich, color=Tissue))+
  geom_smooth(method = "lm") +
  geom_point()
#highest richness in skin 

#otype
ggplot(d, aes(x=Tissue, y=PhenRich, color=orchard.type))+
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

#Which compounds distinguish fruits raised in their respective management systems?----------
#Analysis: NMDS and random forest
###NMDS###
m.NMDS <- metaMDS(d.comp, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS
plot(m.NMDS, type="t")

#Data:     wisconsin(sqrt(d.comp)) 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.1499316 
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
# A, cyanidin_Galactoside, ePicatechin, F, G, I and 3 significant 

#plot 
plot(m1.rf.b.sk)  
getSelectedAttributes(m1.rf.b.sk) #lists all important ones

#[1] "A"                    "F"                    "G"                    "PB1"                 
#[5] "cyanidin_Galactoside" "K"          "U1"                   "I"                   



#performing MANOVA 
d.comp.sk.sel <- data.matrix(d.comp.sk[,getSelectedAttributes(m1.rf.b.sk)])
m1.man.sk <- manova(d.comp.sk.sel ~ d.expl.sk$orchard.type)
summary(m1.man.sk) 
#d.expl.sk$orchard.type   1 0.23681   4.3053      8    111 0.0001535 


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
summary.aov(m1.man.pu)  #follow-up ANOVAs for each individual compound
#0.6069

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
#[1] "E"           "PB2"         "ePicatechin" "U1"          "reynoutrin" 

#Running MANOVAS 
d.comp.se.sel <- data.matrix(d.comp.se[,getSelectedAttributes(m1.rf.b.se)])
m1.man.se <- manova(d.comp.se.sel ~ d.expl.se$orchard.type)
summary(m1.man.se)  
#d.expl.se$orchard.type   1 0.22593   6.5963      5    113 2.013e-05 ***

summary.aov(m1.man.se)  
#E= 0.003836 **
#PB2= 0.001805 **
#epicatechin = 0.0006019 ***
#U1 = 0.004753 **


par(mfrow=c(3,3))
for (i in 1:length(colnames(d.comp.se.sel))){
  d.temp=d.comp.se.sel[,i]
  plot(d.temp ~ d.expl.se$orchard.type, ylab=colnames(d.comp.se.sel)[i])
}
dev.off()




#run GLMS of the signifcant compounds agaisnt latitude 
#Skin 
A.sk <- glmmTMB( A~ orchard.type*Latitude + 
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
ggplot(d.sk, aes(x=Latitude, y=cyanidin_Galactoside/1000, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


#Pulp
G.pu <- glmmTMB( G~ orchard.type*Latitude + 
                      (1|site.code), data=d.pu)
summary(G.pu)
Anova(G.pu)
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


#Which abiotic factors are the most important drivers of fruit chem?------------

#Analysis: principle components analysis followed by PC regression with each variable
#1 for whole fruit, skin, pulp, seed

#can technically use the same first PCA for all 4 different analyses because
#the values of the selected variables are the same

###Principle Components Analysis###
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


###Principle Components Regression###

#Whole Fruit#
pc_clim <- as.data.frame(results$x)
pc_clim <- cbind(pc_clim, c)

#TotalPhen
pc.c.tp <- glmmTMB(TotalPhen ~ PC1 + PC2 + PC3 + PC4 + 
                (1|site.code), data=pc_clim)

summary(pc.c.tp)
#PC2           5415.1     3174.4   1.706    0.088 .  


pc.c.tp.pc2 <- glmmTMB(TotalPhen ~ orchard.type*PC2 + 
                     (1|site.code), data=pc_clim)
summary(pc.c.tp.pc2)
#orchard.typeOrganic     -11490.05    4633.60  -2.480   0.0131 *  
#PC2                       6927.70    4162.34   1.664   0.0960 .   

plot(TotalPhen ~ PC2, data=pc_clim)

results$rotation


#PhenRich
pc.c.pr <- glmmTMB(PhenRich ~ PC1 + PC2 + PC3 + PC4 + 
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


#SKIN#
pc.sk_clim <- cbind(pc_clim, SkinD)

#TotalPhen
pc.sk.tp <- glmmTMB(TotalPhen ~ PC1 + PC2 + PC3 + PC4 + 
                     (1|site.code), data=pc.sk_clim)
summary(pc.sk.tp)
#PC2           5415.1     3174.4   1.706    0.088 .  


pc.sk.tp.pc2 <- glmmTMB(TotalPhen ~ orchard.type*PC2 + 
                         (1|site.code), data=pc.sk_clim)
summary(pc.sk.tp.pc2)
#orchard.typeOrganic     -11490.05    4633.60  -2.480   0.0131 *  
#PC2                       6927.70    4162.34   1.664   0.0960 . 

plot(TotalPhen ~ PC2, data=pc.sk_clim)


#PhenRich
pc.sk.pr <- glmmTMB(PhenRich ~ PC1 + PC2 + PC3 + PC4 + 
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
pc.pu_clim <- cbind(pc_clim, PulpD)

##TotalPhen###
pc.pu.tp <- glmmTMB((TotalPhen/1000000)+0.0001 ~ PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim, family=beta_family (link="logit"))
summary(pc.pu.tp)
#PC2          0.128315   0.055966    2.29   0.0219 *  


pc.pu.tp.pc2 <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type* PC2 + 
                   (1|site.code), data=pc.pu_clim, family=beta_family (link="logit"))
summary(pc.pu.tp.pc2)
#orchard.typeOrganic     -0.26424    0.10203  -2.590   0.0096 ** 


plot(TotalPhen ~ PC2, data=pc.pu_clim)


#PhenRich 
pc.pu.pr <- glmmTMB(PhenRich ~ PC1 + PC2 + PC3 + PC4 + 
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
#nothign 

plot(PhenRich ~ PC1, data=pc.pu_clim)
plot(PhenRich ~ PC2, data=pc.pu_clim)
results$rotation


#SEED#
pc.se_clim <- cbind(pc_clim, SeedD)

#TotalPhen
pc.se.tp <- glmmTMB(TotalPhen ~ PC1 + PC2 + PC3 + PC4 + (1|site.code)
                 , data=pc.se_clim)
summary(pc.se.tp)
#PC2           5415.1     3174.4   1.706    0.088 .  

pc.se.tp.pc2 <- glmmTMB(TotalPhen ~ orchard.type* PC2 + (1|site.code)
                    , data=pc.se_clim)
summary(pc.se.tp.pc2)
#orchard.typeOrganic     -11490.05    4633.60  -2.480   0.0131 *  
#PC2                       6927.70    4162.34   1.664   0.0960 .

plot(TotalPhen ~ PC2, data=pc.se_clim)

results$rotation


#PhenRich 
pc.se.pr <- glmmTMB(PhenRich ~ PC1 + PC2 + PC3 + PC4 +
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


#Which specific management practices are the most important drivers of fruit chem?----

#Total Phenolics
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
mgmt.sk <- glmmTMB(TotalPhen ~ Cultivation + Herbicides + Com_Mul + Mowing +
                      Weed_Mats + Cover_Crops + Acres + (1|site.code), 
                    data=SkinD)
summary(mgmt.sk)
Anova(mgmt.sk)
#Herbicides  3.3487  1    0.06726 .


#Pulp 
mgmt.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ Cultivation + Herbicides + Com_Mul + Mowing +
                      Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                    data=PulpD, family=beta_family(link="logit"))
summary(mgmt.pu)
Anova(mgmt.pu)
#Cover_Crops 3.0087  1    0.08282 .
#Cultivation 4.4598  1    0.03470 *
  

#Seed
mgmt.se <- glmmTMB(TotalPhen ~ Cultivation + Herbicides + Com_Mul + Mowing +
                      Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                    data=SeedD)
summary(mgmt.se)
Anova(mgmt.se)
#Nothing 


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
mgmt.pr.sk <- glmmTMB(PhenRich ~ Cultivation + Herbicides + Com_Mul + Mowing +
                      Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                    data=SkinD)
summary(mgmt.pr.sk)
Anova(mgmt.pr.sk)
#Herbicides  3.8902  1   0.048569 * 


mgmt.pr.sk.h <- glmmTMB(PhenRich ~ orchard.type*Herbicides+ (1|site.code), 
                    data=SkinD)
summary(mgmt.pr.sk.h)
Anova(mgmt.pr.sk.h)
#Herbicides              4.7020  1    0.03013 *


#Pulp 
mgmt.pr.pu <- glmmTMB(PhenRich ~ Cultivation + Herbicides + Com_Mul + Mowing +
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
mgmt.pr.se <- glmmTMB(PhenRich ~ Cultivation + Herbicides + Com_Mul + Mowing +
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

#Which pest or diseases presence has the most significant affect on fruit chem?-----

#first part explained in the PhysicalxClimate Code

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
#Root.Rot       12.4586  1  0.0004161 ***


#How does fruit quality compare to total phenolics and phenolic richness--------
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

