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
library(readr)
library(GGally)
library(emmeans) 
library(ggpubr) #multiplotting 
library(viridis) #color pallete 

#Read Data Organization and Restructuring---------------------------------------
#Read in Data sheets
#Orchard Level Data (24 obs)
Orchard <- read_csv("Ozowara_et_al_OrchardLevelData.csv")
Orchard$orchard.num <- as.factor(as.character(Orchard$orchard.num))

#Tree Level Data (120 obs)
Tree <- read_csv("Ozowara_et_al_TreeLevelData.csv")
View(Tree_Level_Data)
Tree$orchard.num <- as.factor(as.character(Tree$orchard.num))
Tree$Tree <- as.factor(as.character(Tree$Tree))


#Phenolics Data (359 obs)
d <-read_csv("Ozowara_et_al_FruitLevelData.csv")
d$orchard.num <- as.factor(as.character(d$orchard.num))

#Organizing Phenolics Data for Q1#

#creating totalphen and phenrich
d <- d %>%
  mutate(TotalPhen=rowSums(across(8:41)),
         PhenRich=rowSums(across(8:41)!=0))

#tranforms total phenolics for model analysis
d <- d %>%
  mutate(TotalPhenTrans=((TotalPhen/1000000)+0.0001))

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

#creating pest index
Orchard <- Orchard %>%
  mutate(Pest.Index=rowSums(across(21:38))/18)

#categorical latitude data to be used for RF and pest questiosn 
Orchard <- Orchard %>% 
  mutate(lat_cat=cut(Latitude, breaks=c(-Inf, 42, Inf), labels=c("low", "high")))

c <- left_join(Orchard, Tree.sum, by="orchard.num")

#dividing chem data by tissue type to the orchard level  
SkinD <- filter(d, Tissue=="SKIN")
PulpD <- filter(d, Tissue=="PULP")
SeedD <- filter(d, Tissue=="SEED")

#Skin to orchard 
SkinD <- SkinD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich","TotalPhenTrans")
               , mean, na.rm = TRUE)
SkinD <- left_join(SkinD, Tree.sum, by="orchard.num")
SkinD <- left_join(SkinD, Orchard, by="orchard.num")

#Pulp to orchard 
PulpD <- PulpD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich","TotalPhenTrans")
               , mean, na.rm = TRUE)
PulpD <- left_join(PulpD, Tree.sum, by="orchard.num")
PulpD <- left_join(PulpD, Orchard, by="orchard.num")
shapiro.test(PulpD$TotalPhen)

#Seed to orchard 
SeedD <- SeedD %>%
  group_by(orchard.num) %>%
  summarise_at(c("TotalPhen", "PhenRich", "TotalPhenTrans")
               , mean, na.rm = TRUE)
SeedD <- left_join(SeedD, Tree.sum, by="orchard.num")
SeedD <- left_join(SeedD, Orchard, by="orchard.num")


#Q1: How do management systems interact with broad climatic changes across latitude?----
#Analysis: Using GLMMs with tree level data to explore how measured qualities interact
#across latitude 
#Q1-A: Physical Quality------------------------------------------------------------
#binding latitude to the 120 obs data set#
TreeLat <- left_join(Tree, Orchard %>% select(orchard.num, Latitude, lat_cat), 
                     by = c("orchard.num"))

#SSC
Lat1<- glmmTMB(SSC ~ orchard.type*Latitude + maturity.index + (1|site.code/orchard.num), data=TreeLat)
summary(Lat1)
Anova(Lat1) 

plot(SSC ~ maturity.index, data=TreeLat)

#calculate effect size
emmeans(Lat1, ~ Latitude*maturity.index, at = list(Latitude = c(34:48)), type= "reponse")

#Firmness 
Lat2<- glmmTMB(Firmness ~ orchard.type*Latitude + maturity.index + (1|site.code/orchard.num), data=TreeLat)
summary(Lat2)
Anova(Lat2)
#Latitude              34.9002  1   3.47e-09 ***

#calculate effect size
emmeans(Lat2, ~ Latitude, at = list(Latitude = c(34:48)), type= "reponse")

#Average Weight 
Lat3<- glmmTMB(avgwgt ~ orchard.type*Latitude + maturity.index + (1|site.code/orchard.num), data=TreeLat)
summary(Lat3)
Anova(Lat3) 
#orchard.type:Latitude 9.6577  1   0.001886 **

#calculate effect size
emmeans(Lat3, ~ Latitude|orchard.type, at = list(Latitude = c(34:48)), type= "reponse")


###Investigating Latitude and Average weight 
#Splitting latitude in high and low 
Tree_low <- filter(TreeLat, Latitude<42)
Tree_high <- filter(TreeLat, Latitude>42)

#Low Latitude Sites 
avg_low<- glmmTMB(avgwgt ~ orchard.type + maturity.index +(1|site.code/orchard.num), data=Tree_low)
summary(avg_low)
Anova(avg_low) 

#High Latitude Sites 
avg_hgh<- glmmTMB(avgwgt ~ orchard.type + maturity.index + (1|site.code/orchard.num), data=Tree_high)
summary(avg_hgh)
Anova(avg_hgh)

## Calculate the effect size
emmeans(avg_hgh,pairwise~orchard.type, type="response")

#Q1-B: Fruit Chemistry-------------------------------------------------------------
#binding latitude to the 359 obs data set 
ChemLat <- left_join(d, Orchard %>% select(orchard.num, Latitude), 
                     by = c("orchard.num"))
ChemLat <- left_join(ChemLat, TreeLat %>% select(orchard.num, Tree, maturity.index), 
                     by = c("orchard.num", "Tree"))

#Total Phenolics with beta distribution 
tp1 <- glmmTMB(TotalPhenTrans ~ orchard.type*Tissue*Latitude + maturity.index +
                 (1|site.code/orchard.num), data=ChemLat, family=beta_family(link="logit"))
summary(tp1)
Anova(tp1)

#binding latitudinal data to tissue data sets
d.sk <- left_join(d.sk, Orchard[,c(2,4)], by="orchard.num")
d.sk <- left_join(d.sk, TreeLat %>% select(orchard.num, Tree, maturity.index), 
                  by = c("orchard.num", "Tree"))

d.pu <- left_join(d.pu, Orchard[,c(2,4)], by="orchard.num")
d.pu <- left_join(d.pu, TreeLat %>% select(orchard.num, Tree, maturity.index), 
                  by = c("orchard.num", "Tree"))


d.se <- left_join(d.se, Orchard[,c(2,4)], by="orchard.num")
d.se <- left_join(d.se, TreeLat %>% select(orchard.num, Tree, maturity.index), 
                  by = c("orchard.num", "Tree"))

###total phenolics###
#skin
tp.sk <- glmmTMB(TotalPhenTrans ~ orchard.type*Latitude + maturity.index + (1|site.code/orchard.num), 
                 data=d.sk, family=beta_family(link="logit"))
summary(tp.sk)
Anova(tp.sk)
#nothing 

#pulp
tp.pu <- glmmTMB(TotalPhenTrans ~ orchard.type*Latitude + maturity.index + (1|site.code/orchard.num), 
                 data=d.pu, family=beta_family(link="logit"))
summary(tp.pu)
Anova(tp.pu)
#nothing 

#seed
tp.se <- glmmTMB(TotalPhenTrans ~ orchard.type*Latitude + maturity.index +(1|site.code/orchard.num/Tree), 
                 data=d.se, family=beta_family(link="logit"))
summary(tp.se)
Anova(tp.se)
#orchard.type          8.0185  1    0.00463 **
 
#calculate effect size
emmeans(tp.se, ~orchard.type, type="response")


###phenolics richness###
pr1 <- glmmTMB(PhenRich~ orchard.type*Latitude*Tissue + maturity.index + (1|site.code/orchard.num), data=ChemLat, 
               family=poisson(link="log"))
summary(pr1)
Anova(pr1)

#skin
pr.sk <- glmmTMB(PhenRich ~ orchard.type*Latitude+maturity.index+(1|site.code/orchard.num), data=d.sk,
                 family=poisson(link="log"))
summary(pr.sk)
Anova(pr.sk)


#pulp
pr.p <- glmmTMB(PhenRich ~ orchard.type*Latitude+maturity.index+(1|site.code/orchard.num), data=d.pu,
                family=poisson(link="log"))
summary(pr.p)
Anova(pr.p)


#seed
pr.se <- glmmTMB(PhenRich ~ orchard.type*Latitude+maturity.index+(1|site.code/orchard.num), data=d.se,
                 family=poisson(link="log"))
summary(pr.se)
Anova(pr.se)

#maturity.index        10.6801  1   0.001083 **
#orchard.type:Latitude  3.0240  1   0.082042 . 

#calculate effect size
emmeans(pr.se,~maturity.index, type="response")

#Q1-C: NMDS + PERMANOVA---------------------------------------------------------
#left joining the "expl" by latitude  

d.expl <- left_join(d.expl, TreeLat %>% select(orchard.num, Tree, maturity.index,lat_cat), 
          by = c("orchard.num", "Tree"))

d.expl.sk <- left_join(d.expl.sk, TreeLat %>% select(orchard.num, Tree, maturity.index, lat_cat), 
                    by = c("orchard.num", "Tree"))

d.expl.pu <- left_join(d.expl.pu, TreeLat %>% select(orchard.num, Tree, maturity.index,lat_cat), 
                       by = c("orchard.num", "Tree"))

d.expl.se <- left_join(d.expl.se, TreeLat %>% select(orchard.num, Tree, maturity.index, lat_cat), 
                       by = c("orchard.num", "Tree"))

###NMDS for Skin 
m.NMDS.sk <- metaMDS(d.comp.sk, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.sk

plot(m.NMDS.sk, type="t")

#Extract the site scores (coordinates of points in ordination space)
nmds_scores <- vegan::scores(m.NMDS.sk, display = "sites")

# Convert to data frame
data.scores <- as.data.frame(nmds_scores)
data.scores$SampleID <- rownames(data.scores)

# Assuming you have a metadata data frame with 'SampleID' and 'Group' columns
# Example: Creating metadata (replace with your actual metadata)
metadata <- data.frame(SampleID = rownames(d.comp.sk),
                       lat_cat = factor(rep(c("high", "low"), length.out = nrow(d.expl.sk))),
                       orchard.type = factor(rep(c("Organic", "Conventional"), length.out = nrow(d.expl.sk))))

# Merge NMDS scores with metadata
merged_data <- merge(data.scores, metadata, by = "SampleID")

names(merged_data)
# Plot with ggplot2
sk.nmds <- ggplot(merged_data, aes(x = NMDS1, y = NMDS2, color=orchard.type, shape= merged_data$lat_cat)) +
  geom_point(size = 2) +
  theme_classic() +
  scale_color_manual(values = c("#3b0f70", "#de4968"), name = "Management System")+
  scale_shape_manual(values = c(16, 17), name = "Latitudinal Category") +  
  stat_ellipse(aes(group = orchard.type), level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "Skin",
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme(legend.position = "bottom")

print(sk.nmds)


sk.nmds  <- sk.nmds  + guides(color = "none")


#PERMANOVA can test whether the visualized differences are significant
m.perm1 <- adonis2(d.comp.sk~orchard.type*lat_cat.x+maturity.index, data=d.expl.sk)
m.perm1


###NMDS for Pulp 
m.NMDS.pu <- metaMDS(d.comp.pu, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.pu


#Extract the site scores (coordinates of points in ordination space)
nmds_scores <- vegan::scores(m.NMDS.pu, display = "sites")

# Convert to data frame
data.scores <- as.data.frame(nmds_scores)
data.scores$SampleID <- rownames(data.scores)

# Assuming you have a metadata data frame with 'SampleID' and 'Group' columns
# Example: Creating metadata (replace with your actual metadata)
metadata <- data.frame(SampleID = rownames(d.comp.pu),
                       Latitudinal.Category = factor(rep(c("High", "Low"), length.out = nrow(d.expl.pu))),
                       Management.System = factor(rep(c("Organic", " Conventional"), length.out = nrow(d.expl.pu))))

# Merge NMDS scores with metadata
merged_data <- merge(data.scores, metadata, by = "SampleID")

# Plot with ggplot2
pu.nmds <- ggplot(merged_data, aes(x = NMDS1, y = NMDS2, color = Management.System, shape = Latitudinal.Category)) +
  geom_point(size = 2) +
  theme_classic() +
  scale_color_manual(values = c("#3b0f70", "#de4968"), name = "Management System")+
  stat_ellipse(aes(group = Management.System), level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "Pulp",
       x = " ",
       y = " ") +
  theme(legend.position = "bottom")




#PERMANOVA
m.perm2 <- adonis2(d.comp.pu~orchard.type*lat_cat+maturity.index, data=d.expl.pu)
m.perm2


###NMDS for Seed 
m.NMDS.se <- metaMDS(d.comp.se, distance = "bray", trymax=100, autotransform =FALSE)
m.NMDS.se

nmds_scores <- vegan::scores(m.NMDS.se, display = "sites")

# Convert to data frame
data.scores <- as.data.frame(nmds_scores)
data.scores$SampleID <- rownames(data.scores)

# Assuming you have a metadata data frame with 'SampleID' and 'Group' columns
# Example: Creating metadata (replace with your actual metadata)
metadata <- data.frame(SampleID = rownames(d.comp.se),
                       Latitudinal.Category = factor(rep(c("High", "Low"), length.out = nrow(d.expl.se))),
                       Management.System = factor(rep(c("Organic", " Conventional"), length.out = nrow(d.expl.se))))

# Merge NMDS scores with metadata
merged_data <- merge(data.scores, metadata, by = "SampleID")

# Plot with ggplot2
se.nmds <- ggplot(merged_data, aes(x = NMDS1, y = NMDS2, color = Management.System, shape = Latitudinal.Category)) +
  geom_point(size = 2) +
  theme_classic() +
  scale_color_manual(values = c("#3b0f70", "#de4968"), name = "Management System")+
  stat_ellipse(aes(group = Management.System), level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "Seed",
       x = " ",
       y = " ") +
  theme(legend.position = "bottom")


se.nmds  <- se.nmds  + guides(color = "none")


#PERMANOVA
m.perm3 <- adonis2(d.comp.se~orchard.type*lat_cat+maturity.index, data=d.expl.se)
m.perm3



##joining all for multiplot 
ggarrange(sk.nmds, pu.nmds, se.nmds, nrow = 1, ncol = 3, labels = c("A", "B", "C"))
ggsave("fig5.png", width=20, height=24, units="cm", dpi=600)






#Q2: Which abiotic factors are the most important drivers of fruit quality?---------
#Analysis: principle components analysis followed by PC regression with each variable
#create a vector for all of the PCs

p_clim <- dplyr::select(c,c("Prox.Water","elevation", 
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
biplot(results,choices=1:2, cex = 0.001,
       col = c(NA, "BLACK"),
       scale = FALSE, xlabs = rep("*", 24), 
       xlab = "Principal Component 1", 
       ylab = "Principal Component 2", 
       main = "PCA Biplot",
       cex.lab = 1.2, 
       cex.axis = 1.1, 
       cex.main = 1.5,
       xlim = c(-2.1, 2.1), # Adjust these limits to zoom in on the x-axis
       ylim = c(-2.1, 2.1))




#calculate total variance explained by each principal component
summary(results)$importance
summary(results)$importance[2,]

var_explained = results$sdev^2 / sum(results$sdev^2)

df <- data.frame(PC=1:7, var_explained=var_explained)

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
p_clim1 <- dplyr::select(c,c("Prox.Water","elevation", 
                             "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip",
                             "Szn.UVI"))

#The first matrix shows the correlation coefficients between the variables  
#the second matrix shows the corresponding p-values.
rcorr(as.matrix(p_clim1))

#now let's visualize this 
corrplot(cor(p_clim1))

testRes = cor.mtest(p_clim1, conf.level = 0.95)

## add all p-values
corrplot(cor(p_clim1), p.mat = testRes$p, insig = 'p-value', sig.level = -1)

## add significant level stars
corrplot(cor(p_clim1), p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')

#Q2-A: Physical Quality---------------------------------------------------------
##Linear models for PC x quality 
pc_clim_phys <- cbind(pc_clim, c)


###Firmness###
p1 <- glmmTMB(Firmness ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p1)

plot(Firmness ~ PC1, data=pc_clim_phys)
plot(Firmness ~ PC2, data=pc_clim_phys)
plot(Firmness ~ PC3, data=pc_clim_phys)

results$rotation

###SSC###
P2 <- glmmTMB(SSC ~ orchard.type + PC1 + PC2 + PC3 + PC4 + (1|site.code), data= pc_clim_phys)
summary(P2)

plot(SSC ~ PC1, data=pc_clim_phys)
plot(SSC ~ PC2, data=pc_clim_phys)

results$rotation

###AVGWGT
p3 <- glmmTMB(avgwgt ~ orchard.type + PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p3)

plot(avgwgt ~ PC2, data=pc_clim_phys)

results$rotation

###Maturity Index 
p4 <- glmmTMB(maturity.index ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_clim_phys)
summary(p4)

plot(maturity.index ~ PC1, data=pc_clim_phys)

results$rotation

#Q2-B: Fruit Chemistry----------------------------------------------------------
#Linear models for PC X total phenolics and phenolic richness 
#SKIN#
pc.sk_clim <- cbind(pc_clim, SkinD)

#TotalPhen
pc.sk.tp <- glmmTMB(TotalPhenTrans ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.sk_clim, family=beta_family (link="logit"))
summary(pc.sk.tp)

#PhenRich
pc.sk.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.sk_clim)
summary(pc.sk.pr)

#PULP#
pc.pu_clim <- cbind(pc_clim, PulpD)

##TotalPhen###
pc.pu.tp <- glmmTMB(TotalPhenTrans ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim, family=beta_family (link="logit"))
summary(pc.pu.tp)

#PhenRich 
pc.pu.pr <- glmmTMB(PhenRich ~ orchard.type+ PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.pu_clim)
summary(pc.pu.pr)

#SEED#
pc.se_clim <- cbind(pc_clim, SeedD)

#TotalPhen
pc.se.tp <- glmmTMB(TotalPhenTrans ~ orchard.type + PC1 + PC2 + PC3 + PC4 + 
                      (1|site.code), data=pc.se_clim, family=beta_family (link="logit"))
summary(pc.se.tp)

#PhenRich 
pc.se.pr <- glmmTMB(PhenRich ~ orchard.type + PC1 + PC2 + PC3 + PC4 +
                      (1|site.code), data=pc.se_clim)
summary(pc.se.pr)

#Q3 Analysis changes####-----------------------------------------------------
#FAMD-------------------------------------------------------------------------------
library(FactoMineR)
library(factoextra)

#here I'm extracting the management varibales for analysis 

p_mgmt <- dplyr::select(c,c("Com_Mul", "Cover_Crops", "Root.Stock", "Age_Cat", 
                            "Training", "Weed_Mats", "Herbicides", "Acres",
                      "Mowing", "Cultivation"))


#character data needs to be converted into factor data 
p_mgmt$Com_Mul <- as.factor(as.character(p_mgmt$Com_Mul))
p_mgmt$Cover_Crops <- as.factor(as.character(p_mgmt$Cover_Crops))
p_mgmt$Root.Stock <- as.factor(as.character(p_mgmt$Root.Stock))
p_mgmt$Training <- as.factor(as.character(p_mgmt$Training))
p_mgmt$Weed_Mats <- as.factor(as.character(p_mgmt$Weed_Mats))
p_mgmt$Herbicides <- as.factor(as.character(p_mgmt$Herbicides))
p_mgmt$Mowing <- as.factor(as.character(p_mgmt$Mowing))
p_mgmt$Cultivation <- as.factor(as.character(p_mgmt$Cultivation))
p_mgmt$Age_Cat <- as.factor(as.character(p_mgmt$Age_Cat))


#removing unfinished survey from orchard 
p_mgmt <- p_mgmt[-c(17), ]


# Load or prepare your data

# Perform FAMD on the dataset
# Adjust `ncp` to set the number of dimensions/components you want to compute
famd_result <- FAMD(p_mgmt, ncp = 10, graph = TRUE)

FAMD (p_mgmt, ncp = 5, sup.var = NULL, ind.sup = NULL, graph = TRUE)

# View the results summary
summary(famd_result)

#percent contribution by each component
famd_result$var$contrib


# Plot the eigenvalues to determine the explained variance per dimension
fviz_screeplot(famd_result)

# Contributions of variables
# View the contribution of each variable (continuous and categorical) to each dimension
fviz_contrib(famd_result, choice = "var", axes = 1, top = 10) # Contributions to Dim.1

###root stock is contributing largely to the variation 

# Variables plot, showing how each variable contributes to the dimensions
fviz_famd_var(famd_result, repel = TRUE)

# Extract FAMD results for further analysis
# Eigenvalues (variance explained by each dimension)
eigenvalues <- famd_result$eig

famd_scores <- as.data.frame(famd_result$ind$coord)

##GLMMs-------------------------------------------------------------------------
c <- c[-c(17), ]

summary(c$Acres)
hist(c$Acres, breaks = 10, main = "Distribution of Acres", xlab = "Acres")


mean_se(c$Acres)


c2 <- cbind(c, famd_scores)


f1 <- glmmTMB(Firmness ~ orchard.type + Dim.1 + Dim.2 + Dim.3 
              + Dim.4 + Dim.5 + Dim.6 +Dim.7 +
                      (1|site.code), data=c2)
summary(f1)

#Dim.2                -1.7584     0.1421 -12.375  < 2e-16 ***
#Dim.5                -2.9959     0.7310  -4.099 4.16e-05 ***
#Dim.6                 2.8210     0.2732  10.325  < 2e-16 ***
#Dim.7                 0.7361     0.3159   2.330   0.0198 *  

plot(Firmness ~ Dim.2, data=c2)
plot(Firmness ~ Dim.5, data=c2)
plot(Firmness ~ Dim.6, data=c2)
plot(Firmness ~ Dim.7, data=c2)




f2 <- glmmTMB(SSC ~ orchard.type + Dim.1 + Dim.2 + Dim.3 
              + Dim.4 + Dim.5 + Dim.6 +Dim.7 +
                (1|site.code), data=c2)
summary(f2)

#Dim.1                0.61646    0.14393   4.283 1.84e-05 ***
#Dim.2               -0.13967    0.04702  -2.971  0.00297 **

plot(SSC ~ Dim.1, data=c2)
plot(SSC ~ Dim.2, data=c2)


f3 <- glmmTMB(avgwgt ~ orchard.type + Dim.1 + Dim.2 + Dim.3 
              + Dim.4 + Dim.5 + Dim.6 +Dim.7 +
                (1|site.code), data=c2)
summary(f3)

#Dim.2                 -4.569      1.909  -2.393  0.01670 *  
#Dim.3                 10.351      3.224   3.211  0.00132 ** 
#Dim.4                 -8.050      3.248  -2.478  0.01320 *  
#Dim.7                  6.508      3.018   2.156  0.03105 * 


###Chemistry
#SKIN#
#removing incomplete row 
SkinD <- SkinD[-c(9), ]

#binding scores 
Skin1 <- cbind(famd_scores, SkinD)

#TotalPhen
mstp <- glmmTMB(TotalPhenTrans ~ orchard.type+ Dim.1 + Dim.2 + Dim.3 
                    + Dim.4 + Dim.5 + Dim.6 +Dim.7 + 
                      (1|site.code), data=Skin1, family=beta_family (link="logit"))
summary(mstp)

#PhenRich
mspr <- glmmTMB(PhenRich ~ orchard.type + Dim.1 + Dim.2 + Dim.3 
                    + Dim.4 + Dim.5 + Dim.6 +Dim.7 +
                      (1|site.code), data=Skin1)
summary(mspr)

#PULP#
PulpD <- PulpD[-c(9), ]

Pulp1 <- cbind(famd_scores, PulpD)

#TotalPhen
mptp <- glmmTMB(TotalPhenTrans ~ orchard.type+ Dim.1 + Dim.2 + Dim.3 
                + Dim.4 + Dim.5 + Dim.6 +Dim.7 + 
                  (1|site.code), data=Pulp1, family=beta_family (link="logit"))
summary(mptp)

#PhenRich
mppr <- glmmTMB(PhenRich ~ orchard.type + Dim.1 + Dim.2 + Dim.3 
                + Dim.4 + Dim.5 + Dim.6 +Dim.7 +
                  (1|site.code), data=Pulp1)
summary(mppr)


#SEED#
SeedD <- SeedD[-c(9), ]

Seed1 <- cbind(famd_scores, SeedD)

#TotalPhen
msdtp <- glmmTMB(TotalPhenTrans ~ orchard.type+ Dim.1 + Dim.2 + Dim.3 
                + Dim.4 + Dim.5 + Dim.6 +Dim.7 + 
                  (1|site.code), data=Seed1, family=beta_family (link="logit"))
summary(msdtp)

#PhenRich
msdpr <- glmmTMB(PhenRich ~ orchard.type + Dim.1 + Dim.2 + Dim.3 
                + Dim.4 + Dim.5 + Dim.6 +Dim.7 +
                  (1|site.code), data=Seed1)
summary(msdpr)


#Dummy Coded Corroplot----------------------------------------------------------

cheese <- read.csv("Cheese.csv")
names(cheese)


p_mgmt <- dplyr::select(cheese,c("orchard.type", "Com_Mul", "Cover_Crops", "R_Ct", "A_C", 
                            "Tr_Ct", "Weed_Mats", "Herbicides", "Acres",
                            "Mowing", "Cultivation","Prox.Water","elevation", 
                            "Szn.Max.Avg", "Szn.Min.Avg","Szn.Temp.Avg","Szn.Total.Precip","Szn.UVI"))



p_mgmt <- p_mgmt[-c(17), ]

#The first matrix shows the correlation coefficients between the variables  
#the second matrix shows the corresponding p-values.
rcorr(as.matrix(p_mgmt))

#now let's visualize this 
corrplot(cor(p_mgmt))

testRes = cor.mtest(p_mgmt, conf.level = 0.95)

## add all p-values
corrplot(cor(p_mgmt), p.mat = testRes$p, insig = 'p-value', sig.level = -1)

## add significant level stars
corrplot(cor(p_mgmt), p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')


### rootstock has 7 levels so these results shouldnt be trusted too much 
### training system follows the same interpretation 

###Cover Crops: UVI, All temps, acres, "rootstock"
### training x orchard type 
### vultivation x herbicides + com_mulch 
### rootstock x uvi + all temps 






