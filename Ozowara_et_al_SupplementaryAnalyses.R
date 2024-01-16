#Ozowara et al., Supplementary Analysis 
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
#Read Data Organization and Restructuring---------------------------------------
#Read in Data sheets#
#Orchard Level Data (24 obs)
Orchard <- read_csv("Ozowara_et_al_OrchardLevelData.csv")
View(Orchard)
Orchard$orchard.num <- as.factor(as.character(Orchard$orchard.num))

#Tree Level Data (120 obs)
Tree <- read_csv("Ozowara_et_al_TreeLevelData.csv")
View(Tree_Level_Data)
Tree$orchard.num <- as.factor(as.character(Tree$orchard.num))

#Phenolics Data (359 obs)
d <-read_csv("Ozowara_et_al_FruitLevelData.csv")
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

#Q4: Which pest or diseases presence has the most significant affect on fruit quality-----
#visualize pest data in indices 

p <- pivot_longer(data=c, cols=21:38, names_to="pest", values_to="index")

ggplot(p, aes(x=pest, y=index))+
  geom_boxplot()+
  geom_jitter(width=0.2, height=0.1)

#removing everything that doesn't pass 3 on index 
#what were left with: 
#pests: aphids, apple maggots, coddling moth, 
#disease: fireblight, powdery mildew, 

#correlation matrix
p_pest <- dplyr::select(c,c(Aphids, Apple.Maggots, Codling.Moth, 
                            Powdery.mildew, Fire.Blight))

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
#aphinds and apple maggots
#apple maggots and coddling moth 

#Q4-A: Physical Quality---------------------------------------------------------
#Examining Relationships of Pest.Index, orchard.type, and latitude 
#pest index 
#ssc
ssc.index <- glmmTMB(SSC ~ orchard.type*Pest.Index+ (1|site.code), 
                     data=c)
summary(ssc.index)
Anova(ssc.index)

#firmness 
firm.index <- glmmTMB(Firmness ~ orchard.type*Pest.Index+ (1|site.code), 
                      data=c)
summary(firm.index)
Anova(firm.index)
#nothing 

#avg wgt 
wgt.index <- glmmTMB(avgwgt ~ orchard.type*Pest.Index+ (1|site.code), 
                     data=c)
summary(wgt.index)
Anova(wgt.index)

ggplot(c, aes(x=Pest.Index, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

#maturity  
mat.index <- glmmTMB(maturity.index ~ orchard.type*Pest.Index+ (1|site.code), 
                     data=c)
summary(mat.index)
Anova(mat.index)
#Pest.Index              11.4270  1  0.0007238 ***

ggplot(c, aes(x=Pest.Index, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()


###responses to species###
###SSC###
ssc_pest <- glmmTMB(SSC ~ orchard.type + Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Fire.Blight+ 
                      (1|site.code), 
                    data=c)
summary(ssc_pest)
Anova(ssc_pest)
#Fire.Blight    5.3607  1    0.02060 *
#orchard.type   5.6523  1    0.01743 * (organic)


ggplot(c, aes(x=Fire.Blight, y=SSC, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#lower SSC with lower reported fireblight

###avgwgt###
wgt_pest <- glmmTMB(avgwgt ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Fire.Blight+ (1|site.code), 
                    data=c)
summary(wgt_pest)
Anova(wgt_pest)
#Apple.Maggots  4.8605  1    0.02748 *
#Powdery.mildew 3.8439  1    0.04993 *
#Fire.Blight    5.6336  1    0.01762 *


ggplot(c, aes(x=Apple.Maggots, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#weight dcreased as apple maggots increase 

ggplot(c, aes(x=Powdery.mildew, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avgwgt decrease as pow mil increases 

ggplot(c, aes(x=Fire.Blight, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avgwgt decrease as fireblight increases 

ggplot(c, aes(x=Codling.Moth, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avgwgt higher in orchards without codling moth, organic weight higher with high pressure 

###firmness###
frm_pest <- glmmTMB(Firmness ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Fire.Blight+  (1|site.code), 
                    data=c)
summary(frm_pest)
Anova(frm_pest)
#Aphids         4.0144  1    0.04511 *
#Apple.Maggots  3.5262  1    0.06040 .
#Powdery.mildew 3.3710  1    0.06636 .

ggplot(c, aes(x=Aphids, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#frimness inreased as aphids increased in conventional orchards 

ggplot(c, aes(x=Apple.Maggots, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increase in conventional, as apple maggots increased
#decrease in organic 

ggplot(c, aes(x=Powdery.mildew, y=Firmness, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decrease as pow mil increase

###maturity###
mat_pest <- glmmTMB(maturity.index ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                      Powdery.mildew+Fire.Blight+ (1|site.code), 
                    data=c)
summary(mat_pest)
Anova(mat_pest)
#Powdery.mildew 5.1418  1    0.02336 *

ggplot(c, aes(x=Powdery.mildew, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decrease in maturity as pow mil increases 


#Q4-B: Fruit Chemistry ---------------------------------------------------------
#Total Phenolics
#Skin 
pest.index.tp.sk <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Pest.Index +(1|site.code), 
                            data=SkinD, family=beta_family(link="logit"))
summary(pest.index.tp.sk)
Anova(pest.index.tp.sk)
#nothing 

pest.tp.sk <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+ (1|site.code), 
                      data=SkinD, family=beta_family(link="logit"))
summary(pest.tp.sk)
Anova(pest.tp.sk)
#Powdery.mildew 5.0525  1    0.02459 *
#Fire.Blight    3.2125  1    0.07308 .

ggplot(SkinD, aes(x=Powdery.mildew, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher total phen in orchards with higher pow mol 

ggplot(SkinD, aes(x=Fire.Blight, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher total phen in orchards with higher fire.blight


#Pulp
pest.index.tp.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Pest.Index +(1|site.code), 
                            data=PulpD, family=beta_family(link="logit"))
summary(pest.index.tp.pu)
Anova(pest.index.tp.pu)
#nothing 

pest.tp.pu <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+ (1|site.code), 
                      data=PulpD, family=beta_family(link="logit"))
summary(pest.tp.pu)
Anova(pest.tp.pu)
#orchard.type   5.1942  1    0.02266 * 
#Aphids         4.6777  1    0.03056 * 
#Apple.Maggots  4.0006  1    0.04548 * 
#Powdery.mildew 7.5087  1    0.00614 **
#Fire.Blight    6.5499  1    0.01049 * 

ggplot(PulpD, aes(x=Aphids, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#higher total phen in conventional , steady decrease 

ggplot(PulpD, aes(x=Apple.Maggots, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#total phen decrease as applemag increase, steadily higher in conventional 

ggplot(PulpD, aes(x=Powdery.mildew, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increases

ggplot(PulpD, aes(x=Fire.Blight, y=TotalPhen, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increases, higher in conventional 

#Seed 
pest.index.tp.se <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type*Pest.Index +(1|site.code), 
                            data=SeedD, family=beta_family(link="logit"))
summary(pest.index.tp.se)
Anova(pest.index.tp.se)


pest.tp.se <- glmmTMB((TotalPhen/1000000)+0.0001 ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+ (1|site.code), 
                      data=SeedD, family=beta_family(link="logit"))
summary(pest.tp.se)
Anova(pest.tp.se)
#orchard.type   13.5904  1  0.0002273 ***


##PhenRich##
#skin 
pest.index.pr.sk <- glmmTMB(PhenRich~ orchard.type*Pest.Index +(1|site.code), 
                            data=SkinD)
summary(pest.index.pr.sk)
Anova(pest.index.pr.sk)


pest.pr.sk <- glmmTMB(PhenRich ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+ (1|site.code), 
                      data=SkinD)
summary(pest.pr.sk)
Anova(pest.pr.sk)
#Powdery.mildew 14.2281  1  0.0001619 ***

ggplot(SkinD, aes(x=Powdery.mildew, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increase as powmil increases 

#pulp
pest.index.pr.pu <- glmmTMB(PhenRich~ orchard.type*Pest.Index +(1|site.code), 
                            data=PulpD)
summary(pest.index.pr.pu)
Anova(pest.index.pr.pu)


pest.pr.pu <- glmmTMB(PhenRich ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+(1|site.code), 
                      data=PulpD)
summary(pest.pr.pu)
Anova(pest.pr.pu)
#Aphids          7.0735  1   0.007823 ** 
#Apple.Maggots   8.1773  1   0.004242 ** 
#Powdery.mildew 19.5452  1  9.825e-06 ***

ggplot(PulpD, aes(x=Apple.Maggots, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#increase in conventional, decrease in organic 

ggplot(PulpD, aes(x=Aphids, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#o starts low eend shigh, conventional starts high and ends low 

ggplot(PulpD, aes(x=Powdery.mildew, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#sharp increase in o, starts low ends high, both increae 


#Seed
pest.index.pr.se <- glmmTMB(PhenRich~ orchard.type*Pest.Index +(1|site.code), 
                            data=SeedD)
summary(pest.index.pr.se)
Anova(pest.index.pr.se)


pest.pr.se <- glmmTMB(PhenRich ~ orchard.type+Aphids+Apple.Maggots+Codling.Moth+
                        Powdery.mildew+Fire.Blight+(1|site.code), 
                      data=SeedD)
summary(pest.pr.se)
Anova(pest.pr.se)
#orchard.type   13.7879  1  0.0002046 ***
#Aphids          4.4591  1  0.0347154 *  
#Codling.Moth   19.2622  1  1.139e-05 ***


ggplot(SeedD, aes(x=Aphids, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#c starts high ends low, o reverse 

ggplot(SeedD, aes(x=Codling.Moth, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#decrease c,increase o



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
#SSC increase as phen rich decreases 

ggplot(PulpD, aes(x=avgwgt, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#avg wgt increase as phen rich increasses 


ggplot(PulpD, aes(x=maturity.index, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#maturity increase as phen rich decreases 


#pulp avgwgt (figure is not correct)
ggplot(PulpD, aes(x=avgwgt, y=Latitude, color=orchard.type)) + 
  geom_smooth(method = "lm") +
  geom_point() + 
  scale_y_continuous(
    "Latitude", 
    sec.axis = sec_axis(~ . * .25, name = "Phenolic Richness")
  )+
  theme_classic()

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
#SSC increase as phen rich decrease 

ggplot(SeedD, aes(x=avgwgt, y=PhenRich, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()
#AVG wgt increase as phen rich decreases 