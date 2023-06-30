#SRW: This is really moving forward fast, but I think we have a fair number
#of details still to work out for the analyses, and it will take a bit more 
#thought and digesting of the data to really think about the best way to organize
#and present it as clear research questions
#Here are my next suggestions:

#1: I'm curious what happened to the rest of your survey data with the values
#for all the individual pests. I think those should be included in the main dataset
#that you important into R, and potentially also in the PCA. At a minimum you
#import the whole dataset and then calculate the pest index in R so that 
#is reproducible and transparent. But I think for some things it would be much
#better to include all the data. At least separating insects from disease??
#I remember we talked about the pest index, which I think is a good idea, I 
#just think we want, for some things, to also break that up

#2: Go through all your variables when first loading the dataset, and make sure
#they are being treated as the right type of variable. For example, orchard number
#should be a categorical variable, and all the 0/1 variables in your management
#data should also be categorical variables. Right now they are all being treated as 
#numeric, but they are simple yes/no categories. On the other hand, if I remember 
#correctly, the original pest data are on a scale of 0-5, so those should be numeric

#3: Once you get those data all organized, run the PCA + linear models (i.e. the 
#PC regression) for the management practices (parallel to what we have for climate, 
#but feeding all the management data into the PCA)

#4: After thinking more about the scale of your different datasets, and how to
#organize all this in a more strightforward manner, you could consider organizing 
#the results like this:

#Q1: How do management systems (organic vs conventional) shape fruit quality?

#Analyses: one linear mixed models for each fruit quality trait, using "tree level" data
#----and orchard.type as a predictor (basically part 1 of your Q1 below, but
#----with tree level data).  

#OR, alternatively: #Q1: How do management systems (organic vs conventional) 
#interact with broad climatic changes across latitude to shape fruit quality?

#Analyses: same as above but including latitude as a predictor (see my lengthy 
#---comments on this below)



#Q2: Which specific pest pressures and management practices are the most important 
#drivers of fruit quality

#Analyses: PC regression with all the survey data (management practices and pest 
#---pressures feeding into PCA)
#----------here we will need to see what the results look like...I struggle a little
#----------with the fact that "management practices" and "pest pressures" are really
#----------two different types of factors, but they are also quite intertwined


#Q3: Which abiotic factors are the most important drivers of fruit quality?

#Analyses: PC regression with all the climate factors feeding into PCA


#A few other things, like the full correlation matrix for all your variables,
#could always go in the supplement

#this will also depend on what happens with the phenolics data!! I think the
#questions should ultimately be combined so they just say "fruit chemistry and quality"
#and then you include parallel results for chemistry and for quality under
#each question




rm(list=ls()) # clears work space

###install packages-------------------------------------------------------------
library(ggplot2)
library(MASS)
library(reshape2)
library(tidyverse)
library(tidyverse)
library(dplyr)
library(randomForest)
library(glmmTMB) #mixed models
library(car) #mixed models summary tables
library(vegan) #multivariate stats
library(Boruta) #random forest models
library("purrr")

#Read in data-------------------------------------------------------------------
library(readr)
Tree <- read_csv("Tree Level Data.csv")

#SRW: removed blank rows in csv file so it doesn't read in with 999 rows

View(Tree)
names(Tree)
library(readr)
Orchard <- read_csv("Orchard Level Data.csv")
View(Orchard)
names(Orchard)

#SRW: moved this calculation so that you are doing it on the tree dataset
#before summarizing
Tree <- Tree %>%
  mutate(avgwgt=(apple.wgt.3-bag.weight)/3)


#SRW: Also re-naming the summarized version, so you can keep the tree-level
#data for use in other analyses below
Tree.sum <- Tree %>%
  group_by(orchard.num) %>%
  summarise_at(c("avgwgt", "Firmness", "SSC", "maturity.index", "bag.weight"), mean, na.rm = TRUE)

c <- left_join(Tree.sum, Orchard, by="orchard.num")


view(c)

#PCA----------------------------------------------------------------------------
#principle components analysis 
#first we're just going to look at the climatic variables and their affect

#create practice data set 
p <- c
view(p)
names(p)

p_qual <- dplyr::select(p, c("avgwgt", "Firmness", "SSC", "maturity.index"))
p_clim <- dplyr::select(p,c("Prox.Water", "Longitude","Latitude","elevation"))  

##SRW:  why not have the actual  climate data in there (e.g. temp, precip)??? 
#For me those are the most important and interesting


#calculate principal components
results <- prcomp(p_clim, scale = TRUE)

#this gives you the % variance explained
summary(results)  

#display principal components
results$x

results$rotation
#display the first six scores
head(results$x)
#PC2 is the largest (longitude)


#this plots the results of the PCAs into a two dimensional representation 
biplot(results, scale = 0)


#calculate total variance explained by each principal component
summary(results)$importance
summary(results)$importance[2,]

var_explained = results$sdev^2 / sum(results$sdev^2)

df <- data.frame(PC=1:4, var_explained=var_explained)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


##SRW: Here the next step could be to construct linear models that look at how the
#PCs affect the quality variables

pc_dat <- as.data.frame(results$x)
pc_dat <- cbind(pc_dat, c)

hist(c$Firmness)
hist(c$SSC)
hist(c$avgwgt)

m1 <- glmmTMB(Firmness ~ PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_dat)
summary(m1)

#strong positive effects of PC1 + PC2 and strong negative effect of PC3
plot(Firmness ~ PC1, data=pc_dat)
plot(Firmness ~ PC2, data=pc_dat)
plot(Firmness ~ PC3, data=pc_dat)

#So, to interpret this you go back to look at the loadings on each PC axis.
results$rotation
#these tell you how your original variables are correlated with each axis

#For example, along PC1 lat and long are both increasing and elevation and 
#proximity to water are both decreasing
#So, firmness increases as you go north, west, down, and away from water


##One other thing to consider with these models is whether you want to also
#include management system?? If you really want to get at interactions between
#management system and climate this could be an alternative place to do it

m1 <- glmmTMB(Firmness ~ orchard.type*PC1 + orchard.type*PC2 + orchard.type*PC3 + 
                orchard.type*PC4 + (1|site.code), data=pc_dat)
summary(m1)

#Or, more simply could just include orchard.type as another factor without
#assessing interactions (or simplify to this based on no interactions above)
m1 <- glmmTMB(Firmness ~ orchard.type + PC1 + PC2 + PC3 + PC4 + (1|site.code), data=pc_dat)
summary(m1)


###Quality###

#SRW: I don't think you need to run a PCA for the quality data, unless you want
#to use fruit quality to predict something else? Normally we just do a PCA on 
#groups of predictor variables. For the quality data, I would just 
#look at each quality variable separately as a response variable


results <- prcomp(p_qual, scale = TRUE)

#this gives you the % variance explained
summary(results)  

#display principal components
results$x

results$rotation
#display the first six scores
head(results$x)
#PC2 is the largest (longitude)


#this plots the results of the PCAs into a two dimensional representation 
biplot(results, scale = 0)


#calculate total variance explained by each principal component
summary(results)$importance
summary(results)$importance[2,]

var_explained = results$sdev^2 / sum(results$sdev^2)

df <- data.frame(PC=1:4, var_explained=var_explained)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


#SRW: Next I think you do want to run a PCA + linear models on all the management 
#variables, parallel to the procedure with climate variables so you can look 
#at how all of those affect fruit chem 




#Correlation Matrix-------------------------------------------------------------
#we're going to round out looking at the influence of our variables by creating a,
#correlation matrix. This will help to pin point which interactions to look at, 
#during analysis 

#first create a practice data set
#were going to alter this table for fruit quality analysis 
p1 <- c
#filter out variables that we're not going to look at 
p1 <- p1 %>% dplyr::select(-c("orchard.num", "bag.weight", "apple.wgt.3","site.code", "orchard.type"))

view(p1)

library(Hmisc)
#The first matrix shows the correlation coefficients between the variables  
#the second matrix shows the corresponding p-values.
rcorr(as.matrix(p1))
#any value below 0.05 is not a statistically significant relationship 
# need to sort through and organzie these values 

#now let's visualize this 
library(corrplot)
corrplot(cor(p1))
#red = negative, blue = positive 
#big = high, small = low


##SRW: I think this is very helpful to look at to get to know your data, but 
#might not go in the final paper. You could maybe include in the supplement, but
#it ultimately gets at similar questions to the other analyses, but in a less
#robust way because you can't account for colinearity or random effects



#Question 1---------------------------------------------------------------------
#How do management systems interact with broad antibiotic conditions to shape fruit quality
#we'll divide this into two parts
#1: management systems alone
#2: proxy climatic factors (lat, long, elev)


#SRW: For this analysis you will want to use your tree level dataset!!!
#You already have the model coded to account for repeated samples within each
#orchard in the random effects, so just swap out the dataset you are feeding it
#for the tree level data


Tree$orchard.num <- as.character(Tree$orchard.num) #make sure you also followed my comments above so that "Tree" is still the full tree level dataset, should have 120 observations
Orchard$orchard.num <- as.character(Orchard$orchard.num)

#Adding orchard type and site code data to the tree dataset, I would actually
#suggest having those in the csv file from the start. Also adding latitude in here
#so we can play around with it (see below)
Tree2 <- left_join(Tree, Orchard[,c(1:2,3, 13)], by="orchard.num")
 

###part 1### 
m1<- glmmTMB(SSC ~ orchard.type + (1|site.code/orchard.num), data=Tree2)
summary(m1)
Anova(m1) #p=0.294

m2<- glmmTMB(Firmness ~ orchard.type + (1|site.code/orchard.num), data=Tree2)
summary(m2)
Anova(m2) #p=0.9851

m3<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/orchard.num), data=Tree2)
summary(m3)
Anova(m3) #p=0.8712

m4<- glmmTMB(maturity.index ~ orchard.type + (1|site.code/orchard.num), data=Tree2)
summary(m4)
Anova(m4) #p=0.6817


#SRW: ALTERNATIVELY, we could also include latitude in these models. This is
#what I was envisioning initially when we talked about analyses, but I didn't 
#think through then as deeply about the different scales of the data 
m3<- glmmTMB(avgwgt ~ orchard.type*Latitude + (1|site.code/orchard.num), data=Tree2)
summary(m3)
Anova(m3)

#I think this would have to come with the caveat that the latitude data were recorded
#at the orchard scale and are just repeated for each tree, but at least for this
#variable we can say with certainty that the actual variation in latitude from
#tree to tree are miniscule compared to the variation from site to site, so 
#if we had the "true values" for each tree, the results would be highly unlikely
#to be different. That is much harder to argue for some of your other orchard
#level variables, especially something like pest damage, but even microclimates
#can vary a lot over short distances

#If we do this, I would suggest ONLY doing it for latitude (not longitude and
#and elevation), largely for simplicity. To me it makes sense to focus on latitude
#here because your study was designed to capture broad variation across latitude, we
#have a huge latitudinal gradient, and it is the best variable we have to capture
#broad scale variation in climate (just like organic/conventional is the best
#variable we have to capture broad scale variation in management practices)







#aggregate the data so that we have some idea of the differences between orchard types
aggregate(SSC~orchard.type, data=c, FUN=mean)
#Conventional 10.88182
#Organic 11.78462

aggregate(Firmness~orchard.type, data=c, FUN=mean)
#Conventional 22.26473
#Organic 22.44308

aggregate(avgwgt~orchard.type, data=c, FUN=mean)
#Conventional 100.29091
#Organic 98.67692

aggregate(maturity.index~orchard.type, data=c, FUN=mean)
#Conventional 3.781818
#Organic 4.138462

###part 2 using Proxy Climatic Variables
#values not correct!

##SRW: for all these models, if you want to look at the interaction between
#climate and orchard type, we would need to use the orchard level data. So, you
#would not need the orchard number as part of the random effects structure
#(because you only have one sample per orchard)

#I went ahead and did that because I wanted to check the results and to be honest
#since we are not finding much interesting here I am not sure it is worth
#including these interactive models (at least not in the main text). 


###SSC###
s1 <- glmmTMB(SSC ~ Latitude*orchard.type + (1|site.code), data=c)
summary(s1)
Anova(s1)#Latitude p=1.906e-05 

s2<- glmmTMB(SSC ~ elevation*orchard.type + (1|site.code), data=c)
summary(s2)
Anova(s2)

s3<- glmmTMB(SSC ~ Longitude*orchard.type + (1|site.code), data=c)
summary(s3)
Anova(s3)

###Firmness###
f1 <- glmmTMB(Firmness ~ Latitude*orchard.type + (1|site.code), data=c)
summary(f1)
Anova(f1)#Latitude p=4.112e-09 ***


f2<- glmmTMB(Firmness ~ elevation*orchard.type + (1|site.code), data=c)
summary(f2)
Anova(f2)

f3<- glmmTMB(Firmness ~ Longitude*orchard.type + (1|site.code), data=c)
summary(f3)
Anova(f3)

###Average Weight
w1 <- glmmTMB(avgwgt ~ Latitude*orchard.type + (1|site.code), data=c)
summary(w1)
Anova(w1)#Latitude:orchard.type p=0.001822 **

##SRW: Here with this strong interaction, a typical next step would be to explore 
#that interaction. First by plotting with both factors

ggplot(c, aes(x=Latitude, y=avgwgt, color=orchard.type))+
  geom_point() +
  geom_smooth(method="lm")

#then you can split the data into subgroups to explore the interaction. You could
#split by orchard type and look at the effects of latitude
#separately for organic and conventional

w1 <- glmmTMB(avgwgt ~ Latitude, data=c[which(c$orchard.type=='Organic'),])
summary(w1)
Anova(w1)

w1 <- glmmTMB(avgwgt ~ Latitude, data=c[which(c$orchard.type=='Conventional'),])
summary(w1)
Anova(w1)

#Or bin the data by latitude and look at the effects of orchard type separately
#for high and low latitudes

w1 <- glmmTMB(avgwgt ~ orchard.type, data=c[which(c$Latitude<42),])
summary(w1)
Anova(w1)

w1 <- glmmTMB(avgwgt ~ orchard.type, data=c[which(c$Latitude>42),])
summary(w1)
Anova(w1)


#No significant effects when separating the data into groups. Based on the plot
#you can see why there is an interaction, it seems weight decreases with latitude,
#but only for organic. But the pattern is not strong enough to hold up when
#we bin the data by orchard type. This leaves the interpretation kind of murky,
#for me I don't think we have strong enough effects or big enough sample size 
#to know whether something real is going on here. Also I don't have a clear 
#biological hypothesis in mind for why we might expect an interaction here



w2<- glmmTMB(avgwgt ~ elevation*orchard.type + (1|site.code), data=c)
summary(w2)
Anova(w2) #elevation p=0.02534 *

w3<- glmmTMB(avgwgt ~ Longitude*orchard.type + (1|site.code), data=c)
summary(w3)
Anova(w3)


#Maturity  
mt1 <- glmmTMB(maturity.index ~ Latitude*orchard.type + (1|site.code), data=c)
summary(mt1)
Anova(mt1)#Latitude 0.05779 .

mt2<- glmmTMB(maturity.index ~ elevation*orchard.type + (1|site.code), data=c)
summary(mt2)
Anova(mt2)

mt3<- glmmTMB(maturity.index ~ Longitude*orchard.type + (1|site.code), data=c)
summary(mt3)
Anova(mt3)


###traits influence on each other 

##SRW: In linear models, we can only look at one response variable at a time. 
#so you cannot add multiple variables to the response side of this equation
#the way you wrote these, it is just adding together all the values for 
#firmness and SSC and then using that as a response variable, which is nonsensical
#It's the same as doing this:

add <- c$SSC + c$Firmness
add

sxf<- glmmTMB( add ~ orchard.type + (1|site.code/orchard.num), data=c)
summary(sxf)
Anova(sxf)


sxf<- glmmTMB(SSC+Firmness ~ orchard.type + (1|site.code/orchard.num), data=c)
summary(sxf)
Anova(sxf)


sxw<- glmmTMB(SSC+avgwgt~orchard.type + (1|site.code/orchard.num), data=c)
summary(sxw)
Anova(sxw)

sxm<- glmmTMB(SSC+maturity.index~orchard.type + (1|site.code/orchard.num), data=c)
summary(sxm)
Anova(sxm)

wxf<- glmmTMB(avgwgt+Firmness ~orchard.type+ (1|site.code/orchard.num), data=c)
summary(wxf)
Anova(wxf)


wxm<- glmmTMB(avgwgt+maturity.index ~orchard.type+ (1|site.code/orchard.num), data=c)
summary(wxf)
Anova(wxf)


fxm<- glmmTMB(maturity.index+Firmness ~orchard.type+ (1|site.code/orchard.num), data=c)
summary(fxm)
Anova(fxm)



##If you want to explore relationships among the quality traits, I think you 
#can get a sense for that from the correlation matrix you did above. If you 
#have a biological hypothesis about why these traits are related and you
#really want to test that more rigorously you could use a linear model but keep
#in mind that with linear models you need to specify a predictor and a response
#so you need to have some logical reason why you think one trait is affecting
#the other and not vice versa



#Question 1 Figures---------------------------------------------------------------
###mgmt and physical quality alone 
#SRW: For these figures, I would suggest plotting the tree level data (one point 
#for each tree) since that is what you will use in the analysis .   
#Also, the legend is currently not showing anything different from the x-axis, 
#so you could get rid of that. 

b1 <- ggplot(c, aes(x=orchard.type, y=SSC, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("SOLUBLE SUGAE CONTENT") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("Conventional", "Organic"))
b1

b2 <- ggplot(c, aes(x=orchard.type, y
                    =Firmness, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Firmness (N)") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("Conventional", "Organic"))
b2

b3 <- ggplot(c, aes(x=orchard.type, y=avgwgt, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Average Apple Weight (g)") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("Conventional", "Organic"))
b3

b4 <- ggplot(c, aes(x=orchard.type, y=maturity.index, color=orchard.type))+
  theme_classic() +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.2))+
  ylab ("Maturity Index") +
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")+
  scale_x_discrete(labels=c("Conventional", "Organic"))
b4

multiplot(b1,b2,b3,b4, cols=2)

###physical quality, mgmt, and proxy climatic variables 

ssc1 = ggplot(c, aes(x=Latitude, y=SSC, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("SSC") +
  xlab ("Latitude")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
ssc1


ssc2 = ggplot(c, aes(x=elevation, y=SSC, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("SSC") +
  xlab ("Elevation (m)")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
ssc2


avgwgt1 = ggplot(c, aes(x=Latitude, y=avgwgt, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("Average Weight (g)") +
  xlab ("Latitude")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
avgwgt1


##SRW: I do like these plots, but depending on whether we decide to keep the 
#analyses looking at how latitude interacts with management system we might not
#need them




#Question 2---------------------------------------------------------------------
#Which abiotic factors are the most important drivers of fruit quality?

###total precipitation x otype###
pl1<- glmmTMB(SSC ~ Szn.Total.Precip*orchard.type + (1|site.code/orchard.num), data=c)
summary(pl1)
Anova(pl1)
#Szn.Total.Precip              7.3440  1   0.006729 **

pl2<- glmmTMB(Firmness ~ Szn.Total.Precip*orchard.type + (1|site.code/orchard.num), data=c)
summary(pl2)
Anova(pl2)


pl3<- glmmTMB(avgwgt ~ Szn.Total.Precip*orchard.type + (1|site.code/orchard.num), data=c)
summary(pl3)
Anova(pl3) #Szn.Total.Precip p=0.007601 **


pl4<- glmmTMB(maturity.index ~ Szn.Total.Precip*orchard.type + (1|site.code/orchard.num), data=c)
summary(pl4)
Anova(pl4)#Szn.Total.Precip p=0.006861 **


###avgerage seasonal temp x otype### 
atl<- glmmTMB(SSC ~ Szn.Temp.Avg*orchard.type + (1|site.code/orchard.num), data=c)
summary(atl)
Anova(atl)#Szn.Temp.Avg p=2.597e-07 ***


at2<- glmmTMB(Firmness ~ Szn.Temp.Avg*orchard.type + (1|site.code/orchard.num), data=c)
summary(at2)
Anova(at2)#Szn.Temp.Avg p=0.0007365 ***


at3<- glmmTMB(avgwgt ~ Szn.Temp.Avg*orchard.type + (1|site.code/orchard.num), data=c)
summary(at3)
Anova(at3)#Szn.Temp.Avg:orchard.type p=0.002653 **


ggplot(c, aes(x=Szn.Temp.Avg, y=avgwgt, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("Average Weight (g)") +
  xlab ("Avg Temp")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

#SRW: potentially interesting--positive effect of temp on weight but only for organic??


at4<- glmmTMB(maturity.index ~ Szn.Temp.Avg*orchard.type + (1|site.code/orchard.num), data=c)
summary(at4)
Anova(at4)#Szn.Temp.Avg p=0.005154 **

###proximity to water x otype###
pxl<- glmmTMB(SSC ~ Prox.Water*orchard.type + (1|site.code/orchard.num), data=c)
summary(pxl)
Anova(pxl)


px2<- glmmTMB(Firmness ~ Prox.Water*orchard.type + (1|site.code), data=c)
summary(px2)
Anova(px2)


px3<- glmmTMB(avgwgt ~ Prox.Water*orchard.type + (1|site.code/orchard.num), data=c)
summary(px3)
Anova(px3)


px4<- glmmTMB(maturity.index ~ Prox.Water*orchard.type + (1|site.code), data=c)
summary(px4)
Anova(px4)

###season avg max x otype###
smx1<- glmmTMB(SSC ~ Szn.Max.Avg*orchard.type + (1|site.code), data=c)
summary(smx1)
Anova(smx1) #Szn.Max.Avg p=1.375e-05 ***
 

smx2<- glmmTMB(Firmness ~ orchard.type*Szn.Max.Avg + (1|site.code/orchard.num), data=c)
summary(smx2)
Anova(smx2) #Szn.Max.Avg p=0.001009 **
 

smx3<- glmmTMB(avgwgt ~ orchard.type*Szn.Max.Avg + (1|site.code/orchard.num), data=c)
summary(smx3)
Anova(smx3) #orchard.type:Szn.Max.Avg p=0.0008262 ***


ggplot(c, aes(x=Szn.Max.Avg, y=avgwgt, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("Average Weight (g)") +
  xlab ("Max Temp")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")

#SRW: basically same result as avg temp





smx4<- glmmTMB(maturity.index ~ orchard.type*Szn.Max.Avg + (1|site.code/orchard.num), data=c)
summary(smx4)
Anova(smx4)#Szn.Max.Avg p=0.008311 **


###season avg min x otype###
sm1<- glmmTMB(SSC ~ orchard.type*Szn.Min.Avg + (1|site.code/orchard.num), data=c)
summary(sm1)
Anova(sm1)#Szn.Min.Avg p=6.414e-05 ***

sm2<- glmmTMB(Firmness ~ orchard.type*Szn.Min.Avg + (1|site.code/orchard.num), data=c)
summary(sm2)
Anova(sm2)#Szn.Min.Avg p=0.002267 **
 

sm3<- glmmTMB(avgwgt ~ orchard.type*Szn.Min.Avg + (1|site.code), data=c)
summary(sm3)
Anova(sm3) #orchard.type:Szn.Min.Avg p=0.02164 *

sm4<- glmmTMB(maturity.index ~ orchard.type*Szn.Min.Avg + (1|site.code), data=c)
summary(sm4)
Anova(sm4)

#Question 2 Figures ------------------------------------------------------------

#looking at fixed variables with avgwgt. It seems that as temps increase, 
#size increases in all but more with organic apples 
fixed1 = ggplot(c, aes(x=Szn.Max.Avg, y=avgwgt, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("SSC") +
  xlab ("Max Temp Avg")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
fixed1

fixed2 = ggplot(c, aes(x=Szn.Min.Avg, y=avgwgt, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("SSC") +
  xlab ("Min Temp Avg")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
fixed2

fixed3 = ggplot(c, aes(x=Szn.Temp.Avg, y=avgwgt, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("SSC") +
  xlab ("Temp Avg")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
fixed3


#Question 3---------------------------------------------------------------------
#Which specific management practices are the most important drivers of fruit quality?
##pest pressure## 

##SRW: For all of these models, since you are using the orchard-level data, you
#would take orchard ID and tree ID out of your random effects structure

##SRW: It is great to look at all these, but ultimately it would be redundant with
#the PC regression, and the PC regression might be the better approach because
#it looks at the correlated variables together and reduces the number of 
#independent statistical tests you are doing (people start to get uncomfortable)
#interpreting p-values when you present hundreds of them...)


pest1 <- glmmTMB(SSC ~ Pest_Index*orchard.type + (1|site.code), data=c)
summary(pest1)
Anova(pest1)
#Pest_Index:orchard.type 2.9152  1    0.08775 .

pest2 <- glmmTMB(avgwgt ~ Pest_Index*orchard.type + (1|site.code), data=c)
summary(pest2)
Anova(pest2)
# none

pest3 <- glmmTMB(Firmness ~ Pest_Index*orchard.type + (1|site.code), data=c)
summary(pest3)
Anova(pest3)
#none 

pest4 <- glmmTMB(maturity.index ~ Pest_Index*orchard.type + (1|site.code), data=c)
summary(pest4)
Anova(pest4)
#none 

#SRW: strong effect of pest index on maturity index??

ggplot(c, aes(x=Pest_Index, y=maturity.index))+
  geom_point()+
  geom_smooth(method="lm")

##Hmnnn, this is weird because the model coefficients show a positive effect
#but overall across the whole dataset the trend is negative

summary(glmmTMB(maturity.index ~ Pest_Index + (1|site.code), data=c))
#result remains strongly positive when you take out the orchard type

summary(glmmTMB(maturity.index ~ Pest_Index, data=c))
#result disappears when you remove site code as a random effect. this suggests
#the within site comparisons are driving the result

ggplot(c, aes(x=Pest_Index, y=maturity.index, color=site.code))+
  geom_point(aes(shape=orchard.type))+
  geom_smooth(method="lm")
#This shows that although the general trend between pest index and maturity index
#is negative, when you compare within individual site pairs, the ones with
#higher pest index also had higher maturity! This is kind of cool, suggests
#that once you control for broad differences in climate, that plants with
#more pest pressure are ripening more quickly. This could be driven by
#jasmonate-mediated pathways (induced by insects and also control ripening)


##Herbicides## 
herb1 <- glmmTMB(SSC ~ Herbicides*orchard.type + (1|site.code), data=c)
summary(herb1)
Anova(herb1)
#Herbicides              5.4980  1    0.01904 *
#orchard.type            2.7336  1    0.09826 .
#Herbicides:orchard.type 1.2220  1    0.26896  


ggplot(c, aes(x=Herbicides, y=SSC, color=orchard.type))+
  geom_point()+
  geom_smooth(method="lm")

##Plotting makes it obvious that this is treating herbicides as a numeric
#variable when really it should be categorical
#fixing this here, but it should be fixed at the very start of the code when 
#you are loading the datasets

herb1 <- glmmTMB(SSC ~ as.character(Herbicides)*orchard.type + (1|site.code), data=c)
summary(herb1)
Anova(herb1)

ggplot(c, aes(x=as.character(Herbicides), y=SSC, color=orchard.type))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=.2))

##SRW: a lot to unpack here, pretty interesting! The strongest result seems to be
#that sugar content is higher for orchards that use herbicides. What I am trying
#to wrap my head around is how that interacts (or doesn't) with 
#the organic/conventional. There is no significant interaction but the effect
#of herbicides seems to be mainly driven by the conventional. And I guess my 
#bigger question is what herbicides are organic orchards
#using?? I thought organic growers mainly relied on alternative methods, like
#cover crops, cultivation, or mulching. so I'm surprised you have so many 
#organic orchards that say they use herbicides



herb2 <- glmmTMB(avgwgt ~ Herbicides*orchard.type + (1|site.code), data=c)
summary(herb2)
Anova(herb2)
#none

herb3 <- glmmTMB(Firmness ~ Herbicides*orchard.type + (1|site.code), data=c)
summary(herb3)
Anova(herb3)
#none 

herb4 <- glmmTMB(maturity.index ~ Herbicides*orchard.type + (1|site.code), data=c)
summary(herb4)
Anova(herb4)
#none 


##Acres## 
Acres <- glmmTMB(SSC ~ Acres*orchard.type + (1|site.code), data=c)
summary(Acres)
Anova(Acres)
#none

Acres <- glmmTMB(avgwgt ~ Acres*orchard.type + (1|site.code), data=c)
summary(Acres)
Anova(Acres)
#Acres              10.8040  1   0.001013 **

Acres <- glmmTMB(Firmness ~ Acres*orchard.type + (1|site.code), data=c)
summary(Acres)
Anova(Acres)
#none 

Acres <- glmmTMB(maturity.index ~ Acres*orchard.type + (1|site.code), data=c)
summary(Acres)
Anova(Acres)
#none 

ggplot(c, aes(x=as.character(Herbicides), y=maturity.index, color=orchard.type))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=.2))

##Very fascinating...I would like to 


##Mowing## 
Mowing <- lm(SSC ~ Mowing*orchard.type, data=c)
summary(Mowing)
Anova(Mowing)
#none

Mowing <- lm(avgwgt ~ Mowing*orchard.type, data=c)
summary(Mowing)
Anova(Mowing)
#none 

Mowing <- lm(Firmness ~ Mowing*orchard.type, data=c)
summary(Mowing)
Anova(Mowing)
#none 

Mowing <- lm(maturity.index ~ Mowing*orchard.type, data=c)
summary(Mowing)
Anova(Mowing)
#none 

#Question 3 Figures-------------------------------------------------------------
mgmt1 = ggplot(c, aes(x=Pest_Index, y=SSC, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("SSC") +
  xlab ("Pest Pressure")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
mgmt1


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





