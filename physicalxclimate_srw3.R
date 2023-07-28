##SRW: Lots of interesting results here!! Definitely getting closer!!

#Next steps:

#1) Re-do/simplify Q1. See my comments below. There are some issues with the 
#----analyses as they are and I think you need to simplify the question structure

#2) I think the PC regressions for Q2 are looking good. You can start thinking
#------a bit more about how best to summarize and present that in the paper
#------I suggest looking at the literature to see how people present
#------graphs and write up the results for similar analyses
#------we can also talk more about this

#3) For Q3, the PC regressions will not make sense because they are not continuous
#-------variables. I did not quite register that before! But instead you can
#-------just do a multiple regression with the raw management variables I think
#-------that will actually make it easier to interpret. See below

#4) For the pest index part of things, I think there is potentially some very
#------interesting stuff there, but again the PC regression is not working 
#------out that well. So this needs some more thought. I suggested
#------some things below



#fruit quality analysis 
rm(list=ls())# clears work space
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
library(readr)
#Read in data-------------------------------------------------------------------
Tree <- read_csv("Tree Level Data.csv")
Orchard <- read_csv("Orchard Level Data.csv")

#convert columns into factors 

#SRW: No! L56-63 are actually converting what were Y/N characters into numeric!! 
#It looks like you changed the actual data from 0/1 to Y/N, so when you do this
#it will automatically import the data as characters
#the easiest way to see what class each variable is is just look at the dataframe 
#and hover over the column name. try doing that before and after you run one of 
#the lines below and you will see it changes it to a 0/1 numeric
#You can also check the class of any object with class(Orchard$Com_Mul)

#Orchard$Com_Mul <- Orchard$Com_Mul %>% as.factor() %>% as.numeric()
#Orchard$Cultivation <- Orchard$Cultivation %>% as.factor() %>% as.numeric()
#Orchard$Herbicides <- Orchard$Herbicides %>% as.factor() %>% as.numeric()
#Orchard$Fire_Mgmt <- Orchard$Fire_Mgmt %>% as.factor() %>% as.numeric()
#Orchard$Weed_Mats <- Orchard$Weed_Mats %>% as.factor() %>% as.numeric()
#Orchard$Cover_Crops <- Orchard$Cover_Crops %>% as.factor() %>% as.numeric()
#Orchard$Mowing <- Orchard$Mowing %>% as.factor() %>% as.numeric()
#view(Orchard)


#SRW: You actually don't need any of the above lines. The only thing you need to
#really need to worry about before running the analyses 
#below is that orchard number should be a categorical variable

#This is a super important thing to think about and take in. It is the first thing
#you want to ask yourself about any variable before you start any analysis, because
#it affects everything about what analyses apply and how you interpret results. 
#Is it a category or a continuous variable??? 
#In R, categorical variables can be specified as characters or factors. There are
#some subtle differences you should read about, but in most basic analyses they
#will be treated the same. Mostly you need to worry about categorical things
#being specified as numeric when they don't actually represent continuous values


Orchard$orchard.num <- as.factor(as.character(Orchard$orchard.num))


#combine data sheets
View(Tree)  
#SRW: as an alternative to using "view" all the time in your code, just click
#on "Tree" in the environment window to see it when you want to check it out
#Using view is fine too, it just breaks up the flow when you are running through
#the code because it shifts you to the other tab

#SRW: also change orchard number to a factor in the tree level data

Tree$orchard.num <- as.factor(as.character(Tree$orchard.num))

Tree <- Tree %>%
  mutate(avgwgt=(apple.wgt.3-bag.weight)/3) 
  
Tree <- Tree %>% 
  dplyr::select(-c(5:6))

#SRW: Here, you should name the summarized version something different so you
#can also keep the tree level data in your environment. Changed "Tree" to "Tree_Sum"

Tree_Sum <- Tree %>%
  group_by(orchard.num) %>%
  summarise_at(c("avgwgt", "Firmness", "SSC", "maturity.index"), mean, na.rm = TRUE)

#SRW: changed tree to Tree_Sum
c <- left_join(Tree_Sum, Orchard, by="orchard.num")

view(c)

#Question 1 
#How do management systems (organic or conventional) shape fruit quality?-------
#analysis: 1 linear model per trait 

#SRW: Why are you using the orchard level data here instead of the tree level?
#I would suggest the tree level data for sure
#Also, the current models had orchard number in them as a random effect but only 
#one measurement per orchard. In general, you only need something in the model as
#a random effect if you have multiple non-independent measurements from that 
#group 
#Also orchard number was being treated as a numeric variable. 


m1<- glmmTMB(SSC ~ orchard.type + (1|site.code/orchard.num), data=c)
summary(m1)
Anova(m1) #p=0.294

m2<- glmmTMB(Firmness ~ orchard.type + (1|site.code/orchard.num), data=c)
summary(m2)
Anova(m2) #p=0.9851

m3<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/orchard.num), data=c)
summary(m3)
Anova(m3) #p=0.8712

m4<- glmmTMB(maturity.index ~ orchard.type + (1|site.code/orchard.num), data=c)
summary(m4)
Anova(m4) #p=0.6817


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

#How do management systems interact with broad climatic changes across latitude?----
#analysis: 1 linear model per proxy trait (lat, long, elev, prox to water) per fruit trait 

##SRW: Several things with this analysis:
# 1) It is repetitive from the analysis above. Both analyses
#tell you the effect of orchard type. I would suggest dropping the simple version
#above. That is because, considering the very strong effects of latitude, it is very
#likely that looking at the effects of orchard type without accounting for that 
#could be misleading. The analysis that includes both orchard type and 
#latitude is a more fair assessment of the effects of orchard type because 
#it looks at how orchard type affects the traits while also accounting for the 
#effects of latitude

# 2) Considering that there will be way too much here for one paper, and we need
#to simplify, and there are no clear effects of longitude, elevation,
#or proximity to water, I would suggest dropping all of those variables. 

# 3) As above, these models were not correct because they are done on the 
#orchard level data but also include orchard number as a random effect. You would
#only need orchard number in the model if you have multiple measures per 
#orchard. Also, orchard number would need to be a factor

#See my further suggestions below after this section for how to remedy all this


###SSC###
s1 <- glmmTMB(SSC ~ Latitude*orchard.type + (1|site.code/orchard.num), data=c)
summary(s1)
Anova(s1)#Latitude p=1.906e-05 

s2<- glmmTMB(SSC ~ elevation*orchard.type + (1|site.code/orchard.num), data=c)
summary(s2)
Anova(s2)

s3<- glmmTMB(SSC ~ Longitude*orchard.type + (1|site.code/orchard.num), data=c)
summary(s3)
Anova(s3)

s4<- glmmTMB(SSC ~ Prox.Water*orchard.type + (1|site.code/orchard.num), data=c)
summary(s4)
Anova(s4)

###Firmness###
f1 <- glmmTMB(Firmness ~ Latitude*orchard.type + (1|site.code/orchard.num), data=c)
summary(f1)
Anova(f1)#Latitude p=4.112e-09 ***


f2<- glmmTMB(Firmness ~ elevation*orchard.type + (1|site.code/orchard.num), data=c)
summary(f2)
Anova(f2)

f3<- glmmTMB(Firmness ~ Longitude*orchard.type + (1|site.code/orchard.num), data=c)
summary(f3)
Anova(f3)

f4<- glmmTMB(Firmness ~ Prox.Water*orchard.type + (1|site.code), data=c)
summary(f4)
Anova(f4)


###Average Weight
w1 <- glmmTMB(avgwgt ~ Latitude*orchard.type + (1|site.code/orchard.num), data=c)
summary(w1)
Anova(w1)#Latitude:orchard.type p=0.001822 **

w2<- glmmTMB(avgwgt ~ elevation*orchard.type + (1|site.code/orchard.num), data=c)
summary(w2)
Anova(w2) #elevation p=0.02534 *

w3<- glmmTMB(avgwgt ~ Longitude*orchard.type + (1|site.code/orchard.num), data=c)
summary(w3)
Anova(w3)

w4<- glmmTMB(avgwgt ~ Prox.Water*orchard.type + (1|site.code), data=c)
summary(w4)
Anova(w4)


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

mt4<- glmmTMB(maturity.index ~ Prox.Water*orchard.type + (1|site.code), data=c)
summary(mt4)
Anova(mt4)



#SRW: Considering all my comments above, my suggestion for Q1 is to look at 
#management and latitude together, using the tree level data. 

#this requires adding latitude to the tree level dataset, even though it is an
#"orchard-level" variable. I don't 100% love that (it would be better if we had
#precise latitude values for each individual tree), but I think it is the best
#solution and okay to do
#because overall the tree to tree variation in latitude would be virtually
#nothing compared to the overall variation in latitude across the entire dataset
#I don't think that is necessarily the case for other orchard level variables
#like temperature, where actually you can get a lot of microclimatic variation

#Like this:

#----start SRW additions------

#first add latitude to tree-level dataframe
Tree <- left_join(Tree, Orchard[,c(2,4)], by="orchard.num")

m1<- glmmTMB(SSC ~ orchard.type*Latitude + (1|site.code/orchard.num), data=Tree)
summary(m1)
Anova(m1) 

#be sure to check some diagnostics
hist(resid(m1))  #looks great
diagnose(m1)

m2<- glmmTMB(Firmness ~ orchard.type*Latitude + (1|site.code/orchard.num), data=Tree)
summary(m2)
Anova(m2)

#be sure to check some diagnostics
hist(resid(m2))  #looks great
diagnose(m2)

m3<- glmmTMB(avgwgt ~ orchard.type*Latitude + (1|site.code/orchard.num), data=Tree)
summary(m3)
Anova(m3) 

#be sure to check some diagnostics
hist(resid(m3))  #looks great
diagnose(m3)

m4<- glmmTMB(maturity.index ~ orchard.type*Latitude + (1|site.code/orchard.num), data=Tree)
summary(m4)
Anova(m4) 

#be sure to check some diagnostics
hist(resid(m4))  #looks okay
diagnose(m4)  #large coefficients for random effects, could be problematic
#could also try this with a poisson model or negative binomial

m4b<- glmmTMB(maturity.index ~ orchard.type*Latitude + (1|site.code/orchard.num), 
             data=Tree, family="compois")
summary(m4b)
Anova(m4b) 

hist(resid(m4b))  #looks okay
diagnose(m4b)  #similar error

anova(m4, m4b)

#similar conclusions or errors either way, but the poisson model is a better fit
#based on AIC scores, so I would suggest reporting results from that 


#So, basically you have strong effects of latitude on all the quality variables
#and no clear overall effects of orchard type

#however, we do have a strong interaction between latitude and orchard type for 
#weight

#to further explore the interaction between orchard type and latitude, first make
#some plots!
ggplot(Tree, aes(x=Latitude, y=avgwgt, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_point()

#actually latitude effects do not really look linear...
ggplot(Tree, aes(x=Latitude, y=avgwgt, color=orchard.type))+
  geom_smooth()+
  geom_point()

#looks an interaction where organic apples are bigger at low latitude and
#conventional apples are bigger at high latitudes

#Typically when you have an interaction, like this, you want to split the data
#by one of the factors and look at effects of other factor indepedently for
#each group. simplest would be to split by orchard type, but then the focus of
#the analysis is on the different effects of latitude for each type, rather than
#differences between type
#alternatively, you can bin the data into different groups based on latitude, 
#and look at effects of type at different groups of orchard. I think in this 
#case there is a big gap between 40 to 45, so we could just split into high and
#low latitudes  

Tree_low <- filter(Tree, Latitude<42)
Tree_high <- filter(Tree, Latitude>42)

m3_low<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/orchard.num), data=Tree_low)
summary(m3_low)
Anova(m3_low) 

m3_high<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/orchard.num), data=Tree_high)
summary(m3_high)
Anova(m3_high)

#conventional apples marginally larger at high latitudes

#another interesting thing, the latitude effects seem to be almost entirely driven
#by the large differences between high and low latitudes, because if you put latitude
#into the models above, you acutally don't see any effect of latitude within the
#high and low groups





#----end SRW additions--------




#Which abiotic factors are the most important drivers of fruit quality?---------
#Analysis: principle components analysis followed by PC regression with each variable
names(c)
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

##Linear models PC x quality 

pc_clim <- as.data.frame(results$x)
pc_clim <- cbind(pc_clim, c)

hist(c$Firmness)
hist(c$SSC)
hist(c$avgwgt)
hist(c$maturity.index)  #SRW: does not approximate normal, would suggest maybe
    #a poisson model as above for this response variable

###Firmness###

p1 <- glmmTMB(Firmness ~ PC1 + PC2 + PC3 + PC4 + PC5+ (1|site.code), data=pc_clim)
summary(p1)
#PC1  2.26e-10 ***
#PC3  1.18e-06 ***
#PC5  0.0267 *  

#SRW: be sure to check some diagnostics!! 
hist(resid(p1))
diagnose(p1)

#strong positive effects of PC1 + PC3 and strong negative effect of PC5
plot(Firmness ~ PC1, data=pc_clim)
plot(Firmness ~ PC3, data=pc_clim)
plot(Firmness ~ PC5, data=pc_clim)

#So, to interpret this you go back to look at the loading on each PC axis.
results$rotation



###SSC###
P2 <- glmmTMB(SSC ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data= pc_clim)
summary(P2)
#PC1 2e-16 ***
#PC2 0.00407 ** 
#PC3 0.00172 ** 
#PC5 0.00243 ** 

plot(SSC ~ PC1, data=pc_clim)  #WOW that is amazing how much variation is explained by PC1
plot(SSC ~ PC2, data=pc_clim)
plot(SSC ~ PC3, data=pc_clim)
plot(SSC ~ PC5, data=pc_clim)


results$rotation


###AVGWGT
p3 <- glmmTMB(avgwgt ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc_clim)
summary(p3)
#PC2 0.0502 .  
#PC3 0.0463 *  
#PC4 0.0378 *  

plot(avgwgt ~ PC2, data=pc_clim)
plot(avgwgt ~ PC3, data=pc_clim)
plot(avgwgt ~ PC4, data=pc_clim)

results$rotation


###Maturity Index 
p4 <- glmmTMB(maturity.index ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc_clim)
summary(p4)
#PC1 0.000889 ***
#PC4 0.025583 *  
 

plot(maturity.index ~ PC1, data=pc_clim)
plot(maturity.index ~ PC4, data=pc_clim)


results$rotation


###correlation matrix
library(Hmisc)
#The first matrix shows the correlation coefficients between the variables  
#the second matrix shows the corresponding p-values.
rcorr(as.matrix(p_clim))
#any value below 0.05 is not a statistically significant relationship 
# need to sort through and organzie these values 

#now let's visualize this 
library(corrplot)
corrplot(cor(p_clim))
#red = negative, blue = positive 
#big = high, small = low

#SRW: There are lots of ways to add R values or p-values to these plots
#to help organize and present the info, see:
?corrplot
testRes = cor.mtest(p_clim, conf.level = 0.95)

## add all p-values
corrplot(cor(p_clim), p.mat = testRes$p, insig = 'p-value', sig.level = -1)

## add significant level stars
corrplot(cor(p_clim), p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')

#SRW: I think some version of this would be good to include in a supplement, but
#you'd have to play around with it to make it pretty

#Also try...
library(GGally)
ggpairs(p_clim) 



#Which specific management practices are the most important drivers of fruit quality?----
#Analysis: principle components analysis followed by PC regression with each variable

##SRW: This analysis does not work because basically all your management practices
#are categorical variables. You can only feed continuous numeric variables into a 
#PCA. I'm sorry I did not fully think that through before

##I would instead suggest something like...

m <- glmmTMB(SSC ~ Cultivation + Herbicides + Com_Mul + Mowing +
              Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
              data=c)
summary(m)
Anova(m)

#Whoa, crazy effects of herbicides and fire mgmt...

plot(SSC ~ as.factor(Herbicides), data=c)
plot(SSC ~ as.factor(Fire_Mgmt), data=c)

#Okay, so fire management doesn't make sense to include because you only have one
#orchard that does it. 

#Herbicides are weird because they are showing up as such a big effect in the stats
#but it is not clear from the plot. This needs a lot more thought. Also I am 
#really confused about the herbicide variable because it does not correspond
#with organic vs conventional. What herbicides would organic orchards be using??

c %>% count(orchard.type, Herbicides)


#SRW: everything below from here to L682 should be taken out
#since you can't do a PC regression with categorical variables

names(c)
p_mgmt <- dplyr::select(c,c("Cultivation","Herbicides","Com_Mul", "Mowing",
"Weed_Mats", "Cover_Crops","Fire_Mgmt","Acres"))

#calculate principal components
results1 <- prcomp(p_mgmt, scale = TRUE)

#this gives you the % variance explained
summary(results1)  

#display principal components
results1$x

results1$rotation
#display the first six scores
head(results1$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results1,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))


#calculate total variance explained by each principal component
summary(results1)$importance
summary(results1)$importance[2,]

var_explained1 = results1$sdev^2 / sum(results1$sdev^2)

df1 <- data.frame(PC=1:8, var_explained1=var_explained1)

#create scree plot
ggplot(df1, aes(x=PC, y=var_explained1)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


##Linear models PC x quality 

pc_mgmt <- as.data.frame(results1$x)
pc_mgmt <- cbind(pc_mgmt, c)

###Firmness###
p5 <- glmmTMB(Firmness ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 
              + PC8 +(1|site.code), data=pc_mgmt)
summary(p5)
#nothing

###SSC###
P6 <- glmmTMB(SSC ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 
              + PC8 + (1|site.code), data= pc_mgmt)
summary(P6)
#nothing

###AVGWGT
p7 <- glmmTMB(avgwgt ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 
              + PC8 + (1|site.code), data=pc_mgmt)
summary(p7)
#PC3                   -8.477      3.704  -2.289   0.0221 *  
#PC4                   -7.329      4.423  -1.657   0.0975 .  
#PC5                   13.522      5.644   2.396   0.0166 *  
#PC6                  -13.325      5.891  -2.262   0.0237 *  

plot(avgwgt ~ PC3, data=pc_mgmt)
plot(avgwgt ~ PC4, data=pc_mgmt)
plot(avgwgt ~ PC5, data=pc_mgmt)
plot(avgwgt ~ PC6, data=pc_mgmt)

results1$rotation


###Maturity Index 
p8 <- glmmTMB(maturity.index ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 
              + PC8 + (1|site.code), data=pc_mgmt)
summary(p8)
#orchard.typeOrganic  7.54e-06 ***
#PC1 0.0347 *
#PC6 0.0525 .    
#PC7 0.0237 *    
#PC8 0.0671 . 


plot(maturity.index ~ PC1, data=pc_mgmt)
plot(maturity.index ~ PC6, data=pc_mgmt)
plot(maturity.index ~ PC7, data=pc_mgmt)
plot(maturity.index ~ PC8, data=pc_mgmt)

results1$rotation


###correlation matrix 
library(Hmisc)
#The first matrix shows the correlation coefficients between the variables  
#the second matrix shows the corresponding p-values.
rcorr(as.matrix(p_mgmt))
#any value below 0.05 is not a statistically significant relationship 
# need to sort through and organzie these values 

#now let's visualize this 
library(corrplot)
corrplot(cor(p_mgmt))
#red = negative, blue = positive 
#big = high, small = low


#Which pest or diseases presence has the most significant affect on fruit quality-----


p_pest <- dplyr::select(c,c("Anthracnose","European Canker","Bullseye Rot", "Powdery mildew",     
"Apple scab","Root Rot","Fire Blight","Apple Maggots","Codling Moth","Aphids","Tree Borer",
"Cedar Apple Rust","Bitter Rot","Leaf Roller","Trhips",
"Scale"))

p_pest[is.na(p_pest)] <- 1
view(p_pest)

#calculate principal components
results2 <- prcomp(p_pest, scale = TRUE)

#this gives you the % variance explained
summary(results2)  

#display principal components
results2$x

results2$rotation
#display the first six scores
head(results2$x)

#this plots the results of the PCAs into a two dimensional representation 
biplot(results2,
       col = c('darkblue', 'red'),
       scale = TRUE, xlabs = rep("*", 24))

#SRW: oh wow, pretty interesting, it's like fireblight is orthogonal 
#to everything else...


#calculate total variance explained by each principal component
summary(results2)$importance
summary(results2)$importance[2,]

var_explained2 = results2$sdev^2 / sum(results2$sdev^2)

df2 <- data.frame(PC=1:16, var_explained2=var_explained2)

#create scree plot
ggplot(df2, aes(x=PC, y=var_explained2)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

#we'll use 1 through 8 
##Linear models PC x quality 

p_pest <- as.data.frame(results2$x)
p_pest <- cbind(p_pest, c)

###Firmness###
f1 <- glmmTMB(Firmness ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 
              + PC8 +(1|site.code), data=p_pest)
summary(f1)
#orchard.typeOrganic   2.0141     0.4691   4.293 1.76e-05 ***
#PC5                   1.6039     0.4989   3.215   0.0013 ** 
#PC6                  -2.8775     0.5627  -5.114 3.16e-07 ***
#PC8                  -2.4218     0.8705  -2.782   0.0054 ** 


plot(Firmness ~ PC5, data=p_pest)
plot(Firmness ~ PC6, data=p_pest)
plot(Firmness ~ PC8, data=p_pest)

results2$rotation

#SRW: A lot to unpack here!! There are a lot of very strong effects!
#Also, a very weird thing is that all of a sudden when
#you account for all these pest variables the orchard type has a huge effect!!
#I would not have thought initially to include that along with the pest variables
#but that is weird

#One question is whether it really makes sense to include all of these, and whether
#the PC regression is the best approach. Part of that depends on how much variation
#you have in these scales, and how much they are correlated with one another
library(GGally)
ggpairs(c, columns=25:42) 

d <- pivot_longer(data=c, cols=25:42, names_to="pest", values_to="index")

ggplot(d, aes(x=pest, y=index))+
  geom_boxplot()+
  geom_jitter(width=0.2, height=0.1)


#SRW: based on these plots you could consider removing all the ones that are just 
#1s and a few 2s or 3s. Basically, anything where anything above 1 is an outlier
#in the boxplots. That is most of them actually, it would only leave you with
#8 pest/diseases. Then, try the correlation plot again. If they are not strongly 
#correlated, you could just include each individual pest/disease value
#as a predictor, similar to what I suggested for the management data, rather
#than doing the PCA. The advantage of this is that it is easier to interpret than
#a PCA

#SRW: I would also suggest trying to create separate overall indices for 
#insect pressure and disease pressure and trying to look at those. With so much
#data you can't get around needing to look at things a lot of different ways!!







#pest index alone 
f2 <- glmmTMB(Firmness ~ Pest_Index*orchard.type +(1|site.code), data=c)
summary(f2)
Anova(f2)
#nothing 

###SSC###
ssc1 <- glmmTMB(SSC ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 
              + PC8 +(1|site.code), data=p_pest)
summary(ssc1)
#PC5                  0.7325145  0.1795326   4.080  4.5e-05 ***
#PC6                  0.5875009  0.3214876   1.827 0.067633 .  
#PC7                  1.0489073  0.3187258   3.291 0.000999 ***
#PC8                 -0.9180447  0.4629784  -1.983 0.047377 *  


plot(SSC ~ PC5, data=p_pest)
plot(SSC ~ PC6, data=p_pest)
plot(SSC ~ PC7, data=p_pest)
plot(SSC ~ PC8, data=p_pest)

results2$rotation

#pest index alone 
ssc2 <- glmmTMB(SSC ~ Pest_Index*orchard.type +(1|site.code), data=c)
summary(ssc2)
Anova(ssc2)
#Pest_Index:orchard.type 3.0044  1    0.08304 .

###avgwgt###
wgt1 <- glmmTMB(avgwgt ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 
                + PC8 +(1|site.code), data=p_pest)
summary(wgt1)
#PC1                   3.5095     1.8541   1.893 0.058380 . 
#PC3                   7.3827     2.8087   2.629 0.008575 **   
#PC5                  -9.8890     2.6309  -3.759 0.000171 ***
#PC6                 -12.5362     2.7650  -4.534 5.79e-06 ***
#PC8                   9.8685     3.3920   2.909 0.003622 ** 


plot(avgwgt ~ PC1, data=p_pest)
plot(avgwgt ~ PC3, data=p_pest)
plot(avgwgt ~ PC5, data=p_pest)
plot(avgwgt ~ PC6, data=p_pest)
plot(avgwgt ~ PC8, data=p_pest)

results2$rotation

#pest index alone 
wgt2 <- glmmTMB(avgwgt ~ Pest_Index*orchard.type +(1|site.code), data=c)
summary(wgt2)
Anova(wgt2)
#nothing


###maturity.index###
mi1 <- glmmTMB(maturity.index ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 
                + PC8 +(1|site.code), data=p_pest)
summary(mi1)
#PC1                  0.152563   0.075699   2.015  0.04386 *  
#PC2                  0.232106   0.086275   2.690  0.00714 ** 
#PC4                  0.137820   0.064195   2.147  0.03180 * 
#PC8                  0.334125   0.141356   2.364  0.01809 *  


plot(maturity.index ~ PC1, data=p_pest)
plot(maturity.index ~ PC2, data=p_pest)
plot(maturity.index ~ PC4, data=p_pest)
plot(maturity.index ~ PC8, data=p_pest)


results2$rotation

#pest index alone 
mi2 <- glmmTMB(maturity.index ~ Pest_Index*orchard.type +(1|site.code), data=c)
summary(mi2)
Anova(mi2)
#Pest_Index              10.9447  1  0.0009387 ***

#SRW: interesting!! Potentially high pest pressure leads to faster maturity?
#this would make a lot of sense based on phytohormonal responses 
#jasmonic acid controls fruit maturation as well as induced response to insects
#See Whitehead and Poveda 2011
#Or, could just be that this is co-varying with latitude as well...
m <- glmmTMB(Pest_Index ~ Latitude +(1|site.code), data=c)
summary(m)
Anova(m)  
#Does not seem to be that! Worth playing around with this a bit more and looking
#at some plots to really understand it but that is pretty interesting





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

#Figures---------------------------------------------------------------
#How do management systems (organic or conventional) shape fruit quality?
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

#How do management systems interact with broad climatic changes across latitude?

lat1 = ggplot(c, aes(x=Latitude, y=SSC, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("SSC") +
  xlab ("Latitude")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
lat1


lat2 = ggplot(c, aes(x=Latitude, y=Firmness, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("Firmness") +
  xlab ("Latitude)")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
lat2


lat3 = ggplot(c, aes(x=Latitude, y=avgwgt, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("Average Weight (g)") +
  xlab ("Latitude")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
lat3


lat4 = ggplot(c, aes(x=Latitude, y=maturity.index, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("Maturity Index") +
  xlab ("Latitude")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
lat4

multiplot(lat1, lat2, lat3, lat4, cols=2)

elv1 = ggplot(c, aes(x=elevation, y=avgwgt, color=orchard.type)) +
  theme_classic() +
  geom_point() +
  ylab ("Avergae Weight") +
  xlab ("Elevation")+
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
elv1



#Which abiotic factors are the most important drivers of fruit quality?
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




#PCA MGMT traits 
ggplot(c, aes(x = Herbicides, y = maturity.index)) +
  geom_point(aes(color = orchard.type)) +
  geom_smooth(method=glm, se=FALSE)+
  facet_wrap(~orchard.type)+
  scale_color_viridis_d()




