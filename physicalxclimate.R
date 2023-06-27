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
View(Tree)
names(Tree)
library(readr)
Orchard <- read_csv("Orchard Level Data.csv")
View(Orchard)
names(Orchard)

Tree <- Tree %>%
  group_by(orchard.num) %>%
  summarise_at(c("apple.wgt.3", "Firmness", "SSC", "maturity.index", "bag.weight"), mean, na.rm = TRUE)

c <- left_join(Tree, Orchard, by="orchard.num")

c <-c[-25,]

c <- c %>%
  mutate(avgwgt=(c$apple.wgt.3-c$bag.weight)/3)

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


###Quality###

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


#Question 1---------------------------------------------------------------------
#How do management systems interact with broad antibiotic conditions to shape fruit quality
#we'll divide this into two parts
#1: management systems alone
#2: proxy climatic factors (lat, long, elev)

###part 1### 
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

###part 2 using Proxy Climatic Variables
#values not correct!
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


#Maturity  
mt1 <- glmmTMB(maturity.index ~ Latitude*orchard.type + (1|site.code/orchard.num), data=c)
summary(mt1)
Anova(mt1)#Latitude 0.05779 .

mt2<- glmmTMB(maturity.index ~ elevation*orchard.type + (1|site.code/orchard.num), data=c)
summary(mt2)
Anova(mt2)

mt3<- glmmTMB(maturity.index ~ Longitude*orchard.type + (1|site.code/orchard.num), data=c)
summary(mt3)
Anova(mt3)


###traits influence on each other 
sxf<- glmmTMB( SSC+Firmness~ orchard.type + (1|site.code/orchard.num), data=c)
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


#Question 1 Figures---------------------------------------------------------------
###mgmt and physical quality alone 
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
pest1 <- glmmTMB(SSC ~ Pest_Index*orchard.type + (1|site.code/onum/Tree), data=c)
summary(pest1)
Anova(pest1)
#Pest_Index:orchard.type 2.9152  1    0.08775 .

pest2 <- glmmTMB(avgwgt ~ Pest_Index*orchard.type + (1|site.code/onum/Tree), data=c)
summary(pest2)
Anova(pest2)
# none

pest3 <- glmmTMB(Firmness ~ Pest_Index*orchard.type + (1|site.code/onum/Tree), data=c)
summary(pest3)
Anova(pest3)
#none 

pest4 <- glmmTMB(maturity.index ~ Pest_Index*orchard.type + (1|site.code/onum/Tree), data=c)
summary(pest4)
Anova(pest4)
#none 

##Herbicides## 
herb1 <- glmmTMB(SSC ~ Herbicides*orchard.type + (1|site.code/onum/Tree), data=c)
summary(herb1)
Anova(herb1)
#Herbicides              5.4980  1    0.01904 *
#orchard.type            2.7336  1    0.09826 .
#Herbicides:orchard.type 1.2220  1    0.26896  

herb2 <- glmmTMB(avgwgt ~ Herbicides*orchard.type + (1|site.code/onum/Tree), data=c)
summary(herb2)
Anova(herb2)
#none

herb3 <- glmmTMB(Firmness ~ Herbicides*orchard.type + (1|site.code/onum/Tree), data=c)
summary(herb3)
Anova(herb3)
#none 

herb4 <- glmmTMB(maturity.index ~ Herbicides*orchard.type + (1|site.code/onum/Tree), data=c)
summary(herb4)
Anova(herb4)
#none 


##Acres## 
Acres <- glmmTMB(SSC ~ Acres*orchard.type + (1|site.code/onum/Tree), data=c)
summary(Acres)
Anova(Acres)
#none

Acres <- glmmTMB(avgwgt ~ Acres*orchard.type + (1|site.code/onum/Tree), data=c)
summary(Acres)
Anova(Acres)
#Acres              10.8040  1   0.001013 **

Acres <- glmmTMB(Firmness ~ Acres*orchard.type + (1|site.code/onum/Tree), data=c)
summary(Acres)
Anova(Acres)
#none 

Acres <- glmmTMB(maturity.index ~ Acres*orchard.type + (1|site.code/onum/Tree), data=c)
summary(Acres)
Anova(Acres)
#none 

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
