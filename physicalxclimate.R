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

#combine data sheets
View(Tree)
Tree <- Tree %>%
  mutate(avgwgt=(apple.wgt.3-bag.weight)/3) 
  
Tree <- Tree %>% 
  dplyr::select(-c(5:6))

Tree <- Tree %>%
  group_by(orchard.num) %>%
  summarise_at(c("avgwgt", "Firmness", "SSC", "maturity.index"), mean, na.rm = TRUE)

c <- left_join(Tree, Orchard, by="orchard.num")

view(c)


#Question 1 
#How do management systems (organic or conventional) shape fruit quality?-------
#analysis: 1 linear model per trait 
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
hist(c$maturity.index)

###Firmness###

p1 <- glmmTMB(Firmness ~ PC1 + PC2 + PC3 + PC4 + PC5+ (1|site.code), data=pc_clim)
summary(p1)
#PC1  2.26e-10 ***
#PC2  0.7698    
#PC3  1.18e-06 ***
#PC4  0.2002    
#PC5  0.0267 *  

#strong positive effects of PC1 + PC2 and strong negative effect of PC3
plot(Firmness ~ PC1, data=pc_clim)
plot(Firmness ~ PC2, data=pc_clim)
plot(Firmness ~ PC3, data=pc_clim)

#So, to interpret this you go back to look at the loading on each PC axis.
results$rotation



###SSC###
P2 <- glmmTMB(SSC ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data= pc_clim)
summary(P2)
#PC1 2e-16 ***
#PC2 0.00407 ** 
#PC3 0.00172 ** 
#PC4 0.49452    
#PC5 0.00243 ** 

plot(SSC ~ PC1, data=pc_clim)
plot(SSC ~ PC2, data=pc_clim)
plot(SSC ~ PC3, data=pc_clim)
plot(SSC ~ PC5, data=pc_clim)


results$rotation


###AVGWGT
p3 <- glmmTMB(avgwgt ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc_clim)
summary(p3)
#PC1 0.2898    
#PC2 0.0502 .  
#PC3 0.0463 *  
#PC4 0.0378 *  
#PC5 0.6157


plot(avgwgt ~ PC2, data=pc_clim)
plot(avgwgt ~ PC3, data=pc_clim)
plot(avgwgt ~ PC4, data=pc_clim)

results$rotation


###Maturity Index 
p4 <- glmmTMB(maturity.index ~ orchard.type + PC1 + PC2 + PC3 + PC4 + PC5 + (1|site.code), data=pc_clim)
summary(p4)
#PC1 0.000889 ***
#PC2 0.649219    
#PC3 0.268255    
#PC4 0.025583 *  
#PC5 0.350911 

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


#Which specific management practices are the most important drivers of fruit quality?----
#Analysis: principle components analysis followed by PC regression with each variable
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

names(c)

p_pest <- dplyr::select(c,c("Anthracnose","European Canker","Bullseye Rot", "Powdery mildew",     
"Apple scab","Root Rot","Fire Blight","Apple Maggots","Codling Moth","Aphids","Tree Borer",
"Cedar Apple Rust","Bitter Rot","Leaf Roller","Horned Caterpillars","Trhips","Stemble",
"Scale", "Pest_Index"))


#calculate principal components
results2 <- prcomp(p_pest, scale = TRUE)

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
pc_mgmt <- cbind(pc_mgmt, c



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

#Q1 Figures---------------------------------------------------------------
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

#Q2 Figures ------------------------------------------------------------

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

ggplot(c, aes(x = Szn.UVI, y = SSC)) +
  geom_point(aes(color = orchard.type)) +
  geom_smooth(method=glm, se=FALSE)+
  facet_wrap(~orchard.type)+
  scale_color_viridis_d()

#Q3 Figures-------------------------------------------------------------



