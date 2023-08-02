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
library("purrr")
library(Hmisc) #correlation matrix
library(corrplot)#plotting correlation matrix
library(readr)
library(GGally)
#Read in data-------------------------------------------------------------------
Tree <- read_csv("Tree Level Data.csv")
Orchard <- read_csv("Orchard Level Data.csv")

Orchard$orchard.num <- as.factor(as.character(Orchard$orchard.num))
Tree$orchard.num <- as.factor(as.character(Tree$orchard.num))

#combine data sheets
View(Tree)
Tree <- Tree %>%
  mutate(avgwgt=(apple.wgt.3-bag.weight)/3) 
  
Tree <- Tree %>% 
  dplyr::select(-c(5:6))

Tree.Sum <- Tree %>%
  group_by(orchard.num) %>%
  summarise_at(c("avgwgt", "Firmness", "SSC", "maturity.index"), mean, na.rm = TRUE)

c <- left_join(Tree.Sum, Orchard, by="orchard.num")

view(c)

#aggregate the data-------------------------------------------------------------
#so that we have some idea of the differences between orchard types
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
#Analysis: Using GLMMs with tree level data to explore how measured qualities interact
#across latitude 

#creating a new data set to work with 
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
#using poisson distribution 
Lat4<- glmmTMB(maturity.index ~ orchard.type*Latitude + (1|site.code/orchard.num), 
              data=TreeLat, family="compois")
summary(Lat4)
Anova(Lat4)  
#orchard.type          0.0051  1    0.94320  
#Latitude              4.5525  1    0.03287 *
#orchard.type:Latitude 1.2588  1    0.26188 

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


#Which abiotic factors are the most important drivers of fruit quality?---------
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
#PC3  1.18e-06 ***
#PC5  0.0267 *  

#strong positive effects of PC1 + PC2 and strong negative effect of PC3
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

plot(SSC ~ PC1, data=pc_clim)
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


#Which specific management practices are the most important drivers of fruit quality?----
#GLMMs with all management variables and physical traitss 

#SSC
mgmt1 <- glmmTMB(SSC ~ Cultivation + Herbicides + Com_Mul + Mowing +
               Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
             data=c)
summary(mgmt1)
Anova(mgmt1)
#Herbicides  47.9815  1  4.303e-12 ***

ssc_herb <- glmmTMB(SSC ~  orchard.type*Herbicides+ (1|site.code), data=c)
summary(ssc_herb)
Anova(ssc_herb)
#orchard.type            4.0006  1   0.045485 * 
#Herbicides              7.2867  1   0.006947 **


ggplot(c, aes(x=Herbicides, y=SSC, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()

#avgwgt
mgmt2 <- glmmTMB(avgwgt ~ Cultivation + Herbicides + Com_Mul + Mowing +
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
mgmt3 <- glmmTMB(Firmness ~ Cultivation + Herbicides + Com_Mul + Mowing +
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
mgmt4 <- glmmTMB(maturity.index ~ Cultivation + Herbicides + Com_Mul + Mowing +
                   Weed_Mats + Cover_Crops + Fire_Mgmt + Acres + (1|site.code), 
                 data=c)
summary(mgmt4)
Anova(mgmt4)
#Herbicides  4.7331  1    0.02959 *
#Mowing      3.9845  1    0.04592 *
#Weed_Mats   4.9615  1    0.02592 *
#Acres       4.4712  1    0.03447 *

mat_herb <- glmmTMB(maturity.index ~  orchard.type*Herbicides+ (1|site.code), data=c)
summary(mat_herb)
Anova(mat_herb)
#orchard.type:Herbicides  6.6521  1   0.009904 ** 
#organic orchards that used herbicides had higher maturity values 
ggplot(c, aes(x=Herbicides, y=maturity.index))+
  geom_smooth(method = "lm") +
  geom_boxplot()

ggplot(c, aes(x=Herbicides, y=maturity.index, color=orchard.type))+
  geom_smooth(method = "lm") +
  geom_boxplot()




#Which pest or diseases presence has the most significant affect on fruit quality-----

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
p_pest <- dplyr::select(c,c("Powdery mildew",     
"Apple scab","Root Rot","Fire Blight","Apple Maggots","Codling Moth","Aphids","Bitter Rot"))

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

#GLMMs models 
applemaggots= c$`Apple Maggots`
codmoth = c$`Codling Moth`
powmil= c$`Powdery mildew`
bitterrot= c$`Bitter Rot`
applescab=c$`Apple scab`
rootrot= c$`Root Rot`
fireblight= c$`Fire Blight`


###SSC###
ssc_pest <- glmmTMB(SSC ~ Aphids+applemaggots+codmoth+
powmil+bitterrot+applescab+rootrot+fireblight+ (1|site.code), 
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

#Figures
firessc <- ggplot(c, aes(x=fireblight, y=SSC, color=orchard.type)) +
  geom_point() +
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
firessc

rootwgt <- ggplot(c, aes(x=rootrot, y=avgwgt, color=orchard.type)) +
  geom_point() +
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
rootwgt

matpest <- ggplot(c, aes(x=Pest_Index, y=maturity.index, color=orchard.type)) +
  geom_point() +
  geom_smooth(method=glm, se=FALSE)+
  scale_color_manual(values=c("#3EBCD2", "#9A607F"),name="Management System")
matpest



multiplot(firessc, rootwgt,matpest, cols=2)







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
