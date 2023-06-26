rm(list=ls()) # clears work space

#notes from SRW-------------------------------------------------------------------------- 

#STEP 1: re-organize the data into 3 separate spreadsheets for orchard-level,
#tree-level, and fruit-level data. It sounds like your import issues were
#because you had blank rows in the datasheets? But each csv file should only 
#include the variables that were recorded at that level. So, the chemistry data 
#would only be in the "fruit level" datasheet. The "orchard level" datasheet
#should have just one row per orchard. In the "fruit-level" datasheet you do 
#still want to have columns for orchard ID and tree ID, even though those would
#be repeated, because they are categorical variables and still the "true value"
#for each fruit. I hope that makes sense

#STEP 2: double check all the data, I found cases where you have multiple sets 
#of data on fruit quality for some trees (see below), or multiple sets of 
#management data for the same orchard. Doing Step 1 should help with this

#STEP 3: re-do PC regressions with only the predictor variables feeding
#into the PCAs




###install packages-------------------------------------------------------------
library(ggplot2)
library(MASS)
library(reshape2)
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
chap_all_dat <- read_csv("chap_all_dat.csv")
View(chap_all_dat)

#SRW: you seem to have a stray 0 in column 52 in the csv file, so it is adding 
#those to the datatable...should go and delete this in the csv file. For now
#I just changed the columns to 52 below

#were going to alter this table for fruit quality analysis 
library(readr)
c <- read_csv("chap_all_dat.csv")
View(c)
#filter out variables that we're not going to look at 
c <- c %>% dplyr::select(-c(26:52))
view(c)
#remove duplicated rows from the data 
c <- c %>% distinct()
view(c)


##SRW: This should now be a "tree-level" dataset, but you have a lot of trees that
#are still duplicated in multiple rows. This is because the rest of the data
#is not exactly the same. I found cases where you have multiple sets of fruit 
#quality data for the same tree (e.g. tree 3_4)
#and also cases where you have multiple sets of management data for the same
#orchard (e.g. orchard 23)

c$Tree[which(duplicated(c$Tree))]  #list of the ones that are duplicated

View(filter(c, Tree=="23_5"))  #to see the data for particular tree


#SRW: Another note, once you have the data organized into different datasets
#for different levels, you will need to summarize the data at the more detailed
#level to analyze it alongside the data at the more coarse level. So, for example,
#you would take the mean value of fruit firmness for each orchard

phys <- read_csv("chap_all_dat.csv")
phys <- phys[, -c(51:52)]

phys_sum <- phys %>%
  group_by(Orchard.num) %>%
  summarise_at(c("avgwgt", "Firmness", "SSC", "maturity.index"), mean, na.rm = TRUE)

#Then you can combine this new summarized data with your climate data by 
#cross-referencing the orchard number, e.g.  if you had a datatable "clim" that 
#just had the climate data at the orchard level, you would join it with the 
#phys_sum table like this:

left_join(phys_sum, clim, by=Orchard.num)



#PCA----------------------------------------------------------------------------
#principle components analysis 
#first we're just going to look at the climatic variables and their affect

#create practice data set 
p <- read_csv("chap_all_dat.csv")
View(p)

#filter out variables that we're not going to look at 

p <- p %>% dplyr::select(-c(1:13, 26:52 ))
view(p)
#remove duplicated rows from the data 
p <- p %>% distinct()
view(p)

#SRW: here you have a mix of predictor variables and response variables in the 
#PCA. The idea with a PC regression is that you only feed into the PCA
#the predictor variables, and usually only the ones of a certain class that you 
#expect to be strongly correlated with one another. I left all the true climate
#variables along with the "proxies"(lat/long/elev/prox.water), but you could also
#try this with ONLY the true climate variables
#ALSO--the climate variables are recorded at the orchard level, so you should
#be running this stuff with a new dataset that has only one row per orchard!!

#suggest dividing your table into the predictors and responses
p_qual <- select(p, c("avgwgt", "Firmness", "SSC", "maturity.index"))
p_clim <- select(p, 1:8)  
    #SRW: just a note, it is always better coding practice to select by column
    #names instead of numbers because if later your dataset changes slightly and
    #the column indices are different it can introduce a lot of errors into your analysis
    #I did that for the quality variables as an example but was too lazy to
    #type them all out for the climate


#calculate principal components
results <- prcomp(p_clim, scale = TRUE)

#SRW: also see
summary(results)  #note this already gives you the % variance explained
results$x

#reverse the signs   SRW:  not sure why you are doing this???
#results$rotation <- -1*results$rotation

#display principal components
results$rotation #here we see that PC7 has the highest values for all components 

#reverse the signs of the scores  #SRW: not sure why you were doing this???
#this also saves the scores by orchard number 
#results$x <- -1*results$x

#display the first six scores
head(results$x)
#PC2 is the largest (longitude)


#this plots the results of the PCAs into a two dimensional representation 
biplot(results, scale = 0)




#calculate total variance explained by each principal component
results$sdev^2 / sum(results$sdev^2)
# [1] 4.993391e-01 1.590931e-01 1.245701e-01 7.628088e-02 5.318525e-02 3.992590e-02 2.067177e-02 1.141188e-02
#[9] 7.588457e-03 4.096733e-03 2.383892e-03 1.452838e-03 7.578747e-08

#SRW: can also use
summary(results)$importance
summary(results)$importance[2,]


#We can also create a scree plot – a plot that displays the total variance
#explained by each principal component – to visualize the results of PCA:
#calculate total variance explained by each principal component
var_explained = results$sdev^2 / sum(results$sdev^2)


#SRW: changed to ggplot, qplot is deprecated in my updated version of R

df <- data.frame(PC=1:8, var_explained=var_explained)

#create scree plot
ggplot(df, aes(x=PC, y=var_explained)) + 
  geom_line() + 
  geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)



#second we're just going to look at the mgmt variables and their affect

#SRW:  Note, I did not update this section of code, but you could edit this 
#according to the same suggestions I made above for the climate variables.
#Create a separate "predictor" and "response" dataset, re-run

#create practice data set 
p2 <- read_csv("chap_all_dat.csv")
View(p2)

#filter out variables that we're not going to look at 
p2 <- p2 %>% dplyr::select(-c(1:4, 14:21, 26:50 ))
#remove duplicated rows from the data 
p2 <- p2 %>% distinct()
view(p2)

#calculate principal components
results <- prcomp(p2, scale = TRUE)

#reverse the signs
results$rotation <- -1*results$rotation

#display principal components
results$rotation #here we see that PC9 has the highest values for all components 

#reverse the signs of the scores
#this also saves the scores by orchard number 
results$x <- -1*results$x

#display the first six scores
head(results$x)
#PC3 is the largest (herbicides)


#this plots the results of the PCAs into a two dimensional representation 
biplot(results, scale = 0)


#calculate total variance explained by each principal component
results$sdev^2 / sum(results$sdev^2)
# [1] 4.993391e-01 1.590931e-01 1.245701e-01 7.628088e-02 5.318525e-02 3.992590e-02 2.067177e-02 1.141188e-02
#[9] 7.588457e-03 4.096733e-03 2.383892e-03 1.452838e-03 7.578747e-08


#We can also create a scree plot – a plot that displays the total variance
#explained by each principal component – to visualize the results of PCA:
#calculate total variance explained by each principal component
var_explained = results$sdev^2 / sum(results$sdev^2)

#create scree plot
qplot(c(1:14), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


#PCR----------------------------------------------------------------------------


#SRW: I have not used this package before, but I think this model is essentially
#doing an analysis in two steps: 1) run a PCA on the predictor variables (essentially)
#what you did above, and 2) running a linear model with the PC scores as predictors
#and SSC as a response. So, it would essentially be the same as doing
#




#start with a principle components regression with each response variable (ssc, firm, maturity, wgt)
#once we find the most significant interactions we can move on for further analysis
library(pls)

###SSC MODEL###
#make this example reproducible
set.seed(1)
names(p)

#fit PCR model
#the response variable will be the physical trait measured 
#predictor variables will be a combo of our fixed and proxy
model <- pcr(SSC~Latitude+Szn.Temp.Avg+Prox.Water+elevation+Longitutde+
               Szn.Total.Precip+Szn.Min.Avg+Szn.Max.Avg,
             data=p, scale=TRUE, validation="CV")

#intercept is 2.898    
#drops after one component so we should use 2?
#variance is 88.74     with three components?

validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")

#the model fits with two principle components 


##SRW: also see...

model$coefficients
model$scores
model$loadings


###AVGWGT MODEL###
#make this example reproducible
set.seed(1)

#fit PCR model
#the response variable will be the physical trait measured 
#predictor variables will be a combo of our fixed and proxy
model1 <- pcr(avgwgt~Latitude+Szn.Temp.Avg+Prox.Water+elevation+Longitutde+
                Szn.Total.Precip+Szn.Min.Avg+Szn.Max.Avg,
              data=p, scale=TRUE, validation="CV")

summary(model1)
#intercept is 32.49        
#does not drop off until 4 comps 
#variance is 88.74 with three components?

validationplot(model1)
validationplot(model1, val.type="MSEP")
validationplot(model1, val.type="R2")
#the model fits with two principle components



#Correlation Matrix-------------------------------------------------------------
#we're going to round out looking at the influence of our varibales by creating a,
#correlation matrix. This will help to pin point which interactions to look at, 
#during analysis 

#first create a practice data set
#were going to alter this table for fruit quality analysis 
p1 <- read_csv("chap_all_dat.csv")
View(p1)
#filter out variables that we're not going to look at 
p1 <- p1 %>% dplyr::select(-c(1:4, 26:50 ))
#remove duplicated rows from the data 
p1 <- p1 %>% distinct()
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
#values currently not correct!
m1<- glmmTMB(SSC ~ orchard.type + (1|site.code/Orchard.num/Tree), data=c)
summary(m1)
Anova(m1) 

m2<- glmmTMB(Firmness ~ orchard.type + (1|site.code/Orchard.num/Tree), data=c)
summary(m2)
Anova(m2)

m3<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/Orchard.num/Tree), data=c)
summary(m3)
Anova(m3)#orchard.type p=0.8997

m4<- glmmTMB(maturity.index ~ orchard.type + (1|site.code/Orchard.num/Tree), data=c)
summary(m4)
Anova(m4)

#aggregate the data so that we have some idea of the differences between orchard types
aggregate(SSC~orchard.type, data=c, FUN=mean)
# Conventional 10.88182
# Organic 11.78462

aggregate(Firmness~orchard.type, data=c, FUN=mean)
#Conventional 22.26473
#Organic 22.44308

aggregate(avgwgt~orchard.type, data=c, FUN=mean)
#1 Conventional 100.93273
#    Organic 98.68308

aggregate(maturity.index~orchard.type, data=c, FUN=mean)
# Conventional       3.781818
#    Organic       4.138462

###part 2###
#values not correct!
###SSC###
s1 <- glmmTMB(SSC ~ Latitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(s1)
Anova(s1)
#Latitude p=1.906e-05 ***
#orchard.type p=0.4788    
#Latitude:orchard.type p=0.6393  

s2<- glmmTMB(SSC ~ elevation*orchard.type + (1|site.code/onum/Tree), data=c)
summary(s2)
Anova(s2)
#elevation p=0.1148
#orchard.type p=0.2849
#elevation:orchard.type p=0.6719

s3<- glmmTMB(SSC ~ Longitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(s3)
Anova(s3)
#Longitutde              0.5870  1     0.4436
#orchard.type            0.9519  1     0.3292
#Longitutde:orchard.type 0.0287  1     0.8656


###Firmness###
f1 <- glmmTMB(Firmness ~ Latitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(f1)
Anova(f1)
#Latitude p=4.112e-09 ***
#orchard.type p=0.7735    
#Latitude:orchard.type p=0.3366 

f2<- glmmTMB(Firmness ~ elevation*orchard.type + (1|site.code/onum/Tree), data=c)
summary(f2)
Anova(f2)
#elevation p=0.1808
#orchard.type p=0.9686
#elevation:orchard.type p=0.6368

f3<- glmmTMB(Firmness ~ Longitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(f3)
Anova(f3)

#Longitutde              0.2882  1     0.5914
#orchard.type            0.0030  1     0.9563
#Longitutde:orchard.type 1.8458  1     0.1743


###Average Weight (weight of 3 apples minus the bag divided by 3)###
w1 <- glmmTMB(avgwgt ~ Latitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(w1)
Anova(w1)
#Latitude p=0.570994   
#orchard.type p=0.716308   
#Latitude:orchard.type p=0.001643 ** 

w2<- glmmTMB(avgwgt ~ elevation*orchard.type + (1|site.code/onum/Tree), data=c)
summary(w2)
Anova(w2)
#elevation p=0.005707 **
#orchard.type p=0.465530   
#elevation:orchard.type p=0.087329 .

w3<- glmmTMB(avgwgt ~ Longitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(w3)
Anova(w3)
#Longitude p=0.3660
#orchard.type p=0.7945
#Longitude:orchard.type p=0.1693


#Maturity  
mt1 <- glmmTMB(maturity.index ~ Latitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(mt1)
Anova(mt1)
#Latitude p=0.0489 *
#orchard.type p=0.9454  
#Latitude:orchard.type p=0.3272 

mt2<- glmmTMB(maturity.index ~ elevation*orchard.type + (1|site.code/onum/Tree), data=c)
summary(mt2)
Anova(mt2)
#elevation p=0.5632
#orchard.type p=0.8615
#elevation:orchard.type p=0.4549

mt3<- lm(maturity.index ~ Longitude*orchard.type, data=c)
summary(mt3)
Anova(mt3)
#Longitude                6.807   1  4.5859 0.03433 *
#orchard.type             2.193   1  1.4777 0.22661  
#Longitude:orchard.type   0.157   1  0.1055 0.74589 

###traits influence on each other 
sxf<- glmmTMB( SSC*Firmness~ orchard.type + (1|site.code/onum), data=c)
summary(sxf)
Anova(sxf)

sxw<- glmmTMB(SSC*avgwgt~orchard.type + (1|site.code/onum), data=c)
summary(sxw)
Anova(sxw)

sxm<- glmmTMB(SSC*maturity.index~orchard.type + (1|site.code/onum), data=c)
summary(sxm)
Anova(sxm)

wxf<- glmmTMB(avgwgt*Firmness ~orchard.type+ (1|site.code/onum), data=c)
summary(wxf)
Anova(wxf)


wxm<- glmmTMB(avgwgt*maturity.index ~orchard.type+ (1|site.code/onum), data=c)
summary(wxf)
Anova(wxf)


fxm<- glmmTMB(maturity.index*Firmness ~orchard.type+ (1|site.code/onum), data=c)
summary(fxm)
Anova(fxm)


#Question 1 Figures---------------------------------------------------------------
###mgmt and physcial quality alone 
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
#total precipitation x otype 
pl1<- glmmTMB(SSC ~ Szn.Total.Precip*orchard.type + (1|site.code/onum), data=c)
summary(pl1)
Anova(pl1)
#Szn.Total.Precip              7.3440  1   0.006729 **
#orchard.type                  1.5530  1   0.212698   
#Szn.Total.Precip:orchard.type 1.6851  1   0.194244  

pl2<- glmmTMB(Firmness ~ Szn.Total.Precip*orchard.type + (1|site.code/onum), data=c)
summary(pl2)
Anova(pl2)
#Szn.Total.Precip              0.1953  1     0.6586
#orchard.type                  0.0057  1     0.9397
#Szn.Total.Precip:orchard.type 0.6159  1     0.4326

pl3<- glmmTMB(avgwgt ~ Szn.Total.Precip*orchard.type + (1|site.code/onum), data=c)
summary(pl3)
Anova(pl3)
#Szn.Total.Precip              2.8708  1    0.09020 .
#orchard.type                  0.0067  1    0.93480  
#Szn.Total.Precip:orchard.type 3.9157  1    0.04784 *

pl4<- glmmTMB(maturity.index ~ Szn.Total.Precip*orchard.type + (1|site.code/onum), data=c)
summary(pl4)
Anova(pl4)
#Szn.Total.Precip              7.7214  1   0.005457 **
#orchard.type                  0.0916  1   0.762174   
#Szn.Total.Precip:orchard.type 0.2209  1   0.638342  

#avgerage seasonal temp x otype   
atl<- glmmTMB(SSC ~ Szn.Temp.Avg*orchard.type + (1|site.code/onum), data=c)
summary(atl)
Anova(atl)
#Szn.Temp.Avg              26.5282  1  2.597e-07 ***
#orchard.type               0.1087  1     0.7416    
#Szn.Temp.Avg:orchard.type  0.5756  1     0.4481  

at2<- glmmTMB(Firmness ~ Szn.Temp.Avg*orchard.type + (1|site.code/onum), data=c)
summary(at2)
Anova(at2)
#Szn.Temp.Avg              11.3947  1  0.0007365 ***
#orchard.type               0.0010  1  0.9748785    
#Szn.Temp.Avg:orchard.type  0.2247  1  0.6354858  

at3<- glmmTMB(avgwgt ~ Szn.Temp.Avg*orchard.type + (1|site.code/onum), data=c)
summary(at3)
Anova(at3)
#Szn.Temp.Avg              0.3241  1   0.569146   
#orchard.type              0.1482  1   0.700282   
#Szn.Temp.Avg:orchard.type 9.2367  1   0.002372 **

at4<- glmmTMB(maturity.index ~ Szn.Temp.Avg*orchard.type + (1|site.code/onum), data=c)
summary(at4)
Anova(at4)
#Szn.Temp.Avg              6.9050  1   0.008595 **
#orchard.type              0.0423  1   0.837056   
#Szn.Temp.Avg:orchard.type 0.6519  1   0.419438 

#proximity to water x otype 
pxl<- glmmTMB(SSC ~ Prox.Water*orchard.type + (1|site.code/onum), data=c)
summary(pxl)
Anova(pxl)
#Prox.Water              1.4996  1     0.2207
#orchard.type            0.9875  1     0.3204
#Prox.Water:orchard.type 0.4179  1     0.5180

px2<- glmmTMB(Firmness ~ Prox.Water*orchard.type + (1|site.code/onum), data=c)
summary(px2)
Anova(px2)
#Prox.Water              1.1828  1     0.2768
#orchard.type            0.0378  1     0.8457
#Prox.Water:orchard.type 1.8788  1     0.1705

px3<- glmmTMB(avgwgt ~ Prox.Water*orchard.type + (1|site.code/onum), data=c)
summary(px3)
Anova(px3)
#Prox.Water              0.1974  1     0.6569
#orchard.type            0.0106  1     0.9179
#Prox.Water:orchard.type 0.0027  1     0.9585

px4<- glmmTMB(maturity.index ~ Prox.Water*orchard.type + (1|site.code/onum), data=c)
summary(px4)
Anova(px4)
#NaN

#season avg max x otype 
smx1<- glmmTMB(SSC ~ Szn.Max.Avg*orchard.type + (1|site.code/onum), data=c)
summary(smx1)
Anova(smx1)
#Szn.Max.Avg              18.9030  1  1.375e-05 ***
#orchard.type              0.1568  1     0.6921    
#Szn.Max.Avg:orchard.type  0.0080  1     0.9287 

smx2<- glmmTMB(Firmness ~ orchard.type*Szn.Max.Avg + (1|site.code/onum), data=c)
summary(smx2)
Anova(smx2)
#orchard.type              0.0353  1   0.851057   
#Szn.Max.Avg              10.8113  1   0.001009 **
#orchard.type:Szn.Max.Avg  0.2284  1   0.632709   

smx3<- glmmTMB(avgwgt ~ orchard.type*Szn.Max.Avg + (1|site.code/onum), data=c)
summary(smx3)
Anova(smx3)
#orchard.type              0.1452  1  0.7031488    
#Szn.Max.Avg               0.0018  1  0.9663953    
#orchard.type:Szn.Max.Avg 11.4873  1  0.0007007 ***

smx4<- glmmTMB(maturity.index ~ orchard.type*Szn.Max.Avg + (1|site.code/onum), data=c)
summary(smx4)
Anova(smx4)
#orchard.type             0.3166  1   0.573643   
#Szn.Max.Avg              7.7421  1   0.005395 **
#orchard.type:Szn.Max.Avg 0.1453  1   0.703109 

#season avg min x otype 
sm1<- glmmTMB(SSC ~ orchard.type*Szn.Min.Avg + (1|site.code/onum), data=c)
summary(sm1)
Anova(sm1)
#orchard.type              1.0205  1     0.3124    
#Szn.Min.Avg              15.9763  1  6.414e-05 ***
#orchard.type:Szn.Min.Avg  0.4255  1     0.5142

sm2<- glmmTMB(Firmness ~ orchard.type*Szn.Min.Avg + (1|site.code/onum), data=c)
summary(sm2)
Anova(sm2)
#orchard.type             0.0786  1   0.779246   
#Szn.Min.Avg              9.3198  1   0.002267 **
#orchard.type:Szn.Min.Avg 0.0299  1   0.862609  

sm3<- glmmTMB(avgwgt ~ orchard.type*Szn.Min.Avg + (1|site.code/onum), data=c)
summary(sm3)
Anova(sm3)
#orchard.type             0.0461  1    0.82994  
#Szn.Min.Avg              0.8381  1    0.35995  
#orchard.type:Szn.Min.Avg 5.3243  1    0.02103 *

sm4<- glmmTMB(maturity.index ~ orchard.type*Szn.Min.Avg + (1|site.code/onum), data=c)
summary(sm4)
Anova(sm4)
#orchard.type             0.0811  1     0.7757
#Szn.Min.Avg              2.5495  1     0.1103
#orchard.type:Szn.Min.Avg 1.0085  1     0.3153

#Question 2 Figures ------------------------------------------------------------

#Question 3---------------------------------------------------------------------
#Which specific management practices are the most important drivers of fruit quality?
names(c)

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
