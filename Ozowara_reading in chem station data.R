#Analyzing chromatograms using R 
rm(list=ls()) # clears work space
setwd("C:\\Users/xozow/OneDrive/Desktop/HPLC_DATA")
files <- list.files(all.files=TRUE, pattern="*.CSV", recursive=TRUE,include.dirs=TRUE)

#Create data frames for each UV signal and add RTs


s1 <- read.delim("C:\\Users/xozow/OneDrive/Desktop/HPLC_DATA/RL_20230222 2023-02-22 11-09-27/01-T1-SKIN.D/REPORT01.CSV",
             fileEncoding="utf-16", sep=",", header=FALSE)
PAs_280 <- data.frame(RT=s1[2])
names(PAs_280)[1] <- "RT"

s2 <- read.delim("C:\\Users/xozow/OneDrive/Desktop/HPLC_DATA/RL_20230222 2023-02-22 11-09-27/01-T1-SKIN.D/REPORT02.CSV",
             fileEncoding="utf-16", sep=",", header=FALSE)
PAs_525 <- data.frame(RT=s2[2])
names(PAs_525)[1] <- "RT"

s3 <- read.delim("C:\\Users/xozow/OneDrive/Desktop/HPLC_DATA/RL_20230222 2023-02-22 11-09-27/01-T1-SKIN.D/REPORT03.CSV",
             fileEncoding="utf-16", sep=",", header=FALSE)
PAs_320 <- data.frame(RT=s3[2])
names(PAs_320)[1] <- "RT"

s4 <- read.delim("C:\\Users/xozow/OneDrive/Desktop/HPLC_DATA/RL_20230222 2023-02-22 11-09-27/01-T1-SKIN.D/REPORT04.CSV",
             fileEncoding="utf-16", sep=",", header=FALSE)
PAs_365 <- data.frame(RT=s4[2])
names(PAs_365)[1] <- "RT"



####FILE COMPILING FUNCTION-----------------------------------------------------
for (i in 1:length(32)) {
  
  #get files from Signal 1, 280 nm
  if ( grepl("REPORT01.CSV",files[i]) == TRUE) {
    #read data           
    a <- read.delim(files[i],fileEncoding="utf-16", sep=",", 
                    header=FALSE, na.strings="-")          
    b <- a[5]
    #extract sample name from file name 
    names (b) <- gsub(".D/.*$", "",(gsub( "RL.*/M","XO.*/M", "M", files[i] )))
    #add data to dataframe
    PAs_280 <- cbind (PAs_280, b)
  }
  
  
  #get files from Signal 2, 525 nm
  if ( grepl("REPORT02.CSV",files[i]) == TRUE) {
    #read data           
    a <- read.delim(files[i],fileEncoding="utf-16", sep=",", 
                    header=FALSE, na.strings="-")          
    b <- a[5]
    #extract sample name from file name 
    names (b) <- gsub(".D/.*$", "",(gsub( "RL.*/M","XO.*/M", "M", files[i] )))
    #add data to dataframe
    PAs_525 <- cbind (PAs_525, b)

  }
  
  
  #get files from Signal 3, 320 nm
  if ( grepl("REPORT03.CSV",files[i]) == TRUE) {
    #read data           
    a <- read.delim(files[i],fileEncoding="utf-16", sep=",", 
                    header=FALSE, na.strings="-")          
    b <- a[5]
    #extract sample name from file name 
    names (b) <- gsub(".D/.*$", "",(gsub( "RL.*/M","XO.*/M", "M", files[i] )))
    #add data to dataframe
    PAs_320 <- cbind (PAs_320, b)

  }
  
  #get files from Signal 4, 365 nm
  if ( grepl("REPORT04.CSV",files[i]) == TRUE) {
    #read data           
    a <- read.delim(files[i],fileEncoding="utf-16", sep=",", 
                    header=FALSE, na.strings="-")          
    b <- a[5]
    #extract sample name from file name 
    names (b) <- gsub(".D/.*$", "",(gsub( "RL.*/M","XO.*/M", "M", files[i] )))
    #add data to dataframe
    PAs_365 <- cbind (PAs_365, b)
    
  }
  
}

PAs_280
PAs_525
PAs_320
PAs_365

#Combine all data into single data frame

Signal <- c(rep ("280nm", nrow(PAs_280)), rep ("525nm", nrow(PAs_525)), 
            rep ("320nm", nrow(PAs_320)), rep ("365nm", nrow(PAs_365)))
PAs_Md <- cbind(Signal, rbind(PAs_280, PAs_525, PAs_320, PAs_365))

#create another variable with just the sample numbers (to match the way they are listed
#in the bioassay datasheet)
ID2 <- paste("Md", gsub("[^0-9]","",names(PAs_Md)), sep="")
#used gsub to replace everything but numbers with blanks, then paste Md to start

#create another variable with just the plant part (skin, pulp, seeds)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
part <- substrRight (names(PAs_Md), 2)

PAs_Md <- rbind(ID2, part, PAs_Md)
write.csv (PAs_Md, file="PAs_Md.csv")

View(PAs_Md)
