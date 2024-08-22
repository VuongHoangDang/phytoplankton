#THIS IS TITAN ANALYSIS
#To run TITAN, species with less than 4 occurences needs to be removed in input data "HN_PHY_WQ_DATA_log_1.csv'. 
#Data occurence should be transformed into a log10 scale (with cell number + 1) 
#Note: TITAN running is very time consumming, it is better to run TITAn for each environmental seperately

library(TITAN2)
library(dplyr)

#Read dataset
data <- read.csv("HN_PHY_WQ_DATA_log_1.csv")
head(data)
summary(data)
colnames(data)
#Taxa community: 15 - 191
#Create subset community data

taxa <- data[,c(15:191)]
colnames(taxa)

#
#
#
# TN TITAN ANALYSIS (1)

#Creat subset of Environmental varibales: TN (mg/L)
TN = data$TN
summary(TN)

#Creat new data frame
df <-  data.frame(TN, taxa)
df <-  na.omit(df)

TN_taxa  = df
TN = TN_taxa$TN
species <- colnames(TN_taxa)
species
taxa = TN_taxa[,c(2:178)] 

#After remove NA values of environmental data, occurence of spp might be changed, therefore, need to check the occurence again
#Identify spp occured less than 4 times in new dataframe 

for (i in colnames(taxa)){
  occurences <- sum(taxa[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}

# print spp
# [1] "T.Anabaena.bergii"
# [1] "T.Dolichospermum.cf.circinale"

TN_data.TITAN = TN_taxa %>% 
  select(-T.Anabaena.bergii, -T.Dolichospermum.cf.circinale )
#check occurrences again

colnames(TN_data.TITAN)
taxa.TITAN = TN_data.TITAN[,c(2:176)]

#double check spp occured less than 4 times
for (i in colnames(taxa.TITAN)){
  occurences <- sum(taxa.TITAN[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
#Good, no low occurrence species
#Now let's run TITAN

TN.TITAN <- TN_data.TITAN$TN
taxa.TITAN
#TITAN Calculation
glades.titan_TN <- titan(TN.TITAN, taxa.TITAN)
#great, let's wait for the results 5 mins/bootstrap => 500 bootstraps will be done in 2 days

save.image("TITAN_TN_All_spp.RData")
savehistory("TITAN_TN_All_spp.Rhistory")





#
#
#
# TP TITAN ANALYSIS (2)
data <- read.csv("HN_PHY_WQ_DATA_log_1.csv")
#Creat subset of Environmental varibales: TP (mg/L)
TP = data$TP
summary(TP)
taxa <- data[,c(15:191)]

#Creat new data frame
df <-  data.frame(TP, taxa)
df <-  na.omit(df)

TP_taxa = df
TP = TP_taxa$TP
species <- colnames(TP_taxa)
taxa = TP_taxa[,c(2:178)] 

#After remove NA values of environmental data, occurence of spp might be changed, therefore, need to check the occurence again
#Identify spp occured less than 4 times in new dataframe 

for (i in colnames(taxa)){
  occurences <- sum(taxa[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}

# [1] "T.Anabaena.bergii"
# [1] "T.Dolichospermum.cf.circinale"

TP_data.TITAN <- TP_taxa %>% 
  select(-T.Anabaena.bergii, -T.Dolichospermum.cf.circinale )
#check occurrences again

colnames(TP_data.TITAN)
taxa.TITAN = TP_data.TITAN[,c(2:176)]

for (i in colnames(taxa.TITAN)){
  occurences <- sum(taxa.TITAN[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
##Good, no low occurrence species
#Now let's run TITAN

TP.TITAN <- TP_data.TITAN$TP
taxa.TITAN
#TITAN Calculation
glades.titan_TP <- titan(TP.TITAN, taxa.TITAN)

save.image("TITAN_TP_All_spp.RData")
savehistory("TITAN_TP_All_spp.Rhistory")



#
#
#
# TNTP TITAN ANALYSIS (3)
data <- read.csv("HN_PHY_WQ_DATA_log_1.csv")
#Creat subset of Environmental varibales: TNTP
TN_TP_Ratio = data$TN_TP_Ratio
summary(TN_TP_Ratio)

hist(TN_TP_Ratio)
taxa <- data[,c(15:191)]
#Creat new data frame
df <-  data.frame(TN_TP_Ratio, taxa)
df <-  na.omit(df)

TN_TP_Ratio_taxa = df
TN_TP_Ratio = TN_TP_Ratio_taxa$TN_TP_Ratio
species <- colnames(TN_TP_Ratio_taxa)
species
taxa = TN_TP_Ratio_taxa[,c(2:176)] 

#After remove NA values of environmental data, occurence of spp might be changed, therefore, need to check the occurence again
#Identify spp occured less than 4 times in new dataframe 
for (i in colnames(taxa)){
  occurences <- sum(taxa[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}

TN_TP_Ratio_data.TITAN <- TN_TP_Ratio_taxa %>% 
  select(-T.Anabaena.bergii, -T.Dolichospermum.cf.circinale)

#check occurrences again

colnames(TN_TP_Ratio_data.TITAN)
taxa.TITAN = TN_TP_Ratio_data.TITAN[,c(2:176)]

for (i in colnames(taxa.TITAN)){
  occurences <- sum(taxa.TITAN[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
##Good, no low occurrence species
#Now let's run TITAN

TN_TP_Ratio.TITAN <- TN_TP_Ratio_data.TITAN$TN_TP_Ratio
taxa.TITAN
#TITAN Calculation
glades.titan_TN_TP_Ratio <- titan(TN_TP_Ratio.TITAN, taxa.TITAN)

save.image("TITAN_TN_TP_Ratio_All_spp.RData")
savehistory("TITAN_TN_TP_Ratio_All_spp.Rhistory")


#
#
#
# TEMPERATURE TITAN ANALYSIS (4)
data <- read.csv("HN_PHY_WQ_DATA_log_1.csv")
#Creat subset of Environmental varibales: TEMP
TEMP = data$TEMP
summary(TEMP)

hist(TEMP)
taxa <-  data[,c(15:191)]

df <-  data.frame(TEMP, taxa)
df <-  na.omit(df)

TEMP_taxa = df
TEMP = TEMP_taxa$TEMP
taxa = TEMP_taxa[,c(2:178)] 

#After remove NA values of environmental data, occurence of spp might be changed, therefore, need to check the occurence again
#Identify spp occured less than 4 times in new dataframe 

for (i in colnames(taxa)){
  occurences <- sum(taxa[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}

TEMP_data.TITAN <- TEMP_taxa %>% 
  select(-T.Anabaena.bergii, -T.Dolichospermum.cf.circinale)

#check occurrences again

colnames(TEMP_data.TITAN)
taxa.TITAN = TEMP_data.TITAN[,c(2:176)]

for (i in colnames(taxa.TITAN)){
  occurences <- sum(taxa.TITAN[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
##Good, no low occurrence species
#Now let's run TITAN

TEMP.TITAN <- TEMP_data.TITAN$TEMP
taxa.TITAN
#TITAN Calculation
glades.titan_TEMP <- titan(TEMP.TITAN, taxa.TITAN)

save.image("TITAN_TEMP_All_spp.RData")
savehistory("TITAN_TEMP_All_spp.Rhistory")


#
#
#
# AMMONIUM TITAN ANALYSIS (5)

data <- read.csv("HN_PHY_WQ_DATA_log_1.csv")
AMM = data$AMM
summary(AMM)
taxa <-  data[,c(15:191)]
df <-  data.frame(AMM, taxa)
df <-  na.omit(df)

AMM_taxa = df
AMM = AMM_taxa$AMM
taxa = AMM_taxa[,c(2:178)] 

for (i in colnames(taxa)){
  occurences <- sum(taxa[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
# [1] "T.Anabaena.aphanizomenoides"
# [1] 2
# [1] "T.Anabaena.bergii"
# [1] 2
# [1] "T.Dolichospermum.cf.circinale"

AMM_data.TITAN <- AMM_taxa %>% 
  select(-T.Anabaena.aphanizomenoides, -T.Anabaena.bergii, -T.Dolichospermum.cf.circinale )
#check occurrences again

colnames(AMM_data.TITAN)
taxa.TITAN = AMM_data.TITAN[,c(2:175)]

for (i in colnames(taxa.TITAN)){
  occurences <- sum(taxa.TITAN[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
##Good, no low occurrence species
#Now let's run TITAN

AMM.TITAN <- AMM_data.TITAN$AMM
taxa.TITAN
#TITAN Calculation
glades.titan_AMM <- titan(AMM.TITAN, taxa.TITAN)
#great, let's wait for the results 5 mins/bootstrap => 500 bootstraps will be done in 2 days
save.image("TITAN_AMM_All_spp.RData")
savehistory("TITAN_AMM_All_spp.Rhistory")


#
#
#
# NIT TITAN ANALYSIS (6)

data <- read.csv("HN_PHY_WQ_DATA_log_1.csv")
NIT = data$NIT
summary(NIT)
taxa <-  data[,c(15:191)]
df <-  data.frame(NIT, taxa)
df <-  na.omit(df)

NIT_taxa = df
NIT = NIT_taxa$NIT
colnames(NIT_taxa)
taxa = NIT_taxa[,c(2:178)] 


for (i in colnames(taxa)){
  occurences <- sum(taxa[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
# [1] "T.Anabaena.aphanizomenoides"
# [1] 2
# [1] "T.Anabaena.bergii"
# [1] 2
# [1] "T.Dolichospermum.cf.circinale"

NIT_data.TITAN <- NIT_taxa %>% 
  select(-T.Anabaena.aphanizomenoides, -T.Anabaena.bergii, -T.Dolichospermum.cf.circinale )
#check occurrences again

colnames(NIT_data.TITAN)
taxa.TITAN = NIT_data.TITAN[,c(2:175)]

for (i in colnames(taxa.TITAN)){
  occurences <- sum(taxa.TITAN[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
##Good, no low occurrence species
#Now let's run TITAN

NIT.TITAN <- NIT_data.TITAN$NIT

#TITAN Calculation
glades.titan_NIT <- titan(NIT.TITAN, taxa.TITAN)
save.image("TITAN_NIT_All_spp.RData")
savehistory("TITAN_NIT_All_spp.Rhistory")

#
#
#
# FRP TITAN ANALYSIS (7)

data <- read.csv("HN_PHY_WQ_DATA_log_1.csv")
#Create subset community data

taxa <- data[,c(15:191)]
colnames(taxa)
#Creat subset of Environmental varibales: FRP (mg/L)
FRP = data$FRP
summary(FRP)
#Creat new data frame
df <-  data.frame(FRP, taxa)
df <-  na.omit(df)
FRP_taxa = df
FRP = FRP_taxa$FRP
species <- colnames(FRP_taxa)
taxa = FRP_taxa[,c(2:178)] 


for (i in colnames(taxa)){
  occurences <- sum(taxa[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}

# [1] "T.Anabaena.aphanizomenoides"
# [1] 2
# [1] "T.Anabaena.bergii"
# [1] 0
# [1] "T.Aphanizomenonaceae"
# [1] 0
# [1] "T.Dolichospermum.cf.circinale"
# [1] 0
# [1] "T.Dolichospermum.circinale"
# [1] 3
# [1] "T.Radiocystis"
# [1] 1
# [1] "Anagnostidinema"
# [1] 0
# [1] "Haptophyte"
# [1] 0
# [1] "Monoraphidium.cf"
# [1] 3
# [1] "Peridinoid"
# [1] 3
# [1] "Prorocentrum.non.toxic"
# [1] 3
# [1] "Dimorphococcus"
# [1] 2
# [1] "T.Gymnodinoid"

FRP_data.TITAN <- FRP_taxa %>% 
  select(-T.Anabaena.aphanizomenoides, -T.Anabaena.bergii, -T.Aphanizomenonaceae, - T.Dolichospermum.cf.circinale,
         -T.Dolichospermum.circinale, -T.Radiocystis, -Anagnostidinema, -Anagnostidinema, -Haptophyte, -Monoraphidium.cf, 
         -Peridinoid,-Dimorphococcus, -T.Gymnodinoid  , -Prorocentrum.non.toxic)
#check occurrences again

colnames(FRP_data.TITAN)
taxa.TITAN = FRP_data.TITAN[,c(2:165)]

for (i in colnames(taxa.TITAN)){
  occurences <- sum(taxa.TITAN[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
##Good, no low occurrence species
#Now let's run TITAN

FRP.TITAN <- FRP_data.TITAN$FRP
taxa.TITAN
#TITAN Calculation
glades.titan_FRP <- titan(FRP.TITAN, taxa.TITAN)
#great, let's wait for the results 5 mins/bootstrap => 500 bootstraps will be done in 2 days
save.image("TITAN_FRP_All_spp.RData")
savehistory("TITAN_FRP_All_spp.Rhistory")


#
#
#
# SALINITY TITAN ANALYSIS (8)

data <- read.csv("HN_PHY_WQ_DATA_log_1.csv")

taxa <- data[,c(15:191)]
SAL = data$SAL #unit : psu
df <-  data.frame(SAL, taxa)
df <-  na.omit(df)

SAL_taxa = df

SAL = SAL_taxa$SAL

species <- colnames(SAL_taxa)

taxa = SAL_taxa[,c(2:178)]

for (i in colnames(taxa)){
  
  occurences <- sum(taxa[[i]]>0)
  
  if (occurences < 4){
    
    print(occurences)
    
    print(i)
    
  }
  
}



# [1] "T.Anabaena.aphanizomenoides"
# [1] 1
# [1] "T.Anabaena.bergii"
# [1] 0
# [1] "T.Aphanizomenonaceae"
# [1] 1
# [1] "T.Dolichospermum.cf.circinale"
# [1] 1
# [1] "Anagnostidinema"
# [1] 0
# [1] "Monoraphidium.cf"
# [1] 3
# [1] "Dimorphococcus"



SAL_data.TITAN <- SAL_taxa %>%
  
  select(-T.Anabaena.aphanizomenoides, -T.Anabaena.bergii, - T.Aphanizomenonaceae, 
         -T.Dolichospermum.cf.circinale, -Anagnostidinema, -Monoraphidium.cf, -Dimorphococcus)

#check occurrences again



colnames(SAL_data.TITAN)

taxa.SAL.TITAN = SAL_data.TITAN[,c(2:171)]



for (i in colnames(taxa.SAL.TITAN)){
  
  occurences <- sum(taxa.SAL.TITAN[[i]]>0)
  
  if (occurences < 4){
    
    print(occurences)
    
    print(i)
    
  }
  
}

##Good, no low occurrence species

#Now let's run TITAN



SAL.TITAN <- SAL_data.TITAN$SAL

taxa.SAL.TITAN

glades.titan_SAL <- titan(SAL.TITAN, taxa.SAL.TITAN)
save.image("TITAN_SAL_All_spp.RData")
savehistory("TITAN_SAL_All_spp.Rhistory")



#
#
#
# SILICA TITAN ANALYSIS (9)

data <- read.csv("HN_PHY_WQ_DATA_log_1.csv")

SIL = data$SIL
summary(SIL)
taxa <-  data[,c(15:191)]
df <-  data.frame(SIL, taxa)
df <-  na.omit(df)

SIL_taxa = df
SIL = SIL_taxa$SIL
species <- colnames(SIL_taxa)
taxa = SIL_taxa[,c(2:178)] 


for (i in colnames(taxa)){
  occurences <- sum(taxa[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}

# [1] 2
# [1] "T.Anabaena.bergii"
# [1] 0
# [1] "T.Aphanizomenonaceae"
# [1] 0
# [1] "T.Dolichospermum.cf.circinale"
# [1] 0
# [1] "T.Dolichospermum.circinale"
# [1] 2
# [1] "T.Radiocystis"
# [1] 0
# [1] "Anagnostidinema"
# [1] 0
# [1] "Haptophyte"
# [1] 0
# [1] "Monoraphidium.cf"
# [1] 3
# [1] "Westella.cf"
# [1] 2
# [1] "Dimorphococcus"
# [1] 2
# [1] "Quadrigula"
# [1] 0
# [1] "Rhabdogloea"
# > 

SIL_data.TITAN <- SIL_taxa %>% 
  select(-T.Anabaena.bergii, -T.Aphanizomenonaceae,- T.Dolichospermum.cf.circinale,
         -T.Dolichospermum.circinale,-T.Radiocystis,-Anagnostidinema,-Haptophyte,-Monoraphidium.cf,
         -Westella.cf,-Dimorphococcus,-Quadrigula,-Rhabdogloea)
#check occurrences again

colnames(SIL_data.TITAN)
taxa.TITAN = SIL_data.TITAN[,c(2:166)]

for (i in colnames(taxa.TITAN)){
  occurences <- sum(taxa.TITAN[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
##Good, no low occurrence species
#Now let's run TITAN

SIL.TITAN <- SIL_data.TITAN$SIL
taxa.TITAN
#TITAN Calculation
glades.titan_SIL <- titan(SIL.TITAN, taxa.TITAN)
save.image("TITAN_SIL_All_spp.RData")
savehistory("TITAN_SIL_All_spp.Rhistory")


#
#
#
# TURBILITY TITAN ANALYSIS (10)

data <- read.csv("HN_PHY_WQ_DATA_log_1.csv")
TURB = data$TURB
summary(TURB)
taxa <-  data[,c(15:191)]
df <-  data.frame(TURB, taxa)
df <-  na.omit(df)
TURB_taxa = df
TURB = TURB_taxa$TURB
species <- colnames(TURB_taxa)
taxa = TURB_taxa[,c(2:178)] 


for (i in colnames(taxa)){
  occurences <- sum(taxa[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}

TURB_data.TITAN <- TURB_taxa %>% 
  select(-T.Anabaena.bergii, -T.Aphanizomenonaceae, -T.Dolichospermum.cf.circinale)
#check occurrences again

colnames(TURB_data.TITAN)
taxa.TITAN = TURB_data.TITAN[,c(2:175)]

for (i in colnames(taxa.TITAN)){
  occurences <- sum(taxa.TITAN[[i]]>0)
  if (occurences < 4){
    print(occurences)
    print(i)
  }
}
##Good, no low occurrence species
#Now let's run TITAN

TURB.TITAN <- TURB_data.TITAN$TURB
taxa.TITAN
#TITAN Calculation
glades.titan_TURB <- titan(TURB.TITAN, taxa.TITAN)

save.image("TITAN_TURB_All_spp.RData")
savehistory("TITAN_TURB_All_spp.Rhistory")

#DONE TITAN FOR 10 VARIABLES
#SAVE ALL OUTPUT
save.image("Clustering_Data_pre.RData")



library(dplyr)
library(tibble)
library(clustMixType)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(networkD3)
#EXTRACT CHANGE POINT THRESHOLD DATA
load("Clustering_Data_pre.RData")


##TN (1)
spp_data_TN <- glades.titan_TN$sppmax
spp_data_TN <- data.frame(spp_data_TN)

colnames(spp_data_TN)#Only need spp, change point 10 and 90
spp_data_TN <- spp_data_TN[,c(9,11)]
spp_data_TN <- spp_data_TN %>% 
  rownames_to_column(var = "Taxa")
colnames(spp_data_TN)
spp_data_TN <- spp_data_TN %>% 
  rename(Taxa = Taxa,
         TN10 = X10.,
         TN90 = X90.)

colnames(spp_data_TN)

##TP (2)
spp_data_TP <- glades.titan_TP$sppmax
spp_data_TP <- data.frame(spp_data_TP)

colnames(spp_data_TP)#Only need spp, change point 10 and 90
spp_data_TP <- spp_data_TP[,c(9,11)]
spp_data_TP <- spp_data_TP %>% 
  rownames_to_column(var = "Taxa")
colnames(spp_data_TP)
spp_data_TP <- spp_data_TP %>% 
  rename(Taxa = Taxa,
         TP10 = X10.,
         TP90 = X90.)

colnames(spp_data_TP)

##TN_TP_Ratio (3)
spp_data_TN_TP <- glades.titan_TN_TP_Ratio$sppmax
spp_data_TN_TP <- data.frame(spp_data_TN_TP)

colnames(spp_data_TN_TP)#Only need spp, change point 10 and 90
spp_data_TN_TP <- spp_data_TN_TP[,c(9,11)]
spp_data_TN_TP <- spp_data_TN_TP %>% 
  rownames_to_column(var = "Taxa")
colnames(spp_data_TN_TP)
spp_data_TN_TP <- spp_data_TN_TP %>% 
  rename(Taxa = Taxa,
         TN_TP10 = X10.,
         TN_TP90 = X90.)

colnames(spp_data_TN_TP)



#TEMP (4)
spp_data_TEMP <- glades.titan_TEMP$sppmax
spp_data_TEMP <- data.frame(spp_data_TEMP)

colnames(spp_data_TEMP)#Only need spp, change point 10 and 90
spp_data_TEMP <- spp_data_TEMP[,c(9,11)]
spp_data_TEMP <- spp_data_TEMP %>% 
  rownames_to_column(var = "Taxa")
colnames(spp_data_TEMP)
spp_data_TEMP <- spp_data_TEMP %>% 
  rename(Taxa = Taxa,
         TEMP10 = X10.,
         TEMP90 = X90.)

colnames(spp_data_TEMP)

#AMM (5)

spp_data_AMM <- glades.titan_AMM$sppmax
spp_data_AMM <- data.frame(spp_data_AMM)

colnames(spp_data_AMM)#Only need spp, change point 10 and 90
spp_data_AMM <- spp_data_AMM[,c(9,11)]
spp_data_AMM <- spp_data_AMM %>% 
  rownames_to_column(var = "Taxa")
colnames(spp_data_AMM)
spp_data_AMM <- spp_data_AMM %>% 
  rename(Taxa = Taxa,
         AMM10 = X10.,
         AMM90 = X90.)

head(spp_data_AMM)

#NIT (6)

spp_data_NIT <- glades.titan_NIT$sppmax
spp_data_NIT <- data.frame(spp_data_NIT)

colnames(spp_data_NIT)#Only need spp, change point 10 and 90
spp_data_NIT <- spp_data_NIT[,c(9,11)]
spp_data_NIT <- spp_data_NIT %>% 
  rownames_to_column(var = "Taxa")
colnames(spp_data_NIT)
spp_data_NIT <- spp_data_NIT %>% 
  rename(Taxa = Taxa,
         NIT10 = X10.,
         NIT90 = X90.)

head(spp_data_NIT)

##FRP (7)
spp_data_frp <- glades.titan_FRP$sppmax
spp_data_frp <- data.frame(spp_data_frp)

colnames(spp_data_frp)#Only need spp, change point 10 and 90
spp_data_frp <- spp_data_frp[,c(9,11)]
spp_data_frp <- spp_data_frp %>% 
  rownames_to_column(var = "Taxa")
colnames(spp_data_frp)
spp_data_frp <- spp_data_frp %>% 
  rename(Taxa = Taxa,
         FRP10 = X10.,
         FRP90 = X90.)

head(spp_data_frp)

#SAL (8)
spp_data_sal <- glades.titan_SAL$sppmax
spp_data_sal <- data.frame(spp_data_sal)

colnames(spp_data_sal)#Only need spp, change point 10 and 90
spp_data_sal <- spp_data_sal[,c(9,11)]
spp_data_sal <- spp_data_sal %>% 
  rownames_to_column(var = "Taxa")
colnames(spp_data_sal)
spp_data_sal <- spp_data_sal %>% 
  rename(Taxa = Taxa,
         SAL10 = X10.,
         SAL90 = X90.)

head(spp_data_sal)
#TURB (9)
spp_data_turb <- glades.titan_TURB$sppmax
spp_data_turb <- data.frame(spp_data_turb)

colnames(spp_data_turb)#Only need spp, change point 10 and 90
spp_data_turb <- spp_data_turb[,c(9,11)]
spp_data_turb <- spp_data_turb %>% 
  rownames_to_column(var = "Taxa")
colnames(spp_data_turb)
spp_data_turb <- spp_data_turb %>% 
  rename(Taxa = Taxa,
         TURB10 = X10.,
         TURB90 = X90.)

head(spp_data_turb)

#SIL (10)
spp_data_SIL <- glades.titan_SIL$sppmax
spp_data_SIL <- data.frame(spp_data_SIL)

colnames(spp_data_SIL)#Only need spp, change point 10 and 90
spp_data_SIL <- spp_data_SIL[,c(9,11)]
spp_data_SIL <- spp_data_SIL %>% 
  rownames_to_column(var = "Taxa")
colnames(spp_data_SIL)
spp_data_SIL <- spp_data_SIL %>% 
  rename(Taxa = Taxa,
         SIL10 = X10.,
         SIL90 = X90.)

head(spp_data_SIL)


#EXPORT DATA FOR K-PROTOTYPE CLASSIFICATION WITH 4 MAIN FACTORS  (TN, TP, TNTP, and TEMPERATURE AS THE RESULTS OF MCA AND PCA)
merged_data <- spp_data_TP %>%
  full_join(spp_data_TN, by = "Taxa") %>%
  full_join(spp_data_TN_TP, by = "Taxa") %>%
  full_join(spp_data_TEMP, by = "Taxa") %>%
  full_join(spp_data_AMM, by = "Taxa") %>%
  full_join(spp_data_NIT, by = "Taxa") %>%
  full_join(spp_data_FRP, by = "Taxa") %>%
  full_join(spp_data_SIL, by = "Taxa") %>%
  full_join(spp_data_sal, by = "Taxa") %>%
  full_join(spp_data_turb, by = "Taxa") 

  data <- na.omit(merged_data)
 
  clustering_data <- data[,c(1:9)]
  
  write.csv(clustering_data,"Clustering_data_4.csv")
  write.csv(data,"Clustering_data_all.csv")
  
  #DONE TITAN ANALYSIS
  
  
  #PLOT TITAN: The ecological threshold range of phytoplankton taxa along the main environmental gradients (TNTP, TP, TN, and TEMP). 
  
  
  # Load necessary libraries
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  # Read the data
  data_total <- read.csv("Clustering_data_4.csv")
  
  ######TNTP
  data_TNTP <- data_total %>% 
    select(Taxa, TN_TP10, TN_TP90) %>% 
    mutate(mymean = (TN_TP10 + TN_TP90) / 2)
  
  data_TNTP$Taxa <- as.factor(data_TNTP$Taxa)
  data_TNTP$Taxa <- factor(data_TNTP$Taxa, levels = data_TNTP$Taxa[order(data_TNTP$Taxa)])
  
  P_TNTP <- ggplot(data_TNTP) +
    geom_segment(aes(x = Taxa, xend = Taxa, y = TN_TP10, yend = TN_TP90, col = mymean), linewidth = 0.5) +
    geom_point(aes(x = Taxa, y = TN_TP10), size = 1, col = "blue") +
    geom_point(aes(x = Taxa, y = TN_TP90), size = 1, col = "red") +
    coord_flip() +
    xlab("") +
    ylab("TNTP") +
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(colour = "black", size = 5),  # Adjust size as needed
          axis.text.x = element_text(colour = "black"))
  
  ######TP
  data_TP <- data_total %>% 
    select(Taxa, TP10, TP90) %>% 
    mutate(mymean = (TP10 + TP90) / 2)
  
  data_TP$Taxa <- as.factor(data_TP$Taxa)
  data_TP$Taxa <- factor(data_TP$Taxa, levels = data_TP$Taxa[order(data_TP$Taxa)])
  
  P_TP <- ggplot(data_TP) +
    geom_segment(aes(x = Taxa, xend = Taxa, y = TP10, yend = TP90, col = mymean), linewidth = 0.5) +
    geom_point(aes(x = Taxa, y = TP10), size = 1, col = "blue") +
    geom_point(aes(x = Taxa, y = TP90), size = 1, col = "red") +
    coord_flip() +
    xlab("") +
    ylab("TP (mg/L)") +
    theme(legend.position = "none") +
    theme(axis.line.x = element_blank(),
          axis.text.y = element_blank(), 
          axis.text.x = element_text(colour = "black"))
  
  ######TN
  data_TN <- data_total %>% 
    select(Taxa, TN10, TN90) %>% 
    mutate(mymean = (TN10 + TN90) / 2)
  
  data_TN$Taxa <- as.factor(data_TN$Taxa)
  data_TN$Taxa <- factor(data_TN$Taxa, levels = data_TN$Taxa[order(data_TN$Taxa)])
  
  P_TN <- ggplot(data_TN) +
    geom_segment(aes(x = Taxa, xend = Taxa, y = TN10, yend = TN90, col = mymean), linewidth = 0.5) +
    geom_point(aes(x = Taxa, y = TN10), size = 1, col = "blue") +
    geom_point(aes(x = Taxa, y = TN90), size = 1, col = "red") +
    coord_flip() +
    xlab("") +
    ylab("TN (mg/L)") +
    theme(legend.position = "none") +
    theme(axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(colour = "black"))
  
  ######TEMP
  data_TEMP <- data_total %>% 
    select(Taxa, TEMP10, TEMP90) %>% 
    mutate(mymean = (TEMP10 + TEMP90) / 2)
  
  data_TEMP$Taxa <- as.factor(data_TEMP$Taxa)
  data_TEMP$Taxa <- factor(data_TEMP$Taxa, levels = data_TEMP$Taxa[order(data_TEMP$Taxa)])
  
  P_TEMP <- ggplot(data_TEMP) +
    geom_segment(aes(x = Taxa, xend = Taxa, y = TEMP10, yend = TEMP90, col = mymean), linewidth = 0.5) +
    geom_point(aes(x = Taxa, y = TEMP10), size = 1, col = "blue") +
    geom_point(aes(x = Taxa, y = TEMP90), size = 1, col = "red") +
    coord_flip() +
    xlab("") +
    ylab(expression(paste("TEMP ("^"o","C)"))) +
    theme(legend.position = "none") +
    theme(axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(colour = "black"))
  
  # Combine the plots and reduce the gap
  p <- (P_TNTP | P_TP | P_TN | P_TEMP) + plot_layout(ncol = 4, widths = unit(rep(1, 4), "null"), heights = unit(rep(1, 1), "null")) &
    theme(plot.margin = margin(0, 0.4, 0, 0, "cm"))  # Adjust margins to reduce gaps
 # p <- (P_TNTP | P_TP | P_TN | P_TEMP)
  # Export the combined plot to a TIFF file with the specified dimensions and resolution
  tiff("TITAN_threshold1.tiff", width = 183, height = 240, units = "mm", res = 1000)
  print(p)
  dev.off()
  
  
  #DENSITY PLOT THRESHOLD

  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(ggthemes) 
  # Assuming data_total is your data frame with updated column names
  colnames(data_total)[colnames(data_total) == "TN_TP10"] <- "TNTP10"
  colnames(data_total)[colnames(data_total) == "TN_TP90"] <- "TNTP90"
  
  # Select and gather data for TNTP
  dataTNTP <- data_total %>%
    select(TNTP10, TNTP90) %>%
    gather(key = "variables", value = "value", TNTP10, TNTP90)
  
  # Create TNTP density plot
  TNTP <- ggplot(data = dataTNTP, aes(x = value, group = variables, fill = variables)) +
    geom_density(adjust = 1, alpha = 0.4) +
    xlab("TNTP") + ylab("") +
    #theme_minimal() +
    theme(legend.position = c(0.8, 0.8),
          legend.title = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  TNTP
  
  # Select and gather data for TP
  dataTP <- data_total %>%
    select(TP10, TP90) %>%
    gather(key = "variables", value = "value", TP10, TP90)
  
  # Create TP density plot
  TP <- ggplot(data = dataTP, aes(x = value, group = variables, fill = variables)) +
    geom_density(adjust = 1, alpha = 0.4) +
    xlab("TP (mg/L)") + ylab("Density") +
    #theme_minimal() +
    theme(legend.position = c(0.8, 0.8),
          legend.title = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  TP
  # Select and gather data for TN
  dataTN <- data_total %>%
    select(TN10, TN90) %>%
    gather(key = "variables", value = "value", TN10, TN90)
  
  # Create TN density plot
  TN <- ggplot(data = dataTN, aes(x = value, group = variables, fill = variables)) +
    geom_density(adjust = 1, alpha = 0.4) +
    xlab("TN (mg/L)") + ylab("Density") +
    #theme_minimal() +
    theme(legend.position = c(0.8, 0.8),
          legend.title = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  TN
  
  # Select and gather data for TEMP
  dataTEMP <- data_total %>%
    select(TEMP10, TEMP90) %>%
    gather(key = "variables", value = "value", TEMP10, TEMP90)
  
  # Create TEMP density plot
  TEMP <- ggplot(data = dataTEMP, aes(x = value, group = variables, fill = variables)) +
    geom_density(adjust = 1, alpha = 0.4) +
    xlab("TEMP (mg/L)") + ylab("Density") +
    theme(legend.position = c(0.8, 0.8),
          legend.title = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  TEMP
 
  # Export images
  
  TP
  tiff("TP.tiff", width = 89, height = 60,units = "mm", res = 1000)
  print(TP)
  dev.off()
  

  TN
  tiff("TN.tiff", width = 89, height = 60,units = "mm", res = 1000)
  print(TN)
  dev.off()
  
  TNTP
  tiff("TNTP.tiff", width = 89, height = 60,units = "mm", res = 1000)
  print(TNTP)
  dev.off()
  
  TEMP
  tiff("TEMP.tiff", width = 89, height = 60,units = "mm", res = 1000)
  print(TEMP)
  dev.off()

  #THEN USE INKCAPE TO COMBINE THOSE FIGURES   
  