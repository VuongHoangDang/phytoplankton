


# setwd("~/UWA_PhD/Data_analysis/CLustering_taxa/Clustering")
#setwd("D:/UWA_PhD/Data_analysis/CLASSIFICATION/CLustering/CLustering_taxa/Clustering")
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

data <- read.csv("K_prototype_classification_input.csv") 
df <- data.frame(data, row.names = 1)
df$Organism <- as.factor(df$Organism)
df$Movement <- as.factor(df$Movement)
df$Trophic <- as.factor(df$Trophic)
df$N_fixation <- as.factor(df$N_fixation)
df$Silica <- as.factor(df$Silica)
set.seed(7)

# apply k-prototypes

Es <- numeric(10)
for(i in 1:10){
  kpres <- kproto(df, k = i, nstart = 5)
  Es[i] <- kpres$tot.withinss
}
plot(1:10, Es, type = "b", ylab = "Objective Function", xlab = "Clusters", pch = 16)
# ggsave("k_number.png",  dpi=600)
# Based on this plot, k = 5 is the best number of cluster



cluster <- kproto(df, k = 5, lambda = NULL, nstart = 25)
# cluster$withinss
cluster$size
summary(cluster)
# cluster$cluster
df$cluster <- cluster$cluster
df$Clusters <- ifelse(df$cluster == 3, "A", 
                      ifelse(df$cluster == 5, "D",
                             ifelse(df$cluster == 2, "E",
                                    ifelse(df$cluster == 4, "B","C"))))
df$cluster <- as.factor(df$cluster)
df$Clusters <- as.factor(df$Clusters)
Clustering_data <- read_csv("Clustering_data.csv")
df$Group <- Clustering_data$Group
data_all <- read.csv("Change_points_10_90.csv")


# write.csv(df,"~/UWA_PhD/Data_analysis/CLustering_taxa/Clustering/Cluster_kprototye5.csv")
colnames(data_all)
df$AMM10 <- data_all$AMM10
df$AMM90 <- data_all$AMM90
df$NIT10 <- data_all$NIT10
df$NIT90 <- data_all$NIT90
df$SIL10 <- data_all$SIL10
df$SIL90 <- data_all$SIL90
df$FRP10 <- data_all$FRP10
df$FRP90 <- data_all$FRP90

#write.csv(df,"K_prototype_clusters_output.csv")



#COMPARE DIFFERENCE BETWEEN CLUSTERS
library(PMCMRplus)
# Normal distribution test
shapiro.test(df$TP90)
shapiro.test(df$TN90)
shapiro.test(df$TN_TP90)
shapiro.test(df$TEMP90)
shapiro.test(df$AMM90)
shapiro.test(df$NIT90)
shapiro.test(df$FRP90)
shapiro.test(df$SIL90)

shapiro.test(df$TP10)
shapiro.test(df$TN10)
shapiro.test(df$TN_TP10)
shapiro.test(df$TEMP10)
shapiro.test(df$AMM10)
shapiro.test(df$NIT10)
shapiro.test(df$FRP10)
shapiro.test(df$SIL10)

shapiro.test(df$BioVol_mm3)

#p-value of all the shapiro test < 0.05, data was not normal distribution, 
# Therefore, non-parametric test was selected

kruskal.test(TP90 ~ Clusters, data = df)
kruskal.test(TN90 ~ Clusters, data = df) # p-value = 1.033e-07
kruskal.test(TN_TP90 ~ Clusters, data = df) #p-value < 2.2e-16
kruskal.test(NIT90 ~ Clusters, data = df)#p-value = 7.212e-09
kruskal.test(AMM90 ~ Clusters, data = df)#p-value = 0.0006405
kruskal.test(FRP90 ~ Clusters, data = df) #p-value = 0.01578
kruskal.test(SIL90 ~ Clusters, data = df) #p-value = 0.000327
kruskal.test(TEMP90 ~ Clusters, data = df)

kruskal.test(TP10 ~ Clusters, data = df) #p-value = 3.509e-06
kruskal.test(TN10 ~ Clusters, data = df) #p-value = 0.002181
kruskal.test(TN_TP10 ~ Clusters, data = df) #p-value = 0.04084
kruskal.test(NIT10 ~ Clusters, data = df) #p-value = 0.03767
kruskal.test(AMM10 ~ Clusters, data = df) #p-value = 0.0001947
kruskal.test(FRP10 ~ Clusters, data = df)
kruskal.test(SIL10 ~ Clusters, data = df) #p-value = 0.003147
kruskal.test(TEMP10 ~ Clusters, data = df)


kruskal.test(BioVol_mm3 ~ Clusters, data = df) #p-value = 8.601e-05

kwAllPairsConoverTest(BioVol_mm3 ~ Clusters, data = df,
                      p.adjust.method = "holm")
# data: BioVol_mm3 by Clusters
# 
# A       B       C       D      
# B 1.00000 -       -       -      
#   C 1.00000 1.00000 -       -      
#   D 0.00896 0.07470 0.02821 -      
#   E 0.00055 0.01427 0.00768 1.00000

kwAllPairsConoverTest(TN90 ~ Clusters, data = df,
                      p.adjust.method = "holm")
# data: TN90 by Clusters
# 
# A       B       C       D      
# B 0.03781 -       -       -      
#   C 0.00012 0.07174 -       -      
#   D 0.01959 5.5e-05 1.8e-07 -      
#   E 0.18870 0.00240 7.9e-06 0.27728


kwAllPairsConoverTest(NIT90 ~ Clusters, data = df,
                      p.adjust.method = "holm")

# data: NIT90 by Clusters
# 
# A       B       C       D      
# B 0.00012 -       -       -      
#   C 6.2e-07 0.11877 -       -      
#   D 1.00000 0.00012 6.2e-07 -      
#   E 1.00000 0.00012 6.2e-07 1.00000

kwAllPairsConoverTest(TN_TP90 ~ Clusters, data = df,
                      p.adjust.method = "holm")
# data: TN_TP90 by Clusters
# 
# A      B      C      D     
# B <2e-16 -      -      -     
#   C <2e-16 0.0965 -      -     
#   D 0.2137 <2e-16 <2e-16 -     
#   E 0.0028 <2e-16 <2e-16 0.2137

kwAllPairsConoverTest(AMM90 ~ Clusters, data = df,
                      p.adjust.method = "holm")
# data: AMM90 by Clusters
# 
# A      B      C      D     
# B 0.0028 -      -      -     
#   C 0.1733 1.0000 -      -     
#   D 0.0229 1.0000 1.0000 -     
#   E 1.0000 0.0411 0.3871 0.1422

kwAllPairsConoverTest(FRP90 ~ Clusters, data = df,
                      p.adjust.method = "holm")

# data: FRP90 by Clusters
# 
# A     B     C     D    
# B 1.000 -     -     -    
#   C 1.000 1.000 -     -    
#   D 1.000 1.000 1.000 -    
#   E 0.024 0.024 1.000 1.000

kwAllPairsConoverTest(SIL90 ~ Clusters, data = df,
                      p.adjust.method = "holm")

# data: SIL90 by Clusters
# 
# A       B       C       D      
# B 0.30039 -       -       -      
#   C 0.05745 0.50508 -       -      
#   D 0.93008 0.50508 0.11884 -      
#   E 0.05745 0.00192 0.00034 0.11884

#DUNN TEST

library(FSA)
dunnTest(BioVol_mm3 ~ Clusters, data = df,, method="bonferron")

# Comparison          Z     P.unadj      P.adj
# 1       A - B -0.4473593 0.654615642 1.00000000
# 2       A - C  0.6245548 0.532263298 1.00000000
# 3       B - C  0.8650121 0.387032158 1.00000000
# 4       A - D -3.1018635 0.001923066 0.01923066
# 5       B - D -2.2993179 0.021486894 0.21486894
# 6       C - D -2.6795289 0.007372585 0.07372585
# 7       A - E -3.8749423 0.000106650 0.00106650
# 8       B - E -2.9310816 0.003377840 0.03377840
# 9       C - E -3.1772758 0.001486656 0.01486656
# 10      D - E -0.5632949 0.573234124 1.00000000

dunnTest(TN90 ~ Clusters, data = df,, method="bonferron")
# Comparison         Z      P.unadj        P.adj
# 1       A - B -2.322203 2.022201e-02 2.022201e-01
# 2       A - C -3.923035 8.744029e-05 8.744029e-04
# 3       B - C -2.015600 4.384181e-02 4.384181e-01
# 4       A - D  2.587807 9.658898e-03 9.658898e-02
# 5       B - D  4.113586 3.895591e-05 3.895591e-04
# 6       C - D  5.250678 1.515401e-07 1.515401e-06
# 7       A - E  1.487419 1.369042e-01 1.000000e+00
# 8       B - E  3.197970 1.383987e-03 1.383987e-02
# 9       C - E  4.526326 6.001795e-06 6.001795e-05
# 10      D - E -0.963489 3.353022e-01 1.000000e+00

dunnTest(TNTP90 ~ Clusters, data = df,, method="bonferron")
# Comparison          Z      P.unadj        P.adj
# 1       A - B -6.4150649 1.407629e-10 1.407629e-09
# 2       A - C -6.1913572 5.964835e-10 5.964835e-09
# 3       B - C -1.3495537 1.771592e-01 1.000000e+00
# 4       A - D  1.0124818 3.113077e-01 1.000000e+00
# 5       B - D  6.0942853 1.099275e-09 1.099275e-08
# 6       C - D  6.1814065 6.353301e-10 6.353301e-09
# 7       A - E  2.1559045 3.109113e-02 3.109113e-01
# 8       B - E  7.1696056 7.521412e-13 7.521412e-12
# 9       C - E  7.0141296 2.313849e-12 2.313849e-11
# 10      D - E  0.9156989 3.598248e-01 1.000000e+00

dunnTest(NIT90 ~ Clusters, data = df,, method="bonferron")
# Comparison           Z      P.unadj        P.adj
# 1       A - B -3.80611746 1.411654e-04 1.411654e-03
# 2       A - C -4.87767557 1.073433e-06 1.073433e-05
# 3       B - C -1.89418456 5.820052e-02 5.820052e-01
# 4       A - D  0.82701410 4.082291e-01 1.000000e+00
# 5       B - D  3.81012678 1.388955e-04 1.388955e-03
# 6       C - D  4.88995952 1.008567e-06 1.008567e-05
# 7       A - E  0.78473400 4.326095e-01 1.000000e+00
# 8       B - E  3.83239203 1.269033e-04 1.269033e-03
# 9       C - E  4.90887028 9.160254e-07 9.160254e-06
# 10      D - E -0.05382052 9.570782e-01 1.000000e+00
dunnTest(AMM90 ~ Clusters, data = df,, method="bonferron")
# Comparison          Z      P.unadj       P.adj
# 1       A - B -3.5331074 0.0004107057 0.004107057
# 2       A - C -2.0934779 0.0363065167 0.363065167
# 3       B - C  0.4520316 0.6512461916 1.000000000
# 4       A - D -2.9123645 0.0035870382 0.035870382
# 5       B - D  0.3766388 0.7064420165 1.000000000
# 6       C - D -0.1406282 0.8881637030 1.000000000
# 7       A - E -0.2845674 0.7759755933 1.000000000
# 8       B - E  2.6940287 0.0070594081 0.070594081
# 9       C - E  1.6872816 0.0915492127 0.915492127
# 10      D - E  2.2254240 0.0260527895 0.260527895
dunnTest(FRP90 ~ Clusters, data = df,, method="bonferron")
# Comparison          Z     P.unadj      P.adj
# 1       A - B  0.4958435 0.620004808 1.00000000
# 2       A - C -1.2045534 0.228375768 1.00000000
# 3       B - C -1.4239344 0.154465467 1.00000000
# 4       A - D -1.2203792 0.222321147 1.00000000
# 5       B - D -1.4518273 0.146549640 1.00000000
# 6       C - D  0.2311750 0.817178863 1.00000000
# 7       A - E -3.0075804 0.002633364 0.02633364
# 8       B - E -2.9757178 0.002923037 0.02923037
# 9       C - E -0.9494877 0.342372640 1.00000000
# 10      D - E -1.4396253 0.149973450 1.00000000
dunnTest(SIL90 ~ Clusters, data = df,, method="bonferron")
# Comparison           Z      P.unadj        P.adj
# 1       A - B -1.69182296 9.067973e-02 0.9067973440
# 2       A - C -2.57212226 1.010772e-02 0.1010771891
# 3       B - C -1.20878713 2.267446e-01 1.0000000000
# 4       A - D -0.08297574 9.338708e-01 1.0000000000
# 5       B - D  1.30667305 1.913238e-01 1.0000000000
# 6       C - D  2.22282790 2.622741e-02 0.2622740690
# 7       A - E  2.55104546 1.074003e-02 0.1074003173
# 8       B - E  3.58063168 3.427646e-04 0.0034276455
# 9       C - E  4.03175215 5.536256e-05 0.0005536256
# 10      D - E  2.16520289 3.037215e-02 0.3037215097

dunnTest(TP10 ~ Clusters, data = df,, method="bonferron")
# Comparison          Z      P.unadj        P.adj
# 1       A - B  2.5322863 1.133214e-02 1.133214e-01
# 2       A - C  3.5020681 4.616616e-04 4.616616e-03
# 3       B - C  1.4934591 1.353170e-01 1.000000e+00
# 4       A - D -0.3454011 7.297928e-01 1.000000e+00
# 5       B - D -2.3590611 1.832124e-02 1.832124e-01
# 6       C - D -3.3406528 8.358169e-04 8.358169e-03
# 7       A - E -2.3656801 1.799699e-02 1.799699e-01
# 8       B - E -4.1211924 3.769164e-05 3.769164e-04
# 9       C - E -4.7414498 2.121942e-06 2.121942e-05
# 10      D - E -1.6512529 9.868694e-02 9.868694e-01

dunnTest(TN10 ~ Clusters, data = df,, method="bonferron")
# Comparison          Z      P.unadj       P.adj
# 1       A - B -1.2075236 0.2272305800 1.000000000
# 2       A - C -0.4401391 0.6598363340 1.000000000
# 3       B - C  0.4045143 0.6858345875 1.000000000
# 4       A - D  3.3077948 0.0009403369 0.009403369
# 5       B - D  3.8239989 0.0001313045 0.001313045
# 6       C - D  2.6572221 0.0078787508 0.078787508
# 7       A - E  0.6308490 0.5281392352 1.000000000
# 8       B - E  1.5413696 0.1232268643 1.000000000
# 9       C - E  0.8211230 0.4115762268 1.000000000
# 10      D - E -2.2749199 0.0229107340 0.229107340

dunnTest(TN_TP10 ~ Clusters, data = df,, method="bonferron")
# Comparison            Z     P.unadj      P.adj
# 1       A - B -0.566284279 0.571200546 1.00000000
# 2       A - C  0.772285431 0.439945392 1.00000000
# 3       B - C  1.078349032 0.280878023 1.00000000
# 4       A - D  0.991500884 0.321441061 1.00000000
# 5       B - D  1.312653039 0.189299889 1.00000000
# 6       C - D -0.004921893 0.996072913 1.00000000
# 7       A - E  2.730847463 0.006317170 0.06317170
# 8       B - E  2.798384846 0.005135888 0.05135888
# 9       C - E  1.150694356 0.249857999 1.00000000
# 10      D - E  1.405607879 0.159840595 1.00000000

dunnTest(NIT10 ~ Clusters, data = df,, method="bonferron")
# Comparison           Z    P.unadj     P.adj
# 1       A - B -2.78222076 0.00539883 0.0539883
# 2       A - C  0.77381920 0.43903771 1.0000000
# 3       B - C  2.55545690 0.01060485 0.1060485
# 4       A - D -0.09259255 0.92622726 1.0000000
# 5       B - D  2.18650340 0.02877880 0.2877880
# 6       C - D -0.74930574 0.45367295 1.0000000
# 7       A - E -0.83591633 0.40320197 1.0000000
# 8       B - E  1.59996282 0.10960683 1.0000000
# 9       C - E -1.25915989 0.20797259 1.0000000
# 10      D - E -0.60834399 0.54295935 1.0000000

dunnTest(AMM10 ~ Clusters, data = df,, method="bonferron")
# Comparison          Z      P.unadj        P.adj
# 1       A - B -0.7899917 4.295326e-01 1.0000000000
# 2       A - C  1.2600431 2.076538e-01 1.0000000000
# 3       B - C  1.6702082 9.487818e-02 0.9487817639
# 4       A - D -1.4033840 1.605024e-01 1.0000000000
# 5       B - D -0.5617145 5.743105e-01 1.0000000000
# 6       C - D -2.0786515 3.764939e-02 0.3764939113
# 7       A - E -4.0989824 4.149706e-05 0.0004149706
# 8       B - E -2.8372649 4.550184e-03 0.0455018363
# 9       C - E -3.8990752 9.656079e-05 0.0009656079
# 10      D - E -2.1814542 2.914984e-02 0.2914984074
dunnTest(SIL10 ~ Clusters, data = df,, method="bonferron")
# Comparison           Z      P.unadj       P.adj
# 1       A - B  1.79553978 0.0725677406 0.725677406
# 2       A - C  0.04486016 0.9642187764 1.000000000
# 3       B - C -1.15501704 0.2480834499 1.000000000
# 4       A - D  1.55362081 0.1202748922 1.000000000
# 5       B - D -0.12825372 0.8979481904 1.000000000
# 6       C - D  1.02507278 0.3053288469 1.000000000
# 7       A - E  3.78195312 0.0001556027 0.001556027
# 8       B - E  1.73126246 0.0834049640 0.834049640
# 9       C - E  2.51351151 0.0119535878 0.119535878
# 10      D - E  1.79422561 0.0727771654 0.727771654


#TEST TO COMPARE CONTRIBUTION OF TAXONOMIC GROUPS TO FUNCTIONAL GROUPS
# COMPARE WITH CLUSTER COMPOSITION

# Load required libraries
library(ggstatsplot)
library(ggplot2)
library(dplyr)
library(gridExtra)
df <- read.csv("K_prototype_clusters_output.csv")
df <- df %>%
  rename(TNTP90 = TN_TP90)
df <- df %>%
  rename(TNTP10 = TN_TP10)
# Assuming data_all is your dataset
data_all = df
# List of groups to iterate over
groups <- unique(data_all$Clusters)

# Initialize an empty list to store results
results <- data.frame(Group = character(), 
                      TN90_pvalue = numeric(), 
                      TP90_pvalue = numeric(), 
                      TN_TP90_pvalue = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each group
for (group in groups) {
  
  # Filter data for the current group
  data_filter <- data_all %>%
    filter(Clusters == group)
  
  # Perform Kruskal-Wallis tests
  TN90_test <- kruskal.test(TN90 ~ Group, data = data_filter)
  TP90_test <- kruskal.test(TP90 ~ Group, data = data_filter)
  TN_TP90_test <- kruskal.test(TN_TP90 ~ Group, data = data_filter)
  
  # Store p-values in the results dataframe
  results <- rbind(results, data.frame(Group = group,
                                       TN90_pvalue = TN90_test$p.value,
                                       TP90_pvalue = TP90_test$p.value,
                                       TN_TP90_pvalue = TN_TP90_test$p.value))
}

# Export the results to a CSV file
#write.csv(results, "kruskal_test_pvalues.csv", row.names = FALSE)


# Load required libraries
library(ggstatsplot)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Assuming data_all is your dataset
data_all = df
# List of groups to iterate over
groups <- unique(data_all$Clusters)

# Initialize an empty list to store results
results <- data.frame(Group = character(), 
                      TN90_pvalue = numeric(), 
                      TP90_pvalue = numeric(), 
                      TN_TP90_pvalue = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each group
for (group in groups) {
  
  # Filter data for the current group
  data_filter <- data_all %>%
    filter(Clusters == group)
  
  # Perform Kruskal-Wallis tests
  TN90_test <- kruskal.test(TN90 ~ Group, data = data_filter)
  TP90_test <- kruskal.test(TP90 ~ Group, data = data_filter)
  TN_TP90_test <- kruskal.test(TNTP90 ~ Group, data = data_filter)
  
  # Store p-values in the results dataframe
  results <- rbind(results, data.frame(Group = group,
                                       TN90_pvalue = TN90_test$p.value,
                                       TP90_pvalue = TP90_test$p.value,
                                       TN_TP90_pvalue = TN_TP90_test$p.value))
}


# results of p-values for comparing mean of all taxonomic groups contributing to each functional group.
# We expect p-values > 0.05 to indicate there is no significant difference. It mean K-prototype picked up
# species of the taxonomic groups have the similar or the same threshold values to build functional groups

# Group TN90_pvalue TP90_pvalue TN_TP90_pvalue
# 1     B  0.02818403   0.3028439      0.1384170
# 2     D  0.42632184   0.5441839      0.5582994
# 3     A  0.16976591   0.6479533      0.4346965
# 4     E  0.08413258   0.1271829      0.1885804
# 5     C  0.23552423   0.2509048      0.5227535
# Export the results to a CSV file

#write.csv(results, "kruskal_test_pvalues.csv", row.names = FALSE)


#TEST TO COMPARE CONTRIBUTION OF TAXONOMIC GROUPS TO FUNCTIONAL GROUPS
# One taxonomic group contributed to many functional groups. We compare the taxonomic groups of different functional groups to answer
# the question are there any differences between taxonomic groups of different functional groups.
library(car)

data_all <- df
# Get unique groups (clusters)
groups <- unique(data_all$Clusters)
group_taxa <- unique(data_all$Group)
group_taxa

# Initialize an empty dataframe to store results
results_levene <- data.frame(Cluster = character(), 
                             TN90_pvalue = numeric(), 
                             TP90_pvalue = numeric(), 
                             TN_TP90_pvalue = numeric(),
                             stringsAsFactors = FALSE)

# Loop through each group
for (group in groups) {
  
  # Filter data for the current group
  data_filter <- data_all %>%
    filter(Clusters == group) %>% 
    select(TN90, TP90, TNTP90, Clusters)
  
  data_filter_1 <- data_all %>%
    filter(Clusters == group)
  in_group <- unique(data_filter_1$Group)
  
  
  data_filter <- data_filter %>%
    rename(Group = Clusters)
  
  # Extract TN90, TP90, and TN_TP90 values for the current group
  TN90_values <- data_filter$TN90
  TP90_values <- data_filter$TP90
  TN_TP90_values <- data_filter$TNTP90
  
  # Extract TN90, TP90, and TN_TP90 values for all groups
  
  groups_to_exclude <- c("Other", "Cryptophyta", "Chloromonadophyta", "Euglenophyta (Euglenoid)")
  
  data_groups <- data_all %>%
    filter(Group %in% in_group) %>%
    select(TN90, TP90, TNTP90, Group)
  
  
  
  # Combine data_filter and data_groups
  final_data <- rbind(data_filter %>% select(TN90, TP90, TNTP90, Group), data_groups)
  
  # Add a column indicating whether the row is from the cluster or the group
  #final_data$Source <- c(data_filter$Clusters, data_groups$Group)
  # Initialize variables for p-values
  TN90_pvalue <- NA
  TP90_pvalue <- NA
  TN_TP90_pvalue <- NA
  
  # Check if there are more than one unique value in the Source column
  if (length(unique(final_data$Group)) > 1) {
    # Perform Levene's tests if there is more than one unique value
    TN90_levene <- leveneTest(TN90 ~ Group, data = final_data)
    TP90_levene <- leveneTest(TP90 ~ Group, data = final_data)
    TN_TP90_levene <- leveneTest(TNTP90 ~ Group, data = final_data)
    
    # Store p-values
    TN90_pvalue <- TN90_levene$`Pr(>F)`[1]
    TP90_pvalue <- TP90_levene$`Pr(>F)`[1]
    TN_TP90_pvalue <- TN_TP90_levene$`Pr(>F)`[1]
  }
  
  # Append results to the dataframe
  results_levene <- rbind(results_levene, data.frame(Cluster = group,
                                                     TN90_pvalue = TN90_pvalue,
                                                     TP90_pvalue = TP90_pvalue,
                                                     TN_TP90_pvalue = TN_TP90_pvalue))
  
  
  # Print summary of results
  print(results_levene)
    }
#COMPARE CONTRIBUTION OF TAXONOMIC GROUP CONTRIBUTION

groups <- unique(data_all$Group)

# Initialize an empty dataframe to store results
results <- data.frame(Group = character(), 
                      TN90_pvalue = character(), 
                      TP90_pvalue = character(), 
                      TN_TP90_pvalue = character(),
                      stringsAsFactors = FALSE)

# Loop through each group
for (group in groups) {
  
  # Filter data for the current group
  data_filter <- data_all %>%
    filter(Group == group)
  
  # Initialize variables for p-values
  TN90_pvalue <- "samegroup"
  TP90_pvalue <- "samegroup"
  TN_TP90_pvalue <- "samegroup"
  
  # Check if there are more than one unique value in the Clusters column
  if (length(unique(data_filter$Clusters)) > 1) {
    # Perform Kruskal-Wallis tests if there is more than one unique value
    TN90_test <- kruskal.test(TN90 ~ Clusters, data = data_filter)
    TP90_test <- kruskal.test(TP90 ~ Clusters, data = data_filter)
    TN_TP90_test <- kruskal.test(TNTP90 ~ Clusters, data = data_filter)
    
    # Store p-values
    TN90_pvalue <- TN90_test$p.value
    TP90_pvalue <- TP90_test$p.value
    TN_TP90_pvalue <- TN_TP90_test$p.value
  }
  
  # Append results to the dataframe
  results <- rbind(results, data.frame(Group = group,
                                       TN90_pvalue = as.character(TN90_pvalue),
                                       TP90_pvalue = as.character(TP90_pvalue),
                                       TN_TP90_pvalue = as.character(TN_TP90_pvalue)))
}
print(results)


# 
# print(results)
# Group         TN90_pvalue        TP90_pvalue       TN_TP90_pvalue
# 1    Bacilariophyta   0.523512577563989   0.27645433847865 0.000110672194632279
# 2       Chlorophyta 5.9009008486091e-05 0.0541344134173458 2.69266387873459e-07
# 3        Cyanophyta 0.00352652548220963  0.584237866464717   5.876673241324e-05
# 4         Dinophyta   0.113846298006658  0.502819928630726    0.113846298006658
# 5       Cryptophyta           samegroup          samegroup            samegroup
# 6       Chrysophyta   0.144340543306881   0.20329701962372     0.15196525252243
# 7      Euglenophyta           samegroup          samegroup            samegroup
# 8 Chloromonadophyta           samegroup          samegroup            samegroup
# 9     Miscellaneous   0.220671361919847  0.479500122186953    0.220671361919847



#HERE, we expected p-values < 0.05 indicated significant difference between groups within a taxonomic group. It means K-prototype picked up
# species have similar threshold, therefore, even they are the same taxonomic groups, but they were assigned into different functional because 
# significant differerence in threshold values.



# Export the results to a CSV file
write.csv(results, "kruskal_test_pvalues_Taxonomic_contribution.csv", row.names = FALSE)


# The final test is COMPARING THE VARIANCE OF THRESHOLD BETWEEN PAIRS OF FUNCTIONAL GROUPS AND TAXONOMIC GROUPS
# if functional group show less variance, it will be better than taxonomic groups.

library(car)
library(dplyr)
library(ggstatsplot)
library(ggplot2)
library(gridExtra)
# Load the dataset
#df = read.csv("K_prototype_clusters_output.csv")
# df <- df %>%
#   rename(TNTP90 = TN_TP90)
# df <- df %>%
#   rename(TNTP10 = TN_TP10)
data_all <- df%>% 
  filter(Group %in% c("Chlorophyta", "Bacilariophyta", "Cyanophyta", "Dinophyta"))

# Get unique clusters (groups)

clusters <- unique(data_all$Clusters)
group_taxa <- unique(data_all$Group)

# Initialize an empty dataframe to store results
results_levene <- data.frame(Cluster = character(), 
                             Taxa = character(),
                             TN90_pvalue = numeric(), 
                             TP90_pvalue = numeric(), 
                             TN_TP90_pvalue = numeric(),
                             stringsAsFactors = FALSE)

# Loop through each cluster
for (i in clusters) {
  
  # Filter data for the current cluster
  data_cluster <- data_all %>%
    filter(Clusters == i) %>% select(TN90, TP90,TNTP90, Clusters)%>%
    rename(Group = Clusters)  # Rename 'Clusters' column to 'Group'
  
  # Loop through each group taxa
  for (u in group_taxa) {
    
    # Filter data for the current group taxa
    data_taxon <- data_all %>%
      filter(Group == u) %>% select(TN90, TP90,TNTP90, Group)
    
    # Skip if there are not enough data points
    if (nrow(data_taxon) < 2) {
      next
    }
    
    combined_data <- bind_rows(data_cluster, data_taxon)
    
    levene_testTP <- leveneTest(TP90 ~ Group, data = combined_data)
    levene_testTN <- leveneTest(TN90 ~ Group, data = combined_data)
    levene_testTNTP <- leveneTest(TNTP90 ~ Group, data = combined_data)
    TP90_pvalue <- levene_testTP$`Pr(>F)`[1]
    TN90_pvalue <- levene_testTN$`Pr(>F)`[1]
    TNTP90_pvalue <- levene_testTNTP$`Pr(>F)`[1]
    
    
    
    # Store results in the dataframe
    results_levene <- rbind(results_levene, 
                            data.frame(Cluster = i,
                                       Taxa = u,
                                       TN90_pvalue = TN90_pvalue,
                                       TP90_pvalue = TP90_pvalue,
                                       TN_TP90_pvalue = TNTP90_pvalue))
    #THIS IS FOR PLOTTING
    plot1 <- ggbetweenstats(
      data = combined_data,
      x = Group,
      y = TP90,type = "nonparametric",
      title = paste("Group:", i,u, "- TP90")
    )

    plot2 <- ggbetweenstats(
      data = combined_data,
      x = Group,
      y = TN90,
      title = paste("Group:", i,u, "- TN90")
    )

    plot3 <- ggbetweenstats(
      data = combined_data,
      x = Group,
      y = TNTP90,
      title = paste("Group:", i,u, "- TN_TP90")
    )
  } # End of group_taxa loop
  
  
  
  # Combine the plots horizontally
  plot <- grid.arrange(plot1, plot2, plot3, nrow = 1, ncol = 3)

  plot
  
} # End of clusters loop

# Print the final results
print(results_levene)

#Here, we expected p < 0.05 indicate signicant differences between variance, if there was significance, then we compare the variance values of the groups
# 
print(results_levene)
# Cluster           Taxa TN90_pvalue TP90_pvalue TN_TP90_pvalue
# 1        B Bacilariophyta  0.16718194   0.2388185   3.800312e-01
# 2        B    Chlorophyta  0.74830039   0.4722176   6.927306e-02
# 3        B     Cyanophyta  0.45944459   0.8248904   7.078004e-02
# 4        B      Dinophyta  0.49363695   0.5502332   9.463602e-01
# 5        D Bacilariophyta  0.84416701   0.8932269   4.847603e-02
# 6        D    Chlorophyta  0.23974214   0.7359568   1.963953e-02
# 7        D     Cyanophyta  0.15681105   0.5864642   1.773105e-02
# 8        D      Dinophyta  1.00000000   0.9551117   3.177532e-01
# 9        A Bacilariophyta  0.81723378   0.9414160   7.629936e-04
# 10       A    Chlorophyta  0.06673353   0.6181099   6.983299e-05
# 11       A     Cyanophyta  0.02618014   0.4357425   5.303954e-05
# 12       A      Dinophyta  0.98206104   0.9143083   8.600760e-02
# 13       E Bacilariophyta  0.95431082   0.5644312   2.256373e-02
# 14       E    Chlorophyta  0.25346786   0.7858211   1.715371e-02
# 15       E     Cyanophyta  0.17944733   0.9986448   1.493988e-02
# 16       E      Dinophyta  0.94505070   0.8144417   1.856364e-01
# 17       C Bacilariophyta  0.04199100   0.5689177   9.687429e-01
# 18       C    Chlorophyta  0.49846078   0.8126579   3.940519e-01
# 19       C     Cyanophyta  0.80708583   0.9668163   4.103413e-01
# 20       C      Dinophyta  0.28248838   0.7945337   6.779087e-01

write.csv(results_levene, "levene_test.csv")

#SENSITIVITY TEST FOR CLASSIFICATION
# We change the inputs for K-prototype
#Sensitivity Test: This test evaluated how different sets of input variables 
#(including traits and various environmental factors) influenced the classification of functional groups.


# (a) Original test: Includes traits, TN, TP, TNTP, and TEMP derived from MCA and PCA.
# (b) Classification test: Similar to the original test, but TN is replaced with NIT.
# (c) Classification test: Includes traits TN, TP, TNTP, TEMP, NIT, AMM, and SIL.
# (d) Classification test: Includes all variables - TN, TP, TNTP, TEMP, NIT, AMM, FRP, SAL, SIL, and TURB.

#TEST A - Original
data <- read.csv("K_prototype_classification_input.csv") 
df <- data.frame(data, row.names = 1)
df$Organism <- as.factor(df$Organism)
df$Movement <- as.factor(df$Movement)
df$Trophic <- as.factor(df$Trophic)
df$N_fixation <- as.factor(df$N_fixation)
df$Silica <- as.factor(df$Silica)
set.seed(7)

cluster <- kproto(df, k = 5, lambda = NULL, nstart = 25)
# cluster$withinss
cluster$size
df_test <- df
# cluster$cluster
df_test$cluster <- cluster$cluster
df_test$TestA <- ifelse(df_test$cluster == 3, "A", 
                      ifelse(df_test$cluster == 5, "D",
                             ifelse(df_test$cluster == 2, "E",
                                    ifelse(df_test$cluster == 4, "B","C"))))
colnames(df_test)
colnames(df)

#TEST B

data <- read.csv("K_prototype_classification_input.csv") 
df <- data.frame(data, row.names = 1)
df$Organism <- as.factor(df$Organism)
df$Movement <- as.factor(df$Movement)
df$Trophic <- as.factor(df$Trophic)
df$N_fixation <- as.factor(df$N_fixation)
df$Silica <- as.factor(df$Silica)
set.seed(7)

data_all <- read.csv("Change_points_10_90.csv")
df$NIT10 <- data_all$NIT10
df$NIT90 <- data_all$NIT90
df <- df %>% select(-TN10, -TN90)

# Select and reorder the columns
df <- df %>% select(TP10, TP90, NIT10, NIT90, TN_TP10, TN_TP90, TEMP10, TEMP90, BioVol_mm3, Organism, Movement, Trophic, N_fixation, Silica)

# Check the column names to confirm the new order
colnames(df)

clusterb <- kproto(df, k = 5, lambda = NULL, nstart = 25)
# cluster$withinss
cluster$size
clusterb$size
# cluster$cluster
df_test$clusterb <- clusterb$cluster
df_test$TestB <- ifelse(df_test$clusterb == 3, "A", 
                        ifelse(df_test$clusterb == 5, "D",
                               ifelse(df_test$clusterb == 2, "E",
                                      ifelse(df_test$clusterb == 4, "B","C"))))
colnames(df_test)

write.csv(df_test, "AB.csv")
#TEST C: Includes traits TN, TP, TNTP, TEMP, NIT, AMM, and SIL.
data <- read.csv("K_prototype_classification_input.csv") 
df <- data.frame(data, row.names = 1)
df$Organism <- as.factor(df$Organism)
df$Movement <- as.factor(df$Movement)
df$Trophic <- as.factor(df$Trophic)
df$N_fixation <- as.factor(df$N_fixation)
df$Silica <- as.factor(df$Silica)
set.seed(7)
data_all <- read.csv("Change_points_10_90.csv")
df$NIT10 <- data_all$NIT10
df$NIT90 <- data_all$NIT90
df$AMM10 <- data_all$AMM10
df$AMM90 <- data_all$AMM90
df$SIL10 <- data_all$SIL10
df$SIL90 <- data_all$SIL90
colnames(df)
# Select and reorder the columns
df <- df %>% select(TP10, TP90, TN10, TN90, TN_TP10, TN_TP90, TEMP10, TEMP90, NIT10, NIT90, AMM10, AMM90, SIL10, SIL90, BioVol_mm3, Organism, Movement, Trophic, N_fixation, Silica)
clusterc <- kproto(df, k = 5, lambda = NULL, nstart = 25)
# cluster$withinss
cluster$size
clusterb$size
clusterc$size
# cluster$cluster


df_test <- read.csv("AB.csv")

df_test$clusterc <- clusterc$cluster
df_test$TestC <- ifelse(df_test$clusterc == 3, "E", 
                        ifelse(df_test$clusterc == 5, "A",
                               ifelse(df_test$clusterc == 2, "B",
                                      ifelse(df_test$clusterc == 4, "C","D"))))

write.csv(df_test, "ABC.csv")
#TEST D: Includes all variables - TN, TP, TNTP, TEMP, NIT, AMM, FRP, SAL, SIL, and TURB.
data <- read.csv("K_prototype_classification_input.csv") 
df <- data.frame(data, row.names = 1)
df$Organism <- as.factor(df$Organism)
df$Movement <- as.factor(df$Movement)
df$Trophic <- as.factor(df$Trophic)
df$N_fixation <- as.factor(df$N_fixation)
df$Silica <- as.factor(df$Silica)
set.seed(7)

df$NIT10 <- data_all$NIT10
df$NIT90 <- data_all$NIT90
df$AMM10 <- data_all$AMM10
df$AMM90 <- data_all$AMM90
df$SIL10 <- data_all$SIL10
df$SIL90 <- data_all$SIL90
df$FRP10 <- data_all$FRP10
df$FRP90 <- data_all$FRP90
df$SAL10 <- data_all$SAL10
df$SAL90 <- data_all$SAL90
df$TURB10 <- data_all$TURB10
df$TURB90 <- data_all$TURB90

colnames(df)
# Select and reorder the columns
df <- df %>% select(TP10, TP90, TN10, TN90, TN_TP10, TN_TP90, TEMP10, TEMP90, NIT10, NIT90, AMM10, AMM90, SIL10, SIL90,FRP10, FRP90, SAL10, SAL90, TURB10, TURB90, BioVol_mm3, Organism, Movement, Trophic, N_fixation, Silica)
clusterd <- kproto(df, k = 5, lambda = NULL, nstart = 25)
# cluster$withinss
cluster$size
clusterb$size
clusterc$size
clusterd$size
# cluster$cluster
df_test$clusterd <- clusterd$cluster
df_test$TestD <- ifelse(df_test$clusterd == 3, "C", 
                        ifelse(df_test$clusterd == 5, "A",
                               ifelse(df_test$clusterd == 2, "D",
                                      ifelse(df_test$clusterd == 4, "B","E"))))


write.csv(df_test, "ABCD.csv")

final_sensitivity_test <- df_test %>% select(X,TestA, TestB, TestC, TestD)
final_sensitivity_test <- final_sensitivity_test %>%
  rename(Taxa = X)
write.csv(final_sensitivity_test, "Final_sensitivity_test.csv")




#VISUALIZE THE CLUSTERS

#Figure 9
#Export mean of clusters for Rada plot

mean_cluster <- df %>%
  group_by(Clusters) %>%
  summarise(across(c(TP10, TP90, TN10, TN90, TN_TP10, TN_TP90, TEMP10, TEMP90, AMM10, AMM90, NIT10, NIT90, SIL10, SIL90, FRP10, FRP90), mean, na.rm = TRUE))

# Print the result
print(mean_cluster)

rada_data <- mean_cluster %>% select(Clusters,AMM90, NIT90, TEMP90, TP90, TN_TP90,TN90 )
rada <- rada_data %>%
  rename(TNTP90 = TN_TP90)

new_data <- data.frame(
  Clusters = c("1", "2"),
  AMM90 = c(0.4, 0),
  NIT90 = c(3, 0),
  TEMP90 = c(32, 0),
  TP90 = c(0.1, 0),
  TN90 = c(4, 0),
  TNTP90 = c(300, 0)
)

data90 <- bind_rows( new_data, rada)


library(fmsb)
library(scales)
library(RColorBrewer)
library(jpeg)
# Set row names
rownames(data90) <- data90$Clusters
data90 <- data90 %>% select(-Clusters)

# Define colors
coul <- c("#FA746B", "#A4A401", "#4DD2A3","#4DC8F9", "#EE98FA")
colors_in <- alpha(coul, 1)

# Create the radar chart
radarchart(
  data90, axistype = 2, plwd = 2, plty = 2, pcol = colors_in,
  cglcol = "black", cglty = 2, axislabcol = "black", cglwd = 0.5, 
  vlcex = 1.1, title = ""
)

radarchart(
  data90, 
  axistype = 2, 
  plwd = 1, 
  plty = 1, 
  pty = 16,
  pcol = colors_in,
  cglcol = "black", 
  cglty = 2, 
  axislabcol = "black", 
  cglwd = 0.7, 
  vlcex = 1.1, # Adjusted label font size
  title = "", 
  caxislabels = c("10", "20", "30", "40", "50"), # Customize tick labels if needed
  cex.axis = 1 # Adjust the size of the axis tick labels
)


radarchart(data90, axistype=2, plwd=2 , plty=2,pcol = colors_in,
           #custom the grid
           cglcol="grey", cglty=3, axislabcol="black", cglwd=0.5, 
           #custom labels
           vlcex=1.1, title = "")


#legend(x=1.2, y=0.5, legend = rownames(data90[-c(1,2),]),col = colors_in, bty = "n", pch=20 , text.col = "grey", cex=1.2, pt.cex=3)

#SCATTER PLOTS
df <- df %>%
  rename(TNTP90 = TN_TP90)
df <- df %>%
  rename(TNTP10 = TN_TP10)

sc <-  df %>% 
  ggplot( aes(x=NIT90, y=TNTP90, color = Clusters)) +
  xlab("NIT90 (mg/L)") + ylab("TNTP90")+
  geom_point(alpha=0.7, size = 2) + 
  stat_ellipse(show.legend=FALSE)+
  
  # scale_size(range = c(3, 12)) +
  # scale_color_viridis(discrete=TRUE, guide=FALSE) +
  # labs(colour = "Functional groups", size = lab1)+
  # scale_color_manual(values = c("red", "blue", "black", "yellow", "purple"))+
  # theme_ipsum() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme_bw()+
  # guides(color = guide_legend(override.aes = list(size = 6)))+
  # scale_size_continuous(breaks=c(1,50,100,150,200,300))+
  
  
  
  theme(legend.text = element_text(size = 12, colour = "black"))+
  theme(legend.title = element_text(size = 12, color="black"))+
  theme(legend.key.size = unit(1, 'line'),              
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm')) +
  theme(text = element_text(size = 12, color="black"),
        axis.text = element_text(size = 12, color="black")) 
# +
#   theme(legend.position = c(1, 0.2))

sc


#BOX PLOT
###BOX TN
df_TN1 <- df
df_TN1$cluster <- as.factor(df_TN1$cluster)
df_TN1$Clusters <- as.factor(df_TN1$Clusters)

library("reshape2")
df_TN <- df_TN1 %>% 
  select(TN10,TN90,Clusters)
df_TN_long <- melt(df_TN, id = "Clusters") 
n <- ggplot(df_TN_long, aes(x = variable, y = value, color = Clusters)) +  # ggplot function
  geom_boxplot()+ 
  theme(legend.title=element_blank(),legend.position = "None")+
  xlab("")+ scale_y_continuous(limits = c(0, 4), breaks = c(0, 2, 4)) +
  ylab (expression(TN~( mg/L ))) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold.italic"),
    axis.text.x = element_text( color="black", 
                                size=12, angle=0),
    axis.text.y = element_text( color="black", 
                                size=12, angle=0),
    plot.margin = unit(c(0.0,0,0,0), "cm"))
# +
#   annotate("text", x = 0.8, y = 3.5, label = "(C)",size = 5)

n
df_TP <- df_TN1 %>% 
  select(TP10,TP90,Clusters)
df_TP_long <- melt(df_TP, id = "Clusters") 
p <- ggplot(df_TP_long, aes(x = variable, y = value, color = Clusters)) +  # ggplot function
  geom_boxplot()+scale_y_continuous(limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  theme(legend.title=element_blank(),legend.position = "None")+
  xlab("")+ 
  ylab (expression(TP~( mg/L ))) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold.italic"),
    axis.text.x = element_text( color="black", 
                                size=12, angle=0),
    axis.text.y = element_text( color="black", 
                                size=12, angle=0),
    plot.margin = unit(c(0.0,0,0,0), "cm"))
# +
#   annotate("text", x = 0.8, y = 0.085, label = "(E)",size = 5)
p
df_TN_TP <- df_TN1 %>% 
  select(TNTP10,TNTP90,Clusters)
colnames(df_TN_TP) <- c("TNTP10", "TNTP90", "Clusters")
df_TN_TP_long <- melt(df_TN_TP, id = "Clusters") 
r <- ggplot(df_TN_TP_long, aes(x = variable, y = value, color = Clusters)) +  # ggplot function
  geom_boxplot()+ scale_y_continuous(limits = c(0, 400), breaks = c(0, 200, 400)) +
  theme(legend.title=element_blank(),legend.position = "None")+
  xlab("")+ 
  ylab ("TNTP") +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text( color="black", 
                                size=12, angle=0),
    axis.text.y = element_text( color="black", 
                                size=12, angle=0),
    plot.margin = unit(c(0.0,0,0,0), "cm"))
# +
#   annotate("text", x = 0.8, y = 270, label = "(G)",size = 5)
r
df_TEMP <- df_TN1 %>% 
  select(TEMP10,TEMP90,Clusters)
df_TEMP_long <- melt(df_TEMP, id = "Clusters") 
t <- ggplot(df_TEMP_long, aes(x = variable, y = value, color = Clusters)) +  # ggplot function
  geom_boxplot()+ scale_y_continuous(limits = c(0, 40), breaks = c(0, 20, 40)) +
  theme(legend.title=element_blank(),legend.position = "None")+
  xlab("")+ 
  ylab (expression(paste(
    "TEMP (°C)"))) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold.italic"),
    axis.text.x = element_text( color="black", 
                                size=12, angle=0),
    axis.text.y = element_text( color="black", 
                                size=12, angle=0),
    plot.margin = unit(c(0.0,0,0,0), "cm"))
# +
#   annotate("text", x = 0.8, y = 27.5, label = "(H)",size = 5)
t
df_AMM <- df_TN1 %>% 
  select(AMM10,AMM90,Clusters)
df_AMM_long <- melt(df_AMM, id = "Clusters") 
a <- ggplot(df_AMM_long, aes(x = variable, y = value, color = Clusters)) +  # ggplot function
  geom_boxplot()+ scale_y_continuous(limits = c(0, 0.04), breaks = c(0, 0.02, 0.04)) +
  theme(legend.title=element_blank(),legend.position = "None")+
  xlab("")+ 
  ylab (expression(AMM~( mg/L ))) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold.italic"),
    axis.text.x = element_text( color="black", 
                                size=12, angle=0),
    axis.text.y = element_text( color="black", 
                                size=12, angle=0),
    plot.margin = unit(c(0.0,0,0,0), "cm"))
# +
#   annotate("text", x = 0.8, y = 0.035, label = "(D)",size = 5)

a
df_NIT <- df_TN1 %>% 
  select(NIT10,NIT90,Clusters)
df_NIT_long <- melt(df_NIT, id = "Clusters") 
ni <- ggplot(df_NIT_long, aes(x = variable, y = value, color = Clusters)) +  # ggplot function
  geom_boxplot()+theme(legend.title=element_blank(), legend.position = "None")+
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 2, 4)) +
  xlab("")+ 
  ylab (expression(NIT~( mg/L ))) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold.italic"),
    axis.text.x = element_text( color="black", 
                                size=12, angle=0),
    axis.text.y = element_text( color="black", 
                                size=12, angle=0),
    plot.margin = unit(c(0.0,0,0,0), "cm"))
# +
#   annotate("text", x = 0.8, y = 2.6, label = "(B)",size = 5)
ni
df_SIL <- df_TN1 %>% 
  select(SIL10,SIL90,Clusters)
df_SIL_long <- melt(df_SIL, id = "Clusters") 
s <- ggplot(df_SIL_long, aes(x = variable, y = value, color = Clusters)) +  # ggplot function
  geom_boxplot()+theme(legend.title=element_blank(), legend.position = "None")+
  scale_y_continuous(limits = c(0, 50), breaks = c(0, 25, 50)) +
  xlab("")+ 
  ylab (expression(SIL~( mg/L ))) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold.italic"),
    axis.text.x = element_text( color="black", 
                                size=12, angle=0),
    axis.text.y = element_text( color="black", 
                                size=12, angle=0),
    plot.margin = unit(c(0.0,0,0,0), "cm"))
# +
#   annotate("text", x = 0.8, y = 45, label = "(F)",size = 5)

s
#https://aslopubs.onlinelibrary.wiley.com/doi/pdf/10.4319/lo.1969.14.6.0941 silicar correlation with diatom population

df_TN1$BioVol_um3 <- df_TN1$BioVol_mm3*10^9
# b <- ggplot(df_TN1, aes(x = Clusters, y = BioVol_um3, color = Clusters)) +  # ggplot function
#   geom_boxplot()+ylim(c(0,8000))+ theme(legend.title=element_blank(), legend.position = "None")+
#   xlab("")+ 
#   ylab(expression(C_VOL~(µm^3))) +
#   theme(
#     axis.title.x = element_text(size = 12, face = "bold"),
#     axis.title.y = element_text(size = 12, face = "bold.italic"),
#     axis.text.x = element_text( color="black", 
#                                 size=12, angle=0),
#     axis.text.y = element_text( color="black", 
#                                 size=12, angle=0),
#     plot.margin = unit(c(0.5,0,0,0), "cm"))


b <- ggplot(df_TN1, aes(x = Clusters, y = BioVol_um3, color = Clusters)) +  # ggplot function
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 9000), breaks = c(0, 4500, 9000)) + 
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold.italic"),
    axis.text.x = element_text(color = "black", size = 12, angle = 0),
    axis.text.y = element_text(color = "black", size = 12, angle = 0),
    plot.margin = unit(c(0.0, 0, 0, 0), "cm")
  ) +
  xlab("") + 
  ylab(expression(C_VOL ~ (µm^3)))
# +
#   annotate("text", x = 1.2, y = 8000, label = "(A)",size = 5)
b
p1 <- (b/n/p/r)
p2 <- (ni/a/s/t)
f <- p1|p2
f

fi <- f|(sc/sc)
fi

# Save the plot to a TIFF file
tiff("combined_plot.tiff", width = 189, height = 130, units = "mm", res = dpi)

fi

dev.off()

#Then add scatter plot to fi

#Figure 8 SANKY PLOT

library(networkD3)
library(htmlwidgets)
library(webshot2)
library(magick)
library(dplyr)

# Assuming sanky_data and data_all are already loaded
sanky_data$Group <- data_all$Group

# Aggregate the data to count the occurrences
table_data <- sanky_data %>%
  group_by(Group, TestA) %>%
  summarise(value = n()) %>%
  ungroup() %>%
  rename(source = Group, target = TestA)

# Print the resulting table
print(table_data)

# Define the order of TestA
order_TestA <- c("A", "B", "C", "D", "E")

# Ensure TestA follows the desired order
table_data$target <- factor(table_data$target, levels = order_TestA)

data_long <- table_data

# Prepare nodes data frame
nodes <- data.frame(name = c(as.character(data_long$source), as.character(data_long$target)) %>% unique())

# Reformat data for networkD3
data_long$IDsource <- match(data_long$source, nodes$name) - 1 
data_long$IDtarget <- match(data_long$target, nodes$name) - 1

# Prepare colour scale
ColourScal <- 'd3.scaleOrdinal().range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

# Create the Sankey diagram
sankey <- sankeyNetwork(Links = data_long, Nodes = nodes,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name", 
                        sinksRight = FALSE, colourScale = ColourScal, 
                        nodeWidth = 40, fontSize = 13, nodePadding = 20, fontFamily = "Arial")

# Save the Sankey diagram as an HTML file
saveWidget(sankey, "sankey.html")

# Convert the HTML file to a temporary JPEG image
webshot2::webshot("sankey.html", file = "figure7_temp.jpeg", vwidth = 183 * 4, vheight = 120 * 4, delay = 0.2)

# Load the image
img <- image_read("figure7_temp.jpeg")

# Resize the image to the desired resolution
img <- image_resize(img, "7204x4724") # 183 mm x 120 mm at 1000 DPI

# Set the density to 1000 DPI and save as TIFF
image_write(img, path = "Figure_8.tiff", format = "tiff", density = "1000x1000")


colnames(sanky_data)

# SUPPLEMENTARY
# PLOT SENSITIVITY TEST (4 sanky plots)
library(networkD3)
library(dplyr)
library(magrittr)
library(webshot2)
library(magick)

# Loop to generate tables for TestA, TestB, TestC, and TestD
tests <- c("TestA", "TestB", "TestC", "TestD")  # List of tests
table_names <- c("TableA", "TableB", "TableC", "TableD")  # Desired table names

for (i in seq_along(tests)) {
  table_data <- sanky_data %>%
    group_by(Group, !!sym(tests[i])) %>%
    summarise(value = n()) %>%
    ungroup() %>%
    rename(source = Group, target = !!sym(tests[i]))
  
  # Print the table for the current test
  cat(paste(table_names[i], ":\n"))
  print(table_data)
  
  # Ensure TestA follows the desired order (if needed)
  # table_data$target <- factor(table_data$target, levels = order_TestA)
  
  data_long <- table_data
  
  # Prepare nodes data frame
  nodes <- data.frame(name = c(as.character(data_long$source), as.character(data_long$target)) %>% unique())
  
  # Reformat data for networkD3
  data_long$IDsource <- match(data_long$source, nodes$name) - 1 
  data_long$IDtarget <- match(data_long$target, nodes$name) - 1
  
  # Prepare colour scale
  ColourScal <- 'd3.scaleOrdinal().range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'
  
  # Create the Sankey diagram
  sankey <- sankeyNetwork(Links = data_long, Nodes = nodes,
                          Source = "IDsource", Target = "IDtarget",
                          Value = "value", NodeID = "name", 
                          sinksRight = FALSE, colourScale = ColourScal, 
                          nodeWidth = 40, fontSize = 13, nodePadding = 20, fontFamily = "Arial")
  
  # Save the Sankey diagram as an HTML file
  saveWidget(sankey, "sankey.html")
  
  # Convert the HTML file to a temporary JPEG image
  temp_jpeg <- paste("figure_", i, "_temp.jpeg", sep="")
  webshot2::webshot("sankey.html", file = temp_jpeg, vwidth = 183 * 4, vheight = 120 * 4, delay = 0.2)
  
  # Load the image
  img <- image_read(temp_jpeg)
  
  # Resize the image to the desired resolution
  img <- image_resize(img, "7204x4724") # 183 mm x 120 mm at 1000 DPI
  
  # Set the density to 1000 DPI and save as TIFF
  tiff_filename <- paste("Figure_", i, ".tiff", sep="")
  image_write(img, path = tiff_filename, format = "tiff", density = "1000x1000")
  
  # Clean up temporary files
  file.remove("sankey.html", temp_jpeg)
  
  cat(paste("Sankey diagram for", table_names[i], "exported as", tiff_filename, "\n\n"))
}


#

# Load necessary libraries
library(dplyr)
library(ggplot2)

data_input <- read.csv("K_prototype_clusters_output.csv")

# data_input$Group <- replace(data_input$Group, data_input$Group == "Bacilariophyta (Diatom)", "Bacilariophyta")
# data_input$Group <- replace(data_input$Group, data_input$Group == "Chlorophyta (Green algae)", "Chlorophyta")
# data_input$Group <- replace(data_input$Group, data_input$Group == "Cyanophyta (Blue green)", "Cyanophyta")
# data_input$Group <- replace(data_input$Group, data_input$Group == "Chrysophyta (Golden brown)", "Chrysophyta")
# data_input$Group <- replace(data_input$Group, data_input$Group == "Euglenophyta (Euglenoid)", "Euglenophyta")
# data_input$Group <- replace(data_input$Group, data_input$Group == "Cylindrospermopsis", "Cyanophyta")


unique(data_input$Group)

data_TNTP90 <- data_input %>% select(TN_TP90, Group, Clusters)


data_long <- data_TNTP90 %>%
  select(TN_TP90, Clusters = Group) %>%
  bind_rows(data_TNTP90 %>%
              select(TN_TP90, Clusters = Clusters))

# Read and filter the data
data <- data_long %>% 
  filter(Clusters %in% c("A", "B", "C", "D", "E", "Chlorophyta", "Bacilariophyta", "Cyanophyta", "Dinophyta"))

# Convert Clusters column to a factor with specified levels
data$Clusters <- factor(data$Clusters, levels = c("A", "B", "C", "D", "E", "Chlorophyta", "Bacilariophyta", "Cyanophyta", "Dinophyta"))

# Calculate mean and standard deviation
summary_data <- data %>%
  group_by(Clusters) %>%
  summarize(mean_value = mean(TN_TP90, na.rm = TRUE),
            sd_value = sd(TN_TP90, na.rm = TRUE))

# Create the plot with mean and standard deviation and add text annotations
p <- ggplot(summary_data, aes(x = Clusters, y = mean_value)) +
  geom_point(color = "red", size = 3) +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2, color = "blue") +
  geom_text(aes(label = sprintf("Mean: %.2f\nSD: %.2f", mean_value, sd_value)), 
            vjust = -1.5, color = "black", size = 3) +
  ylim(0, 400) +  # Set y-axis limit
  labs(title = "",
       x = "",
       y = "TNTP90") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(size = 10, color = "black", angle = -45),  # X-axis tick settings
        axis.text.y = element_text(size = 10, color = "black"),  # Y-axis tick settings
        axis.title = element_text(size = 12, color = "black"))   # Axis label settings

# Display the plot
print(p)

# Save the plot as a TIFF file with specified resolution, height, and width
ggsave("TNTP_Mean_SD.tiff", plot = p, dpi = 1000, height = 100, width = 189, units = "mm")
library(svglite)
ggsave("TNTP_Mean_SD.svg", plot = p, dpi = 1000, height = 100, width = 189, units = "mm")
