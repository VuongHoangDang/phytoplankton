setwd("D:/UWA_PhD/Paper 1/Script")

library(jpeg)
library(tidyr)
library(dbplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(patchwork)
library(zoo)
library(hrbrthemes)
library(PerformanceAnalytics)
library(corrplot)
library(patchwork)
library(Kendall)
library(date)
library(ggfortify)
library(ggplot2)

#MULTICORRELATION ANALYSIS

df <- read.csv("PHY_WQ_dataset_MCA.csv") %>%
  filter(Site =="N42")
colnames(df)
df$TN_TP <- (df$TN/14)/(df$TP/31) #Calculate TN TP Ratio as molecular
df$Quater_year <- as.yearqtr(df$Date  ,           # Convert dates to quarterly (seasonally)
                             format = "%d-%m-%y")

df <- df %>%  unite(SQY, Site, Quater_year, sep = "", remove = FALSE)


df_mean <- data.frame(aggregate(cbind(TN, TP,TN_TP, SAL, TEMP, AMM, NIT,FRP, SIL, TURB,CHLA, BIOVOL,BACIL,
                                      CHLOROPHYTA, CHLOROMONAD, CHRYS, CRYPT, CYANO, DINO, EUGLE, MIS)~ SQY, data = df, function(x) c(mean = mean(x), sd = sd(x), na.rm = TRUE))) 

d_mean <- data.frame(
  TN  = df_mean$TN[,"mean"],
  TP  = df_mean$TP[,"mean"],
  TNTP = df_mean$TN_TP[,"mean"],
  SAL  = df_mean$SAL[,"mean"],
  TEMP  = df_mean$TEMP[,"mean"],
  AMM  = df_mean$AMM[,"mean"],
  NIT  = df_mean$NIT[,"mean"],
  FRP  = df_mean$FRP[,"mean"],
  TURB  = df_mean$TURB[,"mean"],
  SIL  = df_mean$SIL[,"mean"],
  CHLA  = df_mean$CHLA[,"mean"],
  VOL  = df_mean$BIOVOL[,"mean"])

d_mean <- na.omit(d_mean)


tiff("corrplot_mixed.tiff", width = 120, height = 120, units = "mm", res = 2000)

# Set the text size for axis labels and tick marks
par(cex.lab = 0.5, cex.axis = 0.5) # Adjust these values as needed for your desired size

# Generate the plot
corrplot.mixed(cor(d_mean),
               lower = "number",
               upper = "circle",
               tl.col = "black",
               tl.pos = "lt", 
               mar = c(1, 1, 1, 1),
               # cl.cex = 0.8,
               # cl.ratio = 0.15, 
               number.cex = 0.7,
               # number.font = 1
)

# Close the TIFF device
dev.off()

#PRINCIPLE COMPONENT ANALYSIS


data <- read.csv("PHY_WQ_dataset_PCA.csv")%>% 
  filter(Site.Code != "N3001" & Site.Code != "NB13") 

## Organic Phosphorus (OP) = TP - FRP
## Organic Nitrogen (ON) =  TN - NIT - AMM

data[data == "NB11"] <- "N11"

df <- data %>% select(-X,-BACIL,-CHLOROMONAD,-CHRYS,-CRYPT,-CYANO,-DINO,-EUGLE,-MIS,-CHLOROPHYTA,
                      -Calendar.Date,-Year,-Site.Code,-Cells.mL,-CellVolume.mm3.,-CHLA)

df <- na.omit(df)


d_mean <- data.frame(aggregate(cbind(TN, TP, SAL, TEMP, AMM, NIT,FRP,OP,ON, SIL, TURB,BIOVOL)~ Genus.Name, data = df, function(x) c(mean = mean(x), na.rm = TRUE))) 

df_pca <- data.frame(
  Spp = d_mean$Genus.Name,
  TN  = d_mean$TN[,"mean"],
  TP  = d_mean$TP[,"mean"],
  SAL  = d_mean$SAL[,"mean"],
  TEMP  = d_mean$TEMP[,"mean"],
  AMM  = d_mean$AMM[,"mean"],
  NIT  = d_mean$NIT[,"mean"],
  FRP  = d_mean$FRP[,"mean"],
  OP  = d_mean$OP[,"mean"],
  ON  = d_mean$ON[,"mean"],
  SIL  = d_mean$SIL[,"mean"],
  TURB  = d_mean$TURB[,"mean"]
)         

pca <- df_pca[,2:12] %>% select(-OP,- ON)

res.pca <- prcomp(pca, scale. = TRUE)

# Functions to make circle
gg_circle <- function(r, xc, yc, color="black", fill=NA, ...) {
  x <- xc + r*cos(seq(0, pi, length.out=100))
  ymax <- yc + r*sin(seq(0, pi, length.out=100))
  ymin <- yc + r*sin(seq(0, -pi, length.out=100))
  annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
}


#Plot PCA

PCA_quater <- autoplot(res.pca, data = df_pca, colour = "blue", loadings = TRUE,  loadings.colour = 'black',
                       loadings.label.size = 5,loadings.label.colour = "black",
                       loadings.label = TRUE, ylim =  c(-0.4,0.2), xlim = c(-0.3, 0.4))+ 
  theme_minimal() +
  theme(legend.position = c(0.83, 0.73),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 11),
        legend.key = element_rect(fill = "white", color = NA),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text( color="black", 
                                    size=12, angle=0),
        axis.text.y = element_text( color="black", 
                                    size=12, angle=0),
        
        plot.margin = unit(c(1,1,1,1), "cm"))+
  gg_circle(r=0.1, xc=0.23, yc=-0.06, color="blue", fill="blue", alpha=0.1)+
  gg_circle(r=0.09, xc=-0.01, yc= -0.26, color="red", fill="red", alpha=0.1)+
  gg_circle(r=0.12, xc=-0.15, yc= -0.1, color="green", fill="green", alpha=0.1)


PCA_quater

tiff("PCA.tiff", width = 120, height = 110, units = "mm", res = 2000)
PCA_quater
dev.off()



