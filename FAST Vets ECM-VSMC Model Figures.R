## Make FAST Vets ECM-VSM Model Plots

setwd("C:/Users/pewow/Documents/FAST Vets/ECM - VSMC Models") ## Set working directory
Carotid <- read.csv("ECM-VSMC Model Carotid.csv") ## Load carotid data
cfPWV <- read.csv("ECM-VSMC Model cfPWV.csv") ## Load cfPWV data
library(ggplot2) ## Load ggplot library for graphics

## Convert Sex (0=Male / 1=Female) to a factor variable
Carotid$Sex <- factor(Carotid$Sex) 
cfPWV$Sex <- factor(cfPWV$Sex)

########## Figure 1 ##########

## Define model parameters
betaECM <- c(10.382266)
betaVSMC <- c(5.037054)
k <- c(0.068805)
Dref <- c(7.450829)

## Measured participant data
P_Base <- c(77, 135) 
P_NTG <- c(70, 114)
D_Base <- c(6.735000,
       7.290000)
D_NTG <- c(7.355,
           7.705)

## Set diameters for model curves
Dm <- seq(6 , 8, by=0.01)

## Calculate ECM and VSMC portions of model curves
Pm_ECM <- 80*exp(betaECM*(Dm/Dref -1))
Pm_VSMC <- 80*(k/0.1 * exp(betaVSMC*(Dm/(Dref*(1-k)) -1) ) )

## Total model curve = ECM curve + VSMC curve
Pm_Total <- Pm_ECM + Pm_VSMC

## Create polygons to represent amount of BP that the ECM and VSMC support
D_ECM_poly <- c(6.73, seq(6.73, 7.29, by=0.01), 7.29)
P_ECM_poly <- c(0, 80*exp(betaECM*(seq(6.73, 7.29, by=0.01)/Dref -1)), 0 )

D_VSMC_poly <- c(seq(6.73, 7.29, by=0.01), seq(7.29, 6.73, by=-0.01))
P_VSMC_poly <- c(80*exp(betaECM*(seq(6.73, 7.29, by=0.01)/Dref -1)),  
                80*(exp(betaECM*(seq(7.29, 6.73, by=-0.01)/Dref -1)) + 
                      k/0.1 * exp(betaVSMC*(seq(7.29, 6.73, by=-0.01)/(Dref*(1-k)) -1) )))

## Create Figure with ggplot
Data.Model <- data.frame( 
                         D=D_NTG,
                         P = P_NTG
                         )
## Add plot elements in order from background to foreground
F1 <- ggplot(Data.Model) +
  annotate("polygon", x=D_ECM_poly, y = P_ECM_poly, color="lightblue", fill="lightblue")  +
  annotate("polygon", x=D_VSMC_poly, y = P_VSMC_poly, color="pink", fill="pink") +
  annotate("line", x=Dm, y = Pm_Total, color="black", linewidth=1) +
  annotate("line", x=Dm, y = Pm_ECM, color="blue", linewidth=1) +
  annotate("line", x=Dm, y = Pm_VSMC, color="darkred", linewidth=1) +
  annotate("point", x=D_Base, y = P_Base, shape=21, size=3, color="darkgray", fill="white", stroke=2) +
  annotate("point", x=D_NTG, y = P_NTG, shape=22, size=3, color="darkgray",fill="white", stroke=2) +
   ylim(0, 150) +
  xlab("Diameter (mm)") +
  ylab("Pressure (mmHg)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
        text = element_text(size = 16)) ## Make background blank, axis line black, and fontsize 16

F1

## Make extra plot with elements to make legend in power point. I'm an amateur 
## R user and couldn't figure out how to customize the legend in R

F1_L <- ggplot(Data.Model) +
  annotate("line", x=c(6.5, 6.6), y = c(100, 100), color="black", linewidth=1) +
  annotate("line", x=c(6.5, 6.6), y = c(90, 90), color="blue", linewidth=1) +
  annotate("line", x=c(6.5, 6.6), y = c(80, 80), color="darkred", linewidth=1) +
  annotate("point", x=c(6.55), y = c(120), shape=21, size=3, color="darkgray", fill="white", stroke=2) +
  annotate("point", x=c(6.55), y = c(110), shape=22, size=3, color="darkgray",fill="white", stroke=2) +
  ylim(0, 150) +
  xlim(6, 8) +
  xlab("Diameter (mm)") +
  ylab("Pressure (mmHg)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
      text = element_text(size = 16)) ## Make background blank, axis line black, and fontsize 16

F1_L


########## Figure 2 ##########

ggplot(cfPWV, aes(Age, Beta.ECM..DBP., colour = Sex)) +
  # to create a Scatter plot
  geom_point(aes(shape=Sex, color=Sex), size=2.5) +
  # to fit and overlay a loess trendline
  stat_smooth(formula = y ~ x, method = "loess", se=FALSE, span=0.75) +
  ylim(0, 60)+
  ylab("ECM Stiffness")+
  xlab("Age (Years)") +
  ggtitle("A. cfPWV ECM Stiffness") +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
        text = element_text(size = 20)) ## Get rid of auto legend, make background blank, axis line black, and fontsize 16


ggplot(cfPWV, aes(Age, Beta.VSMC..DBP., colour = Sex)) +
  # to create a scatterplot
  geom_point(aes(shape=Sex, color=Sex), size=2.5) +
  # to fit and overlay a loess trendline
  stat_smooth(formula = y ~ x, method = "loess", se=FALSE, span=0.75) +
  ylim(0, 60)+
  ylab("VSMC Stiffness")+
  xlab("Age (Years)") +
  ggtitle("B. cfPWV VSMC Stiffness") +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
        text = element_text(size = 20)) ## Get rid of auto legend, make background blank, axis line black, and fontsize 16


ggplot(Carotid, aes(Age, BetaECM, colour = Sex)) +
  # to create a scatterplot
  geom_point(aes(shape=Sex, color=Sex), size=2.5) +
  # to fit and overlay a loess trendline
  stat_smooth(formula = y ~ x, method = "loess", se=FALSE, span=0.75) +
  ylim(0, 30)+
  ylab("ECM Stiffness")+
  xlab("Age (Years)") +
  ggtitle("C. Carotid ECM Stiffness") +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
        text = element_text(size = 20)) ## Get rid of auto legend, make background blank, axis line black, and fontsize 16


ggplot(Carotid, aes(Age, betaVSM, colour = Sex)) +
  # to create a scatterplot
  geom_point(aes(shape=Sex, color=Sex), size=2.5) +
  # to fit and overlay a loess trendline
  stat_smooth(formula = y ~ x, method = "loess", se=FALSE, span=0.75) +
  ylim(0, 40)+
  ylab("VSMC Stiffness")+
  xlab("Age (Years)") +
  ggtitle("D. Carotid VSMC Stiffness") +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
        text = element_text(size = 20)) ## Get rid of auto legend, make background blank, axis line black, and fontsize 16


ggplot(Carotid, aes(Age, Dref, colour = Sex)) +
  # to create a scatterplot
  geom_point(aes(shape=Sex, color=Sex), size=2.5) +
  # to fit and overlay a loess trendline
  stat_smooth(formula = y ~ x, method = "loess", se=FALSE, span=0.75) +
  ylim(5, 10)+
  ylab("Reference Diameter (mm)")+
  xlab("Age (Years)") +
  ggtitle("E. Carotid Reference Diameter") +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
        text = element_text(size = 20)) ## Get rid of auto legend, make background blank, axis line black, and fontsize 16


ggplot(Carotid, aes(Age, 100*k, colour = Sex))+
  # to create a scatterplot
  geom_point(aes(shape=Sex, color=Sex), size=2.5) +
  # to fit and overlay a loess trendline
  stat_smooth(formula = y ~ x, method = "loess", se=FALSE, span=0.75) +
  ylim(0, 15)+
  ylab("VSMC Tone (%)")+
  xlab("Age (Years)") +
  ggtitle("F. Carotid VSMC Tone") +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
        text = element_text(size = 20)) ## Get rid of auto legend, make background blank, axis line black, and fontsize 16

## Make an extra plot with a legend so I can modify it in power point. I'm an amateur 
## R user and couldn't figure out how to customize the legend in R
ggplot(Carotid, aes(Age, 100*k, colour = Sex))+
  # to create a scatterplot
  geom_point(aes(shape=Sex, color=Sex), size=2.5) +
  # to fit and overlay a loess trendline
  stat_smooth(formula = y ~ x, method = "loess", se=FALSE, span=0.75) +
  ylim(0, 15)+
  ylab("VSMC Tone (%)")+
  xlab("Age (Years)") +
  ggtitle("F. Carotid VSMC Tone") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
        text = element_text(size = 20)) ## make background blank, axis line black, and fontsize 16


