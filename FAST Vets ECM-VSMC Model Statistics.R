## Make FAST Vets ECM-VSM Model Plots

setwd("C:/Users/pewow/Documents/FAST Vets/ECM - VSMC Models") ## Set working directory
Carotid <- read.csv("ECM-VSMC Model Carotid.csv") ## Load carotid data
cfPWV <- read.csv("ECM-VSMC Model cfPWV.csv") ## Load cfPWV data
library(ggplot2) ## Load ggplot library for graphics
library(rms) ## Load RMS package for multiple linear regression

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

################################################################################
## Carotid Mediation Analyses
library(lavaan) # Load lavaan package for mediation analysis

c_med_data <- data.frame( betaECM = scale(sqrt(Carotid$BetaECM)) ,
                          betaVSM = scale(sqrt(Carotid$betaVSM)),
                          k = scale(sqrt(Carotid$k)),
                          PWV80 = sqrt( 0.5*133.32*80*Carotid$Beta80.Baseline /1050),
                          Age = Carotid$Age,
                          Sex = Carotid$Sex,
                          Smoke = Carotid$FormCurr.Smoke,
                          DM = Carotid$dm,
                          HT = Carotid$Has.HTN,
                          BPMed = Carotid$BP.Med )

c_model <- "
  # Regression paths from co-variates to stiffness parameters
  betaECM ~ a1*Age + a2*Sex + a3*Smoke +a4* DM + a5*HT
  betaVSM ~ b1*Age + b2*Sex + b3*Smoke + b4*DM + b5*HT
  k ~ c1*Age + c2*Sex + c3*Smoke + c4*DM + c5*HT

  # Regression paths from stiffness parameters to PWV80
  PWV80 ~ d1*Age + d2*Sex + d3*Smoke + d4*DM + d5*HT + d6*betaECM + d7*betaVSM + d8*k
  
  # Indirect Effects per 10-years age
  total_Age := 10*(d1 + a1*d6 + b1*d7 + c1*d8)
  direct_Age := 10*d1 
  indirect_ECM := 10*a1*d6
  indirect_VSM := 10*b1*d7
  indirect_k := 10*c1*d8 " 

c_fit <- cfa(c_model, data=c_med_data, se="bootstrap", bootstrap="2000")
parameterEstimates(c_fit)

################################################################################
## cfPWV Mediation Analyses
library(lavaan)

cf_med_data <- data.frame( betaECM = scale(sqrt(cfPWV$Beta.ECM..DBP.)) ,
                           betaVSM = scale(sqrt(cfPWV$Beta.VSMC..DBP.)),
                           PWV80 = sqrt( 0.5*133.32*80*cfPWV$Beta80 /1050),
                           Age = cfPWV$Age,
                           Sex = cfPWV$Sex,
                           Smoke = cfPWV$FormCurr.Smoke,
                           DM = cfPWV$dm,
                           HT = cfPWV$Has.HTN,
                           BPMed = cfPWV$BP.MEd )

cf_model <- "
  # Regression paths from co-variates to stiffness parameters
  betaECM ~ a1*Age + a2*Sex + a3*Smoke +a4* DM + a5*HT
  betaVSM ~ b1*Age + b2*Sex + b3*Smoke + b4*DM + b5*HT

  # Regression paths from stiffness parameters to PWV80
  PWV80 ~ d1*Age + d2*Sex + d3*Smoke + d4*DM + d5*HT + d6*betaECM + d7*betaVSM
  
  # Indirect Effects per 10-years age
  total_Age := 10*(d1 + a1*d6 + b1*d7)
  direct_Age := 10*d1 
  indirect_ECM := 10*a1*d6
  indirect_VSM := 10*b1*d7 "

cf_fit <- cfa(cf_model, data=cf_med_data, se="bootstrap", bootstrap="2000")
parameterEstimates(cf_fit)

################################################################################
## Carotid Mediation Sensitivity Analysis
c_model_BPMed <- "
  # Regression paths from co-variates to stiffness parameters
  betaECM ~ a1*Age + a2*Sex + a3*Smoke +a4* DM + a5*HT + a6*BPMed
  betaVSM ~ b1*Age + b2*Sex + b3*Smoke + b4*DM + b5*HT + b6*BPMed
  k ~ c1*Age + c2*Sex + c3*Smoke + c4*DM + c5*HT + c6*BPMed

  # Regression paths from stiffness parameters to PWV80
  PWV80 ~ d1*Age + d2*Sex + d3*Smoke + d4*DM + d5*HT + d6*betaECM + d7*betaVSM + d8*k + d9*BPMed
  
  # Indirect Effects per 10-years age
  total_Age := 10*(d1 + a1*d6 + b1*d7 + c1*d8)
  direct_Age := 10*d1 
  indirect_ECM := 10*a1*d6
  indirect_VSM := 10*b1*d7
  indirect_k := 10*c1*d8 " 

c_fit_BPMed <- cfa(c_model_BPMed, data=c_med_data, se="bootstrap", bootstrap="2000")
parameterEstimates(c_fit_BPMed)

c_model_noHT <- "
  # Regression paths from co-variates to stiffness parameters
  betaECM ~ a1*Age + a2*Sex + a3*Smoke +a4* DM 
  betaVSM ~ b1*Age + b2*Sex + b3*Smoke + b4*DM 
  k ~ c1*Age + c2*Sex + c3*Smoke + c4*DM 

  # Regression paths from stiffness parameters to PWV80
  PWV80 ~ d1*Age + d2*Sex + d3*Smoke + d4*DM + d6*betaECM + d7*betaVSM + d8*k 
  
  # Indirect Effects per 10-years age
  total_Age := 10*(d1 + a1*d6 + b1*d7 + c1*d8)
  direct_Age := 10*d1 
  indirect_ECM := 10*a1*d6
  indirect_VSM := 10*b1*d7
  indirect_k := 10*c1*d8 " 

c_fit_noHT <- cfa(c_model_noHT, data=c_med_data, se="bootstrap", bootstrap="2000")
parameterEstimates(c_fit_noHT)

################################################################################
## cfPWV Mediation Sensitivity Analysis
cf_model_BPMed <- "
  # Regression paths from co-variates to stiffness parameters
  betaECM ~ a1*Age + a2*Sex + a3*Smoke +a4* DM + a5*HT + a6*BPMed
  betaVSM ~ b1*Age + b2*Sex + b3*Smoke + b4*DM + b5*HT + b6*BPMed

  # Regression paths from stiffness parameters to PWV80
  PWV80 ~ d1*Age + d2*Sex + d3*Smoke + d4*DM + d5*HT + d6*betaECM + d7*betaVSM + d9*BPMed
  
  # Indirect Effects per 10-years age
  total_Age := 10*(d1 + a1*d6 + b1*d7 )
  direct_Age := 10*d1 
  indirect_ECM := 10*a1*d6
  indirect_VSM := 10*b1*d7 "

cf_fit_BPMed <- cfa(cf_model_BPMed, data=cf_med_data, se="bootstrap", bootstrap="2000")
parameterEstimates(cf_fit_BPMed)

cf_model_noHT <- "
  # Regression paths from co-variates to stiffness parameters
  betaECM ~ a1*Age + a2*Sex + a3*Smoke +a4* DM 
  betaVSM ~ b1*Age + b2*Sex + b3*Smoke + b4*DM 

  # Regression paths from stiffness parameters to PWV80
  PWV80 ~ d1*Age + d2*Sex + d3*Smoke + d4*DM + d6*betaECM + d7*betaVSM
  
  # Indirect Effects per 10-years age
  total_Age := 10*(d1 + a1*d6 + b1*d7 )
  direct_Age := 10*d1 
  indirect_ECM := 10*a1*d6
  indirect_VSM := 10*b1*d7 " 

cf_fit_noHT <- cfa(cf_model_noHT, data=cf_med_data, se="bootstrap", bootstrap="2000")
parameterEstimates(cf_fit_noHT)

################################################################################
## Primary Carotid Linear Models with Age, Sex, Smoke, DM, and HTN

## Preliminary Check of Residual normality
dd <- datadist(Carotid)
options(datadist='dd')

a_cECM <- ols(data=Carotid, BetaECM ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )
Carotid$resid <- r <- resid(a_cECM); Carotid$fitted <- fitted(a_cECM)
yl <- ylab('Residuals')
p1 <- ggplot(Carotid, aes(x=fitted, y=resid)) + geom_point() + yl
p2 <- ggplot(Carotid, aes(x=Age, y=resid)) + geom_point()+yl
p3 <- ggplot(Carotid, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, ncol=2)


a_cECM <- ols(data=Carotid, sqrt(BetaECM) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )
Carotid$resid <- r <- resid(a_cECM); Carotid$fitted <- fitted(a_cECM)
yl <- ylab('Residuals')
yl <- ylab('Residuals')
p1 <- ggplot(Carotid, aes(x=fitted, y=resid)) + geom_point() + yl
p2 <- ggplot(Carotid, aes(x=Age, y=resid)) + geom_point()+yl
p3 <- ggplot(Carotid, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, ncol=2)

a_cVSM <- ols(data=Carotid, betaVSM ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )
Carotid$resid <- r <- resid(a_cVSM); Carotid$fitted <- fitted(a_cVSM)
yl <- ylab('Residuals')
p1 <- ggplot(Carotid, aes(x=fitted, y=resid)) + geom_point() + yl
p2 <- ggplot(Carotid, aes(x=Age, y=resid)) + geom_point()+yl
p3 <- ggplot(Carotid, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, ncol=2)


a_cVSM <- ols(data=Carotid, sqrt(betaVSM) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )
Carotid$resid <- r <- resid(a_cVSM); Carotid$fitted <- fitted(a_cVSM)
yl <- ylab('Residuals')
yl <- ylab('Residuals')
p1 <- ggplot(Carotid, aes(x=fitted, y=resid)) + geom_point() + yl
p2 <- ggplot(Carotid, aes(x=Age, y=resid)) + geom_point()+yl
p3 <- ggplot(Carotid, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, ncol=2)


a_k <- ols(data=Carotid, k*100 ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )
Carotid$resid <- r <- resid(a_k); Carotid$fitted <- fitted(a_k)
yl <- ylab('Residuals')
p1 <- ggplot(Carotid, aes(x=fitted, y=resid)) + geom_point() + yl
p2 <- ggplot(Carotid, aes(x=Age, y=resid)) + geom_point()+yl
p3 <- ggplot(Carotid, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, ncol=2)


a_Dref <- ols(data=Carotid, Dref ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )
Carotid$resid <- r <- resid(a_Dref); Carotid$fitted <- fitted(a_Dref)
yl <- ylab('Residuals')
p1 <- ggplot(Carotid, aes(x=fitted, y=resid)) + geom_point() + yl
p2 <- ggplot(Carotid, aes(x=Age, y=resid)) + geom_point()+yl
p3 <- ggplot(Carotid, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, ncol=2)

## Beta ECM and Beta VSMC will be Sqrt() due to residual distribution deviating from normal
a_cECM <- ols(data=Carotid, scale(sqrt(BetaECM)) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )

a_cVSM <- ols(data=Carotid, scale(sqrt(betaVSM)) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )

a_k <- ols(data=Carotid, k*100 ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )

a_Dref <- ols(data=Carotid, Dref ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )

a_cECM
summary(a_cECM, Age=c(60,70) )
a_cVSM
summary(a_cVSM, Age=c(60,70) )
a_k
summary(a_k, Age=c(60,70) )
a_Dref
summary(a_Dref, Age=c(60,70) )

################################################################################
## Primary cfPWV Linear Models with Age, Sex, Smoke, DM, and HTN

## Preliminary Check of Residual normality
dd <- datadist(cfPWV)
options(datadist='dd')

a_cf_ECM <- ols(data=cfPWV, Beta.ECM..DBP. ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )
cfPWV$resid <- r <- resid(a_cf_ECM); cfPWV$fitted <- fitted(a_cf_ECM)
yl <- ylab('Residuals')
p1 <- ggplot(cfPWV, aes(x=fitted, y=resid)) + geom_point() + yl
p2 <- ggplot(cfPWV, aes(x=Age, y=resid)) + geom_point()+yl
p3 <- ggplot(cfPWV, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, ncol=2)


a_cf_ECM <- ols(data=cfPWV, sqrt(Beta.ECM..DBP.) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )
cfPWV$resid <- r <- resid(a_cf_ECM); cfPWV$fitted <- fitted(a_cf_ECM)
yl <- ylab('Residuals')
p1 <- ggplot(cfPWV, aes(x=fitted, y=resid)) + geom_point() + yl
p2 <- ggplot(cfPWV, aes(x=Age, y=resid)) + geom_point()+yl
p3 <- ggplot(cfPWV, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, ncol=2)

a_cf_VSM <- ols(data=cfPWV, Beta.VSMC..DBP. ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )
cfPWV$resid <- r <- resid(a_cf_VSM); cfPWV$fitted <- fitted(a_cf_VSM)
yl <- ylab('Residuals')
p1 <- ggplot(cfPWV, aes(x=fitted, y=resid)) + geom_point() + yl
p2 <- ggplot(cfPWV, aes(x=Age, y=resid)) + geom_point()+yl
p3 <- ggplot(cfPWV, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, ncol=2)


a_cf_VSM <- ols(data=cfPWV, sqrt(Beta.VSMC..DBP.) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )
cfPWV$resid <- r <- resid(a_cf_VSM); cfPWV$fitted <- fitted(a_cf_VSM)
yl <- ylab('Residuals')
p1 <- ggplot(cfPWV, aes(x=fitted, y=resid)) + geom_point() + yl
p2 <- ggplot(cfPWV, aes(x=Age, y=resid)) + geom_point()+yl
p3 <- ggplot(cfPWV, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, ncol=2)

## Beta ECM and Beta VSMC will be Sqrt() due to residual distribution deviating from normal
a_cf_ECM <- ols(data=cfPWV, scale(sqrt(Beta.ECM..DBP.) ) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )

a_cf_VSM <- ols(data=cfPWV, scale(sqrt(Beta.VSMC..DBP.) ) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN )

a_cf_ECM
summary(a_cf_ECM, Age=c(60,70) )
a_cf_VSM
summary(a_cf_VSM, Age=c(60,70) )

################################################################################
## Sensitivity Analysis including BP Meds

dd <- datadist(Carotid)
options(datadist='dd')

a_cECM_BP <- ols(data=Carotid, scale(sqrt(BetaECM)) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN + BP.Med)
a_cVSM_BP <- ols(data=Carotid, scale(sqrt(betaVSM)) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN + BP.Med)
a_k_BP <- ols(data=Carotid, k*100 ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN + BP.Med)
a_Dref_BP <- ols(data=Carotid, Dref ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN + BP.Med)

a_cECM_BP
summary(a_cECM_BP, Age=c(60,70) )
a_cVSM_BP
summary(a_cVSM_BP, Age=c(60,70) )
a_k_BP
summary(a_k_BP, Age=c(60,70) )
a_Dref_BP
summary(a_Dref_BP, Age=c(60,70) )


dd <- datadist(cfPWV)
options(datadist='dd')

a_cf_ECM_BP <- ols(data=cfPWV, scale(sqrt(Beta.ECM..DBP.) ) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN + BP.MEd)
a_cf_VSM_BP <- ols(data=cfPWV, scale(sqrt(Beta.VSMC..DBP.) ) ~ Age + Sex + FormCurr.Smoke + dm + Has.HTN + BP.MEd)

a_cf_ECM_BP
summary(a_cf_ECM_BP, Age=c(60,70) )
a_cf_VSM_BP
summary(a_cf_VSM_BP, Age=c(60,70) )

## Sensitivity Analysis Excluding HTN

dd <- datadist(Carotid)
options(datadist='dd')

a_cECM_noHT <- ols(data=Carotid, scale(sqrt(BetaECM)) ~ Age + Sex + FormCurr.Smoke + dm )
a_cVSM_noHT <- ols(data=Carotid, scale(sqrt(betaVSM)) ~ Age + Sex + FormCurr.Smoke + dm )
a_k_noHT <- ols(data=Carotid, k*100 ~ Age + Sex + FormCurr.Smoke + dm )
a_Dref_noHT <- ols(data=Carotid, Dref ~ Age + Sex + FormCurr.Smoke + dm )

a_cECM_noHT
summary(a_cECM_noHT, Age=c(60,70) )
a_cVSM_noHT
summary(a_cVSM_noHT, Age=c(60,70) )
a_k_noHT
summary(a_k_noHT, Age=c(60,70) )
a_Dref_noHT
summary(a_Dref_noHT, Age=c(60,70) )

dd <- datadist(cfPWV)
options(datadist='dd')

a_cf_ECM_noHT <- ols(data=cfPWV, scale(sqrt(Beta.ECM..DBP.) ) ~ Age + Sex + FormCurr.Smoke + dm )
a_cf_VSM_noHT <- ols(data=cfPWV, scale(sqrt(Beta.VSMC..DBP.) ) ~ Age + Sex + FormCurr.Smoke + dm )

a_cf_ECM_noHT
summary(a_cf_ECM_noHT, Age=c(60,70) )
a_cf_VSM_noHT
summary(a_cf_VSM_noHT, Age=c(60,70) )

