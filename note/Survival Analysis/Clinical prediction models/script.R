# EXAMPLE_SYNTAX paper prof. dr. Y. van der Graaf + Prof. dr. E. Steyerberg
# DATE: 29May2009

library(Design)
library(mice)
library(foreign)

# Importation of SMART-data-file with correct variables and endpoints
# Variables are already truncated to 1th en 99th percentiles
# SMART  <- read.spss('smartdata/SMARTst.sav',use.value.labels=F, to.data.frame=T)

save(SMART, file = "SMART.RData")
# Original SMARTfile not truncated
# SMARTo  <- read.spss('smartdata/SMARTs.sav',use.value.labels=F, to.data.frame=T)
save(SMARTo, file = "SMARTo.RData")
#-------------------------------------------------------------------------------
## PAPER STEP 1: Preliminary#
# See the first row of the data-file; also to get an idea of which variables are included 
SMART[1,]

# Retrieve the dimensions of your file: (How many rows(cases) and how many colums (variables))
dim(SMART)

# Retrieve basic information per variable (e.g: mean, highest value, lowest value, missings etc.)
describe(SMART)				

# Generate an extra window for 1 plot
par(mfrow=c(1,1))

# Selecting 18 variables out of 29 in the loaded SMART-file and immediatly retrieve info again.
SMARTm <- SMART[,Cs(SYSTBP, DIASTBP, SYSTH, DIASTH,
                    CHOL,HDL,LDL,TRIG, HOMOC, GLUT, CREAT, IMT, ALBUMIN, 
                    STENOSIS, DIABETES, SMOKING, PACKYRS, ALCOHOL)]
describe(SMARTm)

# With the following function you determine the fraction of NA's (missing values)
na.patterns <- naclus(SMARTm)
# Plot the 18 variables you just selected in a graph: Fig 23.2. 
plot(na.patterns, ylab="Fraction of NAs in common")

# A different approach to show the NA's is the following:
# First you generate a screen in which two graphs can be shown: Fig 23.1. 
par(mfrow=c(2,2))
naplot(na.patterns, which=c('na per var'))
naplot(na.patterns, which=c('na per obs'))


#-------------------------------------------------------------------------------
## PAPER STEP 2: Coding of predictors. ##
# Start of IMT-variable-example just as in the paper

# Showing the maximum value of IMT in the truncated file.
max(SMART$IMT,na.rm=T)  
# Show how many times this maximum value is exceeded in the original file
sum(SMARTo$IMTO>1.83,na.rm=T)
# Same, but now for the minimum value
min(SMART$IMT,na.rm=T)  
sum(SMARTo$IMTO<.47,na.rm=T)
# show relation between IMT in both the not-truncated file and the truncated file with the cox proportional hazard function
cph(Surv(TEVENT,EVENT) ~  IMT, data=SMART)
cph(Surv(TEVENT,EVENT) ~  IMTO, data=SMARTo)

# again an extra screen
par(mfrow=c(1,2),mar=c(5,5,3,1))
# make an new matrix with the the original and truncated IMT variable. Then name the two colums
IMTdata <- as.matrix(cbind(SMARTo$IMTO, SMART$IMT))
IMTlist <- split(IMTdata, col(IMTdata))
names(IMTlist)  <- Cs(Original, Truncated)
# make a boxplot of this matirx to show the difference of the data when it's truncated 
boxplot(IMTlist, ylab= "IMT (mm)", cex.lab=1.2)  # Fig 23.4

# Now the effect on outcome is shown with the cph-function in three ways 
#to show how the variable which fit the model best: 
#1: Not truncated with restricted cubic spline with 5 knots
dd <- datadist(SMARTo)
options(datadist="dd")
IMTOfit <- cph(Surv(TEVENT,EVENT) ~  rcs(IMTO,5), data=SMARTo)
plot(Predict(IMTOfit), xlab="IMT (mm)", lty=2, lwd=2, conf.int=F)
#2: Truncated with restricted cubic spline with 5 knots 
dd <- datadist(SMART)
options(datadist="dd")
IMTfit <- cph(Surv(TEVENT,EVENT) ~  rcs(IMT,5), data=SMART)
plot(Predict(IMTfit), xlab="", add=T, lty=3, lwd=2,conf.int=F)
#3: Truncated linear
IMTfitlin <- cph(Surv(TEVENT,EVENT) ~  IMT, data=SMART)
plot(Predict(IMTfitlin), xlab="", add=T, lty=1, lwd=2, conf.int=F)
# Putting text in the graphs 
text(x=2,y=.8,"Original,\nspline",adj=0)
text(x=1.8,y=1.1,"Truncated,\nspline",adj=0)
text(x=1.4,y=.8,"Truncated,\nlinear",adj=1)
## End IMT illustration, Fig 23.4

# Start of Age-example just as in the paper
# Example transformation of age to get best fit.
# With cph-function the X2 and dof are measured when the variable in transformed.
dd <- datadist(SMART)
options(datadist="dd")
par(mfrow=c(2,2), mar=c(5,5,3,1))
par(mfrow=c(1,1))
# age related to outcome: age linear --> X2=97
AGEfit <- cph(Surv(TEVENT,EVENT) ~  AGE, data=SMART)  
plot(Predict(AGEfit), xlab="Age (years)", lty=1, lwd=2,conf.int=T, ylim=c(-1.5, 1.5), xlim=c(20,80))
anova(AGEfit)
# age related to outcome: age + age squared --> X2= 125
AGEfit2 <- cph(Surv(TEVENT,EVENT) ~  pol(AGE,2), data=SMART)  # age + age squared
anova(AGEfit2)
plot(Predict(AGEfit2), 
     xlab="Age (years)",
     lty=2, lwd=2,conf.int=T
     )
title(main = paste("Age, 1 df, ", "chi^2"," ", 97,
                       "; Age+","Age^2",", 2 df, ", "chi^2"," ", 125,sep=""), cex.main=1)


# Age related to outcome: 125truncated under 50 / 55 and then linear --> X2= 119
AGEfit3 <- cph(Surv(TEVENT,EVENT) ~  ifelse(AGE>55, (AGE-55),0), data=SMART)  # age linear
anova(AGEfit3)
plot(AGEfit3, xlab="Age (years)", lty=1, lwd=2,conf.int=T, ylim=c(-1.5, 1.5), xlim=c(20,80))
# Age related to outcome: 125truncated under 50 / 55 and then squared --> X2= 130
AGEfit4 <- cph(Surv(TEVENT,EVENT) ~  ifelse(AGE>50, (AGE-50)^2,0), data=SMART)  # age + age squared
anova(AGEfit4)
plot(AGEfit4, xlab="Age (years)", add=T, lty=2, lwd=2,conf.int=T)
title(substitute(paste("Age>55, 1 df, ", chi^2," ", 119,
                       "; ",(Age>50)^2,", 2 df, ", chi^2," ", 130,sep="")), cex.main=1)

# age related to outcome with Spline and 3 df (4 knots) --> X2=125
AGEfit5 <- cph(Surv(TEVENT,EVENT) ~  rcs(AGE,4), data=SMART)  
anovaAge  <- anova(AGEfit5)
plot(AGEfit5, xlab="Age (years)", lty=3, lwd=2,conf.int=T, 
     ylim=c(-1.5, 1.5), xlim=c(20,80))
title(substitute(paste("Spline, 3 df, ", chi^2," ", 125,sep="")), cex.main=1)

# Age related to outcome in quartiles --> X2=93
tapply(SMART$AGE, ifelse(SMART$AGE<50,1,ifelse(SMART$AGE<60,2,ifelse(SMART$AGE<70,3,4))), mean)
SMART$AGEcat  <- as.factor(ifelse(SMART$AGE<50,43,ifelse(SMART$AGE<60,55,ifelse(SMART$AGE<70,65,74))))
AGEfit6 <- cph(Surv(TEVENT,EVENT) ~  as.factor(AGEcat), data=SMART)
anova (AGEfit6)
plot(AGEfit6, xlab="Age (years)", lty=3, lwd=2,conf.int=T, 
     ylim=c(-1.5, 1.5))
title(substitute(paste("Categorized, 3 df, ", chi^2," ", 93,sep="")), cex.main=1)

# Age related to outcome dichotomization --> X2=72
tapply(SMART$AGE, cut2(SMART$AGE,g=2), mean)
AGEfit7 <- cph(Surv(TEVENT,EVENT) ~  ifelse(SMART$AGE<61,51,69), data=SMART) 
anova(AGEfit7)
##End AGE-illustration  Fig 23.5 and Tab 23.4

# Transformations of other continuous vars: RCS with 3 or 4 knots and normal linear is checked
# BMI
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(BMI,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(BMI,3), data=SMART)  
anova(fit) 
# Blood pressure: diastolic
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(DIASTH,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(DIASTH,3), data=SMART)  
anova(fit) 
fit <- cph(Surv(TEVENT,EVENT) ~  DIASTH, data=SMART)  
anova(fit) 
# Blood pressure: Systolic
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(SYSTH,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(SYSTH,3), data=SMART)  
anova(fit) 
fit <- cph(Surv(TEVENT,EVENT) ~  SYSTH, data=SMART)  
anova(fit) 
# Blood pressure: other type of measurement: diastolic
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(DIASTBP,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(DIASTBP,3), data=SMART)  
anova(fit) 
plot(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  DIASTBP, data=SMART)  
anova(fit) 
# Blood pressure: other typeof measurment: systolic
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(SYSTBP,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(SYSTBP,3), data=SMART)  
anova(fit) 
fit <- cph(Surv(TEVENT,EVENT) ~  SYSTBP, data=SMART)  
anova(fit) 
# Lipids: Cholesterol  
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(CHOL,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(CHOL,3), data=SMART)  
anova(fit) 
fit <- cph(Surv(TEVENT,EVENT) ~  CHOL, data=SMART)  
anova(fit) 
# Lipids: HDL 
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(HDL,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(HDL,3), data=SMART)  
anova(fit) 
fit <- cph(Surv(TEVENT,EVENT) ~  HDL, data=SMART)  
anova(fit) 
# Lipids: LDL
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(LDL,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(LDL,3), data=SMART)  
anova(fit) 
fit <- cph(Surv(TEVENT,EVENT) ~  LDL, data=SMART)  
anova(fit) 
# Lipids: Triglyceriden 
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(TRIG,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(TRIG,3), data=SMART)  
anova(fit) 
fit <- cph(Surv(TEVENT,EVENT) ~  TRIG, data=SMART)  
anova(fit) 
# Only HDL seems to have an effect
# Other: Homocysteine 
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(HOMOC,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(HOMOC,3), data=SMART)  
anova(fit) 
fit <- cph(Surv(TEVENT,EVENT) ~  HOMOC, data=SMART)  
anova(fit) 
# Other: glucose 
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(GLUT,4), data=SMART)  
anova(fit)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(GLUT,3), data=SMART)  
anova(fit) 
fit <- cph(Surv(TEVENT,EVENT) ~  GLUT, data=SMART)  
anova(fit) 
# Other: Diabetes  
fit <- cph(Surv(TEVENT,EVENT) ~  DIABETES, data=SMART)  
anova(fit)
# Other: Diabetes + Glucose  
fit <- cph(Surv(TEVENT,EVENT) ~  DIABETES + GLUT, data=SMART)  
anova(fit) 
#Best is to keep diabetes, no glucose

#Other: Creatinine: Also visualiised  in plots, because different transformation is done (logistic)
#Plot 4 graphs for CREAT: Fig 23.6
par(mfrow=c(2,2), mar=c(5,5,3,1))
fit <- cph(Surv(TEVENT,EVENT) ~  CREAT, data=SMART)  
anova(fit) 
plot(fit, xlab="Creatinin (mmol/l)", lty=1, lwd=2,conf.int=T, ylim=c(-.5, 2.5))
title(substitute(paste("Linear, 1 df, ", chi^2," ", 93,sep="")), cex.main=1.2)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(CREAT,4), data=SMART)  
anova(fit)
plot(fit, xlab="Creatinin (mmol/l)", lty=3, lwd=2,conf.int=T, ylim=c(-.5, 2.5))
title(substitute(paste("Spline, 3 df, ", chi^2," ", 116,sep="")), cex.main=1.2)
fit <- cph(Surv(TEVENT,EVENT) ~  rcs(CREAT,3), data=SMART)  
anova(fit) 
plot(fit, xlab="Creatinin (mmol/l)", lty=2, lwd=2,conf.int=T, ylim=c(-.5, 2.5))
title(substitute(paste("Spline, 2 df, ", chi^2," ", 99,sep="")), cex.main=1.2)
fit <- cph(Surv(TEVENT,EVENT) ~  log(CREAT), data=SMART)  
anova(fit) 
plot(fit, xlab="Creatinin (mmol/l)", lty=1, lwd=2,conf.int=T, ylim=c(-.5, 2.5))
title(substitute(paste("Log transformation, 1 df, ", chi^2," ", 131,sep="")), cex.main=1.2)
# Hence: log(CREAT) good choice, Fig 23.6 and Tab 23.4 ##

# Start of History-example just as in the paper
# Score with individual terms, 4 dof
fitsum  <- cph(Surv(TEVENT,EVENT) ~  CEREBRAL+CARDIAC+AAA+PERIPH, data=SMART)
anova (fitsum)
# Score with each affected organ seperately in a score, assuming equal weight 
SMART$HISTCARD  <- SMART$CEREBRAL+SMART$CARDIAC+SMART$AAA+SMART$PERIPH  
fitsum2  <- cph(Surv(TEVENT,EVENT) ~  HISTCARD, data=SMART)
anova(fitsum2)
# Score with each affected organ seperately in a score, AAA doubled 
SMART$HISTCAR2  <- SMART$CEREBRAL+SMART$CARDIAC+2*SMART$AAA+SMART$PERIPH  
describe(SMART$HISTCAR2)
fitsum3  <- cph(Surv(TEVENT,EVENT) ~  HISTCAR2, data=SMART) 
anova(fitsum3)
## End history of cardiovascular disease ##

#-------------------------------------------------------------------------------
## PAPER STEP 3 + 4: Model Specification AND model estimation##
# imputation model with outcome, all potential predictors, and auxilliary variables
SMARTM  <-aregImpute(~I(TEVENT)+EVENT+SEX+I(AGE)+
                       SYSTBP+DIASTBP+SYSTH+DIASTH+
                       DIABETES+CEREBRAL+CARDIAC+AAA+PERIPH+I(HISTCARD)+STENOSIS+
                       I(LENGTH)+I(WEIGHT)+I(BMI)+
                       I(CHOL)+I(HDL)+I(LDL)+I(TRIG)+I(HOMOC)+I(GLUT)+I(CREAT)+I(IMT)+
                       as.factor(ALBUMIN) +as.factor(SMOKING) + I(PACKYRS) + as.factor(ALCOHOL), 
                     n.impute=5, data=SMART)

# Check results of imputation, Fig 23.3
SMARTM
# Use transcan to make 1 singly imputed file.
SMARTc <- SMART
imputed <-impute.transcan(SMARTM, imputation=1, data=SMART, list.out=TRUE, pr=FALSE, check=FALSE)
SMARTc[names(imputed)] <- imputed
SMARTc$HISTCAR2  <- SMARTc$CEREBRAL+SMARTc$CARDIAC+2*SMARTc$AAA+SMARTc$PERIPH  
SMARTc$HISTCAR3  <- ifelse(SMARTc$HISTCARD==4, 3, SMARTc$HISTCARD)         
SMARTc[1:4,]

#Check and investigate your new file, no missings anymore!!!
describe(SMARTc)

# 3 models: full model, LASSO-model, Backward-step model
#1: Full model combined
dd <- datadist(SMARTc)
options(datadist="dd") 
cox01 <-cph(Surv(TEVENT,EVENT) ~  ifelse(AGE>50, (AGE-50)^2,0) + SEX + 
              as.factor(SMOKING) + as.factor(ALCOHOL)+BMI + SYSTH +  
              HDL+DIABETES +HISTCAR2 +HOMOC + log(CREAT)+ as.factor(ALBUMIN)+STENOSIS+IMT,
            data = SMARTc, x=T, y=T, surv=T)
anova(cox01)
summary(cox01)
cox01
# See Table 23.6: slightly different numbers because of new impuation      


# Interactions with sex, 10 df, chi2 14.8
anova(cph(Surv(TEVENT,EVENT) ~  SEX * (ifelse(AGE>55, (AGE-55),0) + 
                                         BMI + HDL+DIABETES +HISTCAR2 + log(CREAT)+ as.factor(ALBUMIN)+STENOSIS+IMT),data = SMARTc))
# specific sex*HISTCAR2, 1 df
anova(cph(Surv(TEVENT,EVENT) ~  ifelse(AGE>55, (AGE-55),0) + 
            BMI + HDL+DIABETES +SEX * HISTCAR2 + log(CREAT)+ as.factor(ALBUMIN)+STENOSIS+IMT,data = SMARTc))
# HISTCAR2 * other predictors, 9 df
anova(cph(Surv(TEVENT,EVENT) ~  SEX + HISTCAR2 * (ifelse(AGE>55, (AGE-55),0) + 
                                                    BMI + HDL+DIABETES + log(CREAT)+ as.factor(ALBUMIN)+STENOSIS+IMT),data = SMARTc))
# specific sex* +CEREBRAL+CARDIAC+AAA+PERIPH
anova(cph(Surv(TEVENT,EVENT) ~  ifelse(AGE>55, (AGE-55),0) + 
            BMI + HDL+DIABETES +SEX*(CEREBRAL+CARDIAC+AAA+PERIPH) + log(CREAT)+ as.factor(ALBUMIN)+STENOSIS+IMT,data = SMARTc))
# Plot interaction history * sex
par(mfrow=c(1,2))
# levels(SMARTc$SEX) <- Cs(Male, Female)
dd <- datadist(SMARTc)
options(datadist="dd")
plot(cph(Surv(TEVENT,EVENT) ~  ifelse(AGE>55, (AGE-55),0) + 
           BMI + HDL+DIABETES +SEX*HISTCAR2 +log(CREAT)+ as.factor(ALBUMIN)+STENOSIS+IMT,data = SMARTc),
     HISTCAR2=NA, SEX=NA, xlim=c(1,5.8), xlab="Sumscore previous symptomatic atherosclerosis")
title('sex*sumscore')
text(x=c(5.2,5.2), y=c(0.93,2.56), Cs(Male, Female), adj=0)

# C-statistic of the full model:
rcorr.cens(-cox01$linear.pred,cox01$y) 


#2: Stepwise selection
backwardstepwise <- fastbw(cox01, rule="aic", type="individual")
backwardstepwise
cox01.step  <- update(cox01, ~ cox01$x[,backwardstepwise$parms.kept])
cox01.step

#3: LASSO-method
## glmpath analysis; 'lasso'
library(glmpath)
memory.size(2048)

## SMART data analysis
SMARTdata <- list(x=cox01$x, time=SMARTc$TEVENT, status=SMARTc$EVENT)
summary(SMARTdata)

## VERY LONG computing time ...
# predictie corrector methode voor berekenen van CPH-model met L1 penalty
SMARTpath <- coxpath(data=SMARTdata)

par(mar=c(5,4,4,6))
plot.coxpath(SMARTpath, cex.lab=1.2)
plot.coxpath(SMARTpath, cex.lab=1.2,type = "aic")
SMARTpath
SMARTpath$aic
round(SMARTpath$b.predictor[SMARTpath$aic==min(SMARTpath$aic),],4)

####################CHECK SMARTHPATH AND CHOOSE ROW IN WHICH HOMOCYST IS DELETED FROM THE MODEL############################ 
SMARTpath$b.predictor
# model ?? is without homocyst
Lasso.coefs <- SMARTpath$b.predictor[#CHOOSE ROW,]
Lasso.coefs <- Lasso.coefs[Lasso.coefs > 0 | Lasso.coefs < 0]
Lasso.coefs

# Alternsative: Goeman's library
library(penalized)
memory.size(2048)

## new SMART data analysis
opt <- optL1(Surv(SMARTc$TEVENT, SMARTc$EVENT), penalized = cox01$x, startbeta=cox01$coefficients, fold = 10, model='cox', maxlambda1=20, standardize=T)
coefficients(opt$fullfit)
cox01.lasso <- penalized(Surv(SMARTc$TEVENT, SMARTc$EVENT) ~ cox01$x, model='cox', standardize=T, lambda1=18)
coef(cox01.lasso)
cox01.lasso.path <- penalized(Surv(SMARTc$TEVENT, SMARTc$EVENT) ~ cox01$x, model='cox', standardize=T, lambda1=1, steps=10)
plotpath(cox01.lasso.path, standardize=T)
title('Standardized coefficients')


#-------------------------------------------------------------------------------
## PAPER STEP 5: Model performance ##

#COX01DEF:
cox01def <-cph(Surv(TEVENT,EVENT) ~  ifelse(AGE>50, (AGE-50)^2,0) +  BMI +
                 HDL+DIABETES +HISTCAR2 + log(CREAT)+ as.factor(ALBUMIN)+STENOSIS+IMT,data = SMARTc, x=T, y=T, surv=T)
dd <- datadist(SMARTc)
options(datadist="dd")
anova(cox01def)
summary(cox01def)

#C-statistic and KM-curves of cox01def
#C-statistic:
rcorr.cens(-cox01def$linear.pred,cox01def$y) 
#KM-curves:
validate(cox01def, B=20, dxy=T)

cox01def$linear.pred
summary(cox01def$linear.predictor)
SMARTc$lp4	<- cut2(as.numeric(cox01def$linear.predictor), g=4)
dd <- datadist(SMARTc)
options(datadist="dd")

table(SMARTc$lp4)
levels(SMARTc$lp4) <- as.character(1:4)

KMrisk4	<- survfit(Surv(TEVENT,EVENT) ~ lp4, data=SMARTc)
KMrisk4      # Fig 23.8

par(las=1)
survplot(KMrisk4, conf="none", n.risk=T, xlab="Follow-up (days)", cex.n.risk=1, adj=.5,
         ylab="Fraction free of cardiovascular event", time.inc=365)

#-------------------------------------------------------------------------------
## PAPER STEP 6: Model performance ##
#Bootstrap:
val.step  <- validate(cox01, bw=T, rule="aic", type="individual")
vcox01def <- validate.cph(cox01def, method='boot', B=200, bw=T,type="individual",dxy=T)
coef.sh 	<- cox01def$coefficients*vcox01def["Slope","index.corrected"]
lp.s		<-cox01def$x%*%coef.sh
rcorr.cens (-lp.s,Surv(SMARTc$TEVENT,SMARTc$EVENT))

# SHRINKAGE FACTORS 
Lasso.coefs <- SMARTpath$b.predictor[SMARTpath$aic==min(SMARTpath$aic),]
Lasso.coefs <- Lasso.coefs[SMARTpath$b.predictor[SMARTpath$aic==min(SMARTpath$aic),] > 0 | 
                             SMARTpath$b.predictor[SMARTpath$aic==min(SMARTpath$aic),] < 0]
Lasso.coefs
cox01def$coef
round(Lasso.coefs/cox01def$coef,3)
mean(Lasso.coefs/cox01def$coef)        # 0.80; stronger shrinkage for weaker effects
# v Houwelingen shrinkage factor, with 17 df --> 0.94
(cox01def$stats[3] -  cox01$stats[4]) / cox01def$stats[3]


#-------------------------------------------------------------------------------
## PAPER STEP 7: Model performance ##

# Nomogram with lasso coefs as in the original paper ( EERST cox01def draaien, daarvoor moet worden geimputeerd).
lassocoefs <- c(AGE=0.0012,BMI=-0.001,HDL=-0.16,DIABETES=0.11,HISTCAR2=0.33,CREAT=0.71,'ALBUMIN=2'=0.13,'ALBUMIN=3'=0.20,STENOSIS=0.16,IMT=0.50)
cox01def$coefficients <- lassocoefs
surv	<- Survival(cox01def)
surv3	<- function(x) surv(3*365.25,lp=x)
surv5	<- function(x) surv(5*365.25,lp=x)
nomogram(cox01def, fun=list(surv3, surv5), lp=F, AGE=c(15,seq(50,85,5)), CREAT=c(60,80,100,150,200,400),
         funlabel=c('3-year survival', '5-year survival'), maxscale=10, lp.at=c(0,1,2,2.4), fun.at=c(.93,.9,.85,.8,.7,.6,.5))


# nomogram ( if you run this before you run the above, you get the nomogram based on the coefficients of the cph (cox01def) fit) 
surv	<- Survival(cox01def)
surv1	<- function(x) surv(1*365.25,lp=x)
surv3	<- function(x) surv(3*365.25,lp=x)
surv5	<- function(x) surv(5*365.25,lp=x)

par(mfrow=c(1,1))
nomogram(cox01def, fun=list(surv3,surv5), lp=F, AGE=c(15,seq(50,85,5)), CREAT=c(60,80,100,150,200,400), funlabel=c('3-year survival', '5-year survival'), maxscale=10,lp.at=c(0,1,2,2.4), fun.at=c(.93,.9,.85,.8,.7,.6,.5))

#-------------------------------------------------------------------------------
## EXTRA
# how I computed the absolute risks per patient

#make baseline database
baseline <- data.frame(AGE=0,BMI=0,HDL=0,DIABETES=0,HISTCAR2=0,CREAT=1,ALBUMIN=1,STENOSIS=0,IMT=0)
summary(survfit(cox01def, newdata=baseline))
plot(survfit(cox01def, newdata=baseline))

#  make new variables in which tranformed variables are transformed.                   
SMARTc_new<-SMARTc
#New columns in SMARTc, columns 35-38:
SMARTc_new$ALBUMIN2<-as.numeric(SMARTc$ALBUMIN==2)
SMARTc_new$ALBUMIN3<-as.numeric(SMARTc$ALBUMIN==3)
SMARTc_new$logCREAT<-log(SMARTc$CREAT)
SMARTc_new$newAGE<-ifelse(SMARTc$AGE>50, (SMARTc$AGE-50)^2,0)

##IMPORTANT: NUMBERS ARE THE NUMBERS OF THE COLOMNS IN ORDER OF THE COEFFICIENTS### THIS CAN CHANGE PER DATASET
SMARTc_new[1,]
cox01def$coefficients
# you can give  the coefficients every value you like it to have : see above ...$coefficients<- .....
linpred<-as.matrix(SMARTc_NEW[,c(38,17,19,5,32,37,35,36,10,25)])%*%cox01def$coefficients # depends wheterh you have overwritten the coeffcients if lasso-beta's are used or not  

#compute risks

base1yr<-min(survfit(cox01def, newdata=baseline)$surv[survfit(cox01def, newdata=baseline)$time<=365.25])
base3yr<-min(survfit(cox01def, newdata=baseline)$surv[survfit(cox01def, newdata=baseline)$time<=1095.75])
base5yr<-min(survfit(cox01def, newdata=baseline)$surv[survfit(cox01def, newdata=baseline)$time<=1826.25])

SMARTc$oneyrsurv <- base1yr^(exp(linpred))
SMARTc$threeyrsurv <- base3yr^(exp(linpred))
SMARTc$fiveyrsurv <- base5yr^(exp(linpred))
SMARTc$linpred <- linpred



## SAVE IMAGE ##
save(list = ls(all=TRUE), file = "#choose name#.RData") 