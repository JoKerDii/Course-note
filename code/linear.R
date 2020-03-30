COPD <- read.csv(file = "survival analysis/COPD_student_dataset.csv", header= TRUE, sep = ",")
head(COPD)
FEV1 <- COPD[,"FEV1"]
SGRQ <- COPD[,"SGRQ"]
SGRQ_FEV1 <- lm(SGRQ~FEV1, data = COPD)
summary(SGRQ_FEV1)

colnames(COPD)
hist(COPD$MWT1Best)
hist(COPD$MWT1Best, main="Histogram of MWT1Best", xlab="MWT1Best", breaks=12) 
subset(COPD, MWT1Best > 650)
subset(COPD, MWT1Best > 600 | MWT1Best < 150) 
hist(COPD$FEV1, main="Histogram of FEV1", xlab="FEV1") 
list("Summary" = summary(COPD$MWT1Best), "Mean" = mean(COPD$MWT1Best, na.rm=TRUE), "Standard Deviation" = sd(COPD$MWT1Best, na.rm=TRUE), "Range" = range(COPD$MWT1Best, na.rm=TRUE), "Inter-Quartile Range" = IQR(COPD$MWT1Best, na.rm=TRUE))
FEV1
list("Summary" = summary(COPD$FEV1), "Mean" = mean(COPD$FEV1, na.rm=TRUE), "Standard Deviation" = sd(COPD$FEV1, na.rm=TRUE), "Range" = range(COPD$FEV1, na.rm=TRUE), "Inter-Quartile Range" = IQR(COPD$FEV1, na.rm=TRUE))
plot(COPD$FEV1, COPD$MWT1Best, xlab = "FEV1", ylab = "MWT1Best") 

cor.test(COPD$FEV1,COPD$MWT1Best, use = "complete.obs", method = "pearson")
cor.test(COPD$FEV1,COPD$MWT1Best, use = "complete.obs", method = "spearman")

hist(COPD$AGE, main="Histogram of AGE", xlab="AGE") 
list("Summary" = summary(COPD$AGE), "Mean" = mean(COPD$AGE, na.rm=TRUE), "Standard Deviation" = sd(COPD$AGE, na.rm=TRUE), "Range" = range(COPD$AGE, na.rm=TRUE), "Inter-Quartile Range" = IQR(COPD$AGE, na.rm=TRUE))
list("Summary" = summary(COPD$MWT1Best), "Mean" = mean(COPD$MWT1Best, na.rm=TRUE), "Standard Deviation" = sd(COPD$MWT1Best, na.rm=TRUE), "Range" = range(COPD$MWT1Best, na.rm=TRUE), "Inter-Quartile Range" = IQR(COPD$MWT1Best, na.rm=TRUE))
plot(COPD$AGE, COPD$MWT1Best, xlab ="AGE", ylab ="MWT1Best")
cor.test(COPD$AGE, COPD$MWT1Best, use="complete.obs", method="pearson")
cor.test(COPD$AGE, COPD$MWT1Best, use="complete.obs", method="spearman")

MWT1Best_FEV1 <- lm(MWT1Best~FEV1, data = COPD)
summary(MWT1Best_FEV1)
confint(MWT1Best_FEV1)
par(mfrow=c(2,2))
plot(MWT1Best_FEV1) 
par(mfrow=c(1,1))
plot(MWT1Best_FEV1) 

MWT1Best_AGE <- lm(MWT1Best~AGE, data = COPD)
summary(MWT1Best_AGE)
confint(MWT1Best_AGE)

predictedVals <- predict(MWT1Best_AGE) # get predicted value for the model
residualVals <- residuals(MWT1Best_AGE) # get residuals between model and data
par(mfrow=c(2,2)) # set plotting format
plot(MWT1Best_AGE)# see residual plot

hist(residualVals, main = "Histogram of residuals", xlab = "Residuals") 


MWT1Best_FEV1_AGE <- lm(MWT1Best~FEV1+AGE, data = COPD)
summary(MWT1Best_FEV1_AGE) 
confint(MWT1Best_FEV1_AGE)

lr1 <- lm(MWT1Best~FVC, data = COPD) # Run the regression, assigning the output to a new variable lr
summary(lr1) # View the output of the regression 
confint(lr1) # View the 95% confidence intervals of the regression 

lr2 <- lm(MWT1Best~AGE, data = COPD)
summary(lr2)
confint(lr2)

lr3 <- lm(MWT1Best~FVC+AGE, data = COPD) 
summary(lr3) 
confint(lr3) 

plot(COPD$AGE, COPD$FVC, xlab ="AGE", ylab ="FVC")
cor.test(COPD$AGE, COPD$FVC, use="complete.obs", method="spearman")  

############# w3
dim(COPD) 
head(COPD) 
class(COPD$AGE) 
summary(COPD$AGE)
hist(COPD$AGE) 
class(COPD$CAT) 
summary(COPD$CAT) 
hist(COPD$CAT) 
class(COPD$COPDSEVERITY) 
table(COPD$COPDSEVERITY, exclude = NULL) 
class(COPD$gender) 
COPD$gender <- as.factor(COPD$gender) 
class(COPD$gender) 
table(COPD$gender, exclude = NULL) 
class(COPD$MWT1Best)
class(COPD$copd)


library(Hmisc)
describe(COPD)
install.packages("gmodels")
library(gmodels)
CrossTable(COPD$copd)
sum(is.na(COPD$copd))
colnames(COPD)
summary(COPD$MWT1Best)
hist(COPD$AGE) 
str(COPD)

mydata <- COPD[,c("AGE", "PackHistory", "FEV1", "FEV1PRED", "FVC", "CAT", "HAD", "SGRQ")]
cor_matrix <- cor(mydata)
round(cor_matrix,2)

pairs(~AGE+PackHistory+FEV1+FEV1PRED+FVC+CAT+HAD+SGRQ,data = COPD)# correlation plot
CrossTable(COPD$hypertension,COPD$IHD)

mlr1 <- lm(MWT1Best~FEV1 + AGE + factor(gender) + factor(COPDSEVERITY) + factor(CAT),data = COPD)
install.packages("mctest")
library(mctest)
imcdiag(model.matrix(mlr1)[,-1],mlr1$model[1], model = "VIF")

COPD$Diabetes <- c(0,1)[as.integer(COPD$Diabetes)]
COPD$AtrialFib <- c(0,1)[as.integer(COPD$AtrialFib)]
DAF <- COPD$Diabetes * COPD$AtrialFib
r1 <- lm(MWT1Best~factor(Diabetes)+factor(AtrialFib)+factor(DAF), data=COPD) 
summary(r1) 
confint(r1)
r2 <- lm(MWT1Best~factor(Diabetes)+factor(AtrialFib)+factor(Diabetes*AtrialFib), data=COPD) 

install.packages("prediction")
library(prediction)
list("Diabetes" = prediction(r2, at = list(Diabetes = c(0,1))),
     "AtrialFib" = prediction(r2, at = list(Diabetes = c(0,1))),
     "Diabetes*AtrialFib" = prediction(r2, at = list(Diabetes = c(0,1))))

r3 <- lm(MWT1Best ~ factor(Diabetes) + factor(IHD), data = COPD )
summary(r3)
confint(r3)
r4 <- lm(MWT1Best ~ factor(Diabetes) + factor(IHD)+ factor(Diabetes*IHD), data = COPD )
summary(r4)
confint(r4)
