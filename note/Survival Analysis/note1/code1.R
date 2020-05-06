### Survival Analysis in R ###

## Import prostate_cancer as prost

# Providing column names

colnames(prost) = c("patient", "treatment", "time", "status", "age", "sh", "size", "index")

head(prost)



### The fundamental Surv function

library(survival)

Surv(prost$time, prost$status)



# Simple Kaplan Meier plot

mykm1 <- survfit(Surv(time, status) ~ 1, data = prost)

plot(mykm1, mark.time = TRUE)



# Stratified KM plot

mykm2 <- survfit(Surv(time, status) ~ treatment, data = prost) ## alternative 2 line chart

plot(mykm2, mark.time = TRUE)



# Advanced plot

library(ggfortify)

autoplot(mykm2)



# Log Rank Test

survival::survdiff(survival::Surv(prost$time, prost$status) ~ prost$treatment)



# Exercise Kaplan Meier

library(SMPracticals)

head(motorette)

plot(motorette$x, motorette$y)



library(survival)

mykm1 <- survfit(Surv(y, cens) ~ 1, data = motorette)

plot(mykm1, mark.time = TRUE)

temp_binary = ifelse(motorette$x < 180, "low", "high"); temp_binary

motorette$temp_binary = temp_binary

survival::survdiff(survival::Surv(motorette$y, motorette$cens) ~ motorette$temp_binary)

mykm2 <- survfit(Surv(y, cens) ~ temp_binary, data = motorette)

plot(mykm2, mark.time = TRUE)



### Cox Proportional Hazards Model

cox <- coxph(Surv(time, status) ~ treatment + age + sh+ size + index,
             
             data = prost)

summary(cox)



# Getting a plot of the model

autoplot(survfit(cox))



# alternative

coxfit = survfit(cox)

autoplot(coxfit)



# Aalens additive regression

aalen <- aareg(Surv(time, status) ~ treatment + age + sh+ size + index,
               
               data = prost)

autoplot(aalen)



# Parametric models

library(flexsurv)

weib_mod = flexsurvreg(Surv(time, status) ~ treatment + age + sh+ size + index,
                       
                       data = prost, dist = "weibull")

exp_mod = flexsurvreg(Surv(time, status) ~ treatment + age + sh+ size + index,
                      
                      data = prost, dist = "exp")

plot(weib_mod)

lines(exp_mod, col="blue")



### Exercise Cox Prop Hazards

# Cox Proportional Hazards

cox <- coxph(Surv(futime, status) ~ trt + stage + bili + riskscore,
             
             data = udca1)

summary(cox)



# Getting a plot of the model

autoplot(survfit(cox))

# Standard Kaplan Meier plot to compare the treatments

kmfit <- survfit(Surv(futime, status) ~ trt, data = udca1)

autoplot(kmfit)

# Can we improve the concordance by eliminating non significant covariates?

cox2 <- coxph(Surv(futime, status) ~ trt,
              
              data = udca1)

summary(cox2)





### Survival Trees

library(ranger)

streefit = ranger(Surv(time, status) ~ treatment + age + sh + size + index,
                  
                  data = prost,
                  
                  importance = "permutation",
                  
                  splitrule = "extratrees",
                  
                  seed = 43)



# Average the survival models

death_times = streefit$unique.death.times

surv_prob = data.frame(streefit$survival)

avg_prob = sapply(surv_prob,mean)



# Plot the survival tree model - average

plot(death_times, avg_prob,
     
     type = "s",
     
     ylim = c(0,1),
     
     col = "red", lwd = 2, bty = "n",
     
     ylab = "Survival Probability", xlab = "Time in Months",
     
     main = "Survival Tree Model\nAverage Survival Curve")



## Comparing the main models we discussed so far

# Set up for ggplot

kmdf <- data.frame(mykm1$time,
                   
                   mykm1$surv,
                   
                   KM = rep("KM", 36))

names(kmdf) <- c("Time in Months","Survival Probability","Model")

coxdf <- data.frame(coxfit$time,
                    
                    coxfit$surv,
                    
                    COX = rep("COX", 36))

names(coxdf) <- c("Time in Months","Survival Probability","Model")

rfdf <- data.frame(death_times,
                   
                   avg_prob,
                   
                   RF = rep("RF", 36))

names(rfdf) <- c("Time in Months","Survival Probability","Model")

plotdf <- rbind(kmdf,coxdf,rfdf)

p <- ggplot(plotdf, aes(x = plotdf[,1] ,
                        
                        y = plotdf[,2] ,
                        
                        color = plotdf[,3]))

p + geom_step(size = 2) +
  
  labs(x = "Time in Months",
       
       y = "Survival Probability",
       
       color = "",
       
       title = "Comparison Plot of\n3 Main Survival Models") +
  
  theme(plot.title = element_text(hjust = 0.5))



## Exercise Missing Data with mgus

summary(mgus)

str(mgus)

# Removing columns

mgusnew = mgus[,-c(5,6)]

# Statistical test for outliers

library(outliers)

grubbs.test(mgusnew$creat, type = 10)

boxplot(mgusnew$creat)

# Identify and replace all values above threshold

x = which(mgusnew$creat > 2.5)

mgusnew$creat[x] = NA

# Check the column maximum

summary(mgusnew)

# NA replacement with mice

library(mice)

mymice <- mice(mgusnew, m = 3, method = "norm")

mymicecomplete <- complete(mymice, 3)

summary(mymicecomplete)
