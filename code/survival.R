# NOTE
############################ w1
g <- read.csv(file = "survival analysis/simulated HF mort data for GMPH (1K) final.csv", header=TRUE, sep=',')
head(g)
dim(g)

install.packages("survival")
install.packages("ggplot")

library(survival) # this is the cornerstone command for survival analysis in R
library(ggplot2) # newer package that does nice plots

gender <- as.factor(g[,"gender"]) # R calls categorical variables factors
fu_time <- g[,"fu_time"] # continuous variable (numeric) 
death <- g[,"death"] # binary variable (numeric) 
age <- g[,"age"]

km_fit <- survfit(Surv(fu_time, death) ~ 1)
plot(km_fit)

summary(km_fit, times = c(1:7,30,60,90*(1:10))) 

km_gender_fit <- survfit(Surv(fu_time, death) ~ gender) 
plot(km_gender_fit)

survdiff(Surv(fu_time, death) ~ gender, rho=0) 


age_65plus <- ifelse(g[,"age"]>=65,1,0) # dichotomise age
table(age_65plus, exclude = NULL) # inspect the numbers - always a good idea
table(age,age_65plus, exclude = NULL) # check - an even better idea...

survdiff(Surv(fu_time, death) ~ age_65plus, rho=0)



############################ w2
cox <- coxph(Surv(fu_time,death) ~ age,data = g)
summary(cox)

install.packages("survminer")
library(survminer)

cox <- coxph(Surv(fu_time, death) ~ ethnicgroup, data = g) # take variables straight from g
summary(cox) # which is wrong as ethnicgroup is not continuous, try this instead:
ethnicgroup <- factor(g[,"ethnicgroup"]) # can also use “as.factor” rather than “factor”
fu_time <- g[,"fu_time"]
death <- g[,"death"]
cox <- coxph(Surv(fu_time, death) ~ ethnicgroup)
summary(cox)

levels(ethnicgroup)<-c(levels(ethnicgroup),"8") # add level 8 to the factor
ethnicgroup[is.na(ethnicgroup)] <- "8" # Change NA to "None"

cox <- coxph(Surv(fu_time, death) ~ ethnicgroup) 
summary(cox) 

############################ w3
summary(age)

t <- table(gender, exclude=NULL)
addmargins(t) # adds the total (a "sum" column)
round(100*prop.table(t),digits=1) # get %s rounded to 1dp

copd <- g[,"copd"]
t <- table(copd, exclude=NULL)
addmargins(t) # adds the total (a "sum" column)
round(100*prop.table(t),digits=1) # get %s rounded to 1dp
copd

prior_dnas <- g[,"prior_dnas"]
t <- table(prior_dnas, exclude=NULL) 
addmargins(t) # adds the total (a "sum" column) 
round(100*prop.table(t),digits=1) # get %s rounded to 1dp 

t <- table(ethnicgroup, exclude=NULL) 
addmargins(t) # adds the total (a "sum" column) 
round(100*prop.table(t),digits=1) # get %s rounded to 1dp 

cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + prior_dnas + 
               ethnicgroup)
summary(cox)

quintile <- as.factor(g[,'quintile'])
cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + quintile + ethnicgroup) 
summary(cox) 

table(quintile, exclude=NULL) 
t <- table(quintile,death) 
t
round(100*prop.table(t,1),digits=1) # row %s 

# Change the reference category
quintile <- as.factor(g[,'quintile'])
quintile <- relevel(quintile, ref = 2) # quintile 1 as the ref cat again
cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + quintile + ethnicgroup)
summary(cox)

# Combine categories
quintile_5groups <- g[,'quintile'] # best start with the original data set, not from “quintile” 
quintile_5groups[quintile_5groups==0] <- 5 
quintile_5groups <- factor(quintile_5groups) 
table(quintile_5groups, exclude=NULL) 

cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + quintile_5groups + 
               ethnicgroup) 
summary(cox) 

# Drop the quintile zero patients
quintile_5groups <- g[,'quintile'] 
quintile_5groups[quintile_5groups==0] <- NA # set the zeroes to missing 
quintile_5groups <- factor(quintile_5groups) 
table(quintile_5groups, exclude=NULL) 

cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + quintile_5groups + 
               ethnicgroup) 
summary(cox) 

# Drop the offending variable
cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + ethnicgroup) 
summary(cox)


############################ w4
# cox.zph(fit, transform="km", global=TRUE)

fit <- coxph(Surv(fu_time, death) ~ gender) # fit the desired model
temp <- cox.zph(fit)# apply the cox.zph function to the desired model
print(temp) # display the results
plot(temp) # plot the curves

fit <- coxph(Surv(fu_time, death) ~ copd,data=Survival) # fit the desired model
temp <- cox.zph(fit)# apply the cox.zph function to the desired model
print(temp,digits = 3) # display the results
plot(temp) # plot the curves

km_fit <- survfit(Surv(fu_time, death) ~ gender) 
autoplot(km_fit)
plot(km_fit, xlab = "time", ylab = "Survival probability") # label the axes 

res.cox <- coxph(Surv(fu_time, death) ~ age) 
ggcoxdiagnostics(res.cox, type = "dfbeta", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 

res.cox <- coxph(Surv(fu_time, death) ~ age) 
ggcoxdiagnostics(res.cox, type = "deviance", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 
ggcoxfunctional(Surv(fu_time, death) ~ age + log(age) + sqrt(age)) 



fit <- coxph(Surv(fu_time, death) ~ gender + tt(gender)) # "tt" is the time -transform function 
summary(fit) 


# make the other covariates 
ihd <- factor(g[,'ihd']) 
valvular <- factor(g[,'valvular_disease']) 
pvd <- factor(g[,'pvd']) 
stroke <- factor(g[,'stroke']) 
pneumonia <- factor(g[,'pneumonia']) 
renal <- factor(g[,'renal_disease']) 
ca <- factor(g[,'cancer']) 
mets <- factor(g[,'metastatic_cancer']) 
mental_health <- factor(g[,'mental_health']) 
ht <- factor(g[,"hypertension"])
cog_imp <- factor(g[,'dementia'])


# Results of the initial model
# run the full model 
cox <- coxph(Surv(fu_time, death) ~ age + gender + ethnicgroup + ihd + 
               valvular + pvd + stroke + copd + pneumonia + ht + renal + 
               ca + mets + mental_health + cog_imp) 
summary(cox) 

# Application of backwards elimination
cox <- coxph(Surv(fu_time, death) ~ age + gender + valvular + pneumonia + mets 
               + cog_imp) 
summary(cox) 

table(cog_imp) 
t <- table(cog_imp,death) 
t
round(100*prop.table(t,1),digits=1) 

# Testing the proportionality assumption on the remaining variables
fit <- coxph(Surv(fu_time, death) ~ age + gender + valvular + pneumonia + 
               mets + cog_imp) # test them all in the same model 
temp <- cox.zph(fit)  
print(temp) 

############################ final
# Install and load relevant packages 

install.packages("survival") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'survival' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##  install.packages("ggplot2") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'ggplot2' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##  install.packages("survminer") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'survminer' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##   

install.packages("ggfortify") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'ggfortify' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 


require(survival)  

## Loading required package: survival 

## Warning: package 'survival' was built under R version 3.5.1 

require(ggplot2)  

## Loading required package: ggplot2 

## Warning: package 'ggplot2' was built under R version 3.5.1 

require(survminer) 

## Loading required package: survminer 

## Warning: package 'survminer' was built under R version 3.5.1 

## Loading required package: ggpubr 

## Warning: package 'ggpubr' was built under R version 3.5.1 

## Loading required package: magrittr 

require(ggfortify) 

## Loading required package: ggfortify 

## Warning: package 'ggfortify' was built under R version 3.5.1 

# Load dataset 
g <- read.csv(file = paste0(getwd(),"/simulated HF mort data for GMPH (1K) final.csv"), header=TRUE, sep=',') 


# Define variables 
gender <- factor(g[,"gender"]) 
fu_time <- g[,"fu_time"]  
death <-  g[,"death"] 
age <- g[,"age"] 
copd <- factor(g[,"copd"]) 
ethnicgroup <- factor(g[,"ethnicgroup"]) 
quintile <- factor(g[,"quintile"]) 
ihd <- factor(g[,'ihd']) 
valvular <- factor(g[,'valvular_disease']) 
pvd <- factor(g[,'pvd']) 
stroke <- factor(g[,'stroke']) 
pneumonia <- factor(g[,'pneumonia']) 
renal <- factor(g[,'renal_disease']) 
ca <- factor(g[,'cancer']) 
mets <- factor(g[,'metastatic_cancer']) 
mental_health <- factor(g[,'mental_health']) 
ht <- factor(g[,"hypertension"]) 
cog_imp <- factor(g[,"senile"]) 
prior_dnas <- g[,"prior_dnas"] 



# Plotting a Kaplan-Meier curve 
###################### 

# 1. Generate the survival curve 
km_fit <- survfit(Surv(fu_time, death) ~ 1) 
# 2b. Alternative plot with ggplot2 
autoplot(km_fit) + theme_bw() # theme_bw() is a predesigned "theme" which makes the plot prettier 
###################### 

# Output the probability of survival at certain times after hospital admission 
summary(km_fit, times = c(1:7,30,60,90*(1:10))) 

## Call: survfit(formula = Surv(fu_time, death) ~ 1) 
##  
##  time n.risk n.event survival std.err lower 95% CI upper 95% CI 
##     1    992      12    0.988 0.00346        0.981        0.995 
##     2    973       7    0.981 0.00435        0.972        0.989 
##     3    963       5    0.976 0.00489        0.966        0.985 
##     4    954       6    0.970 0.00546        0.959        0.980 
##     5    945       5    0.964 0.00590        0.953        0.976 
##     6    938       1    0.963 0.00598        0.952        0.975 
##     7    933       1    0.962 0.00606        0.951        0.974 
##    30    865      39    0.921 0.00865        0.905        0.939 
##    60    809      28    0.891 0.01010        0.871        0.911 
##    90    770      24    0.864 0.01117        0.843        0.887 
##   180    698      43    0.815 0.01282        0.790        0.841 
##   270    653      24    0.787 0.01363        0.760        0.814 
##   360    619      21    0.761 0.01428        0.733        0.789 
##   450    525      44    0.705 0.01554        0.675        0.736 
##   540    429      47    0.639 0.01681        0.607        0.673 
##   630    362      32    0.589 0.01765        0.556        0.625 
##   720    266      43    0.514 0.01876        0.479        0.552 
##   810    190      31    0.448 0.01979        0.411        0.488 
##   900    126      26    0.378 0.02098        0.339        0.421 

# Plotting a Kaplan-Meier curve by gender 
###################### 

# 1. Generate the survival curve 
km_gender_fit <- survfit(Surv(fu_time, death) ~ gender) 

# 2. Plot the curve 
plot(km_gender_fit)

# 2b. Alternative plot with ggplot2 
autoplot(km_gender_fit) + theme_bw() 

###################### 


# Perform log rank test to see whether survival varies by gender 
survdiff(Surv(fu_time, death) ~ gender, rho = 0) 

## Call: 
## survdiff(formula = Surv(fu_time, death) ~ gender, rho = 0) 
##  
##            N Observed Expected (O-E)^2/E (O-E)^2/V 
## gender=1 548      268      271    0.0365     0.082 
## gender=2 452      224      221    0.0448     0.082 
##  
##  Chisq= 0.1  on 1 degrees of freedom, p= 0.8 

# Testing whether those over the age of 65 have different survival to those under it 
###################### 

# 1. Dichotomise age into categorical (binary in this case) variable 
age_65plus <- ifelse(g[,'age']>=65, 1, 0) 

# 2. Perform log rank test 
survdiff(Surv(fu_time, death) ~ age_65plus, rho = 0) 

## Call: 
## survdiff(formula = Surv(fu_time, death) ~ age_65plus, rho = 0) 
##  
##                N Observed Expected (O-E)^2/E (O-E)^2/V 
## age_65plus=0 115       18       67     35.85      41.7 
## age_65plus=1 885      474      425      5.65      41.7 
##  
##  Chisq= 41.7  on 1 degrees of freedom, p= 1e-10 

###################### 

# Plot survival curve by age above or below 65 
###################### 

# 1. Generate survival curve 
km_old_fit <- survfit(Surv(fu_time, death) ~ age_65plus) 

# 2. Plot 
plot(km_old_fit) 
# 2b. Alternative plot in ggplot2 
autoplot(km_old_fit) + theme_bw()
###################### 

# Run Cox regression model with age as predictor (continuous variable) 
###################### 

# 1. Generate model 
cox <- coxph(Surv(fu_time, death) ~ age, data = g) 

# 2. Summarise model 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ age, data = g) 
##  
##   n= 1000, number of events= 492  
##  
##         coef exp(coef) se(coef)     z Pr(>|z|)     
## age 0.056005  1.057602 0.005193 10.78   <2e-16 *** 
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##     exp(coef) exp(-coef) lower .95 upper .95 
## age     1.058     0.9455     1.047     1.068 
##  
## Concordance= 0.651  (se = 0.015 ) 
## Rsquare= 0.129   (max possible= 0.997 ) 
## Likelihood ratio test= 138  on 1 df,   p=<2e-16 
## Wald test            = 116.3  on 1 df,   p=<2e-16 
## Score (logrank) test = 115.7  on 1 df,   p=<2e-16 

###################### 



# Run Cox regression model with quintile as predictor (categorical variable) 
# Changing the reference group to first quintile 
# Removing the zero quintile altogether 
###################### 


# 1. Summarise the variable 
table(quintile, exclude = NULL) 

## quintile 
##    0    1    2    3    4    5 <NA>  
##    4  138  205  211  220  216    6 

# 2. Check levels 
levels(quintile) 

## [1] "0" "1" "2" "3" "4" "5" 

# 3. Generate model 
cox <- coxph(Surv(fu_time, death) ~ quintile) # warning 

## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, : 
## Loglik converged before variable 1,2,3,4,5 ; beta may be infinite. 

# 4. Summarise model 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ quintile) 
##  
##   n= 994, number of events= 489  
##    (6 observations deleted due to missingness) 
##  
##                coef exp(coef)  se(coef)     z Pr(>|z|) 
## quintile1 1.511e+01 3.665e+06 1.202e+03 0.013     0.99 
## quintile2 1.479e+01 2.645e+06 1.202e+03 0.012     0.99 
## quintile3 1.518e+01 3.912e+06 1.202e+03 0.013     0.99 
## quintile4 1.500e+01 3.261e+06 1.202e+03 0.012     0.99 
## quintile5 1.503e+01 3.383e+06 1.202e+03 0.013     0.99 
##  
##           exp(coef) exp(-coef) lower .95 upper .95 
## quintile1   3665107  2.728e-07         0       Inf 
## quintile2   2645345  3.780e-07         0       Inf 
## quintile3   3911559  2.557e-07         0       Inf 
## quintile4   3261463  3.066e-07         0       Inf 
## quintile5   3382505  2.956e-07         0       Inf 
##  
## Concordance= 0.544  (se = 0.015 ) 
## Rsquare= 0.013   (max possible= 0.997 ) 
## Likelihood ratio test= 13.22  on 5 df,   p=0.02 
## Wald test            = 8.39  on 5 df,   p=0.1 
## Score (logrank) test = 10.79  on 5 df,   p=0.06 

# 5. Make the first quintile the reference group 
quintile <- relevel(quintile, ref = "1") 

# 6. Regenerate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ quintile) # warning 

## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, : 
## Loglik converged before variable 1 ; beta may be infinite. 

summary(cox) # still an issue where quintile = 0 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ quintile) 
##  
##   n= 994, number of events= 489  
##    (6 observations deleted due to missingness) 
##  
##                 coef  exp(coef)   se(coef)      z Pr(>|z|)   
## quintile0 -1.511e+01  2.728e-07  1.202e+03 -0.013   0.9900   
## quintile2 -3.261e-01  7.218e-01  1.537e-01 -2.121   0.0339 * 
## quintile3  6.508e-02  1.067e+00  1.495e-01  0.435   0.6633   
## quintile4 -1.167e-01  8.899e-01  1.501e-01 -0.777   0.4369   
## quintile5 -8.024e-02  9.229e-01  1.490e-01 -0.538   0.5903   
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##           exp(coef) exp(-coef) lower .95 upper .95 
## quintile0 2.728e-07  3.665e+06    0.0000       Inf 
## quintile2 7.218e-01  1.385e+00    0.5340    0.9756 
## quintile3 1.067e+00  9.370e-01    0.7962    1.4306 
## quintile4 8.899e-01  1.124e+00    0.6631    1.1942 
## quintile5 9.229e-01  1.084e+00    0.6891    1.2360 
##  
## Concordance= 0.544  (se = 0.015 ) 
## Rsquare= 0.013   (max possible= 0.997 ) 
## Likelihood ratio test= 13.22  on 5 df,   p=0.02 
## Wald test            = 8.39  on 5 df,   p=0.1 
## Score (logrank) test = 10.79  on 5 df,   p=0.06 

# 7. Inspecting quintile variable 
table(quintile, g$death) # Only 4 entries for quintile = 0 and 100% didn't die 

##          
## quintile   0   1 
##        1  60  78 
##        0   4   0 
##        2 111  94 
##        3 105 106 
##        4 116 104 
##        5 109 107 

# 8. Removing quintile = 0 entries as there are only 4 of them 
quintile_5groups <- quintile 
quintile_5groups[quintile_5groups == 0] <- NA # set the zeroes to missing 
quintile_5groups <- factor(quintile_5groups) # this removes 0 as a level as it is an empty category 

# 9. Regenerating the model and summarising 
cox <- coxph(Surv(fu_time, death) ~ quintile_5groups) 
summary(cox) # still an issue where quintile = 0 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ quintile_5groups) 
##  
##   n= 990, number of events= 489  
##    (10 observations deleted due to missingness) 
##  
##                       coef exp(coef) se(coef)      z Pr(>|z|)   
## quintile_5groups2 -0.32606   0.72176  0.15374 -2.121   0.0339 * 
## quintile_5groups3  0.06508   1.06724  0.14949  0.435   0.6633   
## quintile_5groups4 -0.11668   0.88987  0.15010 -0.777   0.4369   
## quintile_5groups5 -0.08024   0.92289  0.14903 -0.538   0.5903   
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##                   exp(coef) exp(-coef) lower .95 upper .95 
## quintile_5groups2    0.7218      1.385    0.5340    0.9756 
## quintile_5groups3    1.0672      0.937    0.7962    1.4306 
## quintile_5groups4    0.8899      1.124    0.6631    1.1942 
## quintile_5groups5    0.9229      1.084    0.6891    1.2360 
##  
## Concordance= 0.542  (se = 0.015 ) 
## Rsquare= 0.009   (max possible= 0.997 ) 
## Likelihood ratio test= 8.62  on 4 df,   p=0.07 
## Wald test            = 8.39  on 4 df,   p=0.08 
## Score (logrank) test = 8.45  on 4 df,   p=0.08 

###################### 

# Run Cox regression model with ethnic group as predictor (categorical variable) 
# Including missing values as another category 
###################### 

# 1. Summarise variable 
table(ethnicgroup, exclude = NULL) 

## ethnicgroup 
##    1    2    3    9 <NA>  
##  889   17   34   17   43 

# 2. Generate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ ethnicgroup) 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ ethnicgroup) 
##  
##   n= 957, number of events= 471  
##    (43 observations deleted due to missingness) 
##  
##                  coef exp(coef) se(coef)      z Pr(>|z|)    
## ethnicgroup2 -0.06428   0.93774  0.32000 -0.201  0.84078    
## ethnicgroup3 -1.19586   0.30244  0.41108 -2.909  0.00362 ** 
## ethnicgroup9  0.07394   1.07674  0.35706  0.207  0.83596    
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##              exp(coef) exp(-coef) lower .95 upper .95 
## ethnicgroup2    0.9377     1.0664    0.5008    1.7558 
## ethnicgroup3    0.3024     3.3064    0.1351    0.6769 
## ethnicgroup9    1.0767     0.9287    0.5348    2.1679 
##  
## Concordance= 0.516  (se = 0.007 ) 
## Rsquare= 0.013   (max possible= 0.997 ) 
## Likelihood ratio test= 12.99  on 3 df,   p=0.005 
## Wald test            = 8.55  on 3 df,   p=0.04 
## Score (logrank) test = 9.61  on 3 df,   p=0.02 

# 3. Add another category (8) 
levels(ethnicgroup) <- c(levels(ethnicgroup),"8")  

# 4. Redefine NA as another group, 8 
ethnicgroup[is.na(ethnicgroup)] <- "8" 

# 5. Regenerate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ ethnicgroup) 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ ethnicgroup) 
##  
##   n= 1000, number of events= 492  
##  
##                  coef exp(coef) se(coef)      z Pr(>|z|)    
## ethnicgroup2 -0.06573   0.93638  0.31999 -0.205  0.83725    
## ethnicgroup3 -1.19368   0.30310  0.41107 -2.904  0.00369 ** 
## ethnicgroup9  0.08160   1.08502  0.35706  0.229  0.81923    
## ethnicgroup8 -0.02353   0.97675  0.22363 -0.105  0.91621    
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##              exp(coef) exp(-coef) lower .95 upper .95 
## ethnicgroup2    0.9364     1.0679    0.5001    1.7532 
## ethnicgroup3    0.3031     3.2992    0.1354    0.6784 
## ethnicgroup9    1.0850     0.9216    0.5389    2.1846 
## ethnicgroup8    0.9767     1.0238    0.6301    1.5140 
##  
## Concordance= 0.518  (se = 0.008 ) 
## Rsquare= 0.013   (max possible= 0.997 ) 
## Likelihood ratio test= 12.95  on 4 df,   p=0.01 
## Wald test            = 8.53  on 4 df,   p=0.07 
## Score (logrank) test = 9.58  on 4 df,   p=0.05 

###################### 



# Investigating our variables in order to best perform a Cox model with multiple predictors 
# Checking for missing values 
# Running a multiple Cox regression 
###################### 

# 1. Summarising age 
summary(g$age) # no NAs 

##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  
##   29.00   73.00   80.00   78.73   87.00  102.00 

hist(g$age) 

# 2. Gender 
gender_table <- table(gender, exclude = NULL) 
addmargins(gender_table) # no NAs 

## gender 
##    1    2  Sum  
##  548  452 1000 

round(100 * prop.table(gender_table), digits = 1) # Percentages rounded to 1 decimal place 

## gender 
##    1    2  
## 54.8 45.2 

# 3. Chronic Obstructive Pulmonary Disease (COPD) 
copd_table <- table(copd, exclude = NULL)  
addmargins(copd_table) # no NAs 

## copd 
##    0    1  Sum  
##  758  242 1000 

round(100 * prop.table(copd_table), digits = 1) # Percentages rounded to 1 decimal place 

## copd 
##    0    1  
## 75.8 24.2 

# 4. Prior OPD appointments missed  
prior_dnas_table <- table(prior_dnas, exclude = NULL)  
addmargins(prior_dnas_table) # no NAs 

## prior_dnas 
##    0    1    2    3    4    5    6    7    8   10  Sum  
##  732  156   50   34   17    3    3    2    1    2 1000 

round(100 * prop.table(prior_dnas_table), digits = 1) # Percentages rounded to 1 decimal place 

## prior_dnas 
##    0    1    2    3    4    5    6    7    8   10  
## 73.2 15.6  5.0  3.4  1.7  0.3  0.3  0.2  0.1  0.2 

# 5. Ethnic group 
ethnicgroup_table <- table(ethnicgroup, exclude = NULL)  
addmargins(ethnicgroup_table) # 4.3% NA 

## ethnicgroup 
##    1    2    3    9    8  Sum  
##  889   17   34   17   43 1000 

round(100 * prop.table(ethnicgroup_table), digits = 1) # Percentages rounded to 1 decimal place 

## ethnicgroup 
##    1    2    3    9    8  
## 88.9  1.7  3.4  1.7  4.3 

# 6. Generate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + prior_dnas + ethnicgroup) 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ age + gender + copd +  
##     prior_dnas + ethnicgroup) 
##  
##   n= 1000, number of events= 492  
##  
##                   coef exp(coef)  se(coef)      z Pr(>|z|)     
## age           0.061999  1.063961  0.005516 11.241  < 2e-16 *** 
## gender2      -0.253460  0.776111  0.094349 -2.686  0.00722 **  
## copd1         0.136649  1.146425  0.103880  1.315  0.18836     
## prior_dnas    0.163461  1.177579  0.039832  4.104 4.07e-05 *** 
## ethnicgroup2 -0.307915  0.734978  0.353009 -0.872  0.38307     
## ethnicgroup3 -0.823643  0.438830  0.414301 -1.988  0.04681 *   
## ethnicgroup9  0.408255  1.504190  0.360737  1.132  0.25775     
## ethnicgroup8 -0.045372  0.955642  0.225204 -0.201  0.84033     
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##              exp(coef) exp(-coef) lower .95 upper .95 
## age             1.0640     0.9399    1.0525    1.0755 
## gender2         0.7761     1.2885    0.6451    0.9338 
## copd1           1.1464     0.8723    0.9352    1.4053 
## prior_dnas      1.1776     0.8492    1.0891    1.2732 
## ethnicgroup2    0.7350     1.3606    0.3680    1.4681 
## ethnicgroup3    0.4388     2.2788    0.1948    0.9884 
## ethnicgroup9    1.5042     0.6648    0.7417    3.0504 
## ethnicgroup8    0.9556     1.0464    0.6146    1.4859 
##  
## Concordance= 0.667  (se = 0.015 ) 
## Rsquare= 0.155   (max possible= 0.997 ) 
## Likelihood ratio test= 168.4  on 8 df,   p=<2e-16 
## Wald test            = 141.7  on 8 df,   p=<2e-16 
## Score (logrank) test = 140  on 8 df,   p=<2e-16 

###################### 


# Investigating whether the assumptions of the Cox model are being broken 
# Testing for proportional hazards assumption (with gender as predictor variable) 
###################### 

# 1. Generate model fit 
fit <- coxph(Surv(fu_time, death) ~ gender) 

# 2. Apply the test to the model 
temp <- cox.zph(fit)     

# 3. Display results 
print(temp) 

##            rho chisq     p 
## gender2 0.0493  1.19 0.275 

# 4. Plot the curves 
plot(temp)     

# 4b. Alternative plot in ggplot 
ggcoxzph(temp) 

###################### 

# Generating other diagnostic plots for Cox Proportional Hazards model 
###################### 

# 1. Define model 
res.cox <- coxph(Surv(fu_time, death) ~ age) 

# Generate diagnostic plots 

# 2. Plotting the estimated changes in the regression coefficients on deleting each patient 
ggcoxdiagnostics(res.cox, type = "dfbeta", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 

# 3. Plotting deviance residuals 
ggcoxdiagnostics(res.cox, type = "deviance", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 

# 4. Plotting Martingale residuals 
fit <- coxph(Surv(fu_time, death) ~ age + log(age) + sqrt(age)) 
ggcoxfunctional(fit, data = g) # note we must specify original dataframe 

###################### 

# Testing proportionality assumption 
# Testing for a statistical relationship between gender and time 
###################### 

# 1. Generate model with time-transform function (tt) 
fit <- coxph(Surv(fu_time, death) ~ gender + tt(gender))  

# 2. Summarise 
summary(fit) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ gender + tt(gender)) 
##  
##   n= 1000, number of events= 492  
##  
##               coef exp(coef) se(coef)      z Pr(>|z|) 
## gender2     0.6405    1.8974   0.9800  0.654    0.513 
## tt(gender) -0.2003    0.8185   0.3182 -0.629    0.529 
##  
##            exp(coef) exp(-coef) lower .95 upper .95 
## gender2       1.8974      0.527    0.2779    12.953 
## tt(gender)    0.8185      1.222    0.4387     1.527 
##  
## Concordance= 0.497  (se = 0.197 ) 
## Rsquare= 0   (max possible= 0.997 ) 
## Likelihood ratio test= 0.49  on 2 df,   p=0.8 
## Wald test            = 0.48  on 2 df,   p=0.8 
## Score (logrank) test = 0.48  on 2 df,   p=0.8 

###################### 

# Backwards elimination to choose predictors for Cox regression 
###################### 



# 1. Run the full model with all of your predictors 
cox <- coxph(Surv(fu_time, death) ~ age + gender + ethnicgroup + ihd + 
               valvular + pvd + stroke + copd + pneumonia + ht + renal + 
               ca + mets + mental_health + cog_imp) 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ age + gender + ethnicgroup +  
##     ihd + valvular + pvd + stroke + copd + pneumonia + ht + renal +  
##     ca + mets + mental_health + cog_imp) 
##  
##   n= 1000, number of events= 492  
##  
##                     coef exp(coef)  se(coef)      z Pr(>|z|)     
## age             0.057681  1.059376  0.005615 10.272  < 2e-16 *** 
## gender2        -0.175075  0.839394  0.096832 -1.808  0.07060 .   
## ethnicgroup2    0.118322  1.125606  0.327954  0.361  0.71826     
## ethnicgroup3   -0.681361  0.505928  0.415695 -1.639  0.10120     
## ethnicgroup9    0.506455  1.659399  0.363311  1.394  0.16332     
## ethnicgroup8   -0.136391  0.872502  0.228592 -0.597  0.55074     
## ihd1            0.182438  1.200139  0.095302  1.914  0.05558 .   
## valvular1       0.227382  1.255309  0.108136  2.103  0.03549 *   
## pvd1           -0.016700  0.983439  0.161517 -0.103  0.91765     
## stroke1         0.059764  1.061586  0.298368  0.200  0.84124     
## copd1           0.066729  1.069005  0.106415  0.627  0.53062     
## pneumonia1      0.434769  1.544606  0.138204  3.146  0.00166 **  
## ht1            -0.054098  0.947339  0.096111 -0.563  0.57352     
## renal1          0.209080  1.232543  0.107546  1.944  0.05188 .   
## ca1             0.271763  1.312276  0.207485  1.310  0.19026     
## mets1           2.176735  8.817473  0.401748  5.418 6.02e-08 *** 
## mental_health1 -0.014738  0.985370  0.180778 -0.082  0.93502     
## cog_imp1        0.282561  1.326522  0.187484  1.507  0.13178     
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##                exp(coef) exp(-coef) lower .95 upper .95 
## age               1.0594     0.9440    1.0478     1.071 
## gender2           0.8394     1.1913    0.6943     1.015 
## ethnicgroup2      1.1256     0.8884    0.5919     2.141 
## ethnicgroup3      0.5059     1.9766    0.2240     1.143 
## ethnicgroup9      1.6594     0.6026    0.8141     3.382 
## ethnicgroup8      0.8725     1.1461    0.5574     1.366 
## ihd1              1.2001     0.8332    0.9957     1.447 
## valvular1         1.2553     0.7966    1.0156     1.552 
## pvd1              0.9834     1.0168    0.7166     1.350 
## stroke1           1.0616     0.9420    0.5915     1.905 
## copd1             1.0690     0.9354    0.8678     1.317 
## pneumonia1        1.5446     0.6474    1.1781     2.025 
## ht1               0.9473     1.0556    0.7847     1.144 
## renal1            1.2325     0.8113    0.9983     1.522 
## ca1               1.3123     0.7620    0.8738     1.971 
## mets1             8.8175     0.1134    4.0121    19.378 
## mental_health1    0.9854     1.0148    0.6914     1.404 
## cog_imp1          1.3265     0.7539    0.9186     1.916 
##  
## Concordance= 0.694  (se = 0.015 ) 
## Rsquare= 0.183   (max possible= 0.997 ) 
## Likelihood ratio test= 202.7  on 18 df,   p=<2e-16 
## Wald test            = 193.7  on 18 df,   p=<2e-16 
## Score (logrank) test = 208.5  on 18 df,   p=<2e-16 

# 2. Run the model with only significant predictors 
cox <- coxph(Surv(fu_time, death) ~ age + gender + valvular + pneumonia + mets + cog_imp) 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ age + gender + valvular +  
##     pneumonia + mets + cog_imp) 
##  
##   n= 1000, number of events= 492  
##  
##                 coef exp(coef)  se(coef)      z Pr(>|z|)     
## age         0.059147  1.060931  0.005441 10.871  < 2e-16 *** 
## gender2    -0.260853  0.770394  0.093497 -2.790  0.00527 **  
## valvular1   0.229659  1.258170  0.107118  2.144  0.03204 *   
## pneumonia1  0.481375  1.618299  0.135192  3.561  0.00037 *** 
## mets1       2.494299 12.113233  0.364373  6.845 7.62e-12 *** 
## cog_imp1    0.241623  1.273314  0.184232  1.312  0.18968     
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##            exp(coef) exp(-coef) lower .95 upper .95 
## age           1.0609    0.94257    1.0497    1.0723 
## gender2       0.7704    1.29804    0.6414    0.9253 
## valvular1     1.2582    0.79480    1.0199    1.5521 
## pneumonia1    1.6183    0.61793    1.2416    2.1093 
## mets1        12.1132    0.08255    5.9307   24.7409 
## cog_imp1      1.2733    0.78535    0.8874    1.8271 
##  
## Concordance= 0.686  (se = 0.015 ) 
## Rsquare= 0.17   (max possible= 0.997 ) 
## Likelihood ratio test= 186.9  on 6 df,   p=<2e-16 
## Wald test            = 179.5  on 6 df,   p=<2e-16 
## Score (logrank) test = 193  on 6 df,   p=<2e-16 

# 3. Test proportionality assumption on these predictors 
cox.zph(cox)  

##                 rho  chisq     p 
## age        -0.04409 0.9098 0.340 
## gender2     0.05474 1.4133 0.235 
## valvular1  -0.04565 1.0202 0.312 
## pneumonia1 -0.04397 0.9540 0.329 
## mets1       0.00615 0.0184 0.892 
## cog_imp1    0.05985 1.8140 0.178 
## GLOBAL           NA 5.4044 0.493 