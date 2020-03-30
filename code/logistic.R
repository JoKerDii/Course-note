g <- read.csv(file = "survival analysis/final diabetes data for R (csv)(2).csv", header= TRUE, sep = ",")
dimnames(g)[[2]]
chol <- g[,"chol"] # cholesterol is continuous, so it’s easy
gender <- as.factor(g[,'gender']) # but gender isn’t.
dm <- as.factor(g[,"dm"]) # neither is dm
t <- table(gender) # store the tabulation for further manipulation
addmargins(t) # this will sum up the gender totals to give an overall total and print the results
round(prop.table(t),digits=3) # get proportions rounded to 3dp
round(100*prop.table(t),digits=1) # get %s rounded to 1dp
dm2 <- factor(dm, exclude=NULL) # make new factor from the old one
table(dm2) # display the counts including the missings (NAs)
summary(chol)
height <- g[,'height']
weight <- g[,'weight']
summary(height)
summary(weight)

height.si <- height*0.0254
weight.si <- weight*0.453592
bmi <- weight.si/height.si^2
summary(bmi)

bmi_categorised <- ifelse(bmi < 18.5, "underweight", 
                          ifelse(bmi >= 18.5 & bmi <= 25, "normal", 
                                 ifelse(bmi > 25 & bmi <= 30, "overweight", 
                                        ifelse(bmi > 30, "obese", NA)))) 
table(bmi_categorised, exclude = NULL) # check that the bmi_categorised variable has worked  

dm_by_bmi_category <- table(bmi_categorised, dm2, exclude = NULL) # frequencies of diabetes by BMI category 
dm_by_bmi_category # check 
round(100 * prop.table(dm_by_bmi_category, margin = 1), digits = 1) 

age <- g[,"age"]
age_categorised <- ifelse(age < 45, "<45", 
                          ifelse(age >= 45 & age <= 64, "45-64", 
                                 ifelse(age >= 65 & age <= 74, "65-74", 
                                        ifelse(age >= 75, ">=75",NA)))) 
table(age_categorised, exclude = NULL) # check that the age_categorised variable has worked  
gender <- factor(g[,'gender'])
gender_by_age_category <- table(age_categorised, gender, exclude = NULL) # frequencies of diabetes by BMI category 
round(100 * prop.table(gender_by_age_category, margin = 2), digits = 1) 




##### Here is the R code to do the cross-tabulations and the resulting output 

# creating "age" variable 
age <- g[,"age"] 

# creating a categorical variable "age_grouped" 
age_grouped <- ifelse(age < 45, "under 45", 
                      ifelse(age >= 45 & age < 65, "45 - 64",  
                             ifelse(age >= 65 & age < 75, "65 - 74",  
                                    ifelse(age >= 75, "75 or over", NA)))) 

# displaying new variable in a table 
table(age_grouped, exclude = NULL) 

## age_grouped 
##    45 - 64    65 - 74 75 or over   under 45  
##        139         41         23        200 

# cross tabulating with gender 
age_group_by_gender <- table(age_grouped, gender, exclude = NULL) 

# display the cross tabulation 
age_group_by_gender 

##             gender 
## age_grouped  female male 
##   45 - 64        75   64 
##   65 - 74        21   20 
##   75 or over     12   11 
##   under 45      126   74 

# display the cross tabulation as proportion of whole sample, converting to percentage and rounding to 1 decimal place 
round(100 * prop.table(age_group_by_gender), digits = 1) 

##             gender 
## age_grouped  female male 
##   45 - 64      18.6 15.9 
##   65 - 74       5.2  5.0 
##   75 or over    3.0  2.7 
##   under 45     31.3 18.4 

# displaying the age frequencies by gender 
round(100 * prop.table(age_group_by_gender, margin = 2), digits = 1) 

##             gender 
## age_grouped  female male 
##   45 - 64      32.1 37.9 
##   65 - 74       9.0 11.8 
##   75 or over    5.1  6.5 
##   under 45     53.8 43.8 

# initialising the age_grouped vector by copying the already existing age vector 
age_grouped <- age 


# below says: if age < 45, then label the value "under 45", if not, then keep it what it already was in age_grouped 
age_grouped <- ifelse(age < 45, "under 45", age_grouped) 

# repeat for the other categories 
age_grouped <- ifelse(age >= 45 & age < 65, "45 - 64", age_grouped) 
age_grouped <- ifelse(age >= 65 & age < 75, "65 - 74", age_grouped)  
age_grouped <- ifelse(age >= 75, "75 or over", age_grouped) 



# check that things make sense 
table(age_grouped, exclude = NULL) 

## age_grouped 
##    45 - 64    65 - 74 75 or over   under 45  
##        139         41         23        200 

# optional extra check for the extra cautious! 
head(cbind(age_grouped, age)) 

##      age_grouped age  
## [1,] "45 - 64"   "46" 
## [2,] "under 45"  "29" 
## [3,] "45 - 64"   "58" 
## [4,] "65 - 74"   "67" 
## [5,] "45 - 64"   "64" 
## [6,] "under 45"  "34" 



glm(dm ~ 1, family=binomial (link=logit))
m <- glm(dm ~ 1, family=binomial (link=logit))
summary(m)
table(m$y)
m <- glm(dm ~ gender, family=binomial (link=logit))
summary(m)
m <- glm(dm ~ age, family=binomial (link=logit))
summary(m)


# create a cross tabulation of age and diabetes status  
dm_by_age <- table(age, dm) 

# output the frequencies of diabetes status by age 
freq_table <- prop.table(dm_by_age, margin = 1) 

# calculate the odds of having diabetes 
odds <- freq_table[, "yes"]/freq_table[, "no"] 

# calculate the log odds 
logodds <- log(odds) 

# plot the ages found in the sample against the log odds of having diabetes 
plot(rownames(freq_table), logodds) 



m <- glm(dm ~ 1, family=binomial (link=logit)) 
summary(m) 
exp(-1.7047)
table(dm)
m <- glm(dm ~ age, family=binomial (link=logit)) 
summary(m) 
m <- glm(dm ~ gender, family=binomial (link=logit)) 
summary(m) 

contrasts(gender)
levels(gender) 
gender <- relevel(gender, ref = "male") 
levels(gender) 

m <- glm(dm ~  gender, family=binomial (link=logit)) 
summary(m) 

m$coefficients # log odds
exp(m$coefficients) # odds ratios

# quiz
location <- factor(g[,'location'], exclude = NULL)
dm_by_loc_category <- table(location, dm2, exclude = NULL) # frequencies of diabetes by BMI category 
round(100 * prop.table(dm_by_loc_category, margin = 1), digits = 1) 
m <- glm(dm ~  location, family=binomial (link=logit)) 
summary(m)
exp(m$coefficients) # odds ratios

# quiz
chol <- g[,"chol"] 
insurance <- g[,'insurance']

dm_by_loc_category <- table(age, dm2, exclude = NULL) # frequencies of diabetes by BMI category 
round(100 * prop.table(dm_by_loc_category, margin = 1), digits = 1) 
m <- glm(dm ~  age + chol + insur_categorised, family=binomial (link=logit)) 
summary(m)
exp(m$coefficients) # odds ratios

insur_categorised <- ifelse(insurance <= 0, "0",  
                           ifelse(insurance >= 1 & insurance < 2, "goverment", 
                                  ifelse(insurance >= 2, "private", NA))) 
m <- glm(dm ~  insur_categorised, family=binomial (link=logit)) 
summary(m)
exp(m$coefficients)

# define the age variable (continuous) 
age <- g[,"age"] 


# create a cross tabulation of age and diabetes status  
dm_by_age <- table(age, dm) # not including NA values because there aren't that many 

# output the frequencies of diabetes status by age 
dm_by_age_prop <- prop.table(dm_by_age, margin = 1) 

# calculate the odds of having diabetes 
odds_age <- dm_by_age_prop[, "yes"]/dm_by_age_prop[, "no"] 

# calculate the log odds 
logodds_age <- log(odds_age) 

# plot the ages found in the sample against the log odds of having diabetes 
plot(rownames(dm_by_age_prop), logodds_age) 

# w4
# design your logistic regression 
full_model <- glm(dm ~ age + chol + insurance, family=binomial (link=logit)) 

# check your model 
summary(full_model) 

# run a null model 
null_model <- glm(dm ~ 1, family=binomial (link=logit)) 

# check 
summary(null_model) 

# calculate McFadden's R-square 
R2 <- 1-logLik(full_model)/logLik(null_model) 
# print it 
R2 

# install a package 
install.packages("DescTools") 
# load package 
require(DescTools) 

# design your logistic regression 
full_model <- glm(dm ~ age + chol + insurance, family=binomial (link=logit)) 

# check your model 
summary(full_model) 

# generate the c-statistic 
Cstat(full_model) 

# install package "ResourceSelection" 
install.packages("ResourceSelection") 
# load package 
require(ResourceSelection) 

# design your logistic regression 
full_model <- glm(dm ~ age + chol + insurance, family = binomial(link = logit)) 
full_model$y

# run Hosmer-Lemeshow test 
HL <- hoslem.test(x = full_model$y, y = fitted(full_model), g = 10) 
HL  

# plot the observed vs expected number of cases for each of the 10 groups 
plot(HL$observed[,"y1"], HL$expected[,"yhat1"]) 

# plot the observed vs expected number of noncases for each of the 10 groups 
plot(HL$observed[,"y0"], HL$expected[,"yhat0"]) 

# plot observed vs. expected prevalence for each of the 10 groups 
plot(x = HL$observed[,"y1"]/(HL$observed[,"y1"]+HL$observed[,"y0"]), 
     y = HL$expected[,"yhat1"]/(HL$expected[,"yhat1"]+HL$expected[,"yhat0"])) 

# install package("generalhoslem") 
install.packages("generalhoslem") 
# load package 
require(generalhoslem) 

# run Hosmer-Lemeshow test 
logitgof(obs = full_model$y, exp = fitted(full_model), g = 10) 

# design your logistic regression 
full_model <- glm(dm ~ age + chol + insurance, family = binomial(link = logit)) 

# analyse table of deviance 
anova(full_model, test = "Chisq") 

##### Make the variables and run the models #####
dm <- as.factor(g[,"dm"]) 
insurance <- as.factor(g[,"insurance"])# let's say 0=none, 1=gov, 2=private 
fh <- as.factor(g[,"fh"]) # 1=FH, 0=no FH 
smoking <- as.factor(g[,"smoking"]) # 1,2,3 
chol <- g[,'chol'] 
hdl <- g[,'hdl'] 
ratio <- g[,'ratio'] 
location <- as.factor(g[,'location']) 
age <- g[,'age'] 
gender <- as.factor(g[,'gender']) 
frame <- as.factor(g[,'frame']) 
systolic <- g[,'bp.1s'] 
diastolic <- g[,'bp.1d'] 

model <- glm(dm ~ age + bmi + chol + hdl + systolic + diastolic, family = 
               binomial(link = logit)) 

summary(model) 

anova(model, test = "Chisq") 

# drop insignificant variables
model1 <- glm(formula = dm ~ age + bmi + chol + hdl, family = binomial(link = logit)) 
summary(model1)

# strange that systolic and diastolic are not significant... 
cor.test(systolic, hdl) # not significant 
cor.test(systolic, bmi) # significant 
cor.test(systolic, chol) # very significant
cor.test(systolic, age) # extremely significant 


# adding gender + location + frame + insurance + smoking 
model <- glm(dm ~ age + bmi + chol + hdl + systolic + diastolic + gender + 
               location + frame + insurance + smoking, family = binomial(link = logit)) 

summary(model) 
