
### Missing Data Handling

# Use FlowerPicks csv - remove ID column

summary(FlowerPicks)

# Vector with complete cases

compcases <- complete.cases(FlowerPicks) == T

# Per default R ignores missing data - na.action argument

boxplot(Score ~ Player,FlowerPicks)



lm(data = FlowerPicks, formula = Score ~ Time)



# Using na.omit to remove al rows with an NA

cleandata <- na.omit(FlowerPicks)

summary(cleandata)

# Data lost several rows - only 5475 observations



# zoo for column wise operations

library(zoo)

x = na.locf(FlowerPicks$Score)

summary(x)





## Machine Learning for NA removal

library(mice)

# Distribution of the missing values

md.pattern(FlowerPicks)

# 1 row has 2 NAs



# Using the mice function

mymice <- mice(FlowerPicks, m = 10,method = "rf")



# The result has class mids

class(mymice)



# Display the calculated data

mymice$imp$Score



# Fill the NAs

mymicecomplete <- complete(mymice, 5)

summary(mymicecomplete)



# Analysis with m variations of the dataset

lmfit <- with(mymice, lm(Score ~ Time))



# Pooling the results

summary(pool(lmfit))









## Outlier Detection: simple method for outlier calculation with ESD



# [mean - t * SD, mean + t * SD]

x = c(rnorm(10), 150); x

t = 3

m = mean(x)

s = sd(x)

b1 = m - s*t; b1 
b2 = m + s*t; b2

y = ifelse(x >= b1 & x <= b2, 0, 1); y 
plot(x, col=y+2)

# simple boxplot

boxplot(x)

boxplot.stats(x)

# package outliers

library(outliers)

dixon.test(x)

grubbs.test(x, type = 11, two.sided = T) # type 11 for 2 sides, opossite

### advanced techniques for multivariate data

library(mvoutlier)

elements = data.frame(Hg = moss$Hg, Fe = moss$Fe, Al = moss$Al, Ni = moss$Ni) 
head(elements)

myout = sign1(elements[,1:4], qqcrit = 0.975); myout

myout = pcout(elements[,1:4])

plot(moss$Fe, moss$Al, col=myout$wfinal01+2)
