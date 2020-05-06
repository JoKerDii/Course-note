

### POSIXt classes in R

x = as.POSIXct("2019-12-25 11:45:34") # nr of seconds

y = as.POSIXlt("2019-12-25 11:45:34")

x; y # it gives the same output, but what is behind it?

unclass(x)

unclass(y)

# what does the number mean?

(50 * 365 * 24 * 60 * 60) - (5.5 * 60 * 60)

y$zone # extracting the elements from POSIXlt

x$zone # not possible since it is simply a number of seconds

# another class based on days

x = as.Date("2019-12-25")

x; class(x)

unclass(x)

50 * 365 - 5 # nr of days since 1970

library(chron)

x = chron("12/25/2019", "23:34:09")

x

class(x)

unclass(x)



### strptime



a = as.character(c("1993-12-30 23:45","1994-11-05 11:43","1992-03-09 21:54"))

class(a)



b = strptime(a, format = "%Y-%m-%d %H:%M")

b; class(b)







### Lets take a look at the package lubridate which has very useful time/date data functions

library(lubridate)

# different ways in how to input dates

ymd(19931123)

dmy(23111993)

mdy(11231993)

# lets use time and date together

mytimepoint <- ymd_hm("1993-11-23 11:23", tz = "Europe/Prague")

mytimepoint

class(mytimepoint)

# extracting the components of it

minute(mytimepoint)

day(mytimepoint)

hour(mytimepoint)

year(mytimepoint)

month(mytimepoint)

# we can even change time values within our object

hour(mytimepoint) <- 14

mytimepoint

# which time zones do we have available

OlsonNames()

# we can take a look at the most common time zones

# but be aware that the time zone recognition also depends on your location and machine



## lets check which day our time point is

wday(mytimepoint)

wday(mytimepoint, label=T, abbr=F) # label to display the name of the day, no abbreviation

# we can calculate which time our timepoint would be in another time zone

with_tz(mytimepoint, tz = "Europe/London")

mytimepoint

# time intervals

time1 = ymd_hm("1993-09-23 11:23", tz = "Europe/Prague")

time2 = ymd_hm("1995-11-02 15:23", tz = "Europe/Prague")

# getting the interval

myinterval = interval(time1, time2); myinterval

class(myinterval) # interval is an object class from lubridate



### Exercise: Creating a Data Frame with lubridate

# lets now build a dataframe with lubridate that contains date and time data

# see the different input formats that are allowed in the ymd function

a = c("1998,11,11", "1983/01/23", "1982:09:04", "1945-05-09", 19821224, "1974.12.03", 19871210)

a = ymd(a, tz = "CET") ;a

# now I am creating a time vector - using different notations of input

b = c("22 4 5", "04;09;45", "11:9:56", "23,15,12", "14 16 34", "8 8 23", "21 16 14")

b = hms(b); b

f = rnorm(7,10); f = round(f, digits = 2); f

date_time_measurement = cbind.data.frame(date = a, time = b, measurement = f)

date_time_measurement





## Calculations with time

minutes(7)

# note that class "Period" needs integers - full numbers

minutes(2.5)

# getting the duration

dminutes(3)

dminutes(3.5)

# how to add minutes and seconds

minutes(2) + seconds(5)

# more calculations

minutes(2) + seconds(75)

# class "duration" to perform addition

as.duration(minutes(2) + seconds(75))

# lubridate has many time classes: period or duration differ!

# which year was a leap year?

leap_year(2009:2014)

ymd(20140101) + years(1)

ymd(20140101) + dyears(1)

# lets do the whole thing with a leap year

leap_year(2016)

ymd(20160101) + years(1)

ymd(20160101) + dyears(1)

# as you see the duration is the one which is always 365 days

# the standard one (the period) makes the year a whole new unit (+1)



## Exercise Lubridate

# create x, with time zone CET and a given time point in 2014 of your choosing

# I use "2014-04-12 23:12"

# the time point consists of year, months, day and hour

# change now the minute of x to 7 and check x in the same line of code

# see which time it would be in London

# create another time point y in 2015 and get the difference between those 2 points



x = ymd_hm(tz = "CET", "2014-04-12 23:12")

minute(x) = 7 ; x

with_tz(x, tz="Europe/London")

y = ymd_hm(tz = "CET", "2015-12-12 09:45")

y-x



### Lubridate: calculating the time column

head(udca)

# simple method for properly formatted columns
x1 =udca$last.dt - udca$entry.dt
class(x1)

# lubridate way is more flexible
library(lubridate)

x2 = as.period(interval(udca$entry.dt, udca$last.dt), 
               unit = "day")
x2 = day(x2)

# add new column to dataframe
udca.new = data.frame(udca, x1, x2); head(udca.new)

str(udca.new)
