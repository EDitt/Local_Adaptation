####  Analysis of Rainfall Data
####  Dataset: "JR daily precip 1983-2015.xls"
####  Obtained from Jasper Ridge Biological Preserve on 9/10/2015 from Nona Chiariello
###         Data were collected from weather station at Jasper Ridge Biological Preserve as well as 
###         Woodside Fire Station from NOAA (used to fill in gaps when precipitation data were missing from JRBP)

library(readxl)
library(tidyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(ggthemes)

##############################################
############### DATA WRANGLING ############### 
##############################################

#read in excel file
precip <- read_excel("JR daily precip 1983-2015.xls")

str(precip)

#combine first two rows as header
names(precip) <- paste(names(precip), precip[1, ], sep = "_")
precip <- precip[-1,]
head(precip)
tail(precip)
length(precip$`1983...2_JRBP`)

#delete last 4 rows
precip <-precip[-c(366:370), ]

#make long format
precip_long <- gather(precip, location, rainfall, "1983...2_JRBP":"2015_JRBP", factor_key = TRUE)

#separate year and location into different columns
precip_long <- separate(precip_long, location, c("Year", "Location"), sep = "_", remove = TRUE)

#remove the "....2" from the year
precip_long$Year <- gsub("\\...*", "", precip_long$Year)

#make wide again for year
precip_wide <- spread(precip_long, Location, rainfall)

str(precip_wide)

#make numbers instead of characters
precip_wide$JRBP<-as.numeric(precip_wide$JRBP)
precip_wide$Wdside<-as.numeric(precip_wide$Wdside)

#### If no data from JRBP (e.g. "NA"), use rainfall data from Woodside weather station:

#label where rainfall data is from first
precip_wide$location = with(precip_wide, ifelse(is.na(precip_wide$JRBP), "Woodside", "JR"))

#Use Woodside data if no JRBP data
precip_wide$rainfall = with(precip_wide, ifelse(is.na(precip_wide$JRBP), precip_wide$Wdside, precip_wide$JRBP))


#make day of year be number
precip_wide$Day <- as.numeric(precip_wide$`...1_day of yr`)

#make Year a factor + number
precip_wide$fYear <- factor(precip_wide$Year)
precip_wide$Year <- as.numeric(precip_wide$Year)

precip_df <- data.frame(precip_wide)

##### Make dates
precip_df$date1 <- as.Date(precip_df$Day - 1, origin = "2014-01-01") #turn day of year into date with arbitrary year
precip_df$date2 <- format(precip_df$date1, format = "%d/%m") #remove arbitrary year
precip_df$month <- format(precip_df$date1, format = "%m") #column for month
precip_df$date3 <- paste(precip_df$date2, precip_df$Year, sep = "/") #add correct year to date
precip_df$date3 <- as.Date(precip_df$date3, "%d/%m/%Y") #date format

aggregate(precip_df$Day, by=list(precip_df$Year), length) #check

colnames(precip_df)[1] <- "Day_of_Year"
write.csv(precip_df, file="JRPrecipitationData.csv",
          quote=FALSE, row.names = FALSE)
# this formatted version uploaded to Dryad

##############################################
######### GROWING SEASON RAINFALL ############
##############################################

### The most relevant "year" for Leptosiphon is from July the previous year (when seeds are dispersed) through June
### (when senescence has occurred). I will therefore convert strict year dates to "grow years"

#total rainfall across years
aggregate(precip_df$rainfall, by=list(precip_df$Year), sum, na.rm=T)

precip_df$month <- as.numeric(precip_df$month)

#rainfall for each month of each year
bymonth <- aggregate(precip_df$rainfall, by=list(precip_df$Year, precip_df$month), sum, na.rm=T)

#Will have "grow years" start on June (from the previous year) until June (of the current year)
#So- if month is 6-12 (June - Dec), grow year will be year + 1
bymonth$grow.year = with(bymonth, ifelse(Group.2 > 6, Group.1 + 1, Group.1))

#I can also define growing season whether pre or post-flowering (e.g. first or second half of growing season)
#Pre-flowering (e.g. first half) is Dec-Feb, post-flowering (e.g. second half) is March-May
bymonth$season2 = with(bymonth, 
                       ifelse((bymonth$Group.2==12 | bymonth$Group.2<=2), "firsthalf",
                              ifelse ((bymonth$Group.2>=3 & bymonth$Group.2<=5), "secondhalf",
                                      "dry")))

aggregate(bymonth$x, by=list(bymonth$season2), length) #check

#Create new dataframe
newdf <- aggregate(bymonth$x, by=list(bymonth$grow.year, bymonth$season2), sum, na.rm=T)

newdf_wide <- spread(newdf, Group.2, x)

head(newdf_wide)

#I can also look at rainfall through the growing season vs. dry season
#what proportion is during growing season (Dec-May)?
#define growing season (Dec-May)

bymonth$season = with(bymonth, 
                      ifelse((bymonth$Group.2==12 | bymonth$Group.2<=5), "grow", "dry"))

newdf2 <- aggregate(bymonth$x, by=list(bymonth$grow.year, bymonth$season), sum, na.rm=T)
newdf2_wide <- spread(newdf2, Group.2, x)

#total overall
newdf3 <- aggregate(bymonth$x, by=list(bymonth$grow.year), sum, na.rm=T)
colnames(newdf3) <- c("Group.1", "Total")

#merge dataframes
Rain_sum1 <- merge(newdf_wide, newdf2_wide, by=c("Group.1", "dry"))
Rain_sum <- merge(newdf3, Rain_sum1, by="Group.1")

#Proportion of rain occurring during the growing season
Rain_sum$Propgrow <- Rain_sum$grow / Rain_sum$Total

#Proportion of growing season rain occurring during the first half
Rain_sum$Propgrow1st <- Rain_sum$firsthalf / Rain_sum$grow
#Proportion of growing season rain occurring during the second half
Rain_sum$Propgrow2nd <- Rain_sum$secondhalf / Rain_sum$grow

#Ratio of first half:second half rain
Rain_sum$FirstrelSecond <- Rain_sum$firsthalf / Rain_sum$secondhalf

#Overall proportion of rain occurring Dec-Feb
Rain_sum$Propgrow1stoverall <- Rain_sum$firsthalf / Rain_sum$Total
#Overall proportion of rain occurring March-May
Rain_sum$Propgrow2ndoverall <- Rain_sum$secondhalf / Rain_sum$Total


##############################################
################ PLOTTING ####################
##############################################

#because only partial data for 2016, will remove
Rain_sum_no2016 <- subset(Rain_sum, Group.1!=2016)

length(Rain_sum_no2016$Total) #33
mean(Rain_sum_no2016$Total) #653.34
sd(Rain_sum_no2016$Total) #238.56
min(Rain_sum_no2016$Total) #252.476
max(Rain_sum_no2016$Total) #1333.502
SE <- sd(Rain_sum_no2016$Total) / (sqrt(length(Rain_sum_no2016$Total)))

#Plot with ggplot
rain<-ggplot(Rain_sum_no2016, aes(x=Group.1, y=Total))
rain + geom_point() + geom_line() + theme_clean()

#just the years I did field work
field_season <- subset (Rain_sum, Group.1>=2012 & Group.1<=2015)

####Plot with base R
par(mar=c(3,4,1,1))
plot(Rain_sum_no2016$Total ~ Rain_sum_no2016$Group.1, ylim=c(200, 1500), type="o", pch=1, axes=T, xlab="", ylab="", cex.axis=1.25, las=1)
points(field_season$Total ~ field_season$Group.1, pch=19, col=c("black"), cex=1) #reciprocal transplant years
abline(h=mean(Rain_sum_no2016$Total), lty=2) #mean rainfall
abline(h=mean(Rain_sum_no2016$Total) - SE, lty=3, col="grey") #mean - 1 SE
abline(h=mean(Rain_sum_no2016$Total) + SE, lty=3, col="grey") #mean + 1 SE

