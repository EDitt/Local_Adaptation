####  Analysis of Field Reciprocal Transplant Data
####  Dataset: ""
####  

library(dplyr)
library(aster)

##############################################
################## DATA SETUP ################ 
##############################################

Field_All <- read.csv("DataFiles/ReciprocalTransplant.csv", header = T)

str(Field_All)

#Use only parental populations in analyses:
Field <- subset(Field_All, Population =="SerpPop" | Population == "SandPop")

Field$Population <- factor(Field$Population) #to use only existing levels

### INDEPENDENT VARIABLES:

levels(Field$Population)

levels(Field$Soil)

Field$Year <- factor(Field$Year) #Year treated as a factor with 4 levels
levels(Field$Year)

levels(Field$Edge) #Edge versus Non-edge

Field$Plot_Rep <- factor(Field$Plot_Rep) #Plot rep is a factor (nested within soil type and year)
levels(Field$Plot_Rep)

### DEPENDENT VARIABLES:

str(Field$Surv_flr) ### Survival to flowering, bernoulli distribution

str(Field$Num_flrs) ### poisson, truncated at 0
hist(Field$Num_flrs)

str(Field$Num_frts) ### poisson
hist(Field$Num_frts)

#The variables that correspond to nodes of the graph are, in the order they are numbered in the graph:
vars<-c("Surv_flr", "Num_flrs", "Num_frts")
#e.g. the 3 fitness components

##############################################
################## ASTER SETUP ###############
##############################################

#Need to remove NA's
fieldnoNA<-subset(Field, Surv_flr>=0 & Num_flrs>=0 & Num_frts>=0) #N=1354 out of 1599

#Very few individuals from Sandstone Population survived to flower on Serpentine Soil:
aggregate(fieldnoNA$Year, by=list(fieldnoNA$Year, fieldnoNA$Soil, fieldnoNA$Population, fieldnoNA$Surv_flr), length)

#Reshape data into long format so there is one line per node of the graph rather than one per individual
redata <- reshape(fieldnoNA, varying = list(vars), direction = "long", timevar = "varb", times = as.factor(vars), v.names = "resp")

names(redata)
#'varb' indicates which original variable the data in 'resp' come from (survival, number of flowers, number of fruits)
#id indicates which original individual the corresponding element of resp came from.

class(redata$varb) #varb must be a factor

#Add a root
redata <- data.frame(redata, root=1) #all values are 1 because non-root nodes ignored

### Graphical Structure of Fitness Components:
### Survival to flowering -> Number of Flowers -> Number of Fruits
pred<-c(0,1,2)

#Distributions to model each fitness component:
### Survival - Bernoulli, Number of Flowers - Truncated Poisson, Number of Fruits - Poisson
fam<-c(1,3,2)
sapply(fam.default(), as.character)[fam]
rbind(vars, fam)

#pred must specify an acyclic graph:
all(pred < seq(along = pred))

## Fitness variable of interest: total number of fruits
fit <- grepl("Num_frts", as.character(redata$varb))
fit <- as.numeric(fit)
redata$fit <- fit

#check
with(redata, sort(unique(as.character(varb)[fit == 0])))
with(redata, sort(unique(as.character(varb)[fit == 1])))

##############################################
################## ASTER MODELS ##############
##############################################

### Fixed-effect models
### Variables of interest are: (Source) Population, Year, Soil type, and their interactions.
### Control for effects of edge position

### Full Model:
aout.full <- aster(resp ~ varb +
                     fit:(Population + Year + SoilType + 
                            Population:SoilType + 
                            Population:Year + 
                            Year:SoilType + 
                            Population:Year:SoilType) +
                     varb:(Edge),
                   pred, fam, varb, id, root, data = redata)

summary(aout.full, show.graph = TRUE) #cannot compute standard errors


#take out 3-way interaction
aout.full2 <- aster(resp ~ varb +
                      fit:(Pop + Year + Soil + 
                             Pop:Soil + Pop:Year + Year:Soil) +
                      varb:(Edge2),
                    pred, fam, varb, id, root, data = redata)
summary(aout.full2)