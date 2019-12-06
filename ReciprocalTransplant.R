####  Analysis of Field Reciprocal Transplant Data
####  Dataset: "ReciprocalTransplant.csv"
####  Fitness data collected from 2 soil types (serpentine and sandstone) from 2012-2015

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
##### ASTER MODELS - FULL DATASET ############
##############################################

### Fixed-effect models
### Variables of interest are: (Source) Population, Year, Soil type, and their interactions.
### Control for effects of edge position

### First, will test significance of 3-way interaction between Population x Soil x Year:

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
names(aout.full$coefficients)

### Why can't it compute standard errors?
info.tol <- sqrt(.Machine$double.eps)
infomat <- aout.full$fisher
fred <- eigen(infomat, symmetric = TRUE)
sally <- fred$values < max(fred$values) * info.tol
apparent <- zapsmall(fred$vectors[ , sally])
rownames(apparent) <- names(aout.full$coefficients)
apparent
### It looks like many of the interactions e.g. Pop x Soil and Year x Soil perfectly predict Pop x Year x Soil
### My guess is that this is due to the very low survival of the sandstone soil on serpentine soil in all four years
### Further, this interaction is n.s. (see below)

aout1 <- aster(resp ~ varb +
                     fit:(Population + Year + SoilType + 
                            Population:SoilType + 
                            Population:Year + 
                            Year:SoilType) + 
                            #Population:Year:SoilType) +
                            varb:(Edge),
                   pred, fam, varb, id, root, data = redata)
summary(aout1)

anova(aout1, aout.full) #ns (p=0.3753)

### Test Pop x Soil interaction on total fitness:

#allow all factors to interact with varb
aout.full2 <- aster(resp ~ varb +
                     varb:(Population + Year + SoilType + 
                            Population:SoilType + 
                            Population:Year + 
                            Year:SoilType + 
                            Population:Year:SoilType +
                     Edge),
                   pred, fam, varb, id, root, data = redata)

### which factors have significant interactions with varb?

aout.Pop <- aster(resp ~ varb + fit:(Population) +
                      varb:(Year + SoilType + 
                              Population:SoilType + 
                              Population:Year + 
                              Year:SoilType + 
                              Population:Year:SoilType +
                              Edge),
                    pred, fam, varb, id, root, data = redata)

anova(aout.Pop, aout.full2) #p<0.0001

aout.Year <- aster(resp ~ varb + fit:(Year) +
                      varb:(Population + SoilType + 
                              Population:SoilType + 
                              Population:Year + 
                              Year:SoilType + 
                              Population:Year:SoilType +
                              Edge),
                    pred, fam, varb, id, root, data = redata)

anova(aout.Year, aout.full2) #p<0.0001

aout.Soil <- aster(resp ~ varb + fit:(SoilType) +
                      varb:(Population + Year +
                              Population:SoilType + 
                              Population:Year + 
                              Year:SoilType + 
                              Population:Year:SoilType +
                              Edge),
                    pred, fam, varb, id, root, data = redata)

anova(aout.Soil, aout.full2) #p<0.0001

aout.PopSoil <- aster(resp ~ varb + fit:(Population:SoilType) +
                      varb:(Population + Year + SoilType + 
                              Population:Year + 
                              Year:SoilType + 
                              Population:Year:SoilType +
                              Edge),
                    pred, fam, varb, id, root, data = redata)

anova(aout.PopSoil, aout.full2) #ns

aout.PopYear <- aster(resp ~ varb + fit:(Population:Year) +
                        varb:(Population + Year + SoilType + 
                                Population:SoilType + 
                                Year:SoilType + 
                                Population:Year:SoilType +
                                Edge),
                      pred, fam, varb, id, root, data = redata)
anova(aout.PopYear, aout.full2) #ns

aout.SoilYear <- aster(resp ~ varb + fit:(Year:SoilType) +
                        varb:(Population + Year + SoilType + 
                                Population:SoilType + 
                                Population:Year + 
                                Population:Year:SoilType +
                                Edge),
                      pred, fam, varb, id, root, data = redata)
anova(aout.SoilYear, aout.full2) #p<0.0001

aout.PopSoilYear <- aster(resp ~ varb + fit:(Population:Year:SoilType) +
                      varb:(Population + Year + SoilType + 
                              Population:SoilType + 
                              Population:Year + 
                              Year:SoilType + 
                              Edge),
                    pred, fam, varb, id, root, data = redata)
anova(aout.PopSoilYear, aout.full2) #p=0.09244

aout.edge <- aster(resp ~ varb + fit:(Edge) +
                      varb:(Population + Year + SoilType + 
                              Population:SoilType + 
                              Population:Year + 
                              Year:SoilType + 
                              Population:Year:SoilType),
                    pred, fam, varb, id, root, data = redata)
anova(aout.edge, aout.full2) #ns

#####
#influence of Population, Soil, and interaction on fitness
aout_PopSoil_full <- aster(resp ~ varb +
                 fit:(Population + SoilType + 
                        Population:SoilType) + 
                 varb:(Year + 
                         Population:Year + 
                         Year:SoilType +
                         Edge),
               pred, fam, varb, id, root, data = redata)

aout_PopSoil_red1 <- aster(resp ~ varb +
                             fit:(Population:SoilType) + 
                             varb:(Population + SoilType + Year + 
                                     Population:Year + 
                                     Year:SoilType +
                                     Edge),
                           pred, fam, varb, id, root, data = redata)
anova(aout_PopSoil_full, aout_PopSoil_red1) #p<0.0001

aout_PopSoil_red2 <- aster(resp ~ varb +
                             fit:(Population + SoilType) + 
                                    #Population:SoilType) + 
                             varb:(Year + 
                                     Population:Year + 
                                     Year:SoilType +
                                     Edge),
                           pred, fam, varb, id, root, data = redata)
anova(aout_PopSoil_red, aout_PopSoil_full, aout.full2) #p<0.0001


aoutsub1 <- aster(resp ~ varb + fit:(SoilType) +
                 varb:(Population + Year + 
                         Population:SoilType + 
                         Population:Year + 
                         Year:SoilType + 
                         #Population:Year:SoilType) +
                         Edge),
               pred, fam, varb, id, root, data = redata)


####
aout_PopSoil <- aster(resp ~ varb +
                 fit:(Population + Year + SoilType + 
                        #Population:SoilType + 
                        Population:Year + 
                        Year:SoilType) + 
                 varb:(Edge),
               pred, fam, varb, id, root, data = redata)

aout_PopYear <- aster(resp ~ varb +
                        fit:(Population + Year + SoilType + 
                               Population:SoilType + 
                               #Population:Year + 
                               Year:SoilType) + 
                        varb:(Edge),
                      pred, fam, varb, id, root, data = redata)

aout_SoilYear <- aster(resp ~ varb +
                        fit:(Population + Year + SoilType + 
                               Population:SoilType + 
                               Population:Year) + 
                               #Year:SoilType) + 
                        varb:(Edge),
                      pred, fam, varb, id, root, data = redata)


anova(aout_PopSoil, aout1) #pop x soil is significant p<0.0001
anova(aout_PopYear, aout1) #pop x year is significant p=0.003263
anova(aout_SoilYear, aout1) #year x soil is significant p=0.004544

# Test main effect of Edge because it is not involved in any significant interactions:
aout_Edge <- aster(resp ~ varb +
                         fit:(Population + Year + SoilType + 
                                Population:SoilType + 
                                Population:Year + 
                         Year:SoilType),
                         #varb:(Edge),
                       pred, fam, varb, id, root, data = redata)
anova(aout_Edge, aout1) #edge is significant p<0.0001

#### NOTE: I did the same model selection using Plot rep (nested within year and soil type) as a random effect using the reaster function
### The only difference in significance was in the Year x Soil interaction which was not significant in these random effect models
### However, because there was only 1 plot per soil type in 2012 and 2013, modelling this as a random effect confounds this interaction

### To incorporate variation across plot replicates, I will test Pop x Soil interactions separately for each year and 
### use random effect models in 2014-2015 to model variation in plot replicate

#put main effects in varb?
aout.test <- aster(resp ~ varb +
                 fit:(Population:SoilType + 
                        Population:Year + 
                        Year:SoilType) + 
                 #Population:Year:SoilType) +
                 varb:(Edge + Population + Year + SoilType),
               pred, fam, varb, id, root, data = redata)
anova(aout1, aout.test) #p<0.0001

aout.big <- aster(resp ~ varb +
                     fit:(Population:SoilType + 
                            Population:Year + 
                            Year:SoilType) + 
                     varb:(Edge + Population + Year + SoilType),
                   pred, fam, varb, id, root, data = redata)

##Biggest model?
aout.full <- aster(resp ~ varb +
                     varb:(Population + Year + SoilType + 
                             Population:SoilType + 
                             Population:Year + 
                             Year:SoilType + 
                             Population:Year:SoilType +
                             Edge),
                   pred, fam, varb, id, root, data = redata)
summary(aout.full) #cannot compute SEs
aout.red1 <- aster(resp ~ varb +
                     fit:(Population:Year:SoilType) +
                     varb:(Population + Year + SoilType + 
                             Population:SoilType + 
                             Population:Year + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.red1, aout.full) #p=0.09244

### Remove 3-way interaction
aout.red2 <- aster(resp ~ varb +
                     varb:(Population + Year + SoilType + 
                             Population:SoilType + 
                             Population:Year + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)

aout.red3 <- aster(resp ~ varb + fit:(Population:SoilType) +
                     varb:(Population + Year + SoilType + 
                             Population:SoilType + 
                             Population:Year + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.red3, aout.red2) #ns - Pop x Soil does not have different effects in different layers of graph

aout.red4 <- aster(resp ~ varb + fit:(Population:Year) +
                     varb:(Population + Year + SoilType + 
                             Population:SoilType + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.red4, aout.red2) #p<0.0001 - Pop x Year

aout.red5 <- aster(resp ~ varb + fit:(Year:SoilType) +
                     varb:(Population + Year + SoilType + 
                             Population:SoilType + 
                             Population:Year + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.red5, aout.red2) #p<0.0001 - Year x Soil

aout.red6 <- aster(resp ~ varb + fit:(Population) +
                     varb:(Year + SoilType + 
                             Population:SoilType + 
                             Population:Year + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.red6, aout.red2) #ns

aout.red7 <- aster(resp ~ varb + fit:(Year) +
                     varb:(Population + SoilType + 
                             Population:SoilType + 
                             Population:Year + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.red7, aout.red2) #ns

aout.red8 <- aster(resp ~ varb + fit:(SoilType) +
                     varb:(Population + Year + 
                             Population:SoilType + 
                             Population:Year + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.red8, aout.red2) #ns

##############################################
########### ASTER MODELS - BY YEAR ###########
##############################################


redata12 <- subset(redata, Year == "2012")
redata13 <- subset(redata, Year == "2013")
redata14 <- subset(redata, Year == "2014")
redata15 <- subset(redata, Year == "2015")

##In 2012-2013, there was no plot rep so will not include this factor

######  2012
aout.full12 <- aster(resp ~ varb +
                       fit:(Population + SoilType + 
                              Population:SoilType) +
                       varb:Edge,
                     pred, fam, varb, id, root, data = redata12)

summary(aout.full12) #cannot compute SEs
aout12_inter <- aster(resp ~ varb +
                        fit:(Population + SoilType) + 
                               #Population:SoilType) +
                        varb:Edge,
                      pred, fam, varb, id, root, data = redata12)

summary(aout12_inter)
anova(aout12_inter, aout.full12) #p=0.0002836 interaction significant

######  2013
aout.full13 <- aster(resp ~ varb +
                       fit:(Population + SoilType + 
                              Population:SoilType) +
                       varb:Edge,
                     pred, fam, varb, id, root, data = redata13)

summary(aout.full13)

aout13_inter <- aster(resp ~ varb +
                        fit:(Population + SoilType) + 
                        #Population:SoilType) +
                        varb:Edge,
                      pred, fam, varb, id, root, data = redata13)

summary(aout13_inter)
anova(aout13_inter, aout.full13) #p=0.001185

##In 2014-2015, there were multiple plots per soil type so will use reaster to test the pop x soil interaction with plot replicate as a random effect

######  2014

rout14.full <- reaster(resp ~ varb +
                         fit:(Population + SoilType + 
                                Population:SoilType) +
                         varb:Edge,
                       list(block = ~ 0 + fit:SoilType:Plot_Rep),
                       pred, fam, varb, id, root, data = redata14)
summary(rout14.full)
rout14.inter <- reaster(resp ~ varb +
                         fit:(Population + SoilType) + 
                                #Population:SoilType) +
                         varb:Edge,
                       list(block = ~ 0 + fit:SoilType:Plot_Rep),
                       pred, fam, varb, id, root, data = redata14)

anova(rout14.inter, rout14.full) #p<0.0001

######  2015

rout15.full <- reaster(resp ~ varb +
                         fit:(Population + SoilType + 
                                Population:SoilType) +
                         varb:Edge,
                       list(block = ~ 0 + fit:SoilType:Plot_Rep),
                       pred, fam, varb, id, root, data = redata15)
summary(rout15.full) #warning message

rout15.inter <- reaster(resp ~ varb +
                          fit:(Population + SoilType) + 
                          #Population:SoilType) +
                          varb:Edge,
                        list(block = ~ 0 + fit:SoilType:Plot_Rep),
                        pred, fam, varb, id, root, data = redata15)
summary(rout15.inter)
anova(rout15.inter, rout15.full) #p<0.0001

##############################################
################ ASTER PREDICT ###############
##############################################

### Predict number of fruits for each population x soil type in each year
### Will use model that includes all pairwise interactions, but not 3-way interaction (can't compute SEs for this model and interaction not significant)

### Dataframe for predictions
newdata <- expand.grid(Population=levels(redata$Population), SoilType=levels(redata$SoilType), Year=levels(redata$Year))

### "Typical Individuals"
newdata$Edge <- "Non-edge"

### Add arbitrary values for response vector
for (v in vars)
  newdata[[v]] <- 1

newdata$root <- 1 #add the root

#now reshape this dataset
renewdata <- reshape(newdata, varying = list(vars), direction = "long", timevar = "varb", times = as.factor(vars), v.names = "resp")

## Fitness variable of interest: total number of fruits
fit <- grepl("Num_frts", as.character(renewdata$varb))
fit <- as.numeric(fit)
renewdata$fit <- fit

#check
with(renewdata, sort(unique(as.character(varb)[fit == 0])))
with(renewdata, sort(unique(as.character(varb)[fit == 1])))

nind <-nrow (newdata) #num rows of df
nnode <- length (vars) #3 variables
amat <-array(0, c(nind, nnode, nind)) #this makes an empty array
dim(amat)

#set the components to be one for variable of interest (number of fruits) and zero for the others
for (i in 1:nind)
  amat[i, grep("Num_frts", vars), i] <- 1


pout <- predict(aout1, newdata = renewdata, varvar = varb, idvar = id, root = root, se.fit = TRUE, amat = amat)
pout$fit


pout <- predict(aout.test, newdata = renewdata, varvar = varb, idvar = id, root = root, se.fit = TRUE, amat = amat)

summary(pout) 
bar <-cbind(pout$fit, pout$se.fit)
dimnames(bar) <-list(as.character(newdata$Pop:newdata$Soil:newdata$Year), c("Estimate", "Std. Error"))
print(bar)

#### Will take out later

ASTER_all_predict <- as.data.frame(bar)

head(ASTER_all_predict)

par(mfrow=c(1,4))

#2012
par(mar=c(3,3,1,0))
plot(c(0.95, 2.2), c(0, 12), col="white", axes=F, frame.plot=F, xlab="", ylab="Lifetime Fitness", cex.lab=1.5, main="")
axis(side=1, at=c(1, 2), labels=c("Sandstone soil", "Serpentine soil"), cex=1.5)
axis(side=2, las=1, cex=1.25)

points(c(1,1), c(ASTER_all_predict[1:2,1]), pch=c(15, 21), cex=1.5, col=c("orchid1", "black"))

points(c(2,2), c(ASTER_all_predict[3:4,1]), pch=c(15, 21), cex=1.5, col=c("orchid1", "black"))

segments(c(1, 1), c(ASTER_all_predict[1:2,1]), c(2, 2), c(ASTER_all_predict[3:4,1]), lty=c(1,2), lwd=1.5, col=c("orchid1", "black"))

arrows(c(1,1), c(ASTER_all_predict[1:2,1])-c(ASTER_all_predict[1:2,2]), c(1,1), c(ASTER_all_predict[1:2,1])+c(ASTER_all_predict[1:2,2]), angle=90, length=0.033, code=3, lwd=1, col=c("orchid1", "black"))

arrows(c(2,2), c(ASTER_all_predict[3:4,1])-c(ASTER_all_predict[3:4,2]), c(2,2), c(ASTER_all_predict[3:4,1])+c(ASTER_all_predict[3:4,2]), angle=90, length=0.033, code=3, lwd=1, col=c("orchid1", "black"))

#2013

par(mar=c(3,3,1,0))
plot(c(0.95, 2.2), c(0, 1.2), col="white", axes=F, frame.plot=F, xlab="", ylab="", cex.lab=1.5, main="")
axis(side=1, at=c(1, 2), labels=c("Sandstone soil", "Serpentine soil"), cex=1.5)
axis(side=2, las=1, cex=1.25)

points(c(1,1), c(ASTER_all_predict[5:6,1]), pch=c(15, 21), cex=1.5, col=c("orchid1", "black"))

points(c(2,2), c(ASTER_all_predict[7:8,1]), pch=c(15, 21), cex=1.5, col=c("orchid1", "black"))

segments(c(1, 1), c(ASTER_all_predict[5:6,1]), c(2, 2), c(ASTER_all_predict[7:8,1]), lty=c(1,2), lwd=1.5, col=c("orchid1", "black"))

arrows(c(1,1), c(ASTER_all_predict[5:6,1])-c(ASTER_all_predict[5:6,2]), c(1,1), c(ASTER_all_predict[5:6,1])+c(ASTER_all_predict[5:6,2]), angle=90, length=0.033, code=3, lwd=1, col=c("orchid1", "black"))

arrows(c(2,2), c(ASTER_all_predict[7:8,1])-c(ASTER_all_predict[7:8,2]), c(2,2), c(ASTER_all_predict[7:8,1])+c(ASTER_all_predict[7:8,2]), angle=90, length=0.033, code=3, lwd=1, col=c("orchid1", "black"))


#2014
par(mar=c(3,3,1,0))
plot(c(0.95, 2.2), c(0, 12), col="white", axes=F, frame.plot=F, xlab="", ylab="", cex.lab=1.5, main="")
axis(side=1, at=c(1, 2), labels=c("Sandstone soil", "Serpentine soil"), cex=1.5)
axis(side=2, las=1, cex=1.25)

points(c(1,1), c(ASTER_all_predict[9:10,1]), pch=c(15, 21), cex=1.5, col=c("orchid1", "black"))

points(c(2,2), c(ASTER_all_predict[11:12,1]), pch=c(15, 21), cex=1.5, col=c("orchid1", "black"))

segments(c(1, 1), c(ASTER_all_predict[9:10,1]), c(2, 2), c(ASTER_all_predict[11:12,1]), lty=c(1,2), lwd=1.5, col=c("orchid1", "black"))

arrows(c(1,1), c(ASTER_all_predict[9:10,1])-c(ASTER_all_predict[9:10,2]), c(1,1), c(ASTER_all_predict[9:10,1])+c(ASTER_all_predict[9:10,2]), angle=90, length=0.033, code=3, lwd=1, col=c("orchid1", "black"))

arrows(c(2,2), c(ASTER_all_predict[11:12,1])-c(ASTER_all_predict[11:12,2]), c(2,2), c(ASTER_all_predict[11:12,1])+c(ASTER_all_predict[11:12,2]), angle=90, length=0.033, code=3, lwd=1, col=c("orchid1", "black"))

#2015
par(mar=c(3,3,1,0))
plot(c(0.95, 2.2), c(0, 3.5), col="white", axes=F, frame.plot=F, xlab="", ylab="", cex.lab=1.5, main="")
axis(side=1, at=c(1, 2), labels=c("Sandstone soil", "Serpentine soil"), cex=1.5)
axis(side=2, las=1, cex=1.25)

points(c(1,1), c(ASTER_all_predict[13:14,1]), pch=c(15, 21), cex=1.5, col=c("orchid1", "black"))

points(c(2,2), c(ASTER_all_predict[15:16,1]), pch=c(15, 21), cex=1.5, col=c("orchid1", "black"))

segments(c(1, 1), c(ASTER_all_predict[13:14,1]), c(2, 2), c(ASTER_all_predict[15:16,1]), lty=c(1,2), lwd=1.5, col=c("orchid1", "black"))

arrows(c(1,1), c(ASTER_all_predict[13:14,1])-c(ASTER_all_predict[13:14,2]), c(1,1), c(ASTER_all_predict[13:14,1])+c(ASTER_all_predict[13:14,2]), angle=90, length=0.033, code=3, lwd=1, col=c("orchid1", "black"))

arrows(c(2,2), c(ASTER_all_predict[15:16,1])-c(ASTER_all_predict[15:16,2]), c(2,2), c(ASTER_all_predict[15:16,1])+c(ASTER_all_predict[15:16,2]), angle=90, length=0.033, code=3, lwd=1, col=c("orchid1", "black"))



