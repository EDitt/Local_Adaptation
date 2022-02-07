####  Model Testing for ASTER predict

library(dplyr)
library(aster)

### Load redata.RData
load("DataFiles/redata.RData")

##############################################
####### MODEL TESTING FOR PREDICT ############
##############################################

### Largest Model without 3-way interaction
aout1 <- aster(resp ~ varb +
                     varb:(Population + Year + SoilType + 
                            Population:SoilType + 
                            Population:Year + 
                            Year:SoilType + 
                            Edge),
                   pred, fam, varb, id, root, data = redata)

#Pop x Soil with 'fit'
aout.popsoil <- aster(resp ~ varb + fit:Population:SoilType + 
                     varb:(Population + Year + SoilType + 
                             Population:Year + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.popsoil, aout1) #p<0.0001

#Year x Soil with 'fit'
aout.yearsoil <- aster(resp ~ varb + fit:Year:SoilType + 
                     varb:(Population + Year + SoilType + 
                             Population:SoilType + 
                             Population:Year + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.yearsoil, aout1) #p<0.0001

#Pop x Year with 'fit'
aout.popyear <- aster(resp ~ varb + fit:Population:Year + 
                     varb:(Population + Year + SoilType + 
                             Population:SoilType + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.popyear, aout1) #p<0.0001

### A model that allows all interactions to vary with varb (different slopes for relationship between each variable and fitness component) fits the best
### Will use this model "aout1" for predict function

# alternative
aout1 <- aster(resp ~ varb + varb:(Year + SoilType + Edge) +
                 fit:(Population + Population:SoilType + 
                         Population:Year + 
                         Year:SoilType),
               pred, fam, varb, id, root, data = redata)

##############################################
############## ASTER PREDICT #################
##############################################

### Predict number of fruits for each population x soil type in each year
### Will use model that includes all pairwise interactions, but not 3-way interaction (can't compute SEs for this model and interaction not significant)
### Model allows for a different slope between covariate and each fitness component

### Dataframe for predictions
newdata <- expand.grid(Population=levels(redata$Population), SoilType=levels(redata$SoilType), Year=levels(redata$Year))

### "Typical Individuals"
newdata$Edge <- "Non-edge"
#newdata$Edge <- "Edge"

### Add arbitrary values for response vector
for (v in vars)
  newdata[[v]] <- 1

newdata$root <- 1 #add the root

#reshape this dataset
renewdata <- reshape(newdata, varying = list(vars), direction = "long", timevar = "varb", times = as.factor(vars), v.names = "resp")

## Fitness variable to predict: total number of fruits
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

#predicted values
pout <- predict(aout1, newdata = renewdata, varvar = varb, idvar = id, root = root, se.fit = TRUE, amat = amat)

summary(pout) 
bar <-cbind(pout$fit, pout$se.fit)
dimnames(bar) <-list(as.character(newdata$Pop:newdata$Soil:newdata$Year), c("Estimate", "Std. Error"))
print(bar)

ASTER_all_predict <- as.data.frame(bar)


##############################################
################## FIGURE ####################
##############################################

#setEPS
#postscript("TotFit_EachYr_Parents_resp.eps", height=3, width=10)

par(mfrow=c(1,4))

#2012
par(mar=c(3,3,1,0))
plot(c(0.95, 2.2), c(0, 12), col="white", axes=F, frame.plot=F, xlab="", ylab="Lifetime Fitness", cex.lab=1.5, main="")
axis(side=1, at=c(1, 2), labels=c("Sandstone", "Serpentine"), cex=1.5)
axis(side=2, las=1, cex=1.25)

points(c(1,1), c(ASTER_all_predict[1:2,1]), pch=c(21, 15), cex=1.5, col=c("black", "orchid1"))

points(c(2,2), c(ASTER_all_predict[3:4,1]), pch=c(21, 15), cex=1.5, col=c("black", "orchid1"))

segments(c(1, 1), c(ASTER_all_predict[1:2,1]), c(2, 2), c(ASTER_all_predict[3:4,1]), lty=c(2,1), lwd=1.5, col=c("black", "orchid1"))

arrows(c(1,1), c(ASTER_all_predict[1:2,1])-c(ASTER_all_predict[1:2,2]), c(1,1), c(ASTER_all_predict[1:2,1])+c(ASTER_all_predict[1:2,2]), angle=90, length=0.033, code=3, lwd=1, col=c("black", "orchid1"))

arrows(c(2,2), c(ASTER_all_predict[3:4,1])-c(ASTER_all_predict[3:4,2]), c(2,2), c(ASTER_all_predict[3:4,1])+c(ASTER_all_predict[3:4,2]), angle=90, length=0.033, code=3, lwd=1, col=c("black", "orchid1"))

#2013

par(mar=c(3,3,1,0))
plot(c(0.95, 2.2), c(0, 1.2), col="white", axes=F, frame.plot=F, xlab="", ylab="", cex.lab=1.5, main="")
axis(side=1, at=c(1, 2), labels=c("Sandstone", "Serpentine"), cex=1.5)
axis(side=2, las=1, cex=1.25)

points(c(1,1), c(ASTER_all_predict[5:6,1]), pch=c(21, 15), cex=1.5, col=c("black", "orchid1"))

points(c(2,2), c(ASTER_all_predict[7:8,1]), pch=c(21, 15), cex=1.5, col=c("black", "orchid1"))

segments(c(1, 1), c(ASTER_all_predict[5:6,1]), c(2, 2), c(ASTER_all_predict[7:8,1]), lty=c(2,1), lwd=1.5, col=c("black", "orchid1"))

arrows(c(1,1), c(ASTER_all_predict[5:6,1])-c(ASTER_all_predict[5:6,2]), c(1,1), c(ASTER_all_predict[5:6,1])+c(ASTER_all_predict[5:6,2]), angle=90, length=0.033, code=3, lwd=1, col=c("black", "orchid1"))

arrows(c(2,2), c(ASTER_all_predict[7:8,1])-c(ASTER_all_predict[7:8,2]), c(2,2), c(ASTER_all_predict[7:8,1])+c(ASTER_all_predict[7:8,2]), angle=90, length=0.033, code=3, lwd=1, col=c("black", "orchid1"))


#2014
par(mar=c(3,3,1,0))
plot(c(0.95, 2.2), c(0, 12), col="white", axes=F, frame.plot=F, xlab="", ylab="", cex.lab=1.5, main="")
axis(side=1, at=c(1, 2), labels=c("Sandstone", "Serpentine"), cex=1.5)
axis(side=2, las=1, cex=1.25)

points(c(1,1), c(ASTER_all_predict[9:10,1]), pch=c(21, 15), cex=1.5, col=c("black", "orchid1"))

points(c(2,2), c(ASTER_all_predict[11:12,1]), pch=c(21, 15), cex=1.5, col=c("black", "orchid1"))

segments(c(1, 1), c(ASTER_all_predict[9:10,1]), c(2, 2), c(ASTER_all_predict[11:12,1]), lty=c(2,1), lwd=1.5, col=c("black", "orchid1"))

arrows(c(1,1), c(ASTER_all_predict[9:10,1])-c(ASTER_all_predict[9:10,2]), c(1,1), c(ASTER_all_predict[9:10,1])+c(ASTER_all_predict[9:10,2]), angle=90, length=0.033, code=3, lwd=1, col=c("black", "orchid1"))

arrows(c(2,2), c(ASTER_all_predict[11:12,1])-c(ASTER_all_predict[11:12,2]), c(2,2), c(ASTER_all_predict[11:12,1])+c(ASTER_all_predict[11:12,2]), angle=90, length=0.033, code=3, lwd=1, col=c("black", "orchid1"))

#2015
par(mar=c(3,3,1,0))
plot(c(0.95, 2.2), c(0, 3.5), col="white", axes=F, frame.plot=F, xlab="", ylab="", cex.lab=1.5, main="")
axis(side=1, at=c(1, 2), labels=c("Sandstone", "Serpentine"), cex=1.5)
axis(side=2, las=1, cex=1.25)

points(c(1,1), c(ASTER_all_predict[13:14,1]), pch=c(21, 15), cex=1.5, col=c("black", "orchid1"))

points(c(2,2), c(ASTER_all_predict[15:16,1]), pch=c(21, 15), cex=1.5, col=c("black", "orchid1"))

segments(c(1, 1), c(ASTER_all_predict[13:14,1]), c(2, 2), c(ASTER_all_predict[15:16,1]), lty=c(2,1), lwd=1.5, col=c("black", "orchid1"))

arrows(c(1,1), c(ASTER_all_predict[13:14,1])-c(ASTER_all_predict[13:14,2]), c(1,1), c(ASTER_all_predict[13:14,1])+c(ASTER_all_predict[13:14,2]), angle=90, length=0.033, code=3, lwd=1, col=c("black", "orchid1"))

arrows(c(2,2), c(ASTER_all_predict[15:16,1])-c(ASTER_all_predict[15:16,2]), c(2,2), c(ASTER_all_predict[15:16,1])+c(ASTER_all_predict[15:16,2]), angle=90, length=0.033, code=3, lwd=1, col=c("black", "orchid1"))

#dev.off()

##############################################
########## SIG PAIRWISE DIFFERENCES ##########
##############################################

## Posthoc significance tests 
## Difference in predicted number of fruits between source populations on each soil type
## P-values should be corrected for 8 tests (x8 for Bonferroni corrected p-values)

## T-test Statistics
SampleSize <- aggregate(fieldnoNA$PositionColumn, by=list(fieldnoNA$Population, fieldnoNA$SoilType, fieldnoNA$Year), length)
Df.Stats <- cbind(SampleSize, bar)
colnames(Df.Stats) <- c("Pop", "Soil", "Year", "N", "Predict", "SE")

Df.Stats$SD <- Df.Stats$SE * sqrt(Df.Stats$N) #standard deviation
Df.Stats$Var <- Df.Stats$SD^2
Df.Stats$VarxN <- Df.Stats$Var * Df.Stats$N

#reshape for t-test
Df.Stats_wide <- reshape(Df.Stats, idvar = c("Soil", "Year"), timevar = "Pop", direction = "wide")

#calculating T-statistics
Df.Stats_wide$Pop_diff <- Df.Stats_wide$Predict.SerpPop - Df.Stats_wide$Predict.SandPop
Df.Stats_wide$CommonVar <- (Df.Stats_wide$VarxN.SerpPop + Df.Stats_wide$VarxN.SandPop) / (Df.Stats_wide$N.SerpPop + Df.Stats_wide$N.SandPop - 2)
Df.Stats_wide$Tstat_denominator <- sqrt((Df.Stats_wide$CommonVar/Df.Stats_wide$N.SerpPop) + (Df.Stats_wide$CommonVar/Df.Stats_wide$N.SandPop))
Df.Stats_wide$Tvalue <- Df.Stats_wide$Pop_diff / Df.Stats_wide$Tstat_denominator
Df.Stats_wide$df <- Df.Stats_wide$N.SerpPop + Df.Stats_wide$N.SandPop - 2

?pt
pt(Df.Stats_wide$Tvalue, Df.Stats_long$df, lower.tail = FALSE)
rownames(Df.Stats_wide)

### taken from: https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha

t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}


Results <- lapply(row.names(Df.Stats_wide), function(x) {t.test2(
  Df.Stats_wide[x,"Predict.SandPop"], Df.Stats_wide[x,"Predict.SerpPop"],
  Df.Stats_wide[x, "SD.SandPop"], Df.Stats_wide[x,"SD.SerpPop"],
  Df.Stats_wide[x,"N.SandPop"], Df.Stats_wide[x,"N.SerpPop"],
  0)
})
names(Results) <- row.names(Df.Stats_wide)


##############

t.test2(Df.Stats_wide[1,"Predict.SandPop"], Df.Stats_wide[1,"Predict.SerpPop"],
        Df.Stats_wide[1, "SD.SandPop"], Df.Stats_wide[1,"SD.SerpPop"],
        Df.Stats_wide[1,"N.SandPop"], Df.Stats_wide[1,"N.SerpPop"],
        0)

rows <- list(1:8)
row.names(Df.Stats_wide)

Results <- t.test2(Df.Stats_wide$Predict.SandPop, Df.Stats_wide$Predict.SerpPop,
        Df.Stats_wide$SD.SandPop, Df.Stats_wide$SD.SerpPop,
        Df.Stats_wide$N.SandPop, Df.Stats_wide$N.SerpPop,
        0)

Results <- apply(Df.Stats_wide, Margin=1,
                 t.test2(Df.Stats_wide$Predict.SandPop, Df.Stats_wide$Predict.SerpPop,
                         Df.Stats_wide$SD.SandPop, Df.Stats_wide$SD.SerpPop,
                         Df.Stats_wide$N.SandPop, Df.Stats_wide$N.SerpPop,
                         0)
)

Results <- apply(Df.Stats_wide, Margin=1, function(x) {
  t.test2(Df.Stats_wide[,"Predict.SandPop"], Df.Stats_wide[,"Predict.SerpPop"],
          Df.Stats_wide[,"SD.SandPop"], Df.Stats_wide[,"SD.SerpPop"],
          Df.Stats_wide[,"N.SandPop"], Df.Stats_wide[,"N.SerpPop"],
                               0)
})





