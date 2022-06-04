####  Model Testing for ASTER predict

library(dplyr)
library(aster)
source("Functions.R")

### Load redata.RData
load("DataFiles/redata.RData")

##############################################
########### MODEL FOR PREDICT ###############
##############################################

aout1 <- aster(resp ~ varb:Edge +
                 fit:(Population + Year + SoilType + 
                         Population:SoilType + 
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
newdata <- expand.grid(Population=levels(redata$Population), SoilType=levels(as.factor(redata$SoilType)), Year=levels(redata$Year))

### "Typical Individuals"
newdata$Edge <- "Non-edge"

### using function
pout <- ASTER_predict(newdata, "Num_frts", aout1)
summary(pout) 
bar <-cbind(pout$fit, pout$se.fit)
dimnames(bar) <-list(as.character(newdata$Pop:newdata$Soil:newdata$Year), c("Estimate", "Std. Error"))
print(bar)

ASTER_all_predict <- as.data.frame(bar)

##############################################
################## FIGURE ####################
##############################################

setEPS
postscript("/Volumes/GoogleDrive/My Drive/Active Projects/Leptosiphon manuscripts/DRAFTS/Part1_LocalAdaptation/Submitting/Revisions/NewFigsTables/TotFit_EachYr_Parents_resp.eps", 
           height=3, width=10)

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

dev.off()

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
Df.Stats_long <- reshape(Df.Stats, idvar = c("Soil", "Year"), timevar = "Pop", direction = "wide")

#calculating T-statistics
Df.Stats_long$Pop_diff <- Df.Stats_long$Predict.SerpPop - Df.Stats_long$Predict.SandPop
Df.Stats_long$CommonVar <- (Df.Stats_long$VarxN.SerpPop + Df.Stats_long$VarxN.SandPop) / (Df.Stats_long$N.SerpPop + Df.Stats_long$N.SandPop - 2)
Df.Stats_long$Tstat_denominator <- sqrt((Df.Stats_long$CommonVar/Df.Stats_long$N.SerpPop) + (Df.Stats_long$CommonVar/Df.Stats_long$N.SandPop))
Df.Stats_long$Tvalue <- Df.Stats_long$Pop_diff / Df.Stats_long$Tstat_denominator
Df.Stats_long$df <- Df.Stats_long$N.SerpPop + Df.Stats_long$N.SandPop - 2


