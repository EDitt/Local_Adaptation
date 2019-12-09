####  Model Testing for ASTER predict

library(dplyr)
library(aster)

### Largest Model without 3-way interaction
aout.full <- aster(resp ~ varb +
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
anova(aout.popsoil, aout.full) #p<0.0001

aout.yearsoil <- aster(resp ~ varb + fit:Year:SoilType + 
                     varb:(Population + Year + SoilType + 
                             Population:SoilType + 
                             Population:Year + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.yearsoil, aout.full) #p<0.0001

aout.popyear <- aster(resp ~ varb + fit:Population:Year + 
                     varb:(Population + Year + SoilType + 
                             Population:SoilType + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.popyear, aout.full) #p<0.0001

aout.SoilType <- aster(resp ~ varb + fit:SoilType +
                     varb:(Population + Year + 
                             Population:SoilType + 
                             Population:Year + 
                             Year:SoilType + 
                             Edge),
                   pred, fam, varb, id, root, data = redata)
anova(aout.SoilType, aout.full) #ns

aout.Population <- aster(resp ~ varb + fit:Population +
                         varb:(SoilType + Year + 
                                 Population:SoilType + 
                                 Population:Year + 
                                 Year:SoilType + 
                                 Edge),
                       pred, fam, varb, id, root, data = redata)
anova(aout.Population, aout.full) #ns

aout.Year <- aster(resp ~ varb + fit:Year + 
                           varb:(SoilType + Population +
                                   Population:SoilType + 
                                   Population:Year + 
                                   Year:SoilType + 
                                   Edge),
                         pred, fam, varb, id, root, data = redata)
anova(aout.Year, aout.full)

aout.Edge <- aster(resp ~ varb + fit:Edge + 
                     varb:(SoilType + Population + Year +
                             Population:SoilType + 
                             Population:Year + 
                             Year:SoilType),
                   pred, fam, varb, id, root, data = redata)
anova(aout.Edge, aout.full) #ns

#should I include main effects interaction with fit?
Mod1 <- aster(resp ~ varb +
                   varb:(Population + Year + SoilType + 
                           Population:SoilType + 
                           Population:Year + 
                           Year:SoilType + 
                           Edge),
                 pred, fam, varb, id, root, data = redata)

Mod2 <- aster(resp ~ varb + fit:(Population + Year + SoilType + Edge) +
                varb:(Population:SoilType + 
                        Population:Year + 
                        Year:SoilType),
              pred, fam, varb, id, root, data = redata)

anova(Mod2, Mod1) #ns

pout <- predict(Mod2, newdata = renewdata, varvar = varb, idvar = id, root = root, se.fit = TRUE, amat = amat)
pout <- predict(aout.full, newdata = renewdata, varvar = varb, idvar = id, root = root, se.fit = TRUE, amat = amat)

summary(pout) 
bar <-cbind(pout$fit, pout$se.fit)
dimnames(bar) <-list(as.character(newdata$Pop:newdata$Soil:newdata$Year), c("Estimate", "Std. Error"))
print(bar)

ASTER_all_predict <- as.data.frame(bar)


####Plot

setEPS
postscript("TotFit_EachYr_Parents_resp.eps", height=3, width=10)

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

