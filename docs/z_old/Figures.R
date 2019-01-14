#Excalibur
frame <- read.csv("ExcaliburGRa.csv", header=T)
head(frame)
modGR1 <- drm(GR ~ Psi, fct=GR501, data=frame, curveid=Perc)#,
              #pmodel=list(~Perc-1, ~1))
summary(modGR1)
plot(modGR1, log="", legendPos=c(-0.9, 1))

modGR2 <- drm(GR ~ Psi, fct=GR503, data=frame, curveid=Perc)
summary(modGR2)
plot(modGR2, log="", legendPos=c(-0.9, 1))
frame$PercF <- factor(frame$Perc)
modGR2b <- drm(GR ~ Psi, fct=GR503, data=frame, curveid=Perc, pmodels=list( ~1, ~PercF-1))
anova(modGR2, modGR2b)
plot(modGR2b, log="", legendPos=c(-0.9, 1))

modGR3 <- drm(y ~ x, fct=GR503, data=frame)
summary(modGR3)

#Festuca arundinacea
y <- c(0, 0, 0, 0, 0.0620803286287862, 0.0952596286952884, 0.124027304016132,
       0.135454008633976, 0.142496108773546, 0.145828499965942, 0.146318563838286,
       0.147403638789466)
x <- c(-2, -1.5, -1.2, -1, -0.8, -0.6, -0.4, -0.25, -0.12, -0.06,
               -0.03, 0)
modGR502 <- drm(y ~ x, fct=GR502)
summary(modGR502)

par(mfrow=c(1,2))
plot(modGR502, log="", xlab="Water potential (MPa)", ylab="GR50", axes=F)
axis(1, at=seq(-2,0,by=0.5), labels=seq(-2,0,by=0.5))
axis(2,at=seq(0, 0.15, by=0.05), labels=seq(0, 0.15, by=0.05))
plot(modGR3, log="", xlab="Water potential (MPa)", ylab="GR50", axes=F, ylim=c(0,1))
axis(1, at=seq(-1.5,0,by=0.5), labels=seq(-1.5,0,by=0.5))
axis(2,at=seq(0, 1, by=0.25), labels=seq(0, 1, by=0.25))
	
lab <- paste(levels(factor(f))[i], "°C", sep="")
	text (1,0.9, labels=lab, cex=1.0, font=3)
	lines(p ~ x, subset=c(f==levels(factor(f))[i]), col="red")
}

mtext("Tempo (d)", las=1, side=1, outer=T, at=c(0.5), padj=2.3, cex=1.0)
mtext("Proporzione di semi germinati", outer=T, side=2, las=0, at=c(0.5), padj=-2.3, cex=1.0)



