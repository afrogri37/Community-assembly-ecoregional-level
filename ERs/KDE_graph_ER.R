library(FD)
library(ks)

source("prep_var_ER.R")

# Use dudi.pco for PCA for being consistent with scaling of scores in analysis 
pco12 <- pcodistrait$li[,1:2]
ll <- pcodistrait$loadings

# Kernel Density Estimation
Ho12 <- Hpi(x=pco12)      # optimal bandwidth estimation
esto12 <- kde(x=pco12, H=Ho12, compute.cont=TRUE)     # kernel density estimation

# Set contour probabilities for drawing contour levels
clo12 <- contourLevels(esto12, prob=c(0.5, 0.05, 0.001), approx=TRUE)

# Setting arrows based on correlations
# pc1:gil,female.max,shr,ae.pas, aq_pas, ae_act
# pc2: gra, pre,isoleg_cem
taball12 <- subset(taball, select = c(gil,female.max,shr,ae.pas, aq_pas, ae_act, gra, pre,isoleg_cem ))
fit12 <- envfit(pco12, taball12)  # use envfit to draw arrows, can be also done using trait loadings
fit122 <- fit12$vectors$arrows*-1 # drawing line segments in arrow opposites direction for pretty layout

# Plot Kernel Density Estimations
pdf("results/KDE_PCoA_biol.pdf", width=5, height=5)
par(mfrow=c(1,1))

# pc1 and pc2
plot(esto12, cont=seq(1,100,by=1), display="filled.contour2", add=FALSE, ylab="PC2", xlab="PC1", cex.axis=0.75, ylim=c(-1, 1), xlim=c(-1, 1),las=1) 
plot(esto12, abs.cont=clo12[1], labels=c(0.5),labcex=0.75, add=TRUE, lwd=0.75, col="grey30")
plot(esto12, abs.cont=clo12[2], labels=c(0.95),labcex=0.75, add=TRUE, lwd=0.5, col="grey60")
plot(esto12, abs.cont=clo12[3], labels=c(0.99),labcex=0.75, add=TRUE, lwd=0.5, col="grey60")
points(pco12[,], pch=16, cex=0.25, col="black") 
plot(fit12, cex=0.90, col=1)
dev.off()
