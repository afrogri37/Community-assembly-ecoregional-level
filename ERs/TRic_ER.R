
source("prep_var_ER.R")

# ER centroid data
centr <- read.csv2(file = "external_files/ecoreg_centroids.csv", header = TRUE, dec = ",")
lat <- centr$lat 
lon <- centr$lon

##TD with 1425 sps   
TRic <- rowSums(tax.mat1)
mod_TRic<-lm(TRic~lat)
plot(TRic~lat, pch=16, ylim=c(150,450), xlim=c(35,70), cex=0.8)
abline(mod_TRic)
r2 <- bquote(italic(r)^2 == .(round(summary(mod_TRic)$r.squared,2)))
p <- bquote(italic(p) == .(round(summary(mod_TRic)$coefficients[2,4])))
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext(p, line = -0.1, adj = 0.9, cex = 1.2, font = 2)
text(TRic~lat, labels = centr$index, pos = 2, cex=0.6)


# model summary and residuals
summary(mod_TRic)
par(mfrow=c(2,2))
plot(mod_TRic)
par(mfrow=c(1,1))

