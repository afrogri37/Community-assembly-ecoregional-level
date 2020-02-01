
library(geiger)
library(picante)

source("prep_var_ER.R")

tax.mat1 <- as.matrix(tax.mat1)
tax.mat2<-as.matrix(tax.mat2)

# ER centroid data
centr <- read.csv2(file = "external_files/ecoreg_centroids.csv", header = TRUE, dec = ",")
lat <- centr$lat 
lon <- centr$lon


# TAXONOMIC ALPHA DIVERSITY
# Calculating taxonomic richness of each site  
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



# PRUNE
prunedphy <- prune.sample(tax.mat1, my.phylo)


# PHYLOGENETIC ALPHA DIVERSITY 

# All sps that exist in phylogeny
tax.mat3 <- match.phylo.comm(my.phylo,tax.mat1)$comm

# Tree only with sps that have traits
phylo.traits <- match.phylo.comm(my.phylo, tax.mat2)$phy                       

# Mean pairwise distance (MPD) within an ecoregion
mpd_ER <- mpd(tax.mat3, cophenetic(my.phylo))


# How closely related is the average pair of species or individuals in a community?
# Phylogenetic community structure

# SES.MPD
phydist <- cophenetic(prunedphy)
ses.mpd.result1 <- ses.mpd(samp=tax.mat3, dis=phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)
ses.mpd.result2 <- ses.mpd(samp=tax.mat3, dis=phydist, null.model = "sample.pool", abundance.weighted = FALSE, runs = 999)
ses.mpd.result3 <- ses.mpd(samp=tax.mat3, dis=phydist, null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
ses.mpd.result4 <- ses.mpd(samp=tax.mat3, dis=phydist, null.model = "independentswap", abundance.weighted = FALSE, runs = 999)
ses.mpd.result1
ses.mpd.result2
ses.mpd.result3
ses.mpd.result4

# mpd.obs.z>0 i mpd.obs.p>0.95 --> phylogenetic evenness
# mpd.obs.z<0 i mpd.obs.p<0.05 --> phylogenetic clustering


tax.lab <- cbind(ses.mpd.result1, lat)
ind.swap <- cbind(ses.mpd.result4, lat)

# Observed Vs Random MPD Graph
plot(mpd.rand.mean~lat, data=tax.lab, ylim=c(1.2,1.4), pch=1, col="violetred", ylab="MPD", yaxt="n", xaxt="n")
axis(1,c(40,50,60))
axis(2, c(1.20,1.25,1.30,1.35,1.40), las=2)
points(mpd.obs~lat, data=tax.lab, col="violetred", pch=16)
abline(mod_mpd, lwd=3, col="violetred")
legend("bottomleft", c("mpd obs", "mpd rand"), col=c("violetred","violetred" ), pch=c(16,1), bty="o", y.intersp=0.5)
r2 <- bquote(italic(r)^2 == .(round(summary(mod_mpd)$r.squared,2)))
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p<0.001", line =-2.3, adj = 0.9, cex = 1.2, font = 2)


plot(mpd.obs.z~lat, data=tax.lab, ylim=c(-4.5,4), pch= 16, col="violetred", 
     ylab="SES.MPD", xlab="latitude", axes=T, las=1)
points(mpd.rand.mean~lat, data=tax.lab, pch=1, col="violetred")
text(mpd.obs.z~lat, data=tax.lab, labels=centr$index, cex=0.6,pos=3)
abline(h=0, lty="dashed")

# ses.mpd~lat correlation
cor.test(tax.lab$lat, tax.lab$mpd.obs.z, method="spearman") # p-value<0.001
cor.test(ind.swap$lat, ind.swap$mpd.obs.z, method="spearman") # p-value<0.001


# COMPARATIVE ANALYSES

# Phylogenetic signal

# Blomberg's K
traitsK <- funspacetr[my.phylo$tip.label,]
K <- multiPhylosignal(traitsK, my.phylo)

# K = 1 --> some degree of conservatism.
# K = 0 --> random pattern of evolution. trait variation has not followed the phylogeny.
# K > 1 --> strong phylogenetic signal and conservatism of traits.


# Pagel's lambda
LambdaPagel <- fitContinuous(phy=prunedphy, dat=funspacetr, model="lambda")

