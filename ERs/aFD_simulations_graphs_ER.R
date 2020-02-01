library(FD)

# ER centroid data
centr <- read.csv2(file = "external_files/ecoreg_centroids.csv", header = TRUE, dec = ",")
lat <- centr$lat 
lon <- centr$lon

# Trying to find precompiled simulations. If not found, run simulations.
generate_simulations <- tryCatch(load("precompiled/FD_indexes.RData"), error=function(c){""})
if(generate_simulations == ""){
  
  library(picante)
  source("external_files/FD_functions.R")
  source("prep_var_ER.R")
  
  # SIMULATIONS TO CALCULATE FD indexes
  
  n.reps <- 999
  n.spp <- nrow(funspace)
  n.sample <- min(rowSums(tax.mat2))  
  
  FRic.out <- matrix(0, ncol=nrow(tax.mat2),nrow=n.reps)
  FDis.out <- matrix(0, ncol=nrow(tax.mat2),nrow=n.reps)
  
  for (i in 0:n.reps) {
    data.temp <- array(0, c(nrow(funspace),ncol(funspace)))
    
    while (min(rowSums(data.temp)) == 0) {
      data.temp <- funspace[sample(n.spp, n.sample, replace = FALSE, prob = NULL),]      
      traits.tax.temp <- data.temp[(intersect(rownames(data.temp), colnames(tax.mat2))),]   
      #community submatrix containing only sps in the traits' sample data.temp:
      tax.mat2.temp <- subset(tax.mat2, select = rownames(traits.tax.temp)) 
    }
    
    # in every iteration a new taxmat has to be created, only with the n.sample sampled sps, in order for FRic to run
    # fills FRic.out by row ([i,]) in every iteration
    FRic.out[i,] <- fric_3d(tax.mat2.temp,data.temp,m=2,prec="Qt")      
    FDis.out[i,] <- fdisp_k_sub(distrait,tax.mat2,tax_sub=colnames(tax.mat2),m=2)$FDis
  }
  
  FRic <- colMeans(FRic.out)
  FDis <- colMeans(FDis.out)
  
  # NULL MODELS FOR FD indexes
  
  FD_Rand <- function(x, m) {
    
    if(m == "FRic"){
      return(fric_3d(randomizeMatrix(tax.mat2, null.model = "independentswap"),x,m=2,prec="Qt"))
    }
    else if (m == "FDis"){
      return(fdisp_k_sub(x,randomizeMatrix(tax.mat2,null.model = "independentswap"),tax_sub=colnames(tax.mat2),m=2)$FDis
)    }
  }
  
  FD_Null <- function(m){
    
    obs.null <- NULL
    
    if(m == "FRic"){
      obs.null <- cbind(FRic, replicate(n.reps, FD_Rand(funspace, "FRic")))
    }
    else if (m == "FDis"){
      obs.null <- cbind(FDis, replicate(n.reps, FD_Rand(distrait, "FDis")))
    }
    
    rownames(obs.null) = rownames(tax.mat2)
    FD.rand.mean <- rowMeans(obs.null[,1:(n.reps+1)])
    ses <- (obs.null[,1] - apply(obs.null, 1, mean)) / apply(obs.null, 1, sd)
    
    #calculate pval for all communities
    pval <- apply(obs.null, MARGIN = 1, rank)[1,]/(n.reps+1) 
    if(m == "FRic"){
      return(cbind(FRic, FRic.rand.mean = FD.rand.mean, FRic.obs.z = ses, FRic.obs.p = pval, lat))
    }
    else if (m == "FDis"){
      return(cbind(FDis, FDis.rand.mean = FD.rand.mean, FDis.obs.z = ses, FDis.obs.p = pval, lat))
    }
    
  }
  
  FRic.null <- FD_Null("FRic")
  FDis.null <- FD_Null("FDis")
  
  
  save.image("precompiled/FD_indexes.RData")
} 

# CALCULATIONS AND PLOTS 

plot(FRic.null[,2]~lat, data=FRic.null, pch=1, col="violetred", ylab="FRic", ylim=c(0.75,1),las=1)
points(FRic~lat, data=FRic.null, col="violetred", pch=16)
text(FRic~lat, label=centr$index, pos=2, cex=0.6)
legend("topright", c("FRic obs","FRic rand" ), col=c("violetred","violetred" ), pch=c(16,1), bty="o", y.intersp=0.7)


# regression models
line.let= 0.5
adj.let=-0.4
cex.let=1.4
lat.seq<-seq(min(lat),max(lat),length.out = 1000)

mod_FRic<-lm(FRic~lat+I(lat^2))

# FRic 
plot(lat,exp(FRic), ylab="FRic ", xlab = "Latitude (degrees)", pch=16, 
     main="FRic~lat", xaxt = "n",yaxt="n", ylim=c(1.9,2.7), col = "purple")
axis(1,c(40,50,60))
axis(2, round(c(exp(0.7),exp(0.8),exp(0.9),exp(1)),2), las=2)
# Fitted values
lines(lat.seq,exp(predict(mod_FRic,data.frame(lat=lat.seq))), lwd=3, col="purple")
r2 <- bquote(italic(r)^2 == .(round(summary(mod_FRic)$r.squared,2)))
p.lat <- bquote(italic(p.lat) == .(round(summary(mod_FRic)$coefficients[2,4],3)))
p.lat2 <- bquote(italic(p.lat2) == .(round(summary(mod_FRic)$coefficients[3,4],3)))

mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext(p.lat, line = -2.3, adj = 0.9, cex = 1.2, font = 2)
mtext(p.lat2, line = -3.5, adj = 0.9, cex = 1.2, font = 2)
legend("topleft", c("FRic obs","FRic rand" ), col=c("purple","purple" ), pch=c(16,1), bty="o", y.intersp=0.7)




plot(FRic.obs.z~lat, data = FRic.null)
text(FRic.obs.z~lat,data = FRic.null, labels=centr$index, pos=1, cex=0.7)
abline(h=0, lty="dashed")



# Statistical test on unimodality 
# anova to examine the significance of the quadratic model
mod_FRic.lin <- lm(FRic~lat)
summary.uni<-summary(mod_FRic)
anova.uni <- anova(mod_FRic.lin,mod_FRic)
anova.uni
summary.uni


# FDis 
plot(FDis~lat)
plot(ses.FDis~lat)
abline(h=0, lty="dashed")
text(ses.FDis~lat, label= centr$index, pos=3, cex=0.6)


plot(FDis~lat, data=FDis.null, pch=16, col="violetred", ylab="FDis", xlab="Latitude")
points(FDis.null[,2]~lat, data=FDis.null, col="violetred", pch=1)
text(FDis~lat, labels=centr$index, cex=0.6, col="grey", pos=4)
legend("topright", c("FDis obs", "FDis rand"), col=c("violetred","violetred" ), pch=c(16,1), bty="o", y.intersp=0.7)

# regression significant
mod_FDis.obs<-lm(FDis~lat)
summary(mod_FDis.obs)

