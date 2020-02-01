
source("prep_var_ER.R")
source("Beta_diversity_calculations_ER.R")


# TAXONOMIC

# The Iberian peninsula is not a refugium despite the high taxonomic richness 

cd1 <- centrdist[,1]
mod_TbNe1<-lm(TbNe$ER1~cd1+I(cd1^2))
mod_TbTu1<-lm(TbTu$ER1~cd1)
summary(mod_TbTu1) 
cd1.seq<-seq(min(cd1),max(cd1),length.out = 1000)


plot(cd1,TbTu$ER1, col="#87deaaff", pch=16, xlab = "Geodist (km)", 
     main="Sorensen starting from \nthe Iberian Peninsula", xaxt = "n", ylim=c(0,0.5))
axis(1,c(1000,2000,3000))
abline(mod_TbTu1, lwd=3, col="#87deaaff")
r2Tu <- bquote(italic(r)^2 == .(round(summary(mod_TbTu1)$r.squared,2)))
mtext(r2Tu, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p<0.01",line = -2.1, adj = 0.9, cex = 1.2, font = 2)


plot(cd1,exp(TbNe$ER1), pch=16, col="#ffdd55ff", ylim=c(exp(0), exp(0.5)), xaxt = "n")
axis(1,c(1000,2000,3000))
# Fitted values
lines(cd1.seq,exp(predict(mod_TbNe1,data.frame(cd1 = cd1.seq))), lwd=3, col="#ffdd55ff")
r2Ne <- bquote(italic(r)^2 == .(round(summary(mod_TbNe1)$r.squared,2)))
mtext(r2Ne, line = -10.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p<0.01",line = -11.1, adj = 0.9, cex = 1.2, font = 2)
text(cd1,exp(TbNe$ER1), labels = centr$index, pos = 3, cex=0.6, col="grey")



# The Pyrenees
plot(TbTu$ER2 ~ centrdist[,2],pch = 16, ylim=c(0,0.6),
     col = "turquoise",  ylab = "Sorensen taxonomic Dissimilarity", xlab = "Geographic distance from the Pyrenees (km)")   #!!!!??
points(TbNe$ER2~centrdist[,2], col = "violetred", pch=16)
legend( "topleft", c("taxonomic beta Turnover", "taxonomic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(TbNe$ER2 ~ centrdist[,2], labels = centr$index, pos = 3, cex=0.6)
text(TbTu$ER2 ~ centrdist[,2], labels = centr$index, pos = 3, cex=0.6)

# Italy mostly turnover IN
plot(TbTu$ER3 ~ centrdist[,3],pch = 16, ylim=c(0,0.5), main = "Sorensen taxonomic dissimilarity\nfrom Italy",
     col = "#87deaaff",  xlab = "Geographic distance from Italy (km)")   #!!!!??
points(TbNe$ER3~centrdist[,3], col = "#ffdd55ff", pch=16)
text(TbNe$ER3 ~ centrdist[,3], labels = centr$index, pos = 3, cex=0.6, col="grey")

# The Alps are a refugium
cd4 <- centrdist[,4]
mod_TbNe4<-lm(TbNe$ER4~cd4+I(cd4^2))
mod_TbTu4<-lm(TbTu$ER4~cd4+I(cd4^2))
summary(mod_TbTu4) #both terms significant, but r2 low

cd4.seq<-seq(min(cd4),max(cd4),length.out = 1000)

plot(cd4,exp(TbNe$ER4), ylab="TbNe", xlab = "Geodist (km)", pch=16, col="#ffdd55ff",
     main="Sorensen dissimilarity starting from the Alps", xaxt = "n", ylim=c(1,exp(0.5)))
points(cd4,exp(TbTu$ER4), col="#87deaaff", pch=16)
axis(1,c(500,1000,1500,2000,3000))
# Fitted values
lines(cd4.seq,exp(predict(mod_TbNe4,data.frame(cd4 = cd4.seq))), lwd=3, col="#ffdd55ff")
lines(cd4.seq,exp(predict(mod_TbTu4,data.frame(cd4 = cd4.seq))), lwd=3, col="#87deaaff")
r2Ne <- bquote(italic(r)^2 == .(round(summary(mod_TbNe4)$r.squared,2)))
r2Tu <- bquote(italic(r)^2 == .(round(summary(mod_TbTu4)$r.squared,2)))
mtext("p<0.01",line = -2.1, adj = 0.9, cex = 1.2, font = 2)
mtext(r2Ne, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext(r2Tu, line = -10.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p<0.01",line = -11.1, adj = 0.9, cex = 1.2, font = 2)
text(cd4,exp(TbNe$ER4), labels = centr$index, pos = 3, cex=0.6, col="grey")



# In Greece strange pattern, similar to the functional beta, but mostly driven by turnover. nestedness only with 4, 8, 9, 10
cd6 <- centrdist[,6]
mod_TbTu6<-lm(TbTu$ER6~cd6) # significant
summary(mod_TbTu6) 

plot(cd6,TbNe$ER6,  xlab = "Geodist (km)", pch=16, col="#ffdd55ff", ylab="TbNe/Tu",
     main="Sorensen dissimilarity starting from\nthe Hellenic-Western Balkan", xaxt = "n", ylim=c(0,0.5))
points(cd6,TbTu$ER6, col="#87deaaff", pch=16)
axis(1,c(0,500,1000,1500,2000,3000))
# Fitted values
abline(mod_TbTu6, lwd=3, col="#87deaaff")
r2Tu <- bquote(italic(r)^2 == .(round(summary(mod_TbTu6)$r.squared,2)))
mtext("p<0.01",line = -2.1, adj = 0.9, cex = 1.2, font = 2)
mtext(r2Tu, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
text(TbNe$ER6 ~ centrdist[,6], labels = centr$index, pos = 3, cex=0.6, col="grey")


#In  Western HIghlands regressions
cd8 <- centrdist[,8]
mod_TbNe8<-lm(TbNe$ER8~cd8)
summary(mod_TbNe8)

plot(cd8,TbNe$ER8, ylab="TbNe/Tu", xlab = "Geodist (km)", pch=16, col="#ffdd55ff",
     main="Sorensen dissimilarity starting from\n the Western Highlands", xaxt = "n", ylim=c(0,0.5))
points(cd8,TbTu$ER8, col="#87deaaff", pch=16)
axis(1,c(0,500,1000,1500,2000,3000))
# Fitted values
abline(mod_TbNe8, lwd=3, col="#ffdd55ff")
r2Ne <- bquote(italic(r)^2 == .(round(summary(mod_TbNe8)$r.squared,2)))
mtext("p<0.01",line = -2.1, adj = 0.8, cex = 1.2, font = 2)
mtext(r2Ne, line = -1.1, adj = 0.8, cex = 1.2, font = 2)
text(TbNe$ER8 ~ centrdist[,8], labels = centr$index, pos = 3, cex=0.6, col="grey")


#IN regressions Central highlands are a refugium
cd9 <- centrdist[,9]
mod_TbNe9<-lm(TbNe$ER9~cd9)
cd9.seq<-seq(min(cd9),max(cd9),length.out = 1000)


plot(cd9,TbNe$ER9, ylab="TbNe/Tu", xlab = "Geodist (km)", pch=16, col="#ffdd55ff",
     main="Sorensen dissimilarity starting from\n the Central Highlands", xaxt = "n", ylim=c(0,0.5))
points(cd9,TbTu$ER9, col="#87deaaff", pch=16)
axis(1,c(0,500,1000,1500,2000,3000))
# Fitted values
abline(mod_TbNe9, lwd=3, col="#ffdd55ff")
r2Ne <- bquote(italic(r)^2 == .(round(summary(mod_TbNe9)$r.squared,2)))
mtext("p<0.01",line = -2.1, adj = 0.9, cex = 1.2, font = 2)
mtext(r2Ne, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
text(TbNe$ER9 ~ centrdist[,9], labels = centr$index, pos = 3, cex=0.6, col="grey")


#ER11
cd11 <- centrdist[,11]
mod_TbNe11<-lm(TbNe$ER11~cd11+I(cd11^2))
summary(mod_TbNe11) 
cd11.seq<-seq(min(cd11),max(cd11),length.out = 1000)


plot(cd11,TbNe$ER11, ylab="TbNe/Tu", xlab = "Geodist (km)", pch=16, col="#ffdd55ff",
     main="Sorensen dissimilarity starting from\n the Central Highlands", xaxt = "n", ylim=c(0,0.5))
points(cd11,TbTu$ER11, col="#87deaaff", pch=16)
axis(1,c(0,500,1000,1500,2000,3000))
# Fitted values
abline(mod_TbNe11, lwd=3, col="#ffdd55ff")
r2Ne <- bquote(italic(r)^2 == .(round(summary(mod_TbNe11)$r.squared,2)))
mtext("p<0.01",line = -2.1, adj = 0.11, cex = 1.2, font = 2)
mtext(r2Ne, line = -1.1, adj = 0.11, cex = 1.2, font = 2)
text(TbNe$ER11 ~ centrdist[,11], labels = centr$index, pos = 3, cex=0.6)



# IN ER20 Regressions
cd20 <- centrdist[,17]
mod_TbNe20<-lm(TbNe$ER20~cd20+I(cd20^2))
mod_TbTu20<-lm(TbTu$ER20~cd20+I(cd20^2))
summary(mod_TbNe20) 

cd20.seq<-seq(min(cd20),max(cd20),length.out = 1000)
plot(cd20,exp(TbNe$ER20), xlab = "Geodist (km)", pch=16, col="#ffdd55ff",
     main="Sorensen dissimilarity starting from\nthe Borealic uplands", xaxt = "n", ylim=c(1,exp(0.5)))
points(cd20,exp(TbTu$ER20), col="#87deaaff", pch=16)
axis(1,c(0,500,1000,1500,2000,2500,3000))
# Fitted values
lines(cd20.seq,exp(predict(mod_TbNe20,data.frame(cd20 = cd20.seq))), lwd=3, col="#ffdd55ff")
lines(cd20.seq,exp(predict(mod_TbTu20,data.frame(cd20 = cd20.seq))), lwd=3, col="#87deaaff")
text(cd20,exp(TbNe$ER20), labels = centr$index, pos = 3, cex=0.6, col="grey")
r2Ne <- bquote(italic(r)^2 == .(round(summary(mod_TbNe20)$r.squared,2)))
r2Tu <- bquote(italic(r)^2 == .(round(summary(mod_TbTu20)$r.squared,2)))
mtext("p<0.01",line = -2.1, adj = 0.9, cex = 1.2, font = 2)
mtext(r2Ne, line = -10.1, adj = 0.9, cex = 1.2, font = 2)
mtext(r2Tu, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p<0.01",line = -11.1, adj = 0.9, cex = 1.2, font = 2)



# Regressions ER22 IN
cd22 <- centrdist[,18]
mod_TbNe22<-lm(TbNe$ER22~cd22+I(cd22^2))
mod_TbTu22<-lm(TbTu$ER22~cd22+I(cd22^2))
summary(mod_TbTu22) 

cd22.seq<-seq(min(cd22),max(cd22),length.out = 1000)

plot(cd22,exp(TbNe$ER22),  xlab = "Geodist (km)", pch=16, col="#ffdd55ff",
     main="Sorensen dissimilarity starting from\nthe Fenno-scandian shield", xaxt = "n", ylim=c(1,exp(0.5)))
points(cd22,exp(TbTu$ER22), col="#87deaaff", pch=16)
axis(1,c(0,500,1000,1500,2200,2500,3000))
# Fitted values
lines(cd22.seq,exp(predict(mod_TbNe22,data.frame(cd22 = cd22.seq))), lwd=3, col="#ffdd55ff")
lines(cd22.seq,exp(predict(mod_TbTu22,data.frame(cd22 = cd22.seq))), lwd=3, col="#87deaaff")
text(cd22,exp(TbNe$ER22), labels = centr$index, pos = 3, cex=0.6, col="grey")
r2Ne <- bquote(italic(r)^2 == .(round(summary(mod_TbNe22)$r.squared,2)))
r2Tu <- bquote(italic(r)^2 == .(round(summary(mod_TbTu22)$r.squared,2)))
mtext("p<0.01",line = -2.1, adj = 0.9, cex = 1.2, font = 2)
mtext(r2Ne, line = -10.1, adj = 0.9, cex = 1.2, font = 2)
mtext(r2Tu, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p<0.01",line = -11.1, adj = 0.9, cex = 1.2, font = 2)





# ER13
plot(TbTu$ER13 ~ centrdist[,12],pch = 16, ylim=c(0,0.6), 
     col = "turquoise", ylab = "Sorensen taxonomic Dissimilarity", xlab = "Geographic distance from the Central highlands (km)")   #!!!!??
points(TbNe$ER13~centrdist[,12], col = "violetred", pch=16)
legend( "topleft", c("taxonomic beta Turnover", "taxonomic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(TbNe$ER13 ~ centrdist[,12], labels = centr$index, pos = 3, cex=0.6)
text(TbTu$ER13 ~ centrdist[,12], labels = centr$index, pos = 3, cex=0.6)


#carpathians don't seem to be a refugium whaaaaat?!
#high nestedness only with Pont and the Tundra
plot(TbTu$ER10 ~ centrdist[,10],pch = 16, ylim=c(0,0.6), 
     col = "turquoise", ylab = "Sorensen taxonomic Dissimilarity", xlab = "Geographic distance from the Carpathians (km)")   #!!!!??
points(TbNe$ER10~centrdist[,10], col = "violetred", pch=16)
legend( "topleft", c("taxonomic beta Turnover", "taxonomic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(TbNe$ER10 ~ centrdist[,10], labels = centr$index, pos = 3, cex=0.6)
text(TbTu$ER10 ~ centrdist[,10], labels = centr$index, pos = 3, cex=0.6)



#ER18

cd18 <- centrdist[,16]
mod_TbTu18<-lm(TbTu$ER18~cd18)
summary(mod_TbTu18) 



plot(TbTu$ER18~cd18, col="#87deaaff", pch=16, ylim=c(0,0.5),xlab = "Geographic distance (km)",  main="Sorensen dissimilarity  from\nEngland", xaxt = "n")
abline(mod_TbTu18, col="#87deaaff", lwd=3)
axis(1,c(0,500,1000,1500,2000,2500,3000))
points(TbNe$ER18~cd18, col="#ffdd55ff", pch=16)
text(cd18,TbNe$ER18, labels = centr$index, pos = 3, cex=0.6, col="grey")
r2Tu <- bquote(italic(r)^2 == .(round(summary(mod_TbTu18)$r.squared,2)))
mtext("p<0.001",line = -2.1, adj = 0.9, cex = 1.2, font = 2)
mtext(r2Tu, line = -1.1, adj = 0.9, cex = 1.2, font = 2)





# in order to know which area is nested to which
#TRic.ord <- order(TRic, decreasing = T)

par(mfrow=c(1,2))
matplot(centrdist, TbNe, type = "p", lty = 1:5, lwd = 1, pch = 16,
        col = "violetred", cex = 0.5, 
        xlab = "Geographic distance", ylab = "Tbeta Tu & Ne", xlim = NULL, ylim = c(0,0.5))
matpoints(centrdist, TbTu, type = "p", lty = 1:5, lwd = 1, pch = 16,
          col = "turquoise", cex = 0.5 )

matplot(centrdist,TbSor, type = "p", lty = 1:5, lwd = 1, pch = 16,
        col = 8, cex = 0.5, 
        xlab = "Geographic distance", ylab = "Tbeta Sorensen dissimilarity", xlim = NULL, ylim = NULL)

mnt.Ts <- mantel.rtest(TbSord, centrdistd, nrepet = 9999)  # significant
mnt.Ttu <- mantel.rtest(TbTud, centrdistd, nrepet = 9999)  # significant
mnt.Tne <- mantel.rtest(TbNed,centrdistd, nrepet = 9999)   # not significant

# Regression models
line.let= 0.5
adj.let=-0.4
cex.let=1.4
centrdist.seq<-seq(min(centrdist.v),max(centrdist.v),length.out = 1000)
mod_TbTu <- lm(TbTu.v~centrdist.v+I(centrdist.v^2))


par(mfrow=c(1,1))
plot(centrdist.v,exp(TbTu.v), ylab="TbTu", ylim=c(1,1.648721),xlab = "Geographic distance (km)", pch=16, 
     main="Taxonomic beta turnover", xaxt = "n")
#can't add nestedness points
axis(1,c(1000,2000,3000))

# Fitted values
lines(centrdist.seq,exp(predict(mod_TbTu,data.frame(centrdist.v=centrdist.seq))), lwd=3, col="violetred")
r2 <- bquote(italic(r)^2 == .(round(summary(mod_TbTu)$r.squared,2)))
p <- bquote(italic(p) == .(round(summary(mod_TbTu)$coefficients[2,4],3)))
p2 <- bquote(italic(p2) == .(round(summary(mod_TbTu)$coefficients[3,4],3)))
mtext(p, line = -2.3, adj = 0.9, cex = 1.2, font = 2)
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext(p2, line = -3.4, adj = 0.9, cex = 1.2, font = 2)
plot(centrdist.v,exp(TbNe.v), ylab="TbNe",col=2, xlab = "Geographic distance (km)", pch=16, 
     main="Taxonomic beta Nestedness", xaxt = "n", ylim=c(1,1.648721))
axis(1,c(1000,2000,3000))


mod_TbSor <- lm(TbSor.v~centrdist.v)

plot(TbSor.v~centrdist.v, ylab="TbSor", ylim=c(0,0.5),xlab = "Geographic distance (km)", pch=16, 
     main="Taxonomic beta Sorensen", xaxt = "n",col=3)
axis(1,c(1000,2000,3000))
abline(mod_TbSor, lwd=3, col="violetred")
r2 <- bquote(italic(r)^2 == .(round(summary(mod_TbSor)$r.squared,2)))
p <- bquote(italic(p) == .(round(summary(mod_TbSor)$coefficients[2,4],3)))
mtext(p, line = -2.3, adj = 0.9, cex = 1.2, font = 2)
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)






# PHYLOGENETIC

par(mfrow=c(2,2))

#The Iberian peninsula is not a refugium despite the high taxonomic richness
plot(PbTu$ER1 ~ centrdist[,1],pch = 16, ylim=c(0,0.45),
     col = "turquoise",  ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from the Iberian peninsula (km)")   
points(PbNe$ER1~centrdist[,1], col = "violetred", pch=16)
text(PbNe$ER1 ~ centrdist[,1], labels = centr$index, pos = 3, cex=0.6)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)


#Italy mostly turnover
plot(PbTu$ER3 ~ centrdist[,3],pch = 16, ylim=c(0,0.45), 
     col = "turquoise", ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from Italy (km)")  
points(PbNe$ER3~centrdist[,3], col = "violetred", pch=16)
text(PbNe$ER3 ~ centrdist[,3], labels = centr$index, pos = 3, cex=0.6)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16,  y.intersp=0.8)

# The Alps
plot(PbTu$ER4 ~ centrdist[,4],pch = 16, ylim=c(0,0.45), 
     col = "turquoise", ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from the Alps (km)")  
points(PbNe$ER4~centrdist[,4], col = "violetred", pch=16)
text(PbNe$ER4 ~ centrdist[,4], labels = centr$index, pos = 3, cex=0.6)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16,  y.intersp=0.8)


#Greece strange pattern, similar to the functional beta, but mostly driven by turnover. nestedness only with 4, 8, 9, 12
plot(PbTu$ER6 ~ centrdist[,6],pch = 16, ylim=c(0,0.45), 
     col = "turquoise", ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from the Hellenic western Balkan (km)")  
points(PbNe$ER6~centrdist[,6], col = "violetred", pch=16)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(PbNe$ER6 ~ centrdist[,6], labels = centr$index, pos = 3, cex=0.6)


#Central highlands
plot(PbTu$ER7 ~ centrdist[,7],pch = 16, ylim=c(0,0.45), 
     col = "turquoise", ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from the Central highlands (km)")  
points(PbNe$ER7~centrdist[,7], col = "violetred", pch=16)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(PbNe$ER7 ~ centrdist[,7], labels = centr$index, pos = 3, cex=0.6)


#Western highlands 
plot(PbTu$ER8 ~ centrdist[,8],pch = 16, ylim=c(0,0.45), 
     col = "turquoise",  ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from the Western highlands (km)")  
points(PbNe$ER8~centrdist[,8], col = "violetred", pch=16)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(PbNe$ER8 ~ centrdist[,8], labels = centr$index, pos = 3, cex=0.6)

#Western highlands 
plot(PbTu$ER9 ~ centrdist[,9],pch = 16, ylim=c(0,0.45), 
     col = "turquoise",  ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from the Western highlands (km)")  
points(PbNe$ER9~centrdist[,9], col = "violetred", pch=16)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(PbNe$ER9 ~ centrdist[,9], labels = centr$index, pos = 3, cex=0.6)

#carpathians 
#high nestedness only with Pont and the Tundra
plot(PbTu$ER10 ~ centrdist[,10],pch = 16, ylim=c(0,0.45), 
     col = "turquoise", ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from the Carpathians (km)")   
points(PbNe$ER10~centrdist[,10], col = "violetred", pch=16)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(PbNe$ER10 ~ centrdist[,10], labels = centr$index, pos = 3, cex=0.6)

#ER 11
plot(PbTu$ER11 ~ centrdist[,11],pch = 16, ylim=c(0,0.45), 
     col = "turquoise", ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from ER11 (km)")   
points(PbNe$ER11~centrdist[,11], col = "violetred", pch=16)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(PbNe$ER11 ~ centrdist[,11], labels = centr$index, pos = 3, cex=0.6)


#England
plot(PbTu$ER18 ~ centrdist[,16],pch = 16, ylim=c(0,0.45), 
     col = "turquoise", ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from England (km)")  
points(PbNe$ER18~centrdist[,16], col = "violetred", pch=16)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(PbNe$ER18 ~ centrdist[,16], labels = centr$index, pos = 3, cex=0.6)


#ER 20
plot(PbTu$ER20 ~ centrdist[,17],pch = 16, ylim=c(0,0.45), 
     col = "turquoise", ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from ER20 (km)")   
points(PbNe$ER20~centrdist[,17], col = "violetred", pch=16)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(PbNe$ER20 ~ centrdist[,17], labels = centr$index, pos = 3, cex=0.6)

# ER 22
plot(PbTu$ER22 ~ centrdist[,18],pch = 16, ylim=c(0,0.45), 
     col = "turquoise", ylab = "Sorensen Phylogenetic Dissimilarity", xlab = "Geographic distance from ER22 (km)")   
points(PbNe$ER22~centrdist[,18], col = "violetred", pch=16)
#legend( "topleft", c("Phylogenetic beta Turnover", "Phylogenetic beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(PbNe$ER22 ~ centrdist[,18], labels = centr$index, pos = 3, cex=0.6)



matplot(centrdist, PbTu, type = "p", lty = 1:5, lwd = 1, pch = 16,
        col = "turquoise", cex = 0.5, 
        xlab = "Geographic distance", ylab = "Pbeta Tu & Ne", xlim = NULL, ylim = c(0,0.5))

matpoints(centrdist,PbNe, type = "p", lty = 1:5, lwd = 1, pch = 16,
        col = "violetred", cex = 0.5, 
        xlab = "Geographic distance", ylab = "Pbeta nestedness", xlim = NULL, ylim = c(0,0.5))


# regression models

par(mfrow=c(1,1))
mod_PbSor <- lm(PbSor.v~centrdist.v)
plot(PbSor.v~centrdist.v, ylab="PbSor", xlab = "Geographic distance (km)", pch=16, 
     main="Phylogenetic beta Sorensen", xaxt = "n", col=3, ylim=c(0,0.5))
axis(1,c(1000,2000,3000))
abline(mod_PbSor, lwd=3, col="violetred")
r2 <- bquote(italic(r)^2 == .(round(summary(mod_PbSor)$r.squared,2)))
p <- bquote(italic(p) == .(round(summary(mod_TbSor)$coefficients[2,4],3)))
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p<0.001", line = -2.2, adj = 0.9, cex = 1.2, font = 2)


mod_PbTu <- lm(PbTu.v~centrdist.v)
plot(PbTu.v~centrdist.v, ylab="PbTu", xlab = "Geographic distance (km)", pch=16, 
     main="Phylogenetic beta turnover", xaxt = "n", ylim=c(0,0.5))
axis(1,c(1000,2000,3000))
# Fitted values
abline(mod_PbTu, lwd=3, col="violetred")
r2 <- bquote(italic(r)^2 == .(round(summary(mod_PbTu)$r.squared,2)))
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p<0.001", line = -2.2, adj = 0.9, cex = 1.2, font = 2)
points(PbNe.v~centrdist.v, col=2, pch=16)


matplot(centrdist,PbSor, type = "p", lty = 1:5, lwd = 1, pch = 16,
        col = 8, cex = 0.5, 
        xlab = "Geographic distance", ylab = "Pbeta Sorensen dissimilarity", xlim = NULL, ylim = c(0,0.5))

mnt.Ps <- mantel.rtest(PbSord, centrdistd, nrepet = 9999)  # signif
mnt.Ptu <- mantel.rtest(PbTud, centrdistd, nrepet = 9999)  # signif
mnt.Pne <- mantel.rtest(PbNed,centrdistd, nrepet = 9999)   



# Functional

#order FRic to check for nestedness
order(FRic, decreasing=T)


#The Iberian peninsula is not a refugium despite the high taxonomic richness
plot(FbTu$ER1 ~ centrdist[,1],pch = 16, 
     col = "turquoise", ylim=c(0,0.12), ylab = "Sorensen Functional Dissimilarity", xlab = "Geographic distance from the Iberian peninsula (km)")   #!!!!??
points(FbNe$ER1~centrdist[,1], col = "violetred", pch=16)
legend( "topright", c("Functional beta Turnover", "Functional beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8, cex=0.6)
text(FbNe$ER1 ~ centrdist[,1], labels = centr$index, pos = 3, cex=0.6)


#Italy could be, but not such a strong relationship 
plot(FbTu$ER3 ~ centrdist[,3],pch = 16, 
     col = "turquoise", ylim=c(0,0.12), ylab = "Sorensen Functional Dissimilarity", xlab = "Geographic distance from Italy (km)")   #!!!!??
points(FbNe$ER3~centrdist[,3], col = "violetred", pch=16)
legend( "topright", c("Functional beta Turnover", "Functional beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)

# The Alps are a refugium
plot(FbTu$ER4 ~ centrdist[,4],pch = 16, 
     col = "turquoise", ylim=c(0,0.12), ylab = "Sorensen Functional Dissimilarity", xlab = "Geographic distance from the Alps (km)")   #!!!!??
points(FbNe$ER4~centrdist[,4], col = "violetred", pch=16)
legend( "topright", c("Functional beta Turnover", "Functional beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8, cex=0.7)
text(FbNe$ER4 ~ centrdist[,4], labels = centr$index, pos = 3, cex=0.6)
#text(FbTu$ER4 ~ centrdist[,4], labels = centr$index, pos = 3, cex=0.6)


#Greece strange pattern
plot(FbTu$ER6 ~ centrdist[,6],pch = 16, 
     col = "turquoise", ylim=c(0,0.12), ylab = "Sorensen Functional Dissimilarity", xlab = "Geographic distance from the Hellenic western Balkan (km)")   #!!!!??
points(FbNe$ER6~centrdist[,6], col = "violetred", pch=16)
legend( "topright", c("Functional beta Turnover", "Functional beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(FbNe$ER6 ~ centrdist[,6], labels = centr$index, pos = 3, cex=0.6)
text(FbTu$ER6 ~ centrdist[,6], labels = centr$index, pos = 3, cex=0.6)

#Western highlands are a refugium
plot(FbTu$ER8 ~ centrdist[,8],pch = 16, 
     col = "turquoise", ylim=c(0,0.12), ylab = "Sorensen Functional Dissimilarity", xlab = "Geographic distance from the Western highlands (km)")   #!!!!??
points(FbNe$ER8~centrdist[,8], col = "violetred", pch=16)
legend( "topright", c("Functional beta Turnover", "Functional beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(FbNe$ER8 ~ centrdist[,8], labels = centr$index, pos = 3, cex=0.6)

# Central highlands, too
plot(FbTu$ER9 ~ centrdist[,9],pch = 16, 
     col = "turquoise", ylim=c(0,0.12), ylab = "Sorensen Functional Dissimilarity", xlab = "Geographic distance from the Central highlands (km)")   #!!!!??
points(FbNe$ER9~centrdist[,9], col = "violetred", pch=16)
legend( "topright", c("Functional beta Turnover", "Functional beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)
text(FbNe$ER9 ~ centrdist[,9], labels = centr$index, pos = 3, cex=0.6)
text(FbTu$ER9 ~ centrdist[,9], labels = centr$index, pos = 3, cex=0.6)


#carpathians?
plot(FbTu$ER10 ~ centrdist[,10],pch = 16, 
     col = "turquoise", ylim=c(0,0.12), ylab = "Sorensen Functional Dissimilarity", xlab = "Geographic distance from the Carpathians (km)")   #!!!!??
points(FbNe$ER10~centrdist[,10], col = "violetred", pch=16)
legend( "topright", c("Functional beta Turnover", "Functional beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)

plot(FbTu$ER14 ~ centrdist[,13],pch = 16, 
     col = "turquoise", ylim=c(0,0.12), ylab = "Sorensen Functional Dissimilarity", xlab = "Geographic distance from the Central plains (km)")   #!!!!??
points(FbNe$ER14~centrdist[,13], col = "violetred", pch=16)
text(FbNe$ER14 ~ centrdist[,13], labels = centr$index, pos = 3, cex=0.6)
legend( "topright", c("Functional beta Turnover", "Functional beta Nestedness"), col= c("turquoise","violetred"), pch = 16, y.intersp=0.8)


# checking if ecoregions with more than 250 species are refugia based on nestedness

par(mfrow=c(1,1))
matplot(centrdist, FbTu, type = "p", lty = 1:5, lwd = 1, pch = 16,
        col = "turquoise", cex = 0.5, 
        xlab = "Geographic distance", ylab = "Fbeta Tu / Ne", xlim = NULL, ylim = c(0, 0.08))
matpoints(centrdist, FbNe, lty = 1:5, lwd = 1, pch = 16, col="violetred", cex=0.5)
matplot(centrdist,FbNe, type = "p", lty = 1:5, lwd = 1, pch = 16,
        col = 6, cex = 0.5, 
        xlab = "Geographic distance", ylab = "Fbeta nestedness", xlim = NULL, ylim = NULL)

matplot(centrdist,FbSor, type = "p", lty = 1:5, lwd = 1, pch = 16,
        col = 8, cex = 0.5, 
        xlab = "Geographic distance", ylab = "Fbeta Sorensen dissimilarity", xlim = NULL, ylim = NULL)

#MANTEL
mnt.Fs <- mantel.rtest(FbSord, centrdistd, nrepet = 9999)  #observation 0.611
mnt.Ftu <- mantel.rtest(FbTud, centrdistd, nrepet = 9999)  #observation 0.50
mnt.Fne <- mantel.rtest(FbNed,centrdistd, nrepet = 9999)    #obs 0.33




# regression models

par(mfrow=c(1,1))
mod_FbSor <- lm(FbSor.v~centrdist.v)   #linear model yields more normal residuals
mod_FbSor2 <- lm(FbTu.v~centrdist.v+I(centrdist.v^2))

plot(FbSor.v~centrdist.v, ylab="FbSor", xlab = "Geographic distance (km)", pch=16, 
     main="Functional beta Sorensen", xaxt = "n", col=3, ylim=c(0,0.08))
axis(1,c(1000,2000,3000))
abline(mod_FbSor, lwd=3, col="violetred")
r2 <- bquote(italic(r)^2 == .(round(summary(mod_FbSor)$r.squared,2)))
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p<0.001", line = -2.3, adj = 0.9, cex = 1.2, font = 2)


mod_FbTu <- lm(FbTu.v~centrdist.v+I(centrdist.v^2))
plot(centrdist.v,exp(FbTu.v), ylab="FbTu", xlab = "Geographic distance (km)", pch=16, 
     main="Functional beta turnover", xaxt = "n", ylim=c(1,1.083287))

axis(1,c(1000,2000,3000))

# Fitted values
lines(centrdist.seq,exp(predict(mod_FbTu,data.frame(centrdist.v=centrdist.seq))), lwd=3, col="violetred")
r2 <- bquote(italic(r)^2 == .(round(summary(mod_FbTu)$r.squared,2)))
p <- bquote(italic(p) == .(round(summary(mod_FbTu)$coefficients[2,4],3)))
mtext(p, line = -2.3, adj = 0.9, cex = 1.2, font = 2)
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p2<0.001", line = -3.4, adj = 0.9, cex = 1.2, font = 2)

mod_FbNe <- lm(FbNe.v~centrdist.v)
plot(FbNe.v~centrdist.v,pch=16,col="violetred", ylim=c(0,0.08),xaxt="n", lwd=3, ylab="FbNe")
axis(1,c(1000,2000,3000))
abline(mod_FbNe)
r2 <- bquote(italic(r)^2 == .(round(summary(mod_FbNe)$r.squared,2)))
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("p<0.001", line = -2.4, adj = 0.9, cex = 1.2, font = 2)



