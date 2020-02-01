source("prep_var_ER.R")
source("Beta_diversity_calculations_ER.R")

# These calculations require matrix format

TbTu <- as.matrix(TbTu)
TbSor <- as.matrix(TbSor)
TbNe <- as.matrix(TbNe)
Cont.TTu <- (mean(TbTu/TbSor, na.rm=T))*100
Cont.TNe <- (mean(TbNe/TbSor, na.rm=T))*100
Tbconts <- rbind(Cont.TNe,Cont.TTu)

PbTu <- as.matrix(PbTu)
PbSor <- as.matrix(PbSor)
PbNe <- as.matrix(PbNe)
Cont.PTu <- (mean(PbTu/PbSor, na.rm=T))*100
Cont.PNe <- (mean(PbNe/PbSor, na.rm=T))*100
Pbconts <- rbind(Cont.PNe,Cont.PTu)

FbTu <- as.matrix(FbTu)
FbSor <- as.matrix(FbSor)
FbNe <- as.matrix(FbNe)
Cont.FTu <- (mean(FbTu/FbSor, na.rm=T))*100
Cont.FNe <- (mean(FbNe/FbSor, na.rm=T))*100
Fbconts <- rbind(Cont.FNe,Cont.FTu)


# barplot contributions
allconts <- data.frame(Taxonomic = Tbconts, Phylogenetic = Pbconts, Functional = Fbconts)
barplot(as.matrix(allconts), ylim=c(0,100),
        main="Nestedness/Turnover Contributions \nin beta Diversities Across Europe", ylab="Sorensen beta Dissimilarity",
        beside=TRUE, 
        col=c("#ffdd55ff","#87deaaff" ))

legend("topright",cex=1, c("Nestedness","Turnover"),  
       fill = c("#ffdd55ff","#87deaaff" ), bty="n")

