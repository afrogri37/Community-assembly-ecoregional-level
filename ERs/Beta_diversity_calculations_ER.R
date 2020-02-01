library(betapart)
library(reshape2)

source("prep_var_ER.R")

# Computing BETA diversity

# Computing 3 distance matrices accounting for the (i) turnover (replacement), (ii) nestedness-resultant component, 
# and (iii) total dissimilarity (i.e. the sum of both components).

# 1. Taxonomic diversity
TbTud <- beta.pair(tax.mat2, index.family = "sorensen")$beta.sim
TbNed <- beta.pair(tax.mat2, index.family = "sorensen")$beta.sne
TbSord <- beta.pair(tax.mat2, index.family = "sorensen")$beta.sor


# conversion to dataframes for graphs
TbSor<- as.data.frame(as.matrix(TbSord))
TbNe<- as.data.frame(as.matrix(TbNed))
TbTu<- as.data.frame(as.matrix(TbTud))


# conversion of dataframes to one column vectors
TbTu.v<- melt(TbTu)$value
TbNe.v<- melt(TbNe)$value
TbSor.v<- melt(TbSor)$value


# 2. Phylogenetic diversity
com498 <- as.matrix(com498)   #community matrix

Pbcore <-phylo.betapart.core(com498, my.phylo)
PbSord <- phylo.beta.pair(Pbcore, my.phylo, index.family="sorensen")$phylo.beta.sor
PbTud <- phylo.beta.pair(Pbcore, my.phylo, index.family="sorensen")$phylo.beta.sim
PbNed <- phylo.beta.pair(Pbcore, my.phylo, index.family="sorensen")$phylo.beta.sne


# conversion to dataframes for graphs
PbSor<-as.data.frame(as.matrix(PbSord))
PbNe<-as.data.frame(as.matrix(PbNed))
PbTu<-as.data.frame(as.matrix(PbTud))


# conversion of dataframes to one column vectors
PbTu.v<- melt(PbTu)$value
PbNe.v<- melt(PbNe)$value
PbSor.v<- melt(PbSor)$value




# 3. Functional diversity
# Trying to find precompiled simulations. If not found, run simulations.
generate_simulations <- tryCatch(load("precompiled/bFD_indexes.RData"), error=function(c){""})
if(generate_simulations == ""){
  

  funspace<-as.matrix(funspace[1:2])
  tax.mat2<-as.matrix(tax.mat2)
  
  # Simulations  
  
  n.reps <- 999
  n.spp <- nrow(funspace)
  n.sample <- min(rowSums(tax.mat2)) 
  
  # initializing matrices where the sums will be stored
  FbSor2 <- matrix(0,nrow=18, ncol=18)  
  FbTu2 <-  matrix(0,nrow=18, ncol=18)
  FbNe2 <-  matrix(0,nrow=18, ncol=18)
  

  for (i in 1:n.reps) {
    
    data.temp <- array(0, c(nrow(funspace),ncol(funspace)))
    
    while (min(rowSums(data.temp)) == 0) {
      data.temp <- funspace[sample(n.spp, n.sample, replace = FALSE, prob = NULL),]
      traits.tax.temp <- data.temp[(intersect(rownames(data.temp), colnames(tax.mat2))),]     
      tax.mat2.temp <- subset(tax.mat2, select = rownames(traits.tax.temp)) #community submatrix containing only sps in the traits' sample "data.temp"
    }
    bcore  <- functional.betapart.core(tax.mat2.temp,data.temp, fbc.step=TRUE, multi=FALSE, warning.time= FALSE)
    f.beta.pair.temp <- functional.beta.pair(bcore,index.family="sorensen")

    eval.parent(substitute(FbSor2 <- FbSor2 + as.matrix(f.beta.pair.temp$funct.beta.sor)))
    eval.parent(substitute(FbTu2 <- FbTu2 + as.matrix(f.beta.pair.temp$funct.beta.sim)))
    eval.parent(substitute(FbNe2 <- FbNe2 + as.matrix(f.beta.pair.temp$funct.beta.sne)))
    
  }
  
  
  # Simulations' output  
  # Average matrices
  FbSor <- FbSor2/n.reps
  FbTu <- FbTu2/n.reps
  FbNe <- FbNe2/n.reps
  
  save.image("bFD_indexes.RData")
}

# conversion to dataframes for graphs
FbTu <- as.data.frame(FbTu)
FbNe <- as.data.frame(FbNe)
FbSor <- as.data.frame(FbSor)

FbSord<-as.dist(FbSor)
FbNed<-as.dist(FbNe)
FbTud<-as.dist(FbTu)

# conversion of dataframes to one column vectors
FbTu.v<- melt(FbTu)$value
FbNe.v<- melt(FbNe)$value
FbSor.v<- melt(FbSor)$value



# 4. Geographic distance matrix
source("external_files/GeoDistanceInMetresMatrix.R")

# ER centroid data
centr <- read.csv2(file = "external_files/ecoreg_centroids.csv", header = TRUE, dec = ",")
lat <- centr$lat 
lon <- centr$lon

# Geographic distance matrix, based on Vincenty formula for distance between two Latitude/Longitude points
coord <- subset(centr, select= c(lat,lon), row.names=centr$index)
centrdist <- GeoDistanceInMetresMatrix(coord)     #distances in km
row.names(centrdist) <- centr$index
colnames(centrdist) <- centr$index
centrdist.v<- melt(centrdist)$value
centrdistd <- as.dist(centrdist)

