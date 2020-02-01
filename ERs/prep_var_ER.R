library(FD)

source("external_files/quality_funct_space_fromdist.R")

# Taxonomic data
tax.mat1 <- read.csv("external_files/ecoregions.csv",header=TRUE, row.names = 1) # Loading site x taxa matrix
tax.mat1[is.na(tax.mat1)] <- 0

# Trait data
datatr <- read.csv2("external_files/biol_traits.csv", header = TRUE, dec = ".", row.names=3) 
traits <- subset(datatr, select = c(teg:ae_act), rownames=3)

# Splitting fuzzy,binary,quantitative
tabF<- traits[, c(5:15, 16:23)]                  #fuzzy
tabFp <- prep.fuzzy(tabF, c(4,7,4,4))
tabB<- traits[, c(1:2)]                          #binary
tabBp <-prep.binary(tabB, 2)
tabQ<- traits[, c(3:4)]                          #quantitative
ktab1 <- ktab.list.df(list(tabFp,tabBp,tabQ))
distrait <- dist.ktab(ktab1, c("F","B","Q"))        #distance matrix
taball<-cbind.data.frame(tabFp,tabBp,tabQ)

### Calculating species x species Gower dissimilarity matrix based on biological traits

# Does this data set have missing values?
if(is.euclid(distrait)==F) distrait <- lingoes(distrait)

# Estimating the optimum number of dimensions.
# mean Squared-Deviation of 0.010 means that average deviation between Euclidean distance and 
# Gower's distance is of (0.01)^0.5=0.1, so it can be seen like an average error of 10% (Villeger)
qual_fs_r <- quality_funct_space_fromdist(distrait, nbdim=10) 
qual_fs_r$meanSD # 2D

# PCoA
pcodistrait <- dudi.pco(distrait, scan = F, nf = 4)      #non negative eigenvalues

# Checking the species in traits matrix are the same of that in community matrix
traits.tax <- taball[(intersect(rownames(taball), colnames(tax.mat1))),]
tax.mat2 <- subset(tax.mat1, select = rownames(traits.tax))
#write.table(tax.mat2, "tax.mat2_197_2911.csv",sep=",")

# Selecting taxa ocurring in the empirical dataset (tax.mat1)
funspace <- pcodistrait$li[(intersect(rownames(pcodistrait$li), colnames(tax.mat1))),]

# Merging species' coordinates in the pcos with traits for further analyses
funspacetr <- cbind.data.frame(funspace[,1:2],traits)

#### sto arxeio me to niche conservatism to 3 na fortwnw to source prep var gia na pairnw to funspacetr
#write.table(funspacetr,"funspacetr_2911.csv", sep=",")

#write.table(funspace, "funspace_2911.csv", sep=",")

#loading phylogeny file
my.phylo <- read.tree("external_files/Ultrametric_Trichoptera.tre")

#All sps that exist in phylogeny
com498 <- tax.mat1[,(intersect(colnames(tax.mat1), my.phylo$tip.label))]


