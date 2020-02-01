library(Hmisc)

source("prep_var_ER.R")

# Spearman rank and Pearson correlations between original trait categories and functional space axes
round(cor(taball,pcodistrait$li, method="spearman"),2)
bind_data<-cbind(taball,pcodistrait$li)
#returns both correlations and their significances
pears <- rcorr(as.matrix(bind_data[1:15]),as.matrix(bind_data[16:20]), type="pearson")

round(cor(taball,pcodistrait$li, method="spearman"),2)

# Explained variance by each axis
round(pcodistrait$eig[1:8]/sum(pcodistrait$eig),2)
# Cumulative explained variance 
sum(pcodistrait$eig[1:5]/sum(pcodistrait$eig))

# Correlation between original traits and pco axes, saving only 2 decimals
cor.table1 <- round(cor(traits_wd, pcodistrait$li, method="spearman", use="pairwise.complete.obs"),2)
cor.table1
write.table(cor.table1,"results/pco_cor_table_bio.txt",sep="\t") # Saving table to interpret PCoA axes
