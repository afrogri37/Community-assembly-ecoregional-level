library(Hmisc)

source("prep_var_ER.R")

# Spearman rank and Pearson correlations between original trait categories and functional space axes
round(cor(taball,pcodistrait$li, method="spearman"),3)
bind_data<-cbind(taball,pcodistrait$li)

# Returns both correlations and their significances
pears <- rcorr(as.matrix(bind_data[1:23]),as.matrix(bind_data[24:27]), type="pearson")
pearson_rp <- round(cbind.data.frame(pears$r[1:23,24:27],pears$P[1:23,24:27]),3)
write.table(pearson_rp,"results/pco_cor_table_bio.txt",sep="\t") # Saving table to interpret PCoA axes



# Explained variance by each axis
round(pcodistrait$eig[1:8]/sum(pcodistrait$eig),2)
# Cumulative explained variance 
sum(pcodistrait$eig[1:5]/sum(pcodistrait$eig))

