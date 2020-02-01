library(FD)
library(adegraphics)

# Supplementary functions
source("external_files/FD_functions.R")
source("external_files/quality_funct_space_fromdist.R")
source("prep_var_ER.R")


# 2D scatter plot associating labels with points.
pdf("results/FS_biol_197.pdf", width=10, height=7)
FS <- s.class(pcodistrait$li, fac = datatr$Fam, plines.col = 1:20, col = T)
dev.off()