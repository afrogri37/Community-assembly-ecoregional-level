# Appendix S5. R code and functions.
#
# Functions to estimate community functional features and to perform null models to assess for 
# non-random patterns along gradients of environmental filters.
#
# Code written by Cayetano Gutierrez-Canovas. Cardiff University (UK). 
# Email: GutierrezCanovasC@cardiff.ac.uk

cent.com<-function(abun,fpc,k){
  fpc[,1:k]->fpc.sel
  
  centroid<-matrix(NA,nrow(abun),k)
  colnames(centroid)<-paste("Axis",1:k)
  
  for (i in 1:nrow(abun)) {
    
    pres <- which(abun[i, ] > 0)
    w <- abun[i, pres]
    data.frame(spp=colnames(abun)[pres])->taxa_names
    fpc_sel<-data.frame(taxa=rownames(fpc),fpc.sel)
    fpc_sel <- sqldf("select * from taxa_names,fpc_sel where taxa_names.spp = fpc_sel.taxa")
    fpc_sel[,1]->row.names(fpc_sel)
    fpc_sel[,c(3:(k+2))]->fpc_sel
    centroid[i,] <- apply(fpc_sel, 2, weighted.mean, w = w)
    
  }
  return(centroid)
}

# fric_3d() estimates the Functional Richness of a set of communties
# This function computes the hypervolume to estimate how each community fills
# the functional space
#
# Inputs:
# taxa: community data
# fpc: functional space
# m: number of axes to select
# prec: convex hull precision ("Qt" or "QJ")
#
# Output:
# a vector with the Functional Richness of each community

fric_3d<-function(taxa,fpc,m,prec=c("Qt","QJ")){
  fric.3d<-rep(NA,nrow(taxa))
  convhulln(fpc[,1:m], c("FA",prec))$vol->fric.3d.max
  specnumber(taxa)->ric
  for (com in 1:nrow(taxa)){
    fpc[which(unlist(rep(taxa[com,]))>0),1:m]->tr.com
    if (ric[com]>=m+1) convhulln(tr.com, c("FA",prec))$vol/fric.3d.max->fric.3d[com] else NA->fric.3d[com]
  }
  return(fric.3d)
}

# fric_1d() estimates the Functional Richness of a set of communties
# This function computes the hypervolume to estimate how each community fills
# the functional space
#
# Inputs:
# taxa: community data
# fpc: functional space
# m: number of axes to select
#
# Output:
# a vector with the Functional Richness of each community

fric_1d<-function(taxa,fpc,m){
  
  fric.1d<-rep(NA,nrow(taxa))
  fric.1d.max<-rep(NA,m)
  
  for (i in 1:m) sum(abs(range(fpc[,i])))->fric.1d.max[i]
  
  specnumber(taxa)->ric
  for (com in 1:nrow(taxa)){
    fpc[which(unlist(rep(taxa[com,]))>0),1:m]->tr.com
    
    if (ric[com]>=1) mean(sum(abs(range(tr.com)))/fric.1d.max)->fric.1d[com] else NA->fric.1d[com]
  }
  return(fric.1d)
}

# calc.FR() estimates the Functional Redundancy (FR) of a set of communties
# This function estimates the taxonomic richness for each taxonomic group,
# estating FR as the ratio between species richness and the number of functional groups
# in a given community.
#
# Inputs:
# taxa: community data
# groups: a grouping vector with the Funtional Groups for each taxon
#
# Outputs:
#
# $nbsp: taxonomic richness
# $ric.fgrs: the taxonomic richness for each Functional Group
# $sd.fgr: the among-Functional-group SD for taxonomic richness for each community
# $FGR: the number of Functional Gruoups for each community
# $FR: the Functional Richness for each community

calc.FR<-function(taxa,groups){
  list()->res
  
  if (ncol(taxa)!=length(groups)) stop("Trait and taxonomic matrix have different number of taxa")
  unique(groups)->gr.names
  specnumber(taxa)->res$nbsp

  res$ab.fgrs<-res$ric.fgrs<-data.frame(matrix(NA,nrow(taxa),length(unique(groups))))
  colnames(res$ab.fgrs)<-colnames(res$ric.fgrs)<-gr.names
  
  j<-0
  
  for (i in gr.names) {
    j<-j+1
    if(is.vector(taxa[,which(groups==i)])==T) decostand(taxa[,which(groups==i)],"pa")->res$ric.fgrs[,j] else specnumber(taxa[,which(groups==i)])->res$ric.fgrs[,j]
    }
  
  j<-0
  
  for (i in gr.names) {
  j<-j+1
    if(is.vector(taxa[,which(groups==i)])==T) taxa[,which(groups==i)]->res$ab.fgrs[,j] else rowSums(taxa[,which(groups==i)])->res$ab.fgrs[,j]
  }
  
  apply(res$ric.fgrs,1,sd)->res$sd.fgr
  specnumber(res$ric.fgrs)->res$FGR # functional group richness estimation
  res$FR<-res$nbsp/res$FGR
  res$FR.ab<-rowSums(taxa)/res$FGR
  return(res)
}

# p.val() estimates p-value and the z-score for a given empirical value and null distribution
#
# Inputs:
# obs: empirical value
# null.dist: null distribution of values
# alternative: how empirical and null values were compared. If "two.sided" is selected a bilateral
# comparison is done. "less" provides the probability of obs to be greater than the null values,
# whilst when "greater" the probability of obs to be lower than the null values is estimated
#
# Outputs:
# $z.score: z-score
# $p.value: p-value
# $alternative: how empirical and null values were compared

p.val<-function (obs, null.dist, alternative = c("two.sided", "less", "greater"))
{ 
  alternative <- match.arg(alternative)
  nsimul=length(null.dist)
  z <- (obs - mean(null.dist))/sd(null.dist) # z-statistic
  pless <- sum(obs <= null.dist, na.rm = TRUE) # null values less or equal than real observation
  pmore <- sum(obs >= null.dist, na.rm = TRUE) # null values greater or equal than real observation
  
  p <- switch(alternative, two.sided = 2 * pmin(pless, pmore), less = pless, greater = pmore)
  p <- pmin(1, (p + 1)/(nsimul + 1)) 
  res<-list(z.score=z,p.value=p,alternative=alternative)
  #return(res)
  return(res)
}


# null.func.dist() is an ad hoc function that estimates the null distribution for 
# a set of functional variables. In this case, the functional variables are those used in this study
#
# Inputs:
# taxa: community data
# tr_inf: a matrix with the functional trait types (effect, response, effect/response traits)
# tr_all: a matrix with the functional traits
# fpc: functional space
# m: functional axes to select (For FRic)
# runs: number of randomizations
#
# Outputs:
#
# $fric.s: functional richness null distribution
# $feve.s: functional evenness null distribution
# $fdis.e.s: functional dispersion of the effect traits null distribution
# $fdis.r.s: response diversity null distribution
# $FR.s: functional redundancy (divisive) null distribution
# $FRa.s: functional redundancy (additive) null distribution

null.func.dist<-function(taxa,tr_inf,tr_all,fpc,m,runs=99){
  
  rownames(tr_all)->spp
  
  e_r<-tr_all[,which(tr_inf[,2]!="R")]
  r_r<-tr_all[,which(tr_inf[,2]!="E")]
  diversity(taxa,index="simpson")->tax.simp_r
  colnames(taxa)<-rownames(e)<-rownames(r)<-rownames(tr_all)<-unlist(spp)
  
fric.s<-feve.s<-fdis.e.s<-fdis.r.s<-FR.s<-FRa.s<-matrix(NA,nrow(taxa),runs)

    for (i in 1:runs){
    
    sample(1:nrow(e_r))->ran
    
    fric_3d(taxa,fpc[ran,],m,1)->fric.tmp
    fric.tmp[which(is.na(fric.tmp))]<-0
    fric.tmp->fric.s[,i]
    
    e_r[ran,]->e.s
    r_r[ran,]->r.s
    rownames(e.s)<-rownames(r.s)<-rownames(e)
    
    gowdis(e.s)->e.dist
    gowdis(r.s)->r.dist
    dbFD(tr_all,as.matrix(taxa),corr="none",calc.FRic = F,calc.CWM = F,calc.FDiv = F,messages = F)->fun.feat
    as.numeric(fun.feat$FEve)->feve.s[,i]
    
    fdisp(e.dist,as.matrix(taxa))$FDis->fdis.e.s[,i]
    fdisp(r.dist,as.matrix(taxa))$FDis->fdis.r.s[,i]
    
    FRa.s[,i]<-tax.simp_r-fdis.e.s[,i]
    
    e.clust <- hclust(e.dist, method = "ward.D")
    e.gr.s <- cutree(e.clust, k = 5)
    
    calc.FR(taxa,e.gr.s)$FR->FR.s[,i]
    }

return(list(fric.s=fric.s,feve.s=feve.s,fdis.e.s=fdis.e.s,fdis.r.s=fdis.r.s,FR.s=FR.s,FRa.s=FRa.s))
}

# null.dist.fr() is an ad hoc function that check the non-randomness of a the relatonship
# between Functional redundancy and a putative environmental filter
#
# In each randomization, Functional Gruoups are assigned randomly to each species
# Inputs:
# taxa: community data
# groups: a grouping vector with the Funtional Groups for each taxon
# runs: number of randomizations
#
# Outputs:
#
# $obs.coef: the empirical model coefficients
# $FR.null: the FR null distribution
# $null.coef: the null model coefficients
# $pval: table containing the p-values and z-scores for the comparison between empirical 
#        and null coefficients


null.dist.fr<-function(taxa,groups,predictor,runs){
  
  res<-list()
  # Groups as vector
  groups<-as.numeric(groups)
  
  # Calculating observed coefficients
  calc.FR(taxa,groups)$FR->FR.obs
  obs.mod<-glm(FR.obs ~ predictor)
  res$obs.coef<-obs.mod$coef
  
  # Creating an empty matrix for null FR values
  res$FR.null<-matrix(rep(NA,nrow(taxa)*runs),ncol=runs)
  colnames(res$FR.null)<-1:runs
  rownames(res$FR.null)<-rownames(taxa)
  # Creating an empty matrix for null model coefficients
  res$null.coef<-matrix(rep(NA,length(obs.mod$coef)*runs),ncol=runs)
  colnames(res$null.coef)<-1:runs
  rownames(res$null.coef)<-names(obs.mod$coef)
  
  for (i in 1:runs){
    sample(groups)->groups
    calc.FR(taxa,groups)$FR->FR.sim->res$FR.null[,i] 
    null.mod<-glm(FR.sim ~ predictor)
    null.mod$coef->res$null.coef[,i]
  }
  
  pval.res<-matrix(rep(NA,2*2),ncol=2); colnames(pval.res)<-c("z-score","pval");rownames(pval.res)<-names(obs.mod$coef)
  pval.res[1,]<-c(p.val(res$obs.coef[1],res$null.coef[1,])$z.score,p.val(res$obs.coef[1],res$null.coef[1,])$p.value)
  pval.res[2,]<-c(p.val(res$obs.coef[2],res$null.coef[2,])$z.score,p.val(res$obs.coef[2],res$null.coef[2,])$p.value)
  
  pval.res->res$pval
  
  # Histogram for intercept and coefficients
  par(mar = c(5, 4.25, 2.5, 2))
  par(mfrow=c(1,2))
  
  hist(res$null.coef[1,],lwd=4,main="Intercept",ylab= "Frequency", xlab="", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,xlim=c(min(res$obs.coef[1],res$null.coef[1,]),max(res$obs.coef[1],res$null.coef[1,])))
  abline(v=res$obs.coef[1], lty=3, lwd=7, col="gray")
  
  hist(res$null.coef[2,],lwd=4,main=names(null.mod$coef)[2],ylab= "Frequency", xlab="", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,xlim=c(min(res$obs.coef[2],res$null.coef[2,]),max(res$obs.coef[2],res$null.coef[2,])))
  abline(v=res$obs.coef[2], lty=3, lwd=7, col="gray")
  

  par(mfrow=c(1,1))
  return(res)
}


fdisp2<-function (d, a, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  if (!is.matrix(a)) 
    stop("'a' must be a matrix.")
  if (ncol(a) != n) 
    stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  if (any(sn.d != sn.a)) 
    stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
         "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    warning("At least one community has zero-sum abundances (no species).", 
            "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    stop("At least one species does not occur in any community (zero total abundance across all communities).", 
         "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  dimnames(vectors) <- list(colnames(a), NULL)
  pos <- eig > 0
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}

fdisp_k<-function (d, a, m, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  if (ncol(a) != n) 
    stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  if (any(sn.d != sn.a)) 
    stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
         "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    warning("At least one community has zero-sum abundances (no species).", 
            "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    warning("At least one species does not occur in any community (zero total abundance across all communities).", 
            "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  dimnames(vectors) <- list(colnames(a), NULL)
  pos <- eig > 0
  if (m>0) pos<-c(pos[1:m],rep(F,length(pos)-m))
  
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}

fdisp_k_sub<-function (d, a, tax_sub=NULL, m, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  #if (ncol(a) != (n-length(tax_sub))) 
  #  stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  #if (any(sn.d != sn.a)) 
  #  stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
  #       "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    warning("At least one community has zero-sum abundances (no species).", 
            "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    warning("At least one species does not occur in any community (zero total abundance across all communities).", 
            "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  rownames(vectors) <- attr(d, "Labels")
  vectors <- vectors[(intersect(rownames(vectors), tax_sub)), 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  pos <- eig > 0
  if (m>0) pos<-c(pos[1:m],rep(F,length(pos)-m))
  
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}

# null.func.dist() is an ad hoc function that estimates the null distribution for 
# a set of functional variables. In this case, the functional variables are those used in this study
#
# taxa: community data
# traits: trait data
# m: number of axes
# obs_dat: observed empirical data (responses, fixed and random factors)
# fun_var: functional variable (only FRic_3d is allowed)
# fixed_var: fixed factors
# random_var: random factors
# runs: number of null  model runs
# num.plots: number of panels
# plot_w: plot width
# plot_h: plot height
# prec: parameter indicating convex hull precision "Qt" or"Qj"

# Inputs:
# taxa: community data
# tr_inf: a matrix with the functional trait types (effect, response, effect/response traits)
# tr_all: a matrix with the functional traits
# fpc: functional space
# m: functional axes to select (For FRic)
# runs: number of randomizations
#
# Outputs:
#
# ses and pvalue foreach model coefficient


null.dist.lmer<-function(taxa,traits,m=2,traits.blo,obs_dat,fun_var,fixed_var,random_var=NULL,runs,num.plots=c(1,1),
                         plot_w=5, plot_h=5, prec="Qt"){
  
  obs.formula<-as.formula(paste(fun_var, paste(c(fixed_var,random_var), collapse=" + "), sep=" ~ "))
  
  res<-list()
  
  obs.mod<-lmer(obs.formula,REML=F,data=obs_dat,na.action="na.omit")
  
  res$obs.coef<-summary(obs.mod)$coefficients[,1]
  
  # Creating an empty matrix for null model coefficients
  res$null.coef<-matrix(rep(NA,length(res$obs.coef)*runs),ncol=runs)
  
  colnames(res$null.coef)<-1:runs
  rownames(res$null.coef)<-names(res$obs.coef)
  
  for (i in 1:runs){
    obs_dat->ran_dat
    sample(1:nrow(traits))->ran
    traits[ran,]->traits.r
    rownames(traits.r)<-rownames(traits)
    
    r.r <- prep.fuzzy(traits.r, traits.blo)
    
    # Combining the traits
    r.ktab.r<-ktab.list.df(list(r.r))
    tr.dist.r <- dist.ktab(r.ktab.r, "F")
    
    dudi.pco(tr.dist.r,scannf = F,nf=10)->tr.pco.r
    
    if (fun_var=="RD.ab") fdisp_k(tr.dist.r,taxa,m=7)$FDis->ran_dat$RD.ab
    if (fun_var=="FRic") fric_3d(taxa, tr.pco.r$li, m=3, prec)->ran_dat$FRic
    if (fun_var=="FRic") min(ran_dat$FRic,na.rm = T)->ran_dat$FRic[which(is.na(ran_dat$FRic)==T)]
      
    ran.formula<-obs.formula
    
    null.mod<-lmer(ran.formula,REML=F,data=ran_dat)
    
    res$null.coef[,i]<-summary(null.mod)$coefficients[,1]
    cat("Randomisation number", i, "\n")
  }
  
  pval.res<-matrix(NA,length(res$obs.coef),ncol=2)
  colnames(pval.res)<-c("ses","pval")
  rownames(pval.res)<-names(res$obs.coef)
  
  for (a in 1:nrow(pval.res)) pval.res[a,]<-c(p.val(res$obs.coef[a],res$null.coef[a,])$z.score,p.val(res$obs.coef[a],res$null.coef[a,])$p.value)
  
  pval.res->res$pval
  
  # Histogram for intercept and coefficients
  pdf(file=paste(fun_var,"null",".pdf",sep="_"),onefile=FALSE,width=plot_w,height=plot_h)
  
  par(mfrow=num.plots, mar = c(5, 4.25, 2.5, 2),cex.axis=1.5,cex.lab=1.5,cex.main=1.75)
  
  for (a in 1:nrow(pval.res)){
    hist(res$null.coef[a,],lwd=4,main=rownames(pval.res)[a],ylab= "Frequency", xlab="",xlim=c(min(res$obs.coef[a],res$null.coef[a,]),max(res$obs.coef[a],res$null.coef[a,])))
    abline(v=res$obs.coef[a], lty=2, lwd=5, col="gray")
  }
  
  dev.off()
  par(mfrow=c(1,1))
  return(res)
}

### Functions to prepare matrices, construct the functional space and to estimate functional diversity components
prep.fuzzy.df<-function (traits, col.blocks) 
{
  if (!is.data.frame(traits)) 
    stop("Data.frame expected")
  if (sum(col.blocks) != ncol(traits)) {
    stop("Non convenient data in col.blocks")
  }
  if (is.null(names(col.blocks))) {
    names(col.blocks) <- paste("FV", as.character(1:length(col.blocks)), sep = "")
  }
  f1 <- function(x) {
    a <- sum(x)
    if (is.na(a)) 
      return(rep(0, length(x)))
    if (a == 0) 
      return(rep(0, length(x)))
    return(x/a)
  }
  k2 <- 0
  col.w <- rep(1, ncol(traits))
  
  for (k in 1:(length(col.blocks))) {
    k1 <- k2 + 1
    if (col.blocks[k]==1) k2<-k1 else k2 <- k2 + col.blocks[k]
    X <- as.matrix(traits[, k1:k2])
    if (col.blocks[k]==1) X[which(X[,1]>0),]<-1 else X <- t(apply(X, 1, f1))
    X.marge <- apply(X, 1, sum)
    X.marge <- X.marge
    X.marge <- X.marge/sum(X.marge)
    X.mean <- apply(X * X.marge, 2, sum)
    nr <- sum(X.marge == 0)
    cat(nr, "missing data found in block", k, "\n")
    traits[, k1:k2] <- X
    col.w[k1:k2] <- X.mean
  }
  attr(traits, "col.blocks") <- col.blocks
  attr(traits, "col.freq") <- col.w
  col.num <- factor(rep((1:length(col.blocks)), col.blocks))
  attr(traits, "col.num") <- col.num
  return(traits)
}

feve_k<-function(fpc,taxa,m,rdt){
  
  # Creating variables
  nrow(taxa)->c
  FEve <- rep(NA, c) ; names(FEve) <- row.names(taxa)
  nb.sp<-specnumber(taxa)
  
  # generating the taxonomic matrix arranged according to the replicated trait values
  tax.pool<-ncol(taxa)
  taxa.rep<-data.frame(matrix(NA,nrow(taxa),tax.pool*rdt))
  spp.list<-c(1:(tax.pool*rdt))
  
  for (spp in 1:tax.pool) {paste(rep("spp",rdt),spp,sep="")->spp.list[((spp-1)*rdt+1):(spp*rdt)]}
  
  colnames(taxa.rep)<-spp.list
  
  for (spp in 1:tax.pool){taxa.rep[,((spp-1)*rdt+1):(spp*rdt)]<-taxa[,spp]/rdt}                             
  
  # Estimating Functional Evenness for each community
  
  for (i in 1:c) {
    sppres <- which(taxa.rep[i, ] > 0)
    # number of species in the community
    S <- length(sppres)
    ab <- as.matrix(taxa.rep[i, sppres])
    # scaling of abundances
    abundrel <- ab / sum(ab)
    
    # selecting the c
    tr <- data.frame(fpc[sppres,1:m ])
    
    if (nb.sp[i] > 2) {
      tr.dist <- dist(tr)
      linkmst <- mst(tr.dist)
      mstvect <- as.dist(linkmst)
      abund2 <- matrix(0, nrow = S, ncol = S)
      for (q in 1:S) for (r in 1:S) abund2[q, r] <- abundrel[q] +  abundrel[r]
      abund2vect <- as.dist(abund2)
      EW <- rep(0, S - 1)
      flag <- 1
      for (mv in 1:((S - 1) * S/2)) {
        if (mstvect[mv] != 0) {
          EW[flag] <- tr.dist[mv]/(abund2vect[mv])
          flag <- flag + 1
        }
      }
      minPEW <- rep(0, S - 1)
      OdSmO <- 1/(S - 1)
      for (l in 1:(S - 1)) minPEW[l] <- min((EW[l]/sum(EW)), 
                                            OdSmO)
      FEve[i] <- ((sum(minPEW)) - OdSmO)/(1 - OdSmO)
    } else FEve[i] <- NA
  }
  return(FEve)
}



############### extra

fdis<-function(taxa,fpc,m,rdt) { ## taxa: fauna, fpc: functional principal components
  if (m>1) combn(m,2)->axis.pair else matrix(c(1,1))->axis.pair ## a matrix containing the different combination of axis pairs
  if (m>1) ncol(axis.pair)->m.ncol else m.ncol<-1
  com.fdis<-rep(0,m.ncol)
  apply(taxa,1,function(com.row) 
  {
    taxon.ric<-sum(decostand(com.row,method="pa"))
    
    # selecting the coordinates for the community of the axis p.row, including all random trials
    com.row[which(com.row>0)]->com.abun ## abundances of the community
    for (p.row in 1:m.ncol){
      fpc[which(com.row>0),(((axis.pair[1,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[1,p.row]-1)+rdt)]->a
      if (m>1) fpc[which(com.row>0),(((axis.pair[2,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[2,p.row]-1)+rdt)]->b
      
      axis.a<-axis.pair[1,p.row]
      if (m>1) axis.b<-axis.pair[2,p.row]
      
      # weighted centroids
      if (taxon.ric<2) mean(a)->cent.a else sum(apply(a,2,function(x){x*t(com.abun/(sum(com.abun)*rdt))}))->cent.a ## weighted by 
      if (m>1) if (taxon.ric<2) mean(b)->cent.b else sum(apply(b,2,function(x){x*t(com.abun/(sum(com.abun)*rdt))}))->cent.b
      
      (a-cent.a)^2->a.temp
      if (m>1) (b-cent.b)^2->b.temp else b.temp<-0
      a.temp+b.temp->sq.dist
      if (taxon.ric<2) mean(sq.dist)->w.fdis else sum(apply(sq.dist,2,function(x){x*t(com.abun/(sum(com.abun)*rdt))}))->w.fdis
      
      
      w.fdis->com.fdis[p.row]
    }
    return(FDis<-sum(com.fdis))
  }
  )
}