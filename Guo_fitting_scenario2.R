
mn <- 1 # parameter for cluster parallel running

##############################################################
###### Linear Model Methods with LASSO - 4 biomarkers
##############################################################

library(MASS)
library(Matrix)
library(glmnet)
library(grpreg)
library(grpregOverlap)
library(vennLasso)
library(Hmisc)
library(Formula)
library(ggplot2)
library(rms)
library(erer)
library(dplyr)
library(boot)
library(caret)


### directories ###
name <- c("datas2")
out.name <- c("fulls2_trt")

data.dir <- data.dir.ext <- paste("//work//users//b//e//beibo//Research//paper2//",name,"//",sep="")
out.dir <- paste("//work//users//b//e//beibo//Research//paper2//",out.name,"//",sep="")

data.dir <- data.dir.ext <- c("C://Users//CSCC//OneDrive - University of North Carolina at Chapel Hill//School//Research//thesis//Paper 2//datas2//")




### parameters ###
n1=400 # sample size 
S <- 1000 # number of datasets
B <- 100 # number of bootstraps 

fold.num.cv <- 5 # number of cvs
alpha <- 0.05  # error

# penalty.type <- c("grLasso")
penalty.type <- c("gel")
# penalty.type <- c("grMCP")
# penalty.type <- c("grSCAD")

family.type <- c("gaussian")

# CV tuning? 1 no 0 yes(cv.min) 2 yes(cv.1se)
lambda.0 <- 0
# if no cv tuning (lambda.0 = 1), what pre-specified param to use?  
lambda.pre <- 0.2

# prevalence cutoff
prev.cutoff.lower <- 0
prev.cutoff.upper <- 2

# fold number to determine tunning parameter in fitting
# fold number for cross-validation to determine cv-adjusted estimators
fold.num <- 10
# set prev to 0 if not passing threshold
trt.max <- 0
# trt cutoff for subgroup selection
trt.cutoff <- 0.5
# effect size cutoff for subgroup selection
eff.cutoff <- 0.2

u.type <- 2
gamma <- 0.5







### cross-validation function ###
lm.cv <- function(n1.lm, cv.index, data.in) {

  
  lm.x.mat <- (x.mat[c(cv.index),])
  lm.x.mat.fit.t <- (x.mat.fit.t[c(cv.index),])
  lm.ron1 <- (data.in[c(cv.index),])
  
  ab.1 <-  rep(NA,n1.lm)
  y1.fit <- lm.ron1[,1]
  
  
  if (lambda.0 == 1) {lambda.1 = lambda.pre} 
  if (lambda.0 == 2) {
    cv.result <- cv.grpregOverlap(X = lm.x.mat, y = y1.fit, group = grp.mat, family = family.type, nfolds = fold.num, penalty = penalty.type)
    # lambda.1se
    lambda.1 <- max(cv.result$lambda[which(cv.result$cve <= cv.result$cve[cv.result$min] + cv.result$cvse[cv.result$min])])
    if (lambda.1 == max(cv.result$lambda)) {
      lambda.1 <- max(cv.result$lambda[which(cv.result$cve <= cv.result$cve[cv.result$min] + cv.result$cvse[cv.result$min] & (cv.result$lambda != lambda.1))])
    }
  } else {
    lambda.1 = cv.grpregOverlap(X = lm.x.mat, y = y1.fit, group = grp.mat, family = family.type, nfolds = fold.num, penalty = penalty.type)$lambda.min
  }
  
  fit.1=grpregOverlap(X = lm.x.mat, y = y1.fit, group = grp.mat, penalty = penalty.type, family = family.type, lambda = lambda.1)
  
  
  
  lm.x.mat <- data.frame(x.mat[c(cv.index),])
  lm.x.mat.fit.t <- data.frame(x.mat.fit.t[c(cv.index),])
  lm.ron1 <- data.frame(data.in[c(cv.index),])
  
  
  
  ab.in <- as.matrix(lm.x.mat.fit.t) %*% as.matrix(fit.1$beta)
  y1.type = fit.1$family
  
  ab.in <- signif(ab.in, digits = 8)
  ab.1 <- ab.in
  
  p.1.a <- p.2.a <- p.3.a <- rep(NA,n1.lm)  
  p.raw <- p.raw.index <- rep(NA,n1.lm)
  
  # actual response to determine subgroup
  # >= threshold
  for (i in 1:n1.lm)   {
    sub.index.1 <- which(ab.1 >= ab.1[i])
    sub.index <- list(sub.index.1)
    p.raw.temp <- rep(NA,1)
    
    # cross-test all potential subgroups in all outcomes
    for (m in 1) {
      # coef.t <- unlist(coef[m])
      sub.index.t <- unlist(sub.index[m])
      g.t.t <- sub.index.t[lm.x.mat$t[sub.index.t] == 1]
      g.c.t <- sub.index.t[lm.x.mat$t[sub.index.t] == 0]
      prev <- (length(g.t.t) + length(g.c.t))/n1.lm
      if((length(g.t.t)>=2) & (length(g.c.t)>=2) & 
  
         ((length(g.t.t) + length(g.c.t))/n1.lm >= prev.cutoff.lower) & 
         ((length(g.t.t) + length(g.c.t))/n1.lm <= prev.cutoff.upper)
  
      ) {

        
        if (u.type == 1) {
          score.adj = 0
          p.1.a[i] <- z.prop(sum(y1.fit[g.c.t]), sum(y1.fit[g.t.t]), length(y1.fit[g.c.t]), length(y1.fit[g.t.t])) / (prev)^score.adj
          
        } else if (u.type == 2) {
          score.adj = gamma
          trt.eff.1 <- mean(y1.fit[g.c.t]) - mean(y1.fit[g.t.t])
          p.1.a[i] <- trt.eff.1*(prev)^score.adj
          
        } 
        
      } else {
        p.1.a[i] <-NA
        
      }

      p.raw.temp[m] <- c(p.1.a[i])[m]
    }
    
    if (sum(is.na(p.raw.temp)) == 1) {
      p.raw.index[i] <- NA
      p.raw[i] <- NA
    } else {

      p.raw.index[i] = which.max(p.raw.temp)   # which subgroup gives the best p value 
      p.raw[i] = p.raw.temp[p.raw.index[i]]  # best p value among the 3 subgroups generated from 3 outcomes, respectively
    }
  }  
  p.raw.index.grt <- p.raw.index
  p.raw.grt <- p.raw
  
  # < threshold   
  for (i in 1:n1.lm)   {
    sub.index.1 <- which(ab.1 < ab.1[i])
    sub.index <- list(sub.index.1)
    p.raw.temp <- rep(NA,1)
    
    # cross-test all potential subgroups in all outcomes
    for (m in 1) {
      # coef.t <- unlist(coef[m])
      sub.index.t <- unlist(sub.index[m])
      g.t.t <- sub.index.t[lm.x.mat$t[sub.index.t] == 1]
      g.c.t <- sub.index.t[lm.x.mat$t[sub.index.t] == 0]
      prev <- (length(g.t.t) + length(g.c.t))/n1.lm
      if((length(g.t.t)>=2) & (length(g.c.t)>=2) & 

         ((length(g.t.t) + length(g.c.t))/n1.lm >= prev.cutoff.lower) & 
         ((length(g.t.t) + length(g.c.t))/n1.lm <= prev.cutoff.upper)

      ) {
 
        
        if (u.type == 1) {
          score.adj = 0
          p.1.a[i] <- z.prop(sum(y1.fit[g.c.t]), sum(y1.fit[g.t.t]), length(y1.fit[g.c.t]), length(y1.fit[g.t.t])) / (prev)^score.adj
          
        } else if (u.type == 2) {
          score.adj = gamma
          trt.eff.1 <- mean(y1.fit[g.c.t]) - mean(y1.fit[g.t.t])
          p.1.a[i] <- trt.eff.1*(prev)^score.adj
          
        }    
        
      } else {
        p.1.a[i] <-NA
        
      }
      p.raw.temp[m] <- c(p.1.a[i])[m]
    }
    
    if (sum(is.na(p.raw.temp)) == 1) {
      p.raw.index[i] <- NA
      p.raw[i] <- NA
    } else {
      p.raw.index[i] = which.max(p.raw.temp)   # which subgroup gives the best p value 
      p.raw[i] = p.raw.temp[p.raw.index[i]]  # best p value among the 3 subgroups generated from 3 outcomes, respectively
    }
  }  
  p.raw.index.les <- p.raw.index
  p.raw.les <- p.raw
  
  p.raw.full <- c(p.raw.grt, p.raw.les)
  
  fit.3lm.un <- fit.1
  ab.un <- c(ab.1,ab.1)
  unique.uti <- split(seq_along(p.raw.full), p.raw.full)
  unique.uti.index <- unlist(lapply(unique.uti, function(x){x[1]}))
  unique.uti.sign <- ifelse(unique.uti.index <= n1.lm, "GE","LS")  # 1 great; 0 less
  
  u.cutoff.full <- u.sign.full <- u.utility <- u.trt.diff <- u.prev <- rep(NA, length(unique.uti.index))
  ind.sub.u.list <- vector(mode = "list", length = length(unique.uti.index))
  
  
  for (u in 1:length(unique.uti.index)) {
    index <- unique.uti.index[u]
    sign <- unique.uti.sign[u] 
    ab.fit.un <- ab.un[index]
    
    un.test <- signif(as.matrix(lm.x.mat.fit.t) %*% as.matrix(fit.1$beta), digits = 8)   
    
    un.cutoff <- signif(ab.fit.un, digits = 8)
    
    if (sign == "GE") {
      ind.sub<-which(ifelse(un.test>= un.cutoff ,1,0) == 1) 
    } else {
      ind.sub<-which(ifelse(un.test< un.cutoff ,1,0) == 1)
    } 
    
    ## drop subgroup if not clinicall reasonable ##
    ron1.a <- lm.ron1[ind.sub,]
    
    ind.sub.u.list[[u]] <- ind.sub 
    u.trt <- aggregate(ron1.a[,1] ~ t, data =  ron1.a , mean)
    u.trt.diff[u] <- u.trt[2,2] - u.trt[1,2]  
    u.prev[u] <- length(ind.sub)/n1.lm
    u.utility[u] <- u.trt.diff[u]*(u.prev[u])^gamma  
    u.cutoff.full[u] <- un.cutoff
    u.sign.full[u] <- sign
  }
  
  u.data <- data.frame(rbind(cbind(u.trt.diff, u.prev, u.prev*n1.lm, u.utility,u.cutoff.full,u.sign.full)))
  fit.beta <- as.matrix(fit.1$beta)
  
  out <- list(u.data, fit.1, ind.sub.u.list)
  return(out)
}










### matrix initialization ###
trt.eff.emp <- prev.emp <- gamma.emp <- bs.trt.sd <- rep(NA,S)
prev.full <-  rep(NA,S)
trt.eff.full.emp <-  trt.eff.full.emp.sd <- gamma.full.emp <-  gamma.full.emp.sd <- rep(NA,S)
trt.eff.full.reg <-  trt.eff.full.reg.sd <- gamma.full.reg <-  gamma.full.reg.sd <- rep(NA,S)
trt.eff.full.aug <-  trt.eff.full.aug.sd <- gamma.full.aug <-  gamma.full.aug.sd <- rep(NA,S)

prev.full.ext <-  rep(NA,S)
trt.eff.ext  <- rep(NA,S)
gamma.ext <- rep(NA,S)

prev.cv <-  rep(NA,S)
trt.eff.cv.emp <-  gamma.cv.emp <-  rep(NA,S)
trt.eff.cv.reg <-  gamma.cv.reg <-  rep(NA,S)
trt.eff.cv.aug <-  gamma.cv.aug <-  rep(NA,S)
gamma.cv.emp.sd <- gamma.cv.reg.sd <- gamma.cv.aug.sd <- rep(NA,S)
trt.cv.emp.sd <- trt.cv.reg.sd <- trt.cv.aug.sd <- rep(NA,S)

r.val <- trt.eff.guo <- trt.sd.guo <- gamma.guo <- gamma.sd.guo <- rep(NA,S)

trt.eff.guo.reg <- gamma.guo.reg <- gamma.sd.guo.reg <- rep(NA,S)
trt.eff.guo.aug <- gamma.guo.aug <- gamma.sd.guo.aug <- rep(NA,S)

trt.guo.reg <- trt.sd.guo.reg <- rep(NA,S)
trt.guo.aug <- trt.sd.guo.aug <- rep(NA,S)


bs.trt.eff <- vector(mode = "list", length = S)




for(s in 1:S)
{
  
  set.seed(s+(mn-1)*S)
  
  datnum <- c(s+(mn-1)*S)
  
  cat("s:",s,"\n")
  
  dat <-  readRDS(file = paste(data.dir, "data_", datnum, ".Rds", sep=""))
  
  t <- dat[,2]
  biom1.x1 <- dat[,3]
  biom1.x2 <- dat[,4]  
  biom1.x3 <- dat[,5]
  biom1.x4 <- dat[,6]
  biom1.x5 <- dat[,7]
  biom1.x6 <- dat[,8]
  biom1.x7 <- dat[,9]
  biom1.x8 <- dat[,10]
  biom1.x9 <- dat[,11] 
  biom1.x0 <- dat[,12]
 
  
  f <- as.formula(y~ t+ biom1.x1+biom1.x2+biom1.x3+biom1.x4+biom1.x5+biom1.x6+biom1.x7+biom1.x8+biom1.x9+biom1.x0+
                    t:biom1.x1+t:biom1.x2+t:biom1.x3+t:biom1.x4+t:biom1.x5+t:biom1.x6+
                    t:biom1.x7+t:biom1.x8+t:biom1.x9+t:biom1.x0)
  
  
  f <-  as.formula(y~(t+biom1.x1+biom1.x2+biom1.x3+biom1.x4+biom1.x5+biom1.x6+biom1.x7+biom1.x8+biom1.x9+biom1.x0)^3) 
  x <- model.matrix(f, dat)[,1:112]  # pairwise interaction - no quadratic
  y <- as.matrix(dat$y, ncol=1)
  
  grp.mat <- c(as.list(as.double(c(1:11))),
               
               Map(c, as.double(c(1)), as.double(c(2:11)), as.double(c(12:21))),
               
               Map(c, as.double(c(2)), as.double(c(3:11)), as.double(c(22:30))),
               Map(c, as.double(c(3)), as.double(c(4:11)), as.double(c(31:38))),
               Map(c, as.double(c(4)), as.double(c(5:11)), as.double(c(39:45))),
               Map(c, as.double(c(5)), as.double(c(6:11)), as.double(c(46:51))),
               Map(c, as.double(c(6)), as.double(c(7:11)), as.double(c(52:56))),
               Map(c, as.double(c(7)), as.double(c(8:11)), as.double(c(57:60))),
               Map(c, as.double(c(8)), as.double(c(9:11)), as.double(c(61:63))),
               Map(c, as.double(c(9)), as.double(c(10:11)), as.double(c(64:65))),
               Map(c, as.double(c(10)), as.double(c(11)), as.double(c(66))),
               
               Map(c, as.double(c(1)), as.double(c(2)), as.double(c(3:11)), as.double(c(22:30)), as.double(c(67:(67+8)))),
               Map(c, as.double(c(1)), as.double(c(3)), as.double(c(4:11)), as.double(c(31:38)), as.double(c((67+8+1)):(67+8+1+7))),
               Map(c, as.double(c(1)), as.double(c(4)), as.double(c(5:11)), as.double(c(39:45)), as.double(c((67+8+1+7+1)):(67+8+1+7+1+6))),
               Map(c, as.double(c(1)), as.double(c(5)), as.double(c(6:11)), as.double(c(46:51)), as.double(c((67+8+1+7+1+6+1)):(67+8+1+7+1+6+1+5))),
               Map(c, as.double(c(1)), as.double(c(6)), as.double(c(7:11)), as.double(c(52:56)), as.double(c((67+8+1+7+1+6+1+5+1)):(67+8+1+7+1+6+1+5+1+4))),
               Map(c, as.double(c(1)), as.double(c(7)), as.double(c(8:11)), as.double(c(57:60)), as.double(c((67+8+1+7+1+6+1+5+1+4+1)):(67+8+1+7+1+6+1+5+1+4+1+3))),
               Map(c, as.double(c(1)), as.double(c(8)), as.double(c(9:11)), as.double(c(61:63)), as.double(c((67+8+1+7+1+6+1+5+1+4+1+3+1)):(67+8+1+7+1+6+1+5+1+4+1+3+1+2))),
               Map(c, as.double(c(1)), as.double(c(9)), as.double(c(10:11)), as.double(c(64:65)), as.double(c((67+8+1+7+1+6+1+5+1+4+1+3+1+2+1)):(67+8+1+7+1+6+1+5+1+4+1+3+1+2+1+1))),
               Map(c, as.double(c(1)), as.double(c(10)), as.double(c(11)), as.double(c(66)), as.double(c((67+8+1+7+1+6+1+5+1+4+1+3+1+2+1+1))+1))
  )
  
  

  x.mat <- x[,-1]
  ron1 <- data.frame(cbind(y,x.mat))
  t.fit <- c(rep(1,n1)) 
  f.fit <-  as.formula(y~(t.fit+biom1.x1+biom1.x2+biom1.x3+biom1.x4+biom1.x5+biom1.x6+biom1.x7+biom1.x8+biom1.x9+biom1.x0)^3) 
  x.mat.fit.t <- model.matrix(f.fit, dat)[,1:112]  # pairwise interaction - no quadratic
  
  c.fit <- c(rep(0,n1)) 
  f.c.fit <-  as.formula(y~(c.fit+biom1.x1+biom1.x2+biom1.x3+biom1.x4+biom1.x5+biom1.x6+biom1.x7+biom1.x8+biom1.x9+biom1.x0)^3) 
  x.mat.fit.c <- model.matrix(f.c.fit, dat)[,1:112]  # pairwise interaction - no quadratic
  
  colnames(ron1)[1] <- c("y")
  
  
  test.full <- full.b <- data.frame(cbind(ron1$y,x.mat[,1],x.mat.fit.t))
  
 

  ##### post-hoc subgroups ######
  # apply to whole sample #
  
  # fit on ron1, old data
  cv.full <- lm.cv(n1, 1:n1, ron1)
  
  
  fit.full <- cv.full[[2]]
  fit.beta.full <- fit.full$beta 
  u.data.full <- cv.full[[1]]
  
  # apply best subgroup def to testing set #
  un.test.full <- signif(as.matrix(test.full[,-c(1:2)]) %*% as.matrix(fit.beta.full), digits = 8)
  full.max.index <-  which.max(u.data.full$u.utility)
  un.cutoff.full <- as.numeric(u.data.full$u.cutoff.full[full.max.index])
  un.sign.full <- u.data.full$u.sign.full[full.max.index]
  if (un.sign.full == "GE") {
    ind.sub.full<-which(ifelse( un.test.full >=  un.cutoff.full ,1,0) == 1) 
  } else {
    ind.sub.full<-which(ifelse( un.test.full <  un.cutoff.full ,1,0) == 1)
  }
  ind.sub.t.full <- ind.sub.full[test.full$V2[ind.sub.full] == 1]
  ind.sub.c.full <- ind.sub.full[test.full$V2[ind.sub.full] == 0]
  
  prev.full[s] <- (length(ind.sub.full))/n1
  
  
  # full-sample empirical estimates #
  trt.eff.full.emp[s] <- mean(test.full$V1[ind.sub.t.full]) - mean(test.full$V1[ind.sub.c.full])
  
  # full-sample regression estimates # 
  un.test.full.c.reg <- signif(as.matrix(x.mat.fit.c) %*% as.matrix(fit.beta.full), digits = 8)
  un.test.full.t.reg <- un.test.full.c.reg + signif(as.matrix(x.mat.fit.t) %*% as.matrix(fit.beta.full), digits = 8)

  trt.eff.full.reg[s] <- (sum(un.test.full.t.reg[ind.sub.full]) - sum(un.test.full.c.reg[ind.sub.full]))/length(ind.sub.full)

  # full-sample augmented estimates #
  trt.eff.full.aug.t <- mean(test.full$V1[ind.sub.t.full]) - 
    sum((un.test.full.t.reg*{(t==1)-(length(ind.sub.t.full)/length(ind.sub.full))})[ind.sub.full]) / length(ind.sub.t.full)
  trt.eff.full.aug.c <- mean(test.full$V1[ind.sub.c.full]) - 
    sum((un.test.full.c.reg*{(t==0)-(length(ind.sub.c.full)/length(ind.sub.full))})[ind.sub.full]) / length(ind.sub.c.full)
  trt.eff.full.aug[s] <- trt.eff.full.aug.t - trt.eff.full.aug.c
  # trt.eff.full.aug.sd[s] <- sqrt(var(trt.eff.full.aug.t)  +  var(trt.eff.full.aug.c))
  
  gamma.full.emp[s] <- trt.eff.full.emp[s]*(prev.full[s])^gamma
  gamma.full.reg[s] <- trt.eff.full.reg[s]*(prev.full[s])^gamma
  gamma.full.aug[s] <- trt.eff.full.aug[s]*(prev.full[s])^gamma
  
  
  
  ## apply to external validating sample for reference estimates ##
  n.obs = 10000
  dat.ext <-  readRDS(file = paste(data.dir.ext, "data_ext_", datnum, ".Rds", sep=""))
  
  y.ext <- dat.ext[,1]
  t.ext <- dat.ext[,2]
  biom1.x1.ext <- dat.ext[,3]
  biom1.x2.ext <- dat.ext[,4]
  biom1.x3.ext <- dat.ext[,5]
  biom1.x4.ext <- dat.ext[,6]
  biom1.x5.ext <- dat.ext[,7]
  biom1.x6.ext <- dat.ext[,8]
  biom1.x7.ext <- dat.ext[,9]
  biom1.x8.ext <- dat.ext[,10]
  biom1.x9.ext <- dat.ext[,11] 
  biom1.x0.ext <- dat.ext[,12]
  
  t.fit.ext <- c(rep(1,n.obs*2)) 
  f.fit.ext <-  as.formula(y.ext~(t.fit.ext+biom1.x1.ext+biom1.x2.ext+biom1.x3.ext+biom1.x4.ext+
                                biom1.x5.ext+biom1.x6.ext+biom1.x7.ext+biom1.x8.ext+biom1.x9.ext+biom1.x0.ext)^3) 
  x.mat.ext.t <- model.matrix(f.fit.ext, dat.ext)[,1:112]
  
  un.test.full.ext <- signif(as.matrix(x.mat.ext.t) %*% as.matrix(fit.beta.full), digits = 8)
  if (un.sign.full == "GE") {
    ind.sub.full.ext<-which(ifelse( un.test.full.ext >=  un.cutoff.full ,1,0) == 1) 
  } else {
    ind.sub.full.ext<-which(ifelse( un.test.full.ext <  un.cutoff.full ,1,0) == 1)
  }
  ind.sub.t.full.ext <- ind.sub.full.ext[dat.ext$t.ext[ind.sub.full.ext] == 1]
  ind.sub.c.full.ext <- ind.sub.full.ext[dat.ext$t.ext[ind.sub.full.ext] == 0]
  
  prev.full.ext[s] <- (length(ind.sub.full.ext))/(2*n.obs)
  trt.eff.ext[s]  <- mean(dat.ext$y[ind.sub.t.full.ext]) - mean(dat.ext$y[ind.sub.c.full.ext])
  gamma.ext[s] <-  trt.eff.ext[s]*(prev.full.ext[s])^gamma

  
  
  ### bootstrap inference ###
  bootstrap.full <- function(data, i, d.beta){
    d <- data[i,]
    
    fit.t.index <- 2:(dim(ron1)[2]+1)
    fit.c.index <- (dim(ron1)[2]+2):(dim(ron1)[2]+1+dim(ron1)[2])
    x.index <- (dim(ron1)[2]+2+dim(ron1)[2]):(dim(d)[2])


    ## Zhang, done separately, missing here just as placeholder ##
    
    fold.index <- createFolds(d[,1], k = fold.num.cv, list = TRUE, returnTrain = FALSE)
    
    gamma.j.zhang.emp <- gamma.j.zhang.reg <- gamma.j.zhang.aug <-rep(NA, fold.num.cv)
    trt.j.zhang.emp <- trt.j.zhang.reg <- trt.j.zhang.aug <- rep(NA, fold.num.cv)
    
    gamma.zhang.emp <- mean(gamma.j.zhang.emp,na.rm = TRUE)
    gamma.zhang.reg <- mean(gamma.j.zhang.reg,na.rm = TRUE)
    gamma.zhang.aug <- mean(gamma.j.zhang.aug,na.rm = TRUE)
    
    trt.zhang.emp <- mean(trt.j.zhang.emp,na.rm = TRUE)
    trt.zhang.reg <- mean(trt.j.zhang.reg,na.rm = TRUE)
    trt.zhang.aug <- mean(trt.j.zhang.aug,na.rm = TRUE)
    
    out.zhang <- c(gamma.zhang.emp,gamma.zhang.reg,gamma.zhang.aug,
                   trt.zhang.emp,trt.zhang.reg,trt.zhang.aug)

    
    ## Guo ##
    
    # empirical
    cv.index <- 1:dim(d)[1]
    n1 <- length(cv.index)
    beta.b <- beta.b.trt <- rep(NA,length(u.data.full$u.utility))
    
    cv.full.reg <- lm.cv(n1, 1:n1, d[,c(1,x.index)]) #refit the model to get new set of regression param estimators.
    fit.beta.full <- cv.full.reg[[2]]$beta # bootstrap param est., apply to subgroups from main data 
    
    un.test.b <- signif(as.matrix(d[,fit.t.index]) %*% as.matrix(fit.beta.full), digits = 8) 
    for (s.b in 1:length(beta.b)) {
      if (u.data.full$u.sign.full[s.b] == "GE") {
        ind.sub.b<-which(ifelse( un.test.b >= as.numeric(u.data.full$u.cutoff.full[s.b]) ,1,0) == 1) 
      } else {
        ind.sub.b<-which(ifelse( un.test.b < as.numeric(u.data.full$u.cutoff.full[s.b]) ,1,0) == 1)
      }
      ind.sub.t.b <- ind.sub.b[d$t[ind.sub.b] == 1]
      ind.sub.c.b <- ind.sub.b[d$t[ind.sub.b] == 0]
      trt.eff.b <- mean(d[ind.sub.t.b,1]) - mean(d[ind.sub.c.b,1])
      prev.b <- (length(ind.sub.b))/n1

        beta.b[s.b] <- trt.eff.b*(prev.b)^gamma

        beta.b.trt[s.b] <- trt.eff.b

    }
    gamma.guo.emp <- max(beta.b+d.c.emp, na.rm = TRUE)
    trt.guo.emp <- max(beta.b.trt+d.c.emp.trt, na.rm = TRUE)
    
    
    # regression
    cv.index <- 1:dim(d)[1]
    n1 <- length(cv.index)

    beta.b <- beta.b.trt <- rep(NA,length(u.data.full$u.utility))
    gamma.g.sub <- vector(mode = "list", length = length(u.data.full$u.utility))
    
    for (u in 1:length(u.data.full$u.utility)) {
      gamma.g.sub[[u]] <-  signif(as.matrix(d[(cv.full[[3]][[u]]),fit.t.index]) %*% as.matrix(fit.beta.full), digits = 8)
    }
    trt.eff.b <-  unlist(lapply(gamma.g.sub, mean))
    prev.b <- as.numeric(u.data.full$u.prev)

      beta.b <- trt.eff.b*(prev.b)^gamma

      beta.b.trt <- trt.eff.b

    gamma.guo.reg <- max(beta.b+d.c.reg, na.rm = TRUE)
    trt.guo.reg <- max(beta.b.trt+d.c.reg.trt, na.rm = TRUE)
    
    
    
    # augmented
    cv.index <- 1:dim(d)[1]
    n1 <- length(cv.index)
    x.mat.fit.t <- d[,fit.t.index]
    x.mat.fit.c <- d[,fit.c.index]
    x.mat <- d[, x.index]
    
    beta.b <- beta.b.trt <- rep(NA,length(u.data.full$u.utility))
    gamma.g.sub.aug.c <-  gamma.g.sub.aug.t <-  vector(mode = "list", length = length(u.data.full$u.utility))
    
    gamma.g.sub.c.full <- signif(as.matrix(d[,(dim(ron1)[2]+2):(dim(ron1)[2]+1+dim(ron1)[2])]) %*% as.matrix(fit.beta.full), digits = 8)
    gamma.g.sub.t.full <- gamma.g.sub.c.full + signif(as.matrix(d[,2:(dim(ron1)[2]+1)]) %*% as.matrix(fit.beta.full), digits = 8)
    
    for (u in 1:length(u.data.full$u.utility)) {
      ## use bootstrap re-fitted parameters to calculate regression estimators using main data
      # augmented estimator
      ind.sub.full <- (cv.full[[3]][[u]])
      ind.sub.t.full <- ind.sub.full[d$t[ind.sub.full] == 1]
      ind.sub.c.full <- ind.sub.full[d$t[ind.sub.full] == 0] 
      gamma.g.sub.aug.t[[u]] <- mean(d[ind.sub.t.full,1]) - 
        sum((gamma.g.sub.t.full*{(d$t==1)-(length(ind.sub.t.full)/length(ind.sub.full))})[ind.sub.full]) / length(ind.sub.t.full)
      gamma.g.sub.aug.c[[u]] <- mean(d[ind.sub.c.full,1]) - 
        sum((gamma.g.sub.c.full*{(d$t==0)-(length(ind.sub.c.full)/length(ind.sub.full))})[ind.sub.full]) / length(ind.sub.c.full)
    }
    trt.eff.aug.t <- unlist(lapply(gamma.g.sub.aug.t, mean))
    trt.eff.aug.c <- unlist(lapply(gamma.g.sub.aug.c, mean))
    trt.eff.aug.diff <- trt.eff.aug.t  -  trt.eff.aug.c
    prev.b <- as.numeric(u.data.full$u.prev)

      beta.b <- trt.eff.aug.diff*(prev.b)^gamma

      beta.b.trt <- trt.eff.aug.diff

    gamma.guo.aug <- max(beta.b+d.c.aug, na.rm = TRUE)
    trt.guo.aug <- max(beta.b.trt+d.c.aug.trt, na.rm = TRUE)
    
    out.guo <- c(gamma.guo.emp,gamma.guo.reg,gamma.guo.aug,
                 trt.guo.emp, trt.guo.reg, trt.guo.aug)
    
    out <- c(out.zhang, out.guo)
    
    return(out)
  }
  

  
  
  full.b.aug <- data.frame(cbind(ron1$y,x.mat.fit.t,x.mat.fit.c,x.mat))

  
  
  

  
  ### CV empirical estimates ###
  set.seed(1000*s+mn*S)
  fold.index <- createFolds(ron1$y, k = fold.num.cv, list = TRUE, returnTrain = FALSE)
  
  

  
  ######################
  # Guo
  ######################
  # r <- 1/3  # tuning parameter
  bootstrap.r <- function(data, i){
    d <- data[i,]
    un.test.b <- signif(as.matrix(d[,-c(1:2)]) %*% as.matrix(fit.beta.full), digits = 8) 
    beta.b <- beta.b.trt <- rep(NA, length(u.data.full$u.cutoff.full))
    for (s.b in 1:length(beta.b)) {
      if (u.data.full$u.sign.full[s.b] == "GE") {
        ind.sub.b<-which(ifelse( un.test.b >= as.numeric(u.data.full$u.cutoff.full[s.b]) ,1,0) == 1) 
      } else {
        ind.sub.b<-which(ifelse( un.test.b < as.numeric(u.data.full$u.cutoff.full[s.b]) ,1,0) == 1)
      }
      ind.sub.t.b <- ind.sub.b[d$V2[ind.sub.b] == 1]
      ind.sub.c.b <- ind.sub.b[d$V2[ind.sub.b] == 0]
      trt.eff.b <- mean(d$V1[ind.sub.t.b]) - mean(d$V1[ind.sub.c.b])
      prev.b <- (length(ind.sub.b))/n1

        beta.b[s.b] <- trt.eff.b*(prev.b)^gamma

        beta.b.trt[s.b] <- trt.eff.b

    }
    out <- max(beta.b+d.c.r, na.rm = TRUE)
    out.trt <- max(beta.b.trt+d.c.r.trt, na.rm = TRUE)
    
    out.full <- c(out,out.trt)
    
    return(out.full)
  }
  
  # tuning parameter
  r.seq <- seq(0.05,0.45,0.05)
  h.i.j.r.min <- rep(NA, length(r.seq))
  gamma.max <- max(as.numeric(u.data.full$u.utility))
  trt.max <- max(as.numeric(u.data.full$u.trt.diff))
  
  for (m in 1:length(r.seq)) {
    
    # cat("r:",m,"\n")
    
    h.i.j <- rep(NA, fold.num.cv*length(u.data.full$u.cutoff.full))
    
    for (j in 1:fold.num.cv) {
      
      test_ind <- fold.index[[j]]
      train_ind <- (1:n1)[-test_ind]
      
      # bias-reduced estimator #
      d.c.r <- (1-n1^(r.seq[m] - 0.5))*(gamma.max - as.numeric(u.data.full$u.utility))
      d.c.r.trt <- (1-n1^(r.seq[m] - 0.5))*(trt.max - as.numeric(u.data.full$u.trt.diff))
      train.data <- full.b[train_ind,]
      bs <- boot(train.data, bootstrap.r, B)
      
      bs.gamma <- bs$t[,1]
      bs.trt <- bs$t[,2]
      
      gamma.max.modified <- bs.gamma
      gamma.diff <- gamma.max.modified - gamma.max
      gamma.reduced <- gamma.max - mean(gamma.diff)
      
      trt.max.modified <- bs.trt
      trt.diff <- trt.max.modified - trt.max
      trt.reduced <- trt.max - mean(trt.diff)
      
      
      # calculations on the testing data #
      test.data <- full.b[test_ind,]
      fold.un.test.b <- signif(as.matrix(test.data[,-c(1:2)]) %*% as.matrix(fit.beta.full), digits = 8)
      
      h.i <- rep(NA, length(u.data.full$u.cutoff.full))
      
      # use the reference data to estimate the effect size of the ith subgroup 
      for (i in 1:length(h.i)) {
        
        if (u.data.full$u.sign.full[i] == "GE") {
          fold.ind.sub.b<-which(ifelse( fold.un.test.b >= as.numeric(u.data.full$u.cutoff.full[i]) ,1,0) == 1)
        } else {
          fold.ind.sub.b<-which(ifelse( fold.un.test.b < as.numeric(u.data.full$u.cutoff.full[i]) ,1,0) == 1)
        }
        fold.ind.sub.t.b <- fold.ind.sub.b[test.data$V2[fold.ind.sub.b] == 1]
        fold.ind.sub.c.b <- fold.ind.sub.b[test.data$V2[fold.ind.sub.b] == 0]
        fold.trt.eff.b <- mean(test.data$V1[fold.ind.sub.t.b]) - mean(test.data$V1[fold.ind.sub.c.b])
        # standard error based on large samples 
        fold.prev.b <- (length(fold.ind.sub.b))/(dim(test.data)[1])
        
        # fold.sd.b is to account for the variation in fold.beta.b
        
        fold.beta.b <- fold.trt.eff.b
        fold.sd.b <- sqrt(var(test.data$V1[fold.ind.sub.t.b]) + var(test.data$V1[fold.ind.sub.c.b]))
        
        # evaluation of accuracy #
        h.i[i] <- (gamma.reduced - fold.beta.b)^2 
      }
      
      h.i.j[((j-1)*length(h.i) + 1):(j*length(h.i))] <- h.i
      
    }
    
    h.i.j.full <- cbind(h.i.j,c(1:length(h.i)))
    h.i.j.sum <- aggregate(h.i.j.full[,1], list(h.i.j.full[,2]), FUN=sum)
    
    h.i.j.r.min[m] <- min(h.i.j.sum[,2], na.rm = TRUE)
  }
  
  # tuning parameter #
  r.val[s] <- r.seq[which.min(h.i.j.r.min)]
  
  r <- r.val[s]
  
  
  
  
  # empirical #
  gamma.max.emp <-  max(as.numeric(u.data.full$u.utility))
  d.c.emp <- (1-n1^(r-0.5))*(gamma.max.emp - as.numeric(u.data.full$u.utility))
  
  trt.max.emp <-  max(as.numeric(u.data.full$u.trt.diff))
  d.c.emp.trt <- (1-n1^(r-0.5))*(trt.max.emp - as.numeric(u.data.full$u.trt.diff))
  
  # regression  ##
  gamma.g.sub <- trt.g.sub <- rep(NA,length(u.data.full$u.utility))
  for (u in 1:length(u.data.full$u.utility)) {
    # gamma.g.sub[[u]] <-  signif(as.matrix(full.b[(cv.full[[3]][[u]]),-c(1:2)]) %*% as.matrix(fit.beta.full), digits = 8)
    gamma.g.sub[u] <-  (mean(signif(as.matrix(full.b[(cv.full[[3]][[u]]),-c(1:2)]) %*% 
                                      as.matrix(fit.beta.full), digits = 8)))*(as.numeric(u.data.full$u.prev[u]))^gamma
    trt.g.sub[u] <- (mean(signif(as.matrix(full.b[(cv.full[[3]][[u]]),-c(1:2)]) %*% 
                                   as.matrix(fit.beta.full), digits = 8)))
  }
  gamma.max.reg <-  max(gamma.g.sub) # max regression estimated trt effect, then resulting utility
  d.c.reg <- (1-n1^(r-0.5))*(gamma.max.reg - as.numeric(u.data.full$u.utility))
  
  trt.max.reg <-  max(trt.g.sub)
  d.c.reg.trt <- (1-n1^(r-0.5))*(trt.max.reg - as.numeric(u.data.full$u.trt.diff))
  
  

  
  # augmented #
  gamma.g.sub <- trt.g.sub <- rep(NA,length(u.data.full$u.utility))
  d <- full.b.aug 
  
  gamma.g.sub.c <- gamma.g.sub.t <- vector(mode = "list", length = length(u.data.full$u.utility))
  gamma.g.sub.aug.c <-  gamma.g.sub.aug.t <-  vector(mode = "list", length = length(u.data.full$u.utility))
  
  gamma.g.sub.c.full <- signif(as.matrix(d[,(dim(ron1)[2]+2):(dim(ron1)[2]+1+dim(ron1)[2])]) %*% as.matrix(fit.beta.full), digits = 8)
  gamma.g.sub.t.full <- gamma.g.sub.c.full + signif(as.matrix(d[,2:(dim(ron1)[2]+1)]) %*% as.matrix(fit.beta.full), digits = 8)
  
  for (u in 1:length(u.data.full$u.utility)) {
    ## use bootstrap re-fitted parameters to calculate regression estimators using main data
    gamma.g.sub.t[[u]] <- gamma.g.sub.t.full
    gamma.g.sub.c[[u]] <- gamma.g.sub.c.full  
    # augmented estimator
    ind.sub.full <- (cv.full[[3]][[u]])
    ind.sub.t.full <- ind.sub.full[d$t[ind.sub.full] == 1]
    ind.sub.c.full <- ind.sub.full[d$t[ind.sub.full] == 0] 
    gamma.g.sub.aug.t[[u]] <- mean(d[ind.sub.t.full,1]) - 
      sum((gamma.g.sub.t[[u]]*{(d$t==1)-(length(ind.sub.t.full)/length(ind.sub.full))})[ind.sub.full]) / length(ind.sub.t.full)
    gamma.g.sub.aug.c[[u]] <- mean(d[ind.sub.c.full,1]) - 
      sum((gamma.g.sub.c[[u]]*{(d$t==0)-(length(ind.sub.c.full)/length(ind.sub.full))})[ind.sub.full]) / length(ind.sub.c.full)
  }
  gamma.g.sub <- (unlist(gamma.g.sub.aug.t) - unlist(gamma.g.sub.aug.c))*(as.numeric(u.data.full$u.prev))^gamma 
  gamma.max.aug <- gamma.max <- max(gamma.g.sub)
  d.c.aug <- (1-n1^(r-0.5))*(gamma.max - as.numeric(u.data.full$u.utility))
  
  trt.g.sub <- (unlist(gamma.g.sub.aug.t) - unlist(gamma.g.sub.aug.c))
  trt.max.aug <-  max(trt.g.sub)
  d.c.aug.trt <- (1-n1^(r-0.5))*(trt.max.aug - as.numeric(u.data.full$u.trt.diff))
  
  
  
  # bootstrap results
  
  bs <- boot(full.b.aug, bootstrap.full, B)

  
  ### Guo ###
  # empirical #
  gamma.max.modified <- bs$t[,7]
  gamma.diff <- gamma.max.modified - gamma.max.emp
  T.b <- sqrt(n1)*(gamma.diff)
  c.alpha <- quantile(T.b, 1-alpha/2)
  gamma.sd.guo[s] <- c.alpha/sqrt(n1)
  gamma.guo[s] <- gamma.max.emp -  mean(gamma.diff)
  
  
  trt.max.modified <- bs$t[,10]
  trt.diff <- trt.max.modified - trt.max.emp
  T.b <- sqrt(n1)*(trt.diff)
  c.alpha <- quantile(T.b, 1-alpha/2)
  trt.sd.guo[s] <- c.alpha/sqrt(n1)
  trt.eff.guo[s] <- trt.max.emp -  mean(trt.diff)

  
  # regression  ##
  gamma.max.modified <-  bs$t[,8]
  gamma.diff <- gamma.max.modified - gamma.max.reg
  T.b <- sqrt(n1)*(gamma.diff)
  c.alpha <- quantile(T.b, 1-alpha/2)
  gamma.sd.guo.reg[s] <- c.alpha/sqrt(n1)
  gamma.guo.reg[s] <- gamma.max.reg -  mean(gamma.diff)
  # trt.eff.guo.reg[s] <- gamma.guo.reg[s]/(prev.full[s]^gamma)
  
  trt.max.modified <- bs$t[,11]
  trt.diff <- trt.max.modified - trt.max.reg
  T.b <- sqrt(n1)*(trt.diff)
  c.alpha <- quantile(T.b, 1-alpha/2)
  trt.sd.guo.reg[s] <- c.alpha/sqrt(n1)
  trt.eff.guo.reg[s] <- trt.max.reg -  mean(trt.diff)
  
  
  
  # augmented #
  gamma.max.modified <- bs$t[,9]
  gamma.diff <- gamma.max.modified - gamma.max.aug
  T.b <- sqrt(n1)*(gamma.diff)
  c.alpha <- quantile(T.b, 1-alpha/2)
  gamma.sd.guo.aug[s] <- c.alpha/sqrt(n1)
  gamma.guo.aug[s] <- gamma.max.aug -  mean(gamma.diff)
  # trt.eff.guo.aug[s] <- gamma.guo.aug[s]/(prev.full[s]^gamma)
  
  trt.max.modified <- bs$t[,12]
  trt.diff <- trt.max.modified - trt.max.aug
  T.b <- sqrt(n1)*(trt.diff)
  c.alpha <- quantile(T.b, 1-alpha/2)
  trt.sd.guo.aug[s] <- c.alpha/sqrt(n1)
  trt.eff.guo.aug[s] <- trt.max.aug -  mean(trt.diff)
  
  
}  



library(erer)
output.full <- list(cbind(prev.full.ext,
                          trt.eff.ext,gamma.ext,
                          
                          # full sample #
                          prev.full,
                          trt.eff.full.emp, gamma.full.emp,
                          trt.eff.full.reg, gamma.full.reg,
                          trt.eff.full.aug, gamma.full.aug,
                          
                          # Zhang, placeholder #
                          prev.cv,
                          trt.eff.cv.emp, trt.cv.emp.sd, gamma.cv.emp, gamma.cv.emp.sd,
                          trt.eff.cv.reg, trt.cv.reg.sd, gamma.cv.reg, gamma.cv.reg.sd,
                          trt.eff.cv.aug, trt.cv.aug.sd, gamma.cv.aug, gamma.cv.aug.sd,
                
                          # Guo #
                          r.val,          
                          trt.eff.guo, trt.sd.guo, gamma.guo, gamma.sd.guo, 
                          trt.eff.guo.reg, trt.sd.guo.reg, gamma.guo.reg, gamma.sd.guo.reg,
                          trt.eff.guo.aug, trt.sd.guo.aug, gamma.guo.aug, gamma.sd.guo.aug                    

))
write.list(output.full,  file = paste(out.dir, out.name, "_", mn, ".csv", sep=""))








