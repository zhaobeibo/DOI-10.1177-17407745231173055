library(MASS)

sce <- 2 # subgroup scenarios
name <- c("datas2")
out.dir <- out.dir.ext <- data.dir.ext <- paste("C:\\Users\\zhaob\\OneDrive - University of North Carolina at Chapel Hill\\School\\Research\\thesis\\Paper 2\\",name,"\\",sep="")

n1=400 # sample size
S <- 1000 # datasets generated



corr.type <- 0 # multiple outcome (if used) correlation type
type <- 0 # multiple outcome (if used) effect size type



### Change point Model: beta0*rn.intercept + t*theta+ t*delta*I(X1>c01 & X2>c02 & X3 + X4 > c03 + c04)
# normal
beta0=0.1 # flat control
theta=0.0 # flat trt
delta=0.3*2 # treatment effect for the best subgroup

# 3 correlated outcomes: y1t, y1c ; y2t, y2c ; y3t, y3c
# correlation types: 0 = independent, 1 = high-correlation
if (corr.type == 0) {
  rho1=0.0  #correlation between y1t and y1c ;  y2t and y2c ;  y3t and y3c 
  rho2=0.0  #correlation between y1t and y2t ;  y1t and y3t ;  y2t and y3t ; y1c and y2c ;  y1c and y3c ;  y2c and y3c 
  rho3 = rho1*rho2 #correlation between y1t and y2c ;  y1t and y3c ;  y2t and y1c ; y2t and y3c ;  y3t and y1c ;  y3t and y2c 
  sig2=sig2t=sig2c=1 # common, known variance  
} else if (corr.type == 1) {
  rho1=0.0  #correlation between y1t and y1c ;  y2t and y2c ;  y3t and y3c 
  rho2=0.8  #correlation between y1t and y2t ;  y1t and y3t ;  y2t and y3t ; y1c and y2c ;  y1c and y3c ;  y2c and y3c 
  rho3 = rho1*rho2 #correlation between y1t and y2c ;  y1t and y3c ;  y2t and y1c ; y2t and y3c ;  y3t and y1c ;  y3t and y2c 
  sig2=sig2t=sig2c=1 # common, known variance 
}
# biomarker randomness
bio.sigma = 0.0
mu0 = 0.0








for(s in 1:S)
{
  ### Generate samples ###
  
  set.seed(s)
  cat("s:",s,"\n")
  
  datnum <- c(s)
  
  # 2Kx2K variance-covariance matrix 
  sig = matrix(c(sig2t,rho1,rho2,rho3,rho2,rho3,
                 rho1,sig2c,rho3,rho2,rho3,rho2,
                 rho2,rho3,sig2t,rho1,rho2,rho3,
                 rho3,rho2,rho1,sig2c,rho3,rho2,
                 rho2,rho3,rho2,rho3,sig2t,rho1,
                 rho3,rho2,rho3,rho2,rho1,sig2c),6,6)
  # treatment indicator matrix
  t <- c(rep(1,n1/2),rep(0,n1/2)) 
  # initialize 4 biomarkers
  biom1.x0.t <- runif(n1/2,0,1)
  biom1.x0.c <- runif(n1/2,0,1)  
  biom1.x1.t <- runif(n1/2,0,1)
  biom1.x1.c <- runif(n1/2,0,1)
  biom1.x2.t <- runif(n1/2,0,1)
  biom1.x2.c <- runif(n1/2,0,1)
  biom1.x3.t <-runif(n1/2,0,1) 
  biom1.x4.t <-runif(n1/2,0,1) 
  biom1.x3.c <-runif(n1/2,0,1) 
  biom1.x4.c <-runif(n1/2,0,1) 
  
  biom1.x0 <- c(biom1.x0.t,biom1.x0.c)
  biom1.x1 <- c(biom1.x1.t,biom1.x1.c)  
  biom1.x2 <- c(biom1.x2.t,biom1.x2.c)  
  biom1.x3 <- c(biom1.x3.t,biom1.x3.c)
  biom1.x4 <- c(biom1.x4.t,biom1.x4.c)
  
  biom1.x5.t <- rbinom(n1/2, size=1, prob=0.5)
  biom1.x5.c <- rbinom(n1/2, size=1, prob=0.5)
  biom1.x6.t <- rbinom(n1/2, size=1, prob=0.5)
  biom1.x6.c <- rbinom(n1/2, size=1, prob=0.5)  
  biom1.x7.t <- rbinom(n1/2, size=1, prob=0.5)
  biom1.x7.c <- rbinom(n1/2, size=1, prob=0.5)
  biom1.x8.t <- rbinom(n1/2, size=1, prob=0.5)
  biom1.x8.c <- rbinom(n1/2, size=1, prob=0.5)
  biom1.x9.t <- rbinom(n1/2, size=1, prob=0.5)
  biom1.x9.c <- rbinom(n1/2, size=1, prob=0.5)  
  
  biom1.x5 <- c(biom1.x5.t,biom1.x5.c)
  biom1.x6 <- c(biom1.x6.t,biom1.x6.c)  
  biom1.x7 <- c(biom1.x7.t,biom1.x7.c)  
  biom1.x8 <- c(biom1.x8.t,biom1.x8.c)
  biom1.x9 <- c(biom1.x9.t,biom1.x9.c)  
  
  #Initialize 3 outcomes 
  y1.t<-y1.c<-y2.t<-y2.c<-y3.t<-y3.c<- mu.t.1<-mu.c.1<-mu.t.2<-mu.c.2<-mu.t.3<-mu.c.3<-rep(NA,n1/2)
  
  mu.t<-mu.c<-rep(NA,n1/2)  
  for(i in 1:(n1/2)){
    
    if (sce == 0) {
      ## null scenario ##
      mu.t[i]<- beta0 + theta
    } else if (sce == 1) {
      ## scenario 1 ##
      mu.t[i]<- beta0 + theta+(delta*ifelse(biom1.x0.t[i]  > 0.5, 1,0))
    } else if (sce == 2) {
      ## scenario 2 ##
      delta.1 <- 0.6
      mu.t[i]<- beta0 + theta+(delta.1*ifelse(biom1.x0.t[i] + biom1.x1.t[i]  > 1, 1,0))
    } else if (sce == 3) {
      # ## scenario 3 ##
      delta.1 <- 0.3
      delta.2 <- 0.3
      mu.t[i]<- beta0 + theta+(delta.1*ifelse(biom1.x0.t[i]  > 0.5 & biom1.x1.t[i] > 0.5, 1,0))+
        (delta.2*ifelse(biom1.x0.t[i] + biom1.x1.t[i] > 1, 1,0))
      
    } else if (sce == 4) {
      ## scenario 4 ##
      delta.1 <- 0.6
      mu.t[i]<- beta0 + theta+(delta.1*ifelse(biom1.x0.t[i] + biom1.x1.t[i]  > 1 & biom1.x5.t[i] +  biom1.x6.t[i] >= 1, 1,0))
    } else if (sce == 10) {
      delta.1 <- 0.6
      mu.t[i]<- beta0 + theta+(delta.1*ifelse(biom1.x0.t[i] + biom1.x1.t[i] + biom1.x2.t[i] + biom1.x3.t[i] + biom1.x4.t[i] > 2.5 & 
                                                biom1.x5.t[i] + biom1.x6.t[i] + biom1.x7.t[i] + biom1.x8.t[i] + biom1.x9.t[i] >= 2, 1,0))
    }
    
    
    mu.c[i]<- beta0 # set means in control group to 0
    
    mu =c(mu.t[i],mu.c[i],mu.t[i],mu.c[i],mu.t[i],mu.c[i])
    if (type == 1) {
      mu =c(mu.t[i],mu.c[i],mu.c[i],mu.c[i],mu.c[i],mu.c[i])
    }
    resp=mvrnorm(1,mu,sig)
    y1.t[i]<-resp[1]
    y1.c[i]<-resp[2]
    y2.t[i]<-resp[3]
    y2.c[i]<-resp[4]
    y3.t[i]<-resp[5]
    y3.c[i]<-resp[6]
  }
  

  dat <- data.frame(y=c(y1.t, y1.c), t, biom1.x1,biom1.x2,biom1.x3,biom1.x4,biom1.x5,biom1.x6,biom1.x7,biom1.x8,biom1.x9,biom1.x0)
  
  saveRDS(dat,  file = paste(out.dir, "data_", datnum, ".Rds", sep=""))
  
  
  ### Generate external validating samples for reference estimates ###
  n.obs = 10000
  # treatment indicator matrix
  t.ext <- c(rep(1,n.obs),rep(0,n.obs))
  # init.extialize 4 biomarkers
  biom1.x0.t.ext <- runif(n.obs,0,1)
  biom1.x0.c.ext <- runif(n.obs,0,1)
  biom1.x1.t.ext <- runif(n.obs,0,1)
  biom1.x1.c.ext <- runif(n.obs,0,1)
  biom1.x2.t.ext <- runif(n.obs,0,1)
  biom1.x2.c.ext <- runif(n.obs,0,1)
  biom1.x3.t.ext <-runif(n.obs,0,1)
  biom1.x4.t.ext <-runif(n.obs,0,1)
  biom1.x3.c.ext <-runif(n.obs,0,1)
  biom1.x4.c.ext <-runif(n.obs,0,1)

  biom1.x0.ext <- c(biom1.x0.t.ext,biom1.x0.c.ext)
  biom1.x1.ext <- c(biom1.x1.t.ext,biom1.x1.c.ext)
  biom1.x2.ext <- c(biom1.x2.t.ext,biom1.x2.c.ext)
  biom1.x3.ext <- c(biom1.x3.t.ext,biom1.x3.c.ext)
  biom1.x4.ext <- c(biom1.x4.t.ext,biom1.x4.c.ext)


  biom1.x5.t.ext <- rbinom(n.obs, size=1, prob=0.5)
  biom1.x5.c.ext <- rbinom(n.obs, size=1, prob=0.5)
  biom1.x6.t.ext <- rbinom(n.obs, size=1, prob=0.5)
  biom1.x6.c.ext <- rbinom(n.obs, size=1, prob=0.5)
  biom1.x7.t.ext <- rbinom(n.obs, size=1, prob=0.5)
  biom1.x7.c.ext <- rbinom(n.obs, size=1, prob=0.5)
  biom1.x8.t.ext <- rbinom(n.obs, size=1, prob=0.5)
  biom1.x8.c.ext <- rbinom(n.obs, size=1, prob=0.5)
  biom1.x9.t.ext <- rbinom(n.obs, size=1, prob=0.5)
  biom1.x9.c.ext <- rbinom(n.obs, size=1, prob=0.5)

  biom1.x5.ext <- c(biom1.x5.t.ext,biom1.x5.c.ext)
  biom1.x6.ext <- c(biom1.x6.t.ext,biom1.x6.c.ext)
  biom1.x7.ext <- c(biom1.x7.t.ext,biom1.x7.c.ext)
  biom1.x8.ext <- c(biom1.x8.t.ext,biom1.x8.c.ext)
  biom1.x9.ext<- c(biom1.x9.t.ext,biom1.x9.c.ext)

  #Initialize 3 outcomes
  y1.t<-y1.c<-y2.t<-y2.c<-y3.t<-y3.c<-rep(NA,n.obs)

  mu.t<-mu.c<-rep(NA,n.obs)
  for(i in 1:(n.obs)){

    if (sce == 0) {
      ## null scenario ##
      mu.t[i]<- beta0 + theta
    } else if (sce == 1) {
      ## scenario 1 ##
      mu.t[i]<- beta0 + theta+(delta*ifelse(biom1.x0.t.ext[i]  > 0.5, 1,0))
    } else if (sce == 2) {
      ## scenario 2 ##
      delta.1 <- 0.6
      mu.t[i]<- beta0 + theta+(delta.1*ifelse(biom1.x0.t.ext[i] + biom1.x1.t.ext[i]  > 1, 1,0))
    } else if (sce == 3) {
      # ## scenario 3 ##
      delta.1 <- 0.3
      delta.2 <- 0.3
      mu.t[i]<- beta0 + theta+(delta.1*ifelse(biom1.x0.t.ext[i]  > 0.5 & biom1.x1.t.ext[i] > 0.5, 1,0))+
        (delta.2*ifelse(biom1.x0.t.ext[i] + biom1.x1.t.ext[i] > 1, 1,0))
    } else if (sce == 4) {
      ## scenario 4 ##
      delta.1 <- 0.6
      mu.t[i]<- beta0 + theta+(delta.1*ifelse(biom1.x0.t.ext[i] + biom1.x1.t.ext[i]  > 1 & biom1.x5.t.ext[i] +  biom1.x6.t.ext[i] >= 1, 1,0))
    } else if (sce == 10) {
      delta.1 <- 0.6
      mu.t[i]<- beta0 + theta+(delta.1*ifelse(biom1.x0.t.ext[i] + biom1.x1.t.ext[i]+ biom1.x2.t.ext[i]+ biom1.x3.t.ext[i]+ biom1.x4.t.ext[i]  > 2.5 &
                                                biom1.x5.t.ext[i]+biom1.x6.t.ext[i]+biom1.x7.t.ext[i]+biom1.x8.t.ext[i]+biom1.x9.t.ext[i] >= 2, 1,0))
    }

    mu.c[i]<- beta0 # set means in control group to 0

    mu =c(mu.t[i],mu.c[i],mu.t[i],mu.c[i],mu.t[i],mu.c[i])
    if (type == 1) {
      mu =c(mu.t[i],mu.c[i],mu.c[i],mu.c[i],mu.c[i],mu.c[i])
    }
    resp=mvrnorm(1,mu,sig)
    y1.t[i]<-resp[1]
    y1.c[i]<-resp[2]
    y2.t[i]<-resp[3]
    y2.c[i]<-resp[4]
    y3.t[i]<-resp[5]
    y3.c[i]<-resp[6]
  }

  dat.ext <- data.frame(y=c(y1.t, y1.c), t.ext, biom1.x1.ext,biom1.x2.ext,biom1.x3.ext,biom1.x4.ext,
                        biom1.x5.ext,biom1.x6.ext,biom1.x7.ext,biom1.x8.ext,biom1.x9.ext,biom1.x0.ext)

  saveRDS(dat.ext,  file = paste(out.dir.ext, "data_ext_", datnum, ".Rds", sep=""))
  
}