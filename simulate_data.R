library(MASS)
library(survsim)
library(expm)
library(truncnorm)

simulate_data <- function(n=200, 
                     ZX=c('year','sindoes','doesW', 'trb1'), 
                     YX=c("year",'year2','sindoes'),
                     Vtime=c(0,1,6,12,18,24,30),
                     endTrial=36,
                     alpha, 
                     beta, 
                     eta,
                     Asso ,
                     lambda, 
                     scale=800, shape=15, sigma=0.5, 
                     ran_p=c(1,1),
                     SIGMA=diag(c(1, 1)) , 
                     cen_par=c(5,1000),
                     check=F,
                     delim_val=1.47,
                     regime=1,
                     freq=T, rate=0.2,
                     contaminated=F,
                     percent=0,
                     outlier_type=c("e-outlier","b-outlier","both"),
                     fixbase = NULL,
                     k_sd=3){  # 1.47 in real data
  
  # BEGIN OF THE FUNCTION
  t <- 1:(endTrial*30)
  maxT <- max(t)        
  nVacc=length(Vtime)  # number of vaccination
  pVacc= diff(c(Vtime,endTrial), lag=1)*30   
  
  ni <- length(t)
  sid <- rep(1:n, each=ni)
  if(is.null(fixbase)){
    smpbase = rnorm(n, 0, 1)
    base <- rep(scale(smpbase), each=ni)
  } else{
    base <- rep(fixbase, each=ni)
  }
  
  perdi <- rep(c(pVacc), c(pVacc))
  perd <- rep(perdi, n)  
  time <- rep(t, n)
  doesti <- c(unlist(lapply(pVacc, function(x){1:x})))
  doest <- rep(doesti,n)
  injectionNOi <- c(rep(1:nVacc, pVacc))
  injectionNO <- rep(injectionNOi, n)
  
  dat <- data.frame(sid, time, injectionNO, doest, perd, base)    
  dat <- dat[order(dat$sid, dat$time), ]
  dat$sindoes <- sin(dat$doest*pi/dat$perd)
  #   dat$sinsin <- sin(dat$doest*pi/dat$perd+sin(dat$doest*pi/dat$perd))
  #   dat$sinsin10 <- dat$sinsin*10
  dat$biweek <- (dat$time-1)/15
  dat$doesM <- dat$doest/30
  dat$doesW <- dat$doest/7
  dat$year <- dat$time/365
  dat$year2 <- dat$year^2
  dat$month <- dat$time/30
  dat$month2 <- dat$month^2
  N <- nrow(dat)
  
  
  ran_effects <- mvrnorm(n, mu=rep(0, sum(ran_p)), Sigma=SIGMA)
  
  
  if(contaminated==T){
    if(outlier_type=="b-outlier"){
      outliers = sample(1:n, round(n*percent))
      topup <- sample(c(-k_sd, k_sd)*sqrt(SIGMA[1,1]), length(outliers), 
                      replace=T)
      ran_effects[outliers, 1] <- ran_effects[outliers, 1]+topup
      # print(cov(ran_effects))
      # hist(ran_effects[,1])
    } else if (outlier_type=="both"){
      outliers = sample(1:n, round(n*percent/2))
      topup <- sample(c(-k_sd, k_sd)*sqrt(SIGMA[1,1]), length(outliers), 
                      replace=T)
      ran_effects[outliers, 1] <- ran_effects[outliers, 1]+topup
      # print(cov(ran_effects))
    }
  }
  
  ran_effects <- as.matrix(ran_effects)%*%diag(1/apply(ran_effects,2,sd))
  #   print(cov(ran_effects))
  #   hist(ran_effects[,1])
  
  # in LME
  dat$trb1 <- rep(ran_effects[,1:ran_p[1]],each=ni) 
  # in GLME
  dat$trb2 <- rep(ran_effects[,(ran_p[1]+1):(sum(ran_p[1:2]))],each=ni)
  
  # (1) generate binomial data
  X <- as.matrix(cbind(1, dat[, ZX], dat$trb2*dat$month))
  logit.z <- X%*%alpha 
  dat$p.z <- exp(logit.z)/(1+exp(logit.z))
  dat$z <- rbinom(N, 1, prob=dat$p.z)  
  
  
  # (2) generate longitudinal data with censoring data
  ## generate censoring indicator
  U <- as.matrix(cbind(1, dat[,YX], dat$trb1))
  #   logit.c <- U%*%eta
  #   dat$p.c <- exp(logit.c)/(1+exp(logit.c))
  #   dat$c <- rbinom(N, 1, prob=dat$p.c)
  #   dat$y <- rtruncnorm(N, a=delim_val, b=Inf, mean=U%*%beta, sd=sigma)
  #   dat$y[dat$c==1] = delim_val
  
  dat$y <- rnorm(N, mean=U%*%beta, sd=sigma)
  dat$c = 0
  dat$c[dat$y<= delim_val] = 1
  dat$y[dat$y<=delim_val] = delim_val
  
  if(contaminated==T){
    if(outlier_type=="e-outlier"){
      outliers = sample(1:N, round(N*percent))
      # print(sd(dat$y))
      topup <- sample(c(k_sd, -k_sd)*sigma, length(outliers), replace=T)
      dat$y[outliers] <- sapply(dat$y[outliers]+topup, function(x) return(max(x, delim_val+0.01)))
    } else if(outlier_type=="both"){
      outliers = sample(1:N, round(N*percent/2))
      topup <- sample(c(k_sd,-k_sd)*sigma, length(outliers), replace=T)
      dat$y[outliers] <- sapply(dat$y[outliers]+topup, function(x) return(max(x, delim_val+0.01)))
    } 
    #     dat$y[dat$y<=delim_val] =delim_val
    #     dat$c[dat$y<=delim_val] = 1
  }
  
  
  
  # generate survival time
  #  scale= 40; shape=20
  Xs <- cbind(1, dat[!duplicated(dat$sid,fromLast=TRUE), c('base')])
  #ran_effects <- as.matrix(ran_effects)%*%1/diag(apply(as.matrix(ran_effects),2,sd))
  sexp <- exp(Xs%*%lambda+ran_effects%*%Asso)
  scale1 <- scale*sexp^{-1/shape}
  S <- rweibull(n, shape=shape, scale = scale1) 
  #  hist(S)
  # generate censoring time
  #  cen_par=c(20, 50)
  cen <- rweibull(n, shape=cen_par[1], scale = cen_par[2])
  #  hist(cen)
  # observed event time
  Scen <- pmin(S, cen)
  Scen[Scen>maxT] <- maxT
  # event indicator  
  di <- as.numeric(S <= cen)*as.numeric(S<=maxT)
  
  surv.dat <- data.frame(sid=c(1:n), 
                         event_time=Scen, 
                         event=di)
  surv.dat$obs_time <- ceiling(surv.dat$event_time)
  # mean(surv.dat$event)
  
  
  ## keep the longitudinal data observed before event time
  keep <- unlist(lapply(surv.dat$obs_time, function(x){
    # rep(c(1,0), c(x*15+1, ni-(x*15+1)))
    rep(c(1,0), c(x, ni-x))
  }))
  dat <- dat[which(keep==1),]
  
  
  ids <- unique(dat$sid)
  subsam <- c()
  
  for(i in 1:n){
    #i=1
    subdat <- subset(dat, sid==ids[i])
    maxt <- max(subdat$time)
    if(freq==T){
      rdays <- seq(16, maxt, by=15)
    } else {
      rdays <- c(16, 48, 198,  378, 558, 738, 918, 1080)
    }
    
    rdays <- rdays + sample(-2:2, length(rdays), replace = T)
    rdays <- c(rdays, 31, 181, 361, 541, 721, 901)
    keept <- c(1, rdays, maxt)
    sub <- subset(subdat, time%in%keept)
    subsam <- rbind(subsam, sub)
  }
  finalDat <- subsam
  
  finalSdat <- cbind(surv.dat, base=dat[!duplicated(dat$sid,fromLast=TRUE), 
                                        c('base')])
  
  dat=finalDat
  if(check==T){
    hist(dat$p.z)       
    cat("mean of z:",mean(dat$z), '\n')
    hist(dat$y)
    cat("mean of c", mean(dat$c),'\n')
    
    #    par(mfrow=c(1,2))
    hist(S, main="true survival time")
    hist(cen, main='censoring time')
    # hist(surv.dat$obs_time, main="observed event time")
    hist(surv.dat$obs_time[surv.dat$event==1], main="observed event time of interest")
    hist(surv.dat$obs_time[surv.dat$event==0], main="observed censoring time ")
    cat("event rate of S:", mean(di), '\n')
  }
  
  
  return(list(dat=finalDat, surv.dat=finalSdat, ran_effects))
  
  # END OF THE FUNCTION
}



