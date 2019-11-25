DATA <- simulate_data(n = 100, 
                      ZX = c('month','sindoes','doesW', 'trb1'), 
                      YX = c("year",'year2','sindoes'),
                      Vtime = c(0,1,6,12,18,24,30), endTrial = 36,
                      alpha = c(-1.65, 0.15, 1.8, -0.05, 0.4, 0.15), 
                      beta = c(2, 1, -0.3, 1.5, 0.5), 
                      eta = c(1,-3, 0.9, -3.5,-1.5),
                      Asso = c( -1.5, 2), lambda = c(0, -0.75), 
                      scale = 800, shape = 15, sigma = 0.5, ran_p = c(1,1),
                      SIGMA = matrix(c(1, 0.5, 0.5,1),2,2), 
                      cen_par = c(5, 1000), check = F, delim_val = 1.47,
                      regime = 1, freq = T, contaminated = F, percent = 0)

long.data <- DATA[[1]]
surv.data <- DATA[[2]]

md2 <- lmer(y ~ 1 + year + year2 + sindoes + (1|sid), data = long.data)
nvB <- data.frame(row.names(ranef(md2)$sid), ranef(md2)$sid)
names(nvB) <- c("sid", "nvb21")
nvB[,2] <- scale(nvB[,2], center=T, scale=T)
mydat <- merge(long.data, nvB, by='sid',all=T)

md3 <- glm(c ~ 1+year+year2+sindoes+ nvb21, family=binomial, data=mydat)
md1 <- glmer(z ~ 1+month+sindoes+doesW+nvb21+(month-1|sid), family="binomial", data=mydat)

Sdata <- surv.data
Sdata$b11 <- ranef(md1)$sid[,1]
Sdata$b21 <- ranef(md2)$sid[,1]
Sdata$b11sd <- scale(Sdata$b11, center=T, scale=T)
Sdata$b21sd <- scale(Sdata$b21, center=T, scale=T)
if(surv.dist=="weibull"){
  fitCOX <- survreg(Surv(obs_time, event) ~ base+b21sd+b11sd, data = Sdata, dist = "weibull")
} else if (surv.dist=="cox"){
  fitCOX <- coxph(Surv(obs_time, event) ~ base+ b21sd+b11sd, data = Sdata)   
}


glmeObject1 <- list(
  response = "z",
  reg_equation = "alpha0 + alpha1*month + alpha2*sindoes +
    alpha3*doesW+ alpha4*b21+ alpha5*month*b11",
  distribution = 'binomial',
  random_effects = c('b11', 'b21'),
  fixed_param = list(names = paste0("alpha", 0:3),
                     start_values = c(-1.65, 0.15, 1.8, -0.05)),
  disp_param = list(names = c('alpha4', 'alpha5'),
                    start_values = c(0.4, 0.15),
                    lower = c(0, 0),
                    upper = c(Inf, Inf)))

glmeObject2 <- list(
  response = "y",
  reg_equation = "beta0 + beta1*year + beta2*year2 + beta3*sindoes + beta4*b21",
  distribution = 'normal',
  random_effects = "b21",
  fixed_param = list(names = paste0("beta", 0:3),
                     start_values = c(2, 1, -0.3, 1.5)),
  disp_param = list(names = "beta4",
                    start_values = 0.5,
                    lower = -Inf,
                    upper = Inf),
  left_censoring = list(indicator = "c",
                        limit_value = 1.47,
                        method = "truncated"), # or "tobit"),
  robust = list(Tpt = 2))

glmeObject3 <- list(
  response = "c",
  reg_equation = "eta0 + eta1*year + eta2*year2 + eta3*sindoes + eta4*b21",
  distribution = 'binomial',
  random_effects = "b21",
  fixed_param = list(names = paste0("eta", 0:3),
                     start_values = c(1, -3, 0.9, -3.5)),
  disp_param = list(names = "eta4",
                    start_values = -1.5,
                    lower = -Inf,
                    upper = Inf))

glmeObjects <- list(glmeObject1, glmeObject2, glmeObject3)

survObject <- list(
  response = "obs_time",
  reg_equation = "lambda0*base + lambda1*b21 + lambda2*b11",
  event = "event",
  distribution = NULL,
  fixed_param = list(names = paste0("lambda", 0:2),
                     start_values = c(-0.75, -1.5, 2)),
  initial_fit = fitCOX)

res <- fitHHJM(glmeObjects, survObject, long.data, surv.data,
               subject_id = "sid",
               randeff_info = list(distribution = "t-dist", degree = 3),
               mcmc_size = 50,
               burn_in = 10, thin = 3)
