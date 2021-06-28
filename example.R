# Set up parameters and generate data
library(lme4)
library(survival)
library(devtools)
source("simulate_data.R")
devtools::load_all()

set.seed(100)
DATA <- simulate_data(n = 100, 
                      alpha = c(-1.65, 0.15, 1.8, -0.05, 0.4, 0.15), 
                      beta = c(2, 1, -0.3, 1.5, 0.5), 
                      eta = c(1, -3, 0.9, -3.5, -1.5),
                      Asso = c(-1.5, 2), lambda = c(0, -0.75), sigma = 0.5,
                      ZX = c('month','sindoes','doesW', 'trb1'), 
                      YX = c("year",'year2','sindoes'),
                      Vtime = c(0,1,6,12,18,24,30), endTrial = 36,
                      scale = 800, shape = 15, ran_p = c(1,1),
                      SIGMA = matrix(c(1, 0.5, 0.5,1),2,2), 
                      cen_par = c(5, 1000), check = F, delim_val = 1.47,
                      regime = 1, freq = T, contaminated = T, percent = 0.05,
                      outlier_type = "e-outlier", k_sd = 5)
long.data <- DATA[[1]]
surv.data <- DATA[[2]]

######### A two-step method #########
## To obtain starting values for the proposed method
# Model 1: LME model for Y
md1 <- lmer(y ~ 1 + year + year2 + sindoes + (1|sid), data = long.data)
nvB <- data.frame(row.names(ranef(md1)$sid), ranef(md1)$sid)
names(nvB) <- c("sid", "nvb21")
nvB[,2] <- scale(nvB[,2], center=T, scale=T)
new_longdat <- merge(long.data, nvB, by='sid',all=T)
# Model 2: GLME model for C
md2 <- glm(c ~ 1+year+year2+sindoes+ nvb21, family=binomial, data=new_longdat)
# Model 3: GLME model for Z
md3 <- glmer(z ~ 1+month+sindoes+doesW+nvb21+(month-1|sid), family="binomial", data=new_longdat)
# Model 4: Cox/Weibull PH model for survival data
Sdata <- surv.data
Sdata$b11 <- ranef(md3)$sid[,1]
Sdata$b21 <- ranef(md1)$sid[,1]
Sdata$b11sd <- scale(Sdata$b11, center=T, scale=T)
Sdata$b21sd <- scale(Sdata$b21, center=T, scale=T)

surv.dist <- "cox"
if(surv.dist=="weibull"){
  fitCOX <- survreg(Surv(obs_time, event) ~ base+b21sd+b11sd, data = Sdata, dist = "weibull")
} else if (surv.dist=="cox"){
  fitCOX <- coxph(Surv(obs_time, event) ~ base+ b21sd+b11sd, data = Sdata)   
}

######### Proposed method: robust joint model ##########
# Defining Z-model 
glmeObject1 <- list(
  response = "z",
  reg_equation = "alpha0 + alpha1*month + alpha2*sindoes +
    alpha3*doesW+ alpha4*b21+ alpha5*month*b11",
  distribution = 'binomial',
  random_effects = c('b11', 'b21'),
  fixed_param = list(names = paste0("alpha", 0:3),
                     start_values = fixef(md3)[1:4]),
  disp_param = list(names = c('alpha4', 'alpha5'),
                    start_values = c(fixef(md3)[5], attr(summary(md3)$varcor$sid, 'stddev')),
                    lower = c(0, 0),
                    upper = c(Inf, Inf)))
# Defining Y-model 
glmeObject2 <- list(
  response = "y",
  reg_equation = "beta0 + beta1*year + beta2*year2 + beta3*sindoes + beta4*b21",
  distribution = 'normal',
  random_effects = "b21",
  fixed_param = list(names = paste0("beta", 0:4),
                     start_values = c(fixef(md1), attr(summary(md1)$varcor$sid, 'stddev'))),
  disp_param = NULL,
  left_censoring = list(indicator = "c",
                        limit_value = 1.47,
                        method = "truncated"),
  robust = list(Tpt = 2))

# Defining C-model 
glmeObject3 <- list(
  response = "c",
  reg_equation = "eta0 + eta1*year + eta2*year2 + eta3*sindoes + eta4*b21",
  distribution = 'binomial',
  random_effects = "b21",
  fixed_param = list(names = paste0("eta", 0:3),
                     start_values = coef(md2)[1:4]),
  disp_param = list(names = "eta4",
                    start_values = coef(md2)[5],
                    lower = -Inf, upper = Inf))

glmeObjects <- list(glmeObject1, glmeObject2, glmeObject3)

# Defining survival-model 
if(surv.dist == "cox"){
  distribution <- NULL
  start_values <- summary(fitCOX)$coeff[,1]
} else if(surv.dist=="weibull"){
  distribution <- "weibull"
  start_values <- - summary(fitCOX)$coeff[-1]/summary(fitCOX)$scale
}
survObject <- list(
  response = "obs_time",
  reg_equation = "lambda0*base + lambda1*b21 + lambda2*b11",
  event = "event",
  distribution = distribution,
  fixed_param = list(names = paste0("lambda", 0:2),
                     start_values = start_values),
  initial_fit = fitCOX)

# Fit a robust joint model
res <- fitHHJM(glmeObjects, survObject, long.data, surv.data,
               subject_id = "sid",
               randeff_info = list(distribution = "t-dist", degree = 3),
               mcmc_size = 50, burn_in = 100, thin = 5)
