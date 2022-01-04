library(mrgsolve)
library(magrittr)
library(dplyr)

##Code for Css simulation located in the folder of Table 1/Model files
## Input PBPK model 
humanPBPK.code    <- readRDS (file = "PBPK_Human.rds")

## Loading pbpk model 
mod.human         <- mcode ("humanpbpk", humanPBPK.code)

## Loading human MCMC data
Human.MCMC        <- readRDS(file = "MCMC_Human.rds")

## loading the theta names
theta.names       <- readRDS(file = "theta.names.rds")
which_sig         <- grep("sig", theta.names)

## Define the prediction function of Steady-state concentration (Css) for human
pred.Css <- function(pars_H) {
  
  ## Define the exposure scenario corresponding to the assumption of reverse dosimetry and IVIVE:
  ## Human were administered 1 ug/kg PFOS orally for 50 YEARs to reach steady state  
  pars_H %<>% lapply(exp)
  names(pars_H) <- names(pars_H)
  pars_H        <- pars_H [-which_sig ]
  
  ## Human body weight
  BW.H = 82.3
  
  ## Exposure scenario
  tinterval.H = 24    ## Time interval
  TDoses.H = 365*50   ## Dose times 50 years
  PDOSEoral.H = 1     ## mg/kg/day; The assumed exposure dose
  DOSEoral.H = PDOSEoral.H* BW.H  ## mg; amount of oral dose
  ex.H <- ev(amt = DOSEoral.H, 
             ii = tinterval.H, 
             addl = TDoses.H-1, cmt="AST", replicate = FALSE)
  
  ## Set up the exposure time
  tsamp.H = tgrid(0,tinterval.H*(TDoses.H-1)+49*365,24)  ## 50 YEARS and simulated for 49*365 hours after dosing
  
  out.H <- 
    mod.human %>% 
    param(pars_H) %>%
    Req (AUC_CL, Liver, Plasma)%>%
    update(atol = 1E-6, maxsteps = 5000) %>%
    mrgsim_d(data = ex.H, tgrid = tsamp.H) %>% mutate(Time = time/24, CA = Plasma,CL=Liver)
  
  
  return("human" = out.H)
  
}
iter <- 1000 #1000 iteration
Random <- sample(1:10000, 1000, replace  = FALSE)
Css.CA <- rep(NA, iter)
Css.CL <- rep(NA, iter)

for(i in 1:iter){

  pars_H <- Human.MCMC$pars[Random[i],]
  MC.Css <- pred.Css(pars_H)
  Css_CA <- MC.Css [MC.Css$Time == 365*50, ]$CA #Css based on plasma concentration
  Css_CL <- MC.Css [MC.Css$Time == 365*50, ]$CL #Css based on liver concentration
  Css.CA [i] = (Css_CA*1000)/500.13 ## Convert to unit of invitro dose (uM); molecular weight of PFOS is 500.13 g/mol
  Css.CL [i] = (Css_CL*1000)/500.13
  cat("iteration = ", i , "\n") # caption the result of each iteration per line
  
}

Css <- cbind.data.frame(Css.CA,Css.CL)
saveRDS(Css, file = "Css simulation.RDS")

summary.serum <- data.frame(median = median(Css.CA),
                            lower = quantile(Css.CA,0.025),
                            upper = quantile(Css.CA,0.975))
summary.liver <- data.frame(median = median(Css.CL),
                            lower = quantile(Css.CL,0.025),
                            upper = quantile(Css.CL,0.975))
