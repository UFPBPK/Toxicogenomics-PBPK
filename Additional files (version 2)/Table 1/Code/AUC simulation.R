library(mrgsolve)
library(magrittr)

# Locate in Table 1/Model files
micePBPK.code     <- readRDS (file = "micePBPK.RDS")
humanPBPK.code    <- readRDS (file = "humanPBPK.RDS")

mod.mouse         <- mcode ("micepbpk", micePBPK.code)
mod.human         <- mcode ("humanpbpk", humanPBPK.code)

Mouse.MCMC        <- readRDS(file = "MCMC_Mice.rds")
Human.MCMC        <- readRDS(file = "MCMC_Human.rds")

theta.names       <- readRDS(file = "theta.names.rds")
which_sig         <- grep("sig", theta.names)

pred.AUC <- function(pars.mouse,pars.human) {
  
  
  ## Define the exposure scenario for mice: oral exposure to 1 mg/kg/day 
  
  pars.mouse %<>% lapply(exp)
  names(pars.mouse) <- names(pars.mouse)
  pars.mouse        <- pars.mouse [-which_sig ]
  
  
  BW.mouse                 = 0.025                                      # mouse body weight
  tinterval.mouse          = 24                                        # Time interval
  TDoses.mouse             = 7                                         # Dose times (total number of dosing)
  PDOSEoral.mouse          = 1                                         # This value is from Table 4-8 (Page 4-14) in the 2016 EPA Report, which was based on Seacat et al. 2003
  DOSEoral.mouse           = PDOSEoral.mouse*BW.mouse                      # amount of oral dose
  ex.mouse                 <- ev (ID = 1, amt= DOSEoral.mouse, 
                                  ii = tinterval.mouse, addl=TDoses.mouse-1, cmt="AST", replicate = FALSE)
  tsamp.mouse.rep          = tgrid(0,tinterval.mouse*(TDoses.mouse-1)+24*1,12)   ## Simulation 24*7 + 24*7 hours (7 days)
  
  ## Prediction of AUC
  out.mouse <- 
    mod.mouse %>% 
    param(pars.mouse) %>%
    Req (Plasma, AUC_CA, AUC_CL)%>%
    update(atol = 1E-5, maxsteps = 50000) %>%
    mrgsim_d (data = ex.mouse, tgrid = tsamp.mouse.rep)
  
  outdf.mouse  <- cbind.data.frame(
    Time   = out.mouse$time/24,
    CA     = out.mouse$Plasma,                                 ## ug/ml,PFOS levels in plasma
    AUC.CA = out.mouse$AUC_CA,
    AUC.CL = out.mouse$AUC_CL)
  
  
  ## Human scenario 
  pars.human %<>% lapply(exp)
  names(pars.human) <- names(pars.human)
  pars.human        <- pars.human [-which_sig ]
  
  ###
  BW.human             = 82.3                                 ## kg
  
  ## 
  tinterval.human      = 24                                   ## Time interval
  TDoses.human.A       = 365*50                              ## Dose time for 50 years to reach the steady state
  PDOSEoral.human.A    = 1                                  ## mg/kg; BW Oral dose
  DOSEoral.human.A = PDOSEoral.human.A* BW.human              ## mg; amount of oral dose
  
  ex.human.A <- ev(amt= DOSEoral.human.A, ii = tinterval.human, addl = TDoses.human.A-1, cmt="AST", replicate = FALSE)
  
  ## set up the exposure time
  tsamp.human.A = tgrid(0,tinterval.human*(TDoses.human.A-1)+49*10,24)  ## 50 years and simulated for 49*365 hours (1 year) after dosing
  
  
  
  out.human.A <- 
    mod.human %>% 
    param(pars.human) %>%
    Req (Plasma, AUC_CA, AUC_CL)%>%
    update(atol = 1E-5, rtol = 1e-4, maxsteps = 50000) %>%
    mrgsim_d(data = ex.human.A, tgrid = tsamp.human.A)
  
  outdf.human.A <- cbind.data.frame(
    Time   = out.human.A$time/24,
    CA     = out.human.A$Plasma,                                 ## ug/ml,PFOS levels in plasma
    AUC.CA = out.human.A$AUC_CA,
    AUC.CL = out.human.A$AUC_CL)
  
  
  return(list("mouse"     = outdf.mouse, # save the outdf.mouse result to the object of mouse
              "human.m" = outdf.human.A))
  
}

iter <- 1000
random <-sample(0:50000,iter,replace = FALSE)
AUC <- matrix(NA, ncol = 4, nrow = iter)

for (i in 1:iter){
  
  pars.mouse    = Mouse.MCMC$pars [random[i],]
  pars.human  = Human.MCMC$pars [random[i],]
  pars.human [2] = log(exp(pars.human [2])/1000) # km of human is ug/L and transfered to unit of mg/L corresponding to the unit of animals
  pars.human [4] = log(exp(pars.human [4])/1000)
  
  MC.HED <- pred.AUC(pars.mouse,pars.human) # Generate the results of 5000 Plasma and AUC_CA. Plasma is included here in case we want to use Cmax for deriving POD
  
  AUC [i,1] <- MC.HED$mouse    [MC.HED$mouse$Time == 7, ]$AUC.CA         ## AUC of plasma in mice based on NOAEL 
  AUC [i,2] <- MC.HED$human.m[MC.HED$human.m$Time == 365*50, ]$AUC.CA     ## AUC of plasma in humans based on rat NOAEL
  AUC [i,3] <- MC.HED$mouse    [MC.HED$mouse$Time == 7, ]$AUC.CL         ## AUC of liver in mice based on NOAEL 
  AUC [i,4] <- MC.HED$human.m[MC.HED$human.m$Time == 365*50, ]$AUC.CL     ## AUC of liver in humans based on mouse NOAEL
  cat("iteration = ", i , "\n") # caption the result of each iteration per line
}

colnames(AUC) = c("AUC.CA.mouse","AUC.CA.human.m",
                  "AUC.CL.mouse","AUC.CL.human.m")

AUC.idv  = as.data.frame(AUC) #make a copy of AUC

## Estimated the average serum concentraiton (ASC) 
AUC.idv$Avg.CA.mouse    = AUC.idv$AUC.CA.mouse/(7*24)       # mouse ASC: AUC normalized by exposure duration (7 days or 7*24 hours) 
AUC.idv$Avg.CA.human.m   = AUC.idv$AUC.CA.human.m/(365*50*24)  # Human ASC: AUC normalized by the duration reaching the steady state (Based on the simulation results, it reached steady state at ~25 years.)   

## Estimated the average liver concentration (ALC)
AUC.idv$Avg.CL.mouse     = AUC.idv$AUC.CL.mouse/(7*24)          
AUC.idv$Avg.CL.human.m   = AUC.idv$AUC.CL.human.m/(365*50*24) 

saveRDS(AUC.idv,file = "AUC simulation.RDS")

AUC.summary.mouse.CA <- c(lower = quantile(AUC.idv$Avg.CA.mouse,0.025),
                          median = median(AUC.idv$Avg.CA.mouse),
                          upper = quantile(AUC.idv$Avg.CA.mouse,0.975))

AUC.summary.mouse.CL <- c(lower = quantile(AUC.idv$Avg.CL.mouse,0.025),
                          median = median(AUC.idv$Avg.CL.mouse),
                          upper = quantile(AUC.idv$Avg.CL.mouse,0.975))

AUC.summary.human.CA <- c(lower = quantile(AUC.idv$Avg.CA.human.m,0.025),
                          median = median(AUC.idv$Avg.CA.human.m),
                          upper = quantile(AUC.idv$Avg.CA.human.m,0.975))

AUC.summary.human.CL <- c(lower = quantile(AUC.idv$Avg.CL.human.m,0.025),
                          median = median(AUC.idv$Avg.CL.human.m),
                          upper = quantile(AUC.idv$Avg.CL.human.m,0.975))

