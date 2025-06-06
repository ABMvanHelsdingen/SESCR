# Fit SCR and SESCR to pregenerated data, using TMB
# Last Updated: 28 May 2025


S <- 1 # Set number for the simulations
source1 <- "PR" # name of folder for storing results (OU,PR or SCR)
source2 <- "SESCR" # name of folder where simulations stored (OUSCR, SESCR or SCR)
mask_size <- 4000; runs <- 20

## TMB template
if (!"SCLFunction" %in% getLoadedDLLs()){
  TMB::compile("SCLFunction.cpp")
  dyn.load(TMB::dynlib("SCLFunction"))
}

## Load necessary functions
source("FitFunctions.R")
library(secr)
library(stelfi)
library(TMB)

## Start script
set.seed(S)
file_path = paste("SimStudy-",source1,"-SCL/",sep="")
dir.create(file_path)


pars <- read.csv(paste(source2,"Sims/Pars_",S,".csv",sep=""))
n_sims <- nrow(pars)
camera_locations = as.matrix(read.csv(paste(source2,"Sims/Cameras_",S,".csv",sep="")))[,2:3]


# Output Data frame
output <- as.data.frame(matrix(0, nrow = n_sims, ncol = 23))
names(output) <- c("n_obs", "m", "C_obs",
                   "N_scr", "N_scr_se", "sigma_scr", "sigma_scr_se", "lambda0_scr",
                   "N", "N_se", "sigma", "sigma_se", "lambda0", "beta", "d",
                   "scr_ran", "N_better", "sigma_better", "mask_size", "runs","runtime",
                   "NLL_SESCR", "NLL_SCR")


trap_file = paste("traps_Sims_",S,".csv",sep="")
capt_file = paste("capthist_Sims_",S,".csv",sep="")


for(i in 1:n_sims){
  # Read in raw data
  identifier = paste(S,"_",i,sep="")
  bounds = c(pars$bounds1[i], pars$bounds2[i], pars$bounds3[i], pars$bounds4[i])
  obs = read.csv(paste(source2,"Sims/Sim_Data_",identifier,".csv",sep=""))
  
  # Capture histories
  times <- obs$times
  ids <- obs$ids
  cameras <- obs$cameras
  
  
  m <- max(ids)
  n <- length(times)
  if(n == 0){
    next
  }
  
  # Summaries of the Data
  output$n_obs[i] <- n
  output$m[i] <- m
  output$C_obs[i] <- length(unique(cameras))
  
  # Run SCR, using secr package
  out_scr <- tryCatch(
    {
      write.table(1000*camera_locations, trap_file, sep=",", col.names=FALSE)
      session <- rep("a", n)
      capthist <- data.frame(session, ids, 1, cameras)
      write.table(capthist, capt_file, sep=",", row.names = FALSE, col.names=FALSE)
      secrdata <- secr::read.capthist(capt_file, trap_file, detector = "count", noccasions = 1)
      
      xs <- seq(1000*bounds[1], 1000*bounds[2], length.out=66)
      ys <- seq(1000*bounds[3], 1000*bounds[4], length.out=66)
      coord <- data.frame(x = rep(xs,each=66), y = rep(ys, times=66))
      mask <- secr::read.mask(data=coord)
      
      fit <- secr::secr.fit(secrdata, mask = mask, trace = FALSE, detectfn = "HHN")
      out <- secr::region.N(fit)
      
      # Output from secr
      output$N_scr[i] <- out[2,1] # realised, not expected population
      output$N_scr_se[i] <- out[2,2] # se of population
      out1 <- predict(fit)
      output$sigma_scr[i] <- (out1[3,2] / 1000) # sigma
      output$sigma_scr_se[i] <- (out1[3,3] / 1000)
      output$lambda0_scr[i] <- out1[2,2] / pars$t[i]
      if (!is.na(output$N_scr[i])){
        output$scr_ran[i] <- 1 # success flag
      }
    }, error = function(cond){
      print(i)
      print(cond)
    }
  )

  # FIT IN TMB
  # Pre-Processing function
  out <- prepare_SESCR_NIMBLE(camera_locs = camera_locations, times = times, cameras = cameras,
                              ids = ids, bounds = bounds, nrand=mask_size)
  
  
  # Calculate the last camera where the individual currently detected was seen at
  events <- matrix(0, nrow = n + 1, ncol = 1)
  last_cameras <- numeric(m)
  for(j in 2:n){
    last_cameras[ids[(j-1)]] = cameras[(j-1)]
    events[j, 1] <- last_cameras[ids[j]]
  }
  
  # TMB/C++ indexes from 0, therefore 1 is subtracted from several vectors
  data <- list(J = out$J, m = out$m, K = n, area = out$area,
               camera_locations = camera_locations, times = c(times, pars$t[i]),
               cameras = c(cameras, 0) - 1, animals = c(ids, 0) - 1, last_cameras = events[,1] - 1,
               nrand = out$nrand, mask = out$rpts)
  
  NLL <- Inf
  t0 <- Sys.time()
  for(run in 1:runs){
    
    out_tmb <- tryCatch(
      {
        # random parameter starting points
        param <- list(log_N = log(runif(1,0.01,2)*m), log_sigma = log(runif(1,0.5,2) * output$sigma_scr[i]),
                      logit_d = runif(1,-4,1), log_beta = runif(1,-1,1), log_lambda0 = runif(1,-10,-4))
        
        obj <- TMB::MakeADFun(data = data, parameters = param, 
                              DLL = "SCLFunction",  silent = TRUE)
        opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                             control = list(trace = FALSE))
        
        
        if (opt$objective < NLL){
          coefsr = stelfi::get_coefs(obj)
          # Occasionally the MLEs diverge to nonsensical results
          if ((coefsr[1,1] < 1000) & (coefsr[4,1] < 1)){
            coefs <- coefsr
            NLL = opt$objective
          }
        }
      }, error = function(cond){
        print(cond)
        print(run)
      }
    )
  }
  t1 <- Sys.time()
  
  # Save NLL of both SESCR and SCR
  output$NLL_SESCR[i] <- NLL
  output$NLL_SCR[i] <- obj$fn(c(log(output[i,4]),log(output[i,6]),0,10,
                                log(output[i,8])))
  
  # Output from TMB
  out_results <- tryCatch(
    {
      output$N[i] <- coefs[1,1]
      output$N_se[i] <- coefs[1,2]
      output$sigma[i] <- coefs[4,1]
      output$sigma_se[i] <- coefs[4,2]
      output$lambda0[i] <- coefs[3,1]
      output$beta[i] <- coefs[2,1]
      output$d[i] <- coefs[5,1]
    }, error = function(cond){
      print("No successful run")
    }
  )
  
  output$runtime[i] <- as.numeric(difftime(t1,t0,units="mins"))
  output$runs[i] <- runs
  output$mask_size[i] <- mask_size
  
  
  # binary indicators comparing SCR and SESCR
  if (!is.na(output$N[i]) && !is.na(output$N_scr[i])){
    if (abs(output$N[i] - pars$N[i]) < abs(output$N_scr[i] - pars$N[i])){
      output$N_better[i] <- 1
    }
  }
  
  if (!is.na(output$sigma[i]) && !is.na(output$sigma_scr[i])){
    if (abs(output$sigma[i] - pars$sigma[i]) < abs(output$sigma_scr[i] - pars$sigma[i])){
      output$sigma_better[i] <- 1
    }
  }
  
  # Save after each iteration
  write.csv(output, paste(file_path, "Results_",S,".csv",sep=""))
  print(i)
}
