# Fit SCR and SESCR to pregenerated data, using TMB and MLE
# Last Updated: 26 February 2025

# This block is designed to be run on the NeSI server
args <- commandArgs(trailingOnly = TRUE)
N <- as.numeric(args[1])

# To be adjusted as appropriate
n_sims <- 8
mask_size <- 4000; runs <- 20
source1 <- "PR" # name of folder for storing results (OU,PR or SCR)
source2 <- "SESCR" # name of folder where simulations stored (OUSCR, SESCR or SCR)


library(coda)

## TMB templates
if (!"MLEFunction" %in% getLoadedDLLs()){
  TMB::compile("MLEFunction.cpp")
  dyn.load(TMB::dynlib("MLEFunction"))
}

## START OF SCRIPT ##
set.seed(2+10*N)
source("FitFunctions.R")

file_path = paste("SimStudy-",source1,"-NIMBLE/",sep="")
dir.create(file_path)


pars <- read.csv(paste(source2,"Sims/Pars_",N,".csv",sep=""))
camera_locations = as.matrix(read.csv(paste(source2,"Sims/Cameras_",N,".csv",sep="")))[,2:3]


# Output Data frame
output <- as.data.frame(matrix(0, nrow = nrow(pars), ncol = 28))
names(output) <- c("n_cameras", "N_true", "t", "sigma_true", 
                   "mu0_true", "beta_true", "d_true", "n_obs", "m", "C_obs",
                   "N_scr", "N_scr_se", "sigma_scr", "sigma_scr_se", "lambda0_scr",
                   "N", "N_se", "sigma", "sigma_se", "lambda0", "beta", "d",
                   "scr_ran", "N_better", "sigma_better", "mask_size", "runs","runtime")

if (source1 == "OU"){
  names(output)[5:7] <- c("epsilon_true", "beta_true", "NULL")
} else if (source1 == "SCR"){
  names(output)[5:7] <- c("mu0_true", "NULL1", "NULL2")
}

trap_file = paste("traps_Sims_",N,".csv",sep="")
capt_file = paste("capthist_Sims_",N,".csv",sep="")


for(i in 1:nrow(pars)){
  identifier = paste(N,"_",i,sep="")
  bounds = c(pars$bounds1[i], pars$bounds2[i], pars$bounds3[i], pars$bounds4[i])
  
  
  obs = read.csv(paste(source2,"Sims/Sim_Data_",identifier,".csv",sep=""))
  pars <- read.csv(paste(source2,"Sims/Pars_",N,".csv",sep=""))
  

  times <- obs$times
  ids <- obs$ids
  cameras <- obs$cameras
  m <- max(ids)

  n <- length(times); print(i); print(n); print(m)
  if(n == 0){
    next
  }
  
  # Summary of the Data
  output$n_cameras[i] <- nrow(camera_locations)
  output$N_true[i] <- pars$N[i]
  output$t[i] <- pars$t[i]
  output$sigma_true[i] <- pars$sigma[i]
  
  if (source1 == "OU"){
    output$epsilon_true[i] <- pars$epsilon[i]
    output$beta_true[i] <- 1/pars$tau[i]
  } else if (source1 == "PR"){
    output$mu0_true[i] <- pars$mu0[i]
    output$beta_true[i] <- pars$beta[i]
    output$d_true[i] <- pars$dratio[i]
  } else if (source1 == "SCR"){
    output$mu0_true[i] <- pars$mu0[i]
  }

  
  output$n_obs[i] <- n
  output$m[i] <- m
  output$C_obs[i] <- length(unique(cameras))
  
  
  out_scr <- tryCatch(
    {
      
      write.table(1000*camera_locations, trap_file, sep=",", col.names=FALSE)
      session <- rep("ij", length(times))
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

  # Fit in TMB
  
  # Pre-Processing function
  out <- prepare_SESCR_NIMBLE(camera_locs = camera_locations, times = times, cameras = cameras,
                              ids = ids, bounds = bounds, nrand=mask_size)
  
  # Data as list
  
  events <- matrix(0, nrow = length(times) + 1, ncol = 4)
  events[, 1] <- c(times, pars$t[i])
  events[, 2] <- c(cameras, 0)
  events[, 3] <- c(ids, 0)
  last_cameras <- numeric(m)
  for(j in 2:length(times)){
    last_cameras[ids[(j-1)]] = cameras[(j-1)]
    events[j, 4] <- last_cameras[ids[j]]
    
  }
  
  data <- list(J = out$J, m = out$m, K = length(times) - 1, area = out$area,
               camera_locations = camera_locations, times = events[,1],
               cameras = events[,2] - 1, animals = events[,3] - 1, last_cameras = events[,4] - 1,
               nrand = out$nrand, mask = out$rpts)
  
  # parameter starting points
  NLL <- Inf
  t0 <- Sys.time()
  for(run in 1:runs){
    
    out_tmb <- tryCatch(
      {
        # parameter starting points
        param <- list(log_N = log(runif(1,0.01,2)*m), log_sigma = log(runif(1,0.5,2) * output$sigma_scr[i]),
                      logit_d = runif(1,-4,1), log_beta = runif(1,-1,1), log_lambda0 = runif(1,-10,-4))
        
        obj <- TMB::MakeADFun(data = data, parameters = param, 
                              DLL = "MLEFunction",  silent = TRUE)
        opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                             control = list(trace = FALSE))
        
        
        if (opt$objective < NLL){
          coefsr = stelfi::get_coefs(obj)
          if (coefsr[1,1] < 1000){
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
  output$runtime[i] <- as.numeric(difftime(t1,t0,units="mins"))
  output$runs[i] <- runs
  output$mask_size[i] <- mask_size
  
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
  
  
  # binary indicators comparing SCR and SESCR
  if (!is.na(output$N[i]) && !is.na(output$N_scr[i])){
    
    if (abs(output$N[i] - output$N_true[i]) < abs(output$N_scr[i] - output$N_true[i])){
      output$N_better[i] <- 1
    }
  }
  
  if (!is.na(output$sigma[i]) && !is.na(output$sigma_scr[i])){
    if (abs(output$sigma[i] - output$sigma_true[i]) < abs(output$sigma_scr[i] - output$sigma_true[i])){
      output$sigma_better[i] <- 1
    }
  }
  
  # OUTPUT after each iteration
  write.csv(output, paste(file_path, "Results_",N,".csv",sep=""))
}
