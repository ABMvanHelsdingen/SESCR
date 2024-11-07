# Parameter Recovery Simulation Study
# Last Updated: 8 November 2024

args <- commandArgs(trailingOnly = TRUE)
N <- as.numeric(args[1])
n_sims <- 8

library(coda)

## START OF SCRIPT ##
set.seed(2+10*N)
source("FitFunctions.R")


file_path = "SimStudy-SESCR-TMB/"
dir.create(file_path)


pars <- read.csv(paste("SESCRSims/Pars_",N,".csv",sep=""))
camera_locations = as.matrix(read.csv(paste("SESCRSims/Cameras_",N,".csv",sep="")))[,2:3]

## TMB templates
if (!"MLEFunction" %in% getLoadedDLLs()){
  TMB::compile("MLEFunction.cpp")
  dyn.load(TMB::dynlib("MLEFunction"))
}


# Output Data frame

output <- as.data.frame(matrix(0, nrow = nrow(pars), ncol = 30))
names(output) <- c("n_cameras", "M_true", "t", "sigma_true", "mu0_true", "beta_true", "d_true",
                   "n_obs", "m", "C_obs", "m_1obs", "m_maxobs", "m_1C", "m_maxC",
                   "M_scr", "M_scr_se", "sigma_scr", "sigma_scr_se",
                   "M", "M_se", "sigma", "sigma_se", "g0", "beta", "d",
                   "scr_ran", "M_better", "sigma_better", "run", "runtime")

trap_file = paste("traps_Sims_",N,".csv",sep="")
capt_file = paste("capthist_Sims_",N,".csv",sep="")


for(i in 1:nrow(pars)){
  identifier = paste(N,"_",i,sep="")
  bounds = c(pars$bounds1[i], pars$bounds2[i], pars$bounds3[i], pars$bounds4[i])
  
  obs = read.csv(paste("SESCRSims/Sim_Data_",identifier,".csv",sep=""))
  pars <- read.csv(paste("SESCRSims/Pars_",N,".csv",sep=""))
  
  
  times <- obs$times
  ids <- obs$ids
  cameras <- obs$cameras
  m <- max(ids)
  
  n <- length(times); print(i); print(n); print(m)
  if(n == 0){
    next
  }
  
  # Summary of the Data
  output[i, 1] <- nrow(camera_locations)
  output[i, 2] <- pars$M[i]
  output[i, 3] <- pars$t[i]
  output[i, 4] <- pars$sigma[i]
  output[i, 5] <- pars$mu0[i]
  output[i, 6] <- pars$beta[i]
  output[i, 7] <- pars$dratio[i]
  
  output[i, 8] <- n
  output[i, 9] <- m
  output[i, 10] <- length(unique(cameras))
  
  tab <- table(ids)
  output[i, 11] <- sum(tab == 1)
  output[i, 12] <- max(tab)
  
  TWTable = table(ids, cameras) > 0
  CPI = rowSums(TWTable) # cameras per individual
  
  output[i, 13] <- sum(CPI == 1)
  output[i, 14] <- max(CPI)
  
  out_scr <- tryCatch(
    {
      
      write.table(1000*camera_locations, trap_file, sep=",", col.names=FALSE)
      session <- rep("ij", length(times))
      occs <- pars$t[i]
      capthist <- data.frame(session, ids, occs, cameras)
      write.table(capthist, capt_file, sep=",", row.names = FALSE, col.names=FALSE)
      secrdata <- secr::read.capthist(capt_file, trap_file, detector = "count", noccasions = pars$t[i])
      
      
      xs <- seq(1000*bounds[1], 1000*bounds[2], length.out=66)
      ys <- seq(1000*bounds[3], 1000*bounds[4], length.out=66)
      coord <- data.frame(x = rep(xs,each=66), y = rep(ys, times=66))
      mask <- secr::read.mask(data=coord)
      
      fit <- secr::secr.fit(secrdata, mask = mask, trace = FALSE)
      out <- secr::region.N(fit)
      
      # Output from secr
      output[i, 15] <- out[2,1] # realised, not expected population
      output[i, 16] <- out[2,2] # se of population
      out1 <- predict(fit)
      output[i, 17] <- (out1[3,2] / 1000) # sigma
      output[i, 18] <- (out1[3,3] / 1000)
      if (!is.na(output[i, 15])){
        output[i, 26] <- 1 # success flag
      }
      
      
    }, error = function(cond){
      print(i)
      print(cond)
    }
  )
  
  # Fit in TMB
  
  # Pre-Processing function
  out <- prepare_SESCR_NIMBLE(camera_locs = camera_locations, times = times, cameras = cameras,
                              ids = ids, bounds = bounds, nrand=4000)
  
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
  
  data <- list(N = out$N, m = out$m, Nobs = length(times) - 1, area = out$area,
               camera_locations = camera_locations, times = events[,1],
               cameras = events[,2] - 1, animals = events[,3] - 1, last_cameras = events[,4] - 1,
               nrand = out$nrand, mask = out$rpts)
  
  # parameter starting points
  NLL <- Inf
  t0 <- Sys.time()
  for(run in 1:20){
    
    out_tmb <- tryCatch(
      {
      # parameter starting points
      param <- list(log_M = log(runif(1,0.01,2)*m), log_sigma = log(runif(1,0.5,2) * output[i,17]),
                logit_d = runif(1,-4,1), log_beta = runif(1,-1,1), log_g0 = runif(1,-10,-4))
                
      obj <- TMB::MakeADFun(data = data, parameters = param, 
                        DLL = "SESCR_SCLMR",  silent = TRUE)
      opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                         control = list(trace = FALSE))
    
    
      if (opt$objective < NLL){
        coefsp = stelfi::get_coefs(obj)
        if (coefsp[1,1] < 1000){
          coefs <- coefsp
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
  output[i, 30] <- as.numeric(difftime(t1,t0,units="mins"))
  output[i, 29] <- 4000.20
  
  out_results <- tryCatch(
    {
      output[i, 19] <- coefs[1,1]
      output[i, 20] <- coefs[1,2]
      output[i, 21] <- coefs[4,1]
      output[i, 22] <- coefs[4,2]
      output[i, 23] <- coefs[3,1]
      output[i, 24] <- coefs[2,1]
      output[i, 25] <- coefs[5,1]
    }, error = function(cond){
      print("No successful run")
    }
      
  )

  
  # binary indicators comparing SCR and SESCR
  if (!is.na(output[i, 19]) && !is.na(output[i, 15])){
    
    if (abs(output[i, 19] - output[i, 2]) < abs(output[i, 15] - output[i, 2])){
      output[i, 27] <- 1
    }
  }
  
  if (!is.na(output[i, 21]) && !is.na(output[i, 17])){
    if (abs(output[i, 21] - output[i, 4]) < abs(output[i, 17] - output[i, 4])){
      output[i, 28] <- 1
    }
  }
  
  # OUTPUT after each iteration
  write.csv(output, paste(file_path, "Results_",N,".csv",sep=""))
}
  
