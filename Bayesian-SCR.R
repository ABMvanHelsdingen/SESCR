# SCR Simulation Study
# Last Updated: 8 November 2024

# This block is designed to be run on the NeSI server
args <- commandArgs(trailingOnly = TRUE)
N <- as.numeric(args[1])
n_sims <- 8

library(coda)

## START OF SCRIPT ##
set.seed(2+10*N)
source("FitFunctions.R")

file_path = "SimStudy-SCR-NIMBLE/"
dir.create(file_path)


pars <- read.csv(paste("SCRSims/Pars_",N,".csv",sep=""))
camera_locations = as.matrix(read.csv(paste("SCRSims/Cameras_",N,".csv",sep="")))[,2:3]


# Output Data frame

output <- as.data.frame(matrix(0, nrow = nrow(pars), ncol = 28))
names(output) <- c("n_cameras", "M_true", "t", "sigma_true", "mu0_true",
                   "n_obs", "m", "C_obs", "m_1obs", "m_maxobs", "m_1C", "m_maxC",
                   "M_scr", "M_scr_se", "sigma_scr", "sigma_scr_se",
                   "M", "M_se", "sigma", "sigma_se", "g0", "beta", "d",
                   "scr_ran", "M_better", "sigma_better", "ess_sigma", "runtime")

trap_file = paste("traps_Sims_",N,".csv",sep="")
capt_file = paste("capthist_Sims_",N,".csv",sep="")


for(i in 1:nrow(pars)){
  identifier = paste(N,"_",i,sep="")
  bounds = c(pars$bounds1[i], pars$bounds2[i], pars$bounds3[i], pars$bounds4[i])
 
  
  obs = read.csv(paste("SCRSims/Sim_Data_",identifier,".csv",sep=""))
  pars <- read.csv(paste("SCRSims/Pars_",N,".csv",sep=""))
  

  times <- obs$times
  ids <- obs$ids
  cameras <- obs$cameras
  m <- max(ids)

  n <- length(times); print(i); print(n); print(m)
  if(n == 0){
    next
  }
  
  # Summary of the data
  output[i, 1] <- nrow(camera_locations)
  output[i, 2] <- pars$M[i]
  output[i, 3] <- pars$t[i]
  output[i, 4] <- pars$sigma[i]
  output[i, 5] <- pars$mu0[i]
  
  output[i, 6] <- n
  output[i, 7] <- m
  output[i, 8] <- length(unique(cameras))
  
  tab <- table(ids)
  output[i, 9] <- sum(tab == 1)
  output[i, 10] <- max(tab)
  
  TWTable = table(ids, cameras) > 0
  CPI = rowSums(TWTable) # cameras per individual
  
  output[i, 11] <- sum(CPI == 1)
  output[i, 12] <- max(CPI)
  
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
      output[i, 13] <- out[2,1] # realised, not expected population
      output[i, 14] <- out[2,2] # se of population
      out1 <- predict(fit)
      output[i, 15] <- (out1[3,2] / 1000) # sigma
      output[i, 16] <- (out1[3,3] / 1000)
      if (!is.na(output[i, 15])){
        output[i, 24] <- 1 # success flag
      }
      
      
    }, error = function(cond){
      print(i)
      print(cond)
    }
  )

  # FIT WITH NIMBLE
  out <- prepare_SESCR_NIMBLE(camera_locs = camera_locations, times = times, cameras = cameras,
                              ids = ids, bounds = bounds, nrand=2000)
  M = 10 * m
  
  initA = matrix(0, nrow = M, ncol = 2)
  initA[,1] = runif(M,bounds[1],bounds[2]); initA[,2] = runif(M,bounds[3],bounds[4])
  initA[1:m, ] = out$obsA
  
  
  constants <- list(N = out$N, m = out$m, Nobs = length(times), area = out$area,
                    bounds = bounds, camera_locations = camera_locations, M = M)
  
  
  initsList <- list(g0 = runif(1,0.009,0.027), beta = runif(1,0.5,2),  
                    sigma = 3, Dratio = 0.3, a = initA, psi = 0.2, 
                    z = c(rep(1,m),rep(0,M - 2*m)))
  
  
  events <- matrix(0, nrow = length(times) + 1, ncol = 4)
  events[, 1] <- c(times, pars$t[i])
  events[, 2] <- c(cameras, 0)
  events[, 3] <- c(ids, 0)
  last_cameras <- numeric(m)
  for(j in 2:length(times)){
    last_cameras[ids[(j-1)]] = cameras[(j-1)]
    events[j, 4] <- last_cameras[ids[j]]
    
  }
  data <- list(events = events)
  
  t0 <- Sys.time()
  
  Hmodel <- nimbleModel(mcSESCR_DA, constants = constants, 
                        data = data, inits = initsList)
  
  compiled <- compileNimble(Hmodel)
  configured <- configureMCMC(compiled, enableWAIC = TRUE, monitors = c("Nind", "Sigma", "Beta", "g0", "Dratio"))
  configured$removeSamplers('a', print = FALSE)
  for(ani in 1:M){
    configured$addSampler(target = paste('a[',ani,',1:2]',sep=""), type = 'RW_block',
                          control = list(adaptScaleOnly = TRUE),
                          silent = TRUE)
  }
  
  
  built <- buildMCMC(configured)
  mcmc <- compileNimble(built)
  mcmc.out <- runMCMC(mcmc = mcmc, niter = 16000, nchains = 1, nburnin = 6000, thin =4, 
                      summary = TRUE, WAIC = TRUE)
  t1 <- Sys.time()
  output[i, 28] <- as.numeric(difftime(t1,t0,units="mins"))
  ess <- effectiveSize(mcmc.out$samples)
  output[i, 27] <- ess[3]
  
  # OUTPUT from NIMBLE
  output[i, 17] <- mcmc.out$summary[3,1] # N
  output[i, 18] <- mcmc.out$summary[3,3] # N_se
  output[i, 19] <- mcmc.out$summary[4,1] # sigma
  output[i, 20] <- mcmc.out$summary[4,3] # sigma_se
  output[i, 21] <- mcmc.out$summary[5,1] # g0
  output[i, 22] <- mcmc.out$summary[1,1] # beta
  output[i, 23] <- mcmc.out$summary[2,1] # d
  
  # binary indicators comparing SCR and SESCR
  if (!is.na(output[i, 17]) && !is.na(output[i, 13])){
    
    if (abs(output[i, 17] - output[i, 2]) < abs(output[i, 13] - output[i, 2])){
      output[i, 25] <- 1
    }
  }
  
  if (!is.na(output[i, 19]) && !is.na(output[i, 15])){
    if (abs(output[i, 19] - output[i, 4]) < abs(output[i, 15] - output[i, 4])){
      output[i, 26] <- 1
    }
  }
   

  # OUTPUT after each iteration
  write.csv(mcmc.out$summary, paste(file_path, "Summary_",identifier,".csv",sep=""))
  write.csv(mcmc.out$samples, paste(file_path, "Samples_",identifier,".csv",sep=""))
  
  write.csv(output, paste(file_path, "Results_",N,".csv",sep=""))

}