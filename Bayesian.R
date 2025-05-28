# Fit SESCR to pregenerated data, in a Bayesian framework with NIMBLE
# Last Updated: 28 May 2025

S <- 1 # Set number for the simulations
source1 <- "PR" # name of folder for storing results (OU,PR or SCR)
source2 <- "SESCR" # name of folder where simulations stored (OUSCR, SESCR or SCR)

## Load Necessary functions
source("FitFunctions.R")
library(coda)
library(secr)

## START OF SCRIPT ##
set.seed(S)
file_path = paste("SimStudy-",source1,"-DA/",sep="")
dir.create(file_path)


pars <- read.csv(paste(source2,"Sims/Pars_",S,".csv",sep=""))
n_sims <- nrow(pars)
camera_locations = as.matrix(read.csv(paste(source2,"Sims/Cameras_",S,".csv",sep="")))[,2:3]


# Output Data frame
output <- as.data.frame(matrix(0, nrow = n_sims, ncol = 20))
names(output) <- c("n_obs", "m", "C_obs", 
                   "N_scr", "N_scr_se", "sigma_scr", "sigma_scr_se", "lambda0_scr",
                   "N", "N_se", "sigma", "sigma_se", "lambda0", "beta", "d",
                   "scr_ran", "N_better", "sigma_better", "ess_sigma", "runtime")


trap_file = paste("traps_Sims_",S,".csv",sep="")
capt_file = paste("capthist_Sims_",S,".csv",sep="")


for(i in 1:n_sims){
  # Read in raw data
  identifier = paste(S,"_",i,sep="")
  bounds = c(pars$bounds1[i], pars$bounds2[i], pars$bounds3[i], pars$bounds4[i])
  obs = read.csv(paste(source2,"Sims/Sim_Data_",identifier,".csv",sep=""))
  
  # Capture Histories
  times <- obs$times
  ids <- obs$ids
  cameras <- obs$cameras
  
  m <- max(ids)
  n <- length(times)
  if(n == 0){
    next
  }
  
  # Summaries of the data
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

  # FIT WITH NIMBLE
  out <- prepare_SESCR_NIMBLE(camera_locs = camera_locations, times = times, cameras = cameras,
                              ids = ids, bounds = bounds)
  
  # Superpopulation
  M = 10 * m
  
  # Constants
  constants <- list(J = out$J, m = out$m, K = n, area = out$area,
                    bounds = bounds, camera_locations = camera_locations, M = M)
  
  # Starting ACs are the empirical centroids for observed animals, random for the others.
  initAC = matrix(0, nrow = M, ncol = 2)
  initAC[,1] = runif(M,bounds[1],bounds[2]); initAC[,2] = runif(M,bounds[3],bounds[4])
  initAC[1:m, ] = out$obsAC
  
  # We start with random lambda0 and beta values. d and sigma fixed. 
  # Starting true population is twice the number of observed animals
  initsList <- list(lambda0 = runif(1,0.009,0.027), beta = runif(1,0.2,1),  
                    sigma = 3, Dratio = 0.5, s = initAC, psi = 0.2, 
                    z = c(rep(1,m),rep(0,M - 2*m)))
  
  
  # The NIMBLE function requires:
  # 1. event times and survey duration # 2. Camera IDs
  # 3. Individual IDs # 4. last camera where individual captured was seen
  events <- matrix(0, nrow = n + 1, ncol = 4)
  events[, 1] <- c(times, pars$t[i])
  events[, 2] <- c(cameras, 0)
  events[, 3] <- c(ids, 0)
  
  last_cameras <- numeric(m)
  for(j in 2:n){
    last_cameras[ids[(j-1)]] = cameras[(j-1)]
    events[j, 4] <- last_cameras[ids[j]]
  }
  
  data <- list(events = events)
  
  # Setting up and configuring NIMBLE model
  t0 <- Sys.time()
  Hmodel <- nimbleModel(mcSESCR_DA, constants = constants, 
                        data = data, inits = initsList)
  
  compiled <- compileNimble(Hmodel)
  configured <- configureMCMC(compiled, monitors = c("Nind", "Sigma", "Beta", "lambda0", "Dratio"))
  configured$removeSamplers('s', print = FALSE)
  for(ani in 1:M){
    configured$addSampler(target = paste('s[',ani,',1:2]',sep=""), type = 'RW_block',
                          control = list(adaptScaleOnly = TRUE),
                          silent = TRUE)
  }
  
  
  built <- buildMCMC(configured)
  mcmc <- compileNimble(built)
  
  # Run MCMC Chain in NIMBLE
  mcmc.out <- runMCMC(mcmc = mcmc, niter = 16000, nchains = 1, nburnin = 6000, thin =4, 
                      summary = TRUE)
  t1 <- Sys.time()
  
  # Output from NIMBLE
  output$N[i] <- mcmc.out$summary[3,1]
  output$N_se[i] <- mcmc.out$summary[3,3]
  output$sigma[i] <- mcmc.out$summary[4,1]
  output$sigma_se[i] <- mcmc.out$summary[4,3]
  output$lambda0[i] <- mcmc.out$summary[5,1]
  output$beta[i] <- mcmc.out$summary[1,1]
  output$d[i] <- mcmc.out$summary[2,1]
  
  output$runtime[i] <- as.numeric(difftime(t1,t0,units="mins"))
  output$ess_sigma[i] <- effectiveSize(mcmc.out$samples)[4]
  
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
  write.csv(mcmc.out$summary, paste(file_path, "Summary_",identifier,".csv",sep=""))
  write.csv(mcmc.out$samples, paste(file_path, "Samples_",identifier,".csv",sep=""))
  
  write.csv(output, paste(file_path, "Results_",S,".csv",sep=""))
  print(i)

}
