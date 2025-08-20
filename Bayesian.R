# Fit SESCR to pregenerated data, in a Bayesian framework with NIMBLE
# Last Updated: 17 August 2025

S <- 1 # Set number for the simulations
source <- "SESCR" # name of folder where simulations stored (OUSCR, SESCR or SCR)
size <- 64 # number of points across in the habitat mask

## Load Necessary functions
source("FitFunctions.R")
library(coda)
library(secr)

## START OF SCRIPT ##
set.seed(S)
file_path = paste("SimStudy-",source,"-DA/",sep="")
dir.create(file_path)


pars <- read.csv(paste(source,"Sims/Pars_",S,".csv",sep=""))
n_sims <- nrow(pars)
camera_locations = as.matrix(read.csv(paste(source,"Sims/Cameras_",S,".csv",sep="")))[,2:3]


# Output Data frame
output <- as.data.frame(matrix(0, nrow = n_sims, ncol = 21))
names(output) <- c("n_obs", "m", "C_obs",
                   "N_scr", "N_scr_se", "sigma_scr", "sigma_scr_se", "lambda0_scr",
                   "scr_ran", "runtime_scr",
                   "N", "N_se", "sigma", "sigma_se", "lambda0", "beta", "d",
                   "runtime", "ess_sigma",
                   "N_better", "sigma_better")


trap_file = paste("traps_Sims_",S,".csv",sep="")
capt_file = paste("capthist_Sims_",S,".csv",sep="")


for(i in 1:n_sims){
  # Read in raw data
  identifier = paste(S,"_",i,sep="")
  bounds = c(pars$bounds1[i], pars$bounds2[i], pars$bounds3[i], pars$bounds4[i])
  obs = read.csv(paste(source,"Sims/Sim_Data_",identifier,".csv",sep=""))
  
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
  
  # Make mask
  gap <- 0.5/size # Points must be grid centroids, thus offset from boundaries
  xs <- seq(gap*bounds[2] + (1-gap)*bounds[1], gap*bounds[1] + (1-gap)*bounds[2], length.out=size)
  ys <- seq(gap*bounds[4] + (1-gap)*bounds[3], gap*bounds[3] + (1-gap)*bounds[4], length.out=size)
  coord <- data.frame(x = rep(1000*xs,each=size), y = rep(1000*ys, times=size))
  mask <- secr::read.mask(data=coord)
  
  # Run SCR, using secr package
  t0 <- Sys.time()
  out_scr <- tryCatch(
    {
      write.table(1000*camera_locations, trap_file, sep=",", col.names=FALSE)
      session <- rep("a", n)
      capthist <- data.frame(session, ids, 1, cameras)
      write.table(capthist, capt_file, sep=",", row.names = FALSE, col.names=FALSE)
      secrdata <- secr::read.capthist(capt_file, trap_file, detector = "count", noccasions = 1)
      
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
  
  # Record runtime
  t1 <- Sys.time()
  output$runtime_scr[i] <- as.numeric(difftime(t1,t0,units="mins"))

  # FIT WITH NIMBLE
  out <- prepare_SESCR(camera_locs = camera_locations, times = times, cameras = cameras,
                              ids = ids)
  
  
  # Superpopulation size
  M = 8*m

  
  
  # Constants
  constants <- list(J = out$J, m = out$m, K = n, area = 0.01*maskarea(mask),
                    bounds = bounds, camera_locations = camera_locations, M = M)
  
  # Starting ACs are the empirical centroids for observed animals, random for the others.
  initAC = matrix(0, nrow = M, ncol = 2)
  initAC[,1] = runif(M,bounds[1],bounds[2]); initAC[,2] = runif(M,bounds[3],bounds[4])
  initAC[1:m, ] = out$obsAC
  
  # Starting true population is twice the number of observed animals
  initsList <- list(lambda0 = 0.05, beta = 1,  
                    sigma = 3, Dratio = 0.5, s = initAC, psi = 2*m/M, 
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
  configured <- configureMCMC(compiled, monitors = c("Nind", "Sigma", "Beta", "lambda0", "Dratio", "s"))
  configured$removeSamplers('s', print = FALSE)
  for(ani in 1:M){
    configured$addSampler(target = paste('s[',ani,',1:2]',sep=""), type = 'RW_block',
                          control = list(adaptScaleOnly = TRUE),
                          silent = TRUE)
  }
  
  
  built <- buildMCMC(configured)
  mcmc <- compileNimble(built)
  
  # Run MCMC Chain in NIMBLE
  mcmc.out <- runMCMC(mcmc = mcmc, niter = 20000, nchains = 1, nburnin = 10000, thin =4, 
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
  output$N_better[i] <- compare_estimates(pars$N[i], output$N_scr[i],
                                          output$N[i])
  output$sigma_better[i] <- compare_estimates(pars$sigma[i], output$sigma_scr[i],
                                              output$sigma[i])
   

  # Save after each iteration
  write.csv(mcmc.out$summary, paste(file_path, "Summary_",identifier,".csv",sep=""))
  write.csv(mcmc.out$samples, paste(file_path, "Samples_",identifier,".csv",sep=""))
  
  write.csv(output, paste(file_path, "Results_",S,".csv",sep=""))
  print(i)

}
