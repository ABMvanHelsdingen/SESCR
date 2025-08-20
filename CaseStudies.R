# Loads data from Case Studies and fits secr, frequentist and DA SESCR
# N = 1 is the American Martens data, 2 the Sitka deer, 3 the Nepal leopards
# Last Updated: 20 August 2025

# Outputs for each case study are:
# M: secr mask
# A: SCR population estimates
# B: SCR parameter estimates
# C: SESCR frequentist parameter estimates
# D: metadata about SCR and SESCR fitting
# E: Summaries of MCMC chains for SESCR
# F: full MCMC chains for SESCR
# G: metadata about SESCR Data Augmentation fitting

N <- 3
library(secr)
source("FitFunctions.R")

duration = as.matrix(read.csv("CS/durations.csv", header = FALSE))[N,1]
times = as.vector(read.csv(paste("CS/times_",N,".csv",sep=""), header = TRUE)[,2])
cameras = as.vector(read.csv(paste("CS/cameras_",N,".csv",sep=""), header = TRUE)[,2])
ids = as.vector(read.csv(paste("CS/ids_",N,".csv",sep=""), header = TRUE)[,2])
camera_locations = as.matrix(read.csv(paste("CS/camera_locations_",N,".csv",sep=""), 
                                      header = TRUE)[,2:3])

m <- max(ids)
n <- length(times)

# scaling factor from processing, buffer and spacing
scales = as.matrix(read.csv("CS/scales.csv", header = TRUE))[N,]
buffer = scales[2]; spacing = scales[3]

# MCMC fitting parameters
thinf = 4; niter = 50000; nburnin = 10000

# TMB fitting parameters
runs <- 20

# Define trap locations and masks
traps = data.frame(x = 1000*camera_locations[,1], y = 1000*camera_locations[,2])
trap = make.poly(x=traps$x, y=traps$y)
mask = make.mask(trap,buffer=buffer,spacing=spacing,type="trapbuffer")

# Fit in secr
t0 <- Sys.time()
trap_file = paste("traps_Sims",N,".csv",sep="")
capt_file = paste("capthist_Sims",N,".csv",sep="")
write.table(1000*camera_locations, trap_file, sep=",", col.names=FALSE)
session <- rep("a", n)
capthist <- data.frame(session, ids, 1, cameras)
write.table(capthist, capt_file, sep=",", row.names = FALSE, col.names=FALSE)
secrdata <- secr::read.capthist(capt_file, trap_file, detector = "count", noccasions = 1)


traps <- traps(secrdata)

write.csv(mask, paste("CS/outputM",N,".csv",sep=""))

fit <- secr::secr.fit (secrdata, mask = mask, trace = FALSE, detectfn = "HHN")
out <- secr::region.N(fit)

write.csv(out, paste("CS/outputA",N,".csv",sep=""))

out_secr <- predict(fit)

write.csv(out_secr, paste("CS/outputB",N,".csv",sep=""))

t1 <- Sys.time()
outputD = data.frame(runtime_scr = as.numeric(difftime(t1,t0,units="mins")),
                     mask = nrow(mask), runs = runs, area = 0.01*maskarea(mask))

# SESCR

# TMB Template
if (!"SCLFunction" %in% getLoadedDLLs()){
  TMB::compile("SCLFunction.cpp")
  dyn.load(TMB::dynlib("SCLFunction"))
}

out <- prepare_SESCR(camera_locs = camera_locations, times = times, cameras = cameras,
                            ids = ids)

# The NIMBLE and TMB functions require:
# 1. event times and survey duration # 2. Camera IDs
# 3. Individual IDs # 4. last camera where individual captured was seen
events <- matrix(0, nrow = n + 1, ncol = 4)
events[, 1] <- c(times, duration)
events[, 2] <- c(cameras, 0)
events[, 3] <- c(ids, 0)

last_cameras <- numeric(m)
for(j in 2:n){
  last_cameras[ids[(j-1)]] = cameras[(j-1)]
  events[j, 4] <- last_cameras[ids[j]]
}


data <- list(J = out$J, m = m, K = n, area = 0.01*maskarea(mask),
             camera_locations = camera_locations, times = events[,1],
             cameras = events[,2] - 1, animals = events[,3] - 1, last_cameras = events[,4] - 1,
             nrand = nrow(mask), mask = as.matrix(0.001*mask))

NLL <- Inf
t0 <- Sys.time()
for(run in 1:runs){
  
  out_tmb <- tryCatch(
    {
      # parameter starting points
      param <- list(log_N = log(runif(1,0.01,2)*m), log_sigma = log(0.001 * runif(1,0.5,2) * out_secr[3,2]),
                    logit_d = runif(1,-4,1), log_beta = runif(1,-1,1), log_lambda0 = runif(1,-10,-4))
      
      obj <- TMB::MakeADFun(data = data, parameters = param, 
                            DLL = "SCLFunction",  silent = TRUE)
      opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                           control = list(trace = FALSE))
      
      
      if (opt$objective < NLL){
        coefsr <- summary(TMB::sdreport(obj), "report") 
        # Occasionally the MLEs diverge to nonsensical results
        if (coefsr[1,1] < 1000){
          coefs <- coefsr
          NLL <- opt$objective
        }
      }
    }, error = function(cond){
      print(cond)
      print(run)
    }
  )
}
t1 <- Sys.time()

out_results <- tryCatch(
  {
    write.csv(coefs, paste("CS/outputC",N,".csv",sep=""))
  }, error = function(cond){
    print("No successful run")
  }
  
)

outputD$runtime <- as.numeric(difftime(t1,t0,units="mins"))
write.csv(outputD, paste("CS/outputD",N,".csv",sep=""))

# BAYESIAN FRAMEWORK

# Prepare data, constants and initial values of parameters
data <- list(events = events)
# Superpopulation size
M = 8*m

# Prepare mask, which is not identical to secr mask
maskv = perimeterMask(mask)
mask_bay <- as.matrix(0.001*maskv$perimeter)

a = min(mask_bay[,1]); b = max(mask_bay[,1])
c = min(mask_bay[,2]); d = max(mask_bay[,2])
bounds = c(a,b,c,d)

constants <- list(J = out$J, m = out$m, K = length(times), area = 1e-6 * maskv$area,
                  bounds = bounds, camera_locations = camera_locations, M = M,
                  mask = mask_bay)


initAC = matrix(0, nrow = M, ncol = 2)
initAC[,1] = runif(M,a,b); initAC[,2] = runif(M,c,d)
initAC[1:m, ] = out$obsAC

# Use MLEs as starting values for lambda0,beta and sigma.
initsList <- list(lambda0 = coefs[3,1], beta = 1/coefs[2,1],  
                  sigma = log(coefs[4,1]^-2 / 2), Dratio = 0.5, s = initAC, psi = 2*m/M, 
                  z = c(rep(1,m),rep(0,M - 2*m)))


# Setting up and configuring NIMBLE model
t0 <- Sys.time()
Hmodel <- nimbleModel(mcSESCR_DA_Irregular, constants = constants, 
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
mcmc.out <- runMCMC(mcmc = mcmc, niter = niter, nchains = 1, nburnin = nburnin, thin =thinf, 
                    summary = TRUE)
t1 = Sys.time()

# Save output from NIMBLE
write.csv(mcmc.out$summary, paste("CS/outputE",N,".csv",sep=""))
write.csv(mcmc.out$samples, paste("CS/outputF",N,".csv",sep=""))


# Save metadata about fitting
outputG = data.frame(runtime = as.numeric(difftime(t1,t0,units="mins")),
                     niter = niter, nburnin = nburnin, thinf = thinf, area = constants$area)
write.csv(outputG, paste("CS/outputG",N,".csv",sep=""))
