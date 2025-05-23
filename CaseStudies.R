# Loads data from Case Studies and fits secr, SCL and DA SESCR
# N = 1 is the American Martens data, 2 the Sitka deer and 3 the Belarus lynx
# Last Updated: 23 May 2025

N <- 2
library(secr)
source("FitFunctions.R")

duration = as.matrix(read.csv("CS/durations.csv", header = FALSE))[N,1]
times = as.vector(read.csv(paste("CS/times_",N,".csv",sep=""), header = TRUE)[,2])
cameras = as.vector(read.csv(paste("CS/cameras_",N,".csv",sep=""), header = TRUE)[,2])
ids = as.vector(read.csv(paste("CS/ids_",N,".csv",sep=""), header = TRUE)[,2])
camera_locations = as.matrix(read.csv(paste("CS/camera_locations_",N,".csv",sep=""), 
                                      header = TRUE)[,2:3])

# scaling factor from processing, buffer and spacing
scales = as.matrix(read.csv("CS/scales.csv", header = FALSE))[N,]
scale <- scales[1]; buffer = scales[2]; spacing = scales[3]

# MCMC parameters
thinf = 4; niter = 50000; nburnin = 10000

## create the mask over the region; different buffers were tested 
## until the estimated parameters remained stable
traps = data.frame(x = 1000*camera_locations[,1], y = 1000*camera_locations[,2])
trap = make.poly(x=traps$x, y=traps$y)
mask = make.mask(trap,buffer=buffer/scale,spacing=spacing/scale,type="trapbuffer")

# Fit in secr
trap_file = paste("traps_Sims",N,".csv",sep="")
capt_file = paste("capthist_Sims",N,".csv",sep="")
write.table(1000*camera_locations, trap_file, sep=",", col.names=FALSE)
session <- rep("ij", length(times))
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


# SESCR

# TMB Template
if (!"MLEFunction" %in% getLoadedDLLs()){
  TMB::compile("MLEFunction.cpp")
  dyn.load(TMB::dynlib("MLEFunction"))
}

runs <- 20
m <- max(ids)
n <- length(times)
out <- prepare_SESCR_NIMBLE(camera_locs = camera_locations, times = times, cameras = cameras,
                            ids = ids)

# Calculate the last camera where the individual currently detected was seen at
events <- matrix(0, nrow = n + 1, ncol = 1)
last_cameras <- numeric(m)
for(j in 2:n){
  last_cameras[ids[(j-1)]] = cameras[(j-1)]
  events[j, 1] <- last_cameras[ids[j]]
}


data <- list(J = out$J, m = out$m, K = n, area = out$area,
             camera_locations = camera_locations, times = c(times, pars$t[i]),
             cameras = c(cameras, 0) - 1, animals = c(ids, 0) - 1, last_cameras = events[,1] - 1,
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
                            DLL = "MLEFunction",  silent = TRUE)
      opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                           control = list(trace = FALSE))
      
      
      if (opt$objective < NLL){
        coefsr = stelfi::get_coefs(obj)
        # Sometimes the MLEs diverge to nonsensical results
        if (coefsr[1,1] < 10000){
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

out_results <- tryCatch(
  {
    write.csv(coefs, paste("CS/outputC",N,".csv",sep=""))
  }, error = function(cond){
    print("No successful run")
  }
  
)

outputD = data.frame(runtime = as.numeric(difftime(t1,t0,units="mins")),
                     mask = nrow(mask), runs = runs)

write.csv(outputD, paste("CS/outputD",N,".csv",sep=""))

# 22 November
# Bayesian framework with 1/beta ~ U(0,1000)
mask = 0.001*mask
a = min(mask[,1]); b = max(mask[,1])
c = min(mask[,2]); d = max(mask[,2])
bounds = c(a,b,c,d)
m <- max(ids)
M = 10 * m


initAC = matrix(0, nrow = M, ncol = 2)
initAC[,1] = runif(M,a,b); initAC[,2] = runif(M,c,d)
initAC[1:m, ] = out$obsAC

# Prepare data, constants and initial values of parameters
data <- list(events = events)

constants <- list(J = out$J, m = out$m, K = length(times), area = out$area,
                  bounds = bounds, camera_locations = camera_locations, M = M,
                  mask = as.matrix(mask))


initsList <- list(lambda0 = runif(1,0.0009,0.0027), beta = runif(1,1,5),  
                  sigma = 3, Dratio = 0.5, s = initAC, psi = 0.2, 
                  z = c(rep(1,m),rep(0,M - 2*m)))


t0 <- Sys.time()
Hmodel <- nimbleModel(mcSESCR_DA_Irregular, constants = constants, 
                      data = data, inits = initsList)

compiled <- compileNimble(Hmodel)
configured <- configureMCMC(compiled, monitors = c("Nind", "Sigma", "Beta", "lambda0", "Dratio", "s",
                                                                      "z"))
configured$removeSamplers('s', print = FALSE)
for(ani in 1:M){
  configured$addSampler(target = paste('s[',ani,',1:2]',sep=""), type = 'RW_block',
                        control = list(adaptScaleOnly = TRUE),
                        silent = TRUE)
}


built <- buildMCMC(configured)
mcmc <- compileNimble(built)
mcmc.out <- runMCMC(mcmc = mcmc, niter = niter, nchains = 1, nburnin = nburnin, thin =thinf, 
                    summary = TRUE)

t1 = Sys.time()

write.csv(mcmc.out$summary, paste("CS/outputE",N,".csv",sep=""))

# OUTPUT from NIMBLE
write.csv(mcmc.out$samples, paste("CS/outputF",N,".csv",sep=""))

# NIMBLE metadata
outputG = data.frame(runtime = as.numeric(difftime(t1,t0,units="mins")),
                     niter = niter, nburnin = nburnin, thinf = thinf)

write.csv(outputG, paste("CS/outputG",N,".csv",sep=""))