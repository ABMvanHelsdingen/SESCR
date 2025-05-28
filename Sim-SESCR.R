# Simulate from SESCR Model and save output
# Last Updated: 23 May 2025

n_sims <- 1000
S <- 1 # study number

## START OF SCRIPT ##
set.seed(S)
source("SimulateFunctions.R")
dir.create("SESCRSims")

N <- 50; t <- 10; C <- 25
bounds <- c(0,1,0,1)
mu0 <- runif(n_sims, 0.005,0.01)
beta <- runif(n_sims,0.2,1)
sigma <- runif(n_sims,0.12,0.18)
dratio <- runif(n_sims,0.2,0.5)
  
pars <- data.frame(N = N, t = t, C = C, mu0 = mu0, beta = beta, sigma = sigma,
                     dratio = dratio, bounds1 = bounds[1], bounds2 = bounds[2], 
                     bounds3 = bounds[3], bounds4 = bounds[4])

locs <- seq(0.25,0.75,length.out=sqrt(C))
camera_locations <- matrix(unlist(expand.grid(locs, locs)), nrow = C, ncol = 2)
write.csv(camera_locations, paste("SESCRSims/Cameras_",S,".csv",sep=""))
write.csv(pars, paste("SESCRSims/Pars_",S,".csv",sep=""))


for(i in 1:n_sims){
  identifier = paste(S,"_",i,sep="")
  
  centers <- matrix(0, nrow = pars$N[i], ncol = 2)
  centers[,1] <- runif(pars$N[i], bounds[1], bounds[2])
  centers[,2] <- runif(pars$N[i], bounds[3], bounds[4])
  Sigma = matrix(c(pars$sigma[i] ^ 2, 0,0,pars$sigma[i] ^2), nrow = 2)
  d = (pars$dratio[i] ^2) * Sigma
  obs <- simSESCR(N = pars$N[i], mu0 = pars$mu0[i], beta = pars$beta[i], 
                        d = d, sigma = Sigma, camera_locations = camera_locations, t = pars$t[i],
                        s = centers)
  # Save simulation
  sim_data <- data.frame(times = obs$times, ids = obs$ids, cameras = obs$cameras)
  write.csv(sim_data, paste("SESCRSims/Sim_Data_",identifier,".csv",sep=""))
}


