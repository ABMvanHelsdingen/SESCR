# Simulate from SESCR Model and save output
# Last Updated: 8 November 2024
args <- commandArgs(trailingOnly = TRUE)
N <- as.numeric(args[1])
n_sims <- 8

## START OF SCRIPT ##
set.seed(2+10*N)
source("SimulateFunctions.R")
dir.create("SESCRSims")

Ns <- 50; ts <- 10
Cs <- 25; mu0s <- runif(n_sims, 0.005,0.01)
betas <- runif(n_sims,0.2,1); sigmas <- runif(n_sims,0.12,0.18); dratios <- runif(n_sims,0.2,0.5)
bounds1 <- 0; bounds2 <- 1; bounds3 <- 0; bounds4 <- 1
  
pars <- data.frame(N = Ns, t = ts, C = Cs, mu0 = mu0s, beta = betas, sigma = sigmas,
                     dratio = dratios, bounds1 = bounds1, bounds2 = bounds2, 
                     bounds3 = bounds3, bounds4 = bounds4)
  
traps <- secr::make.grid(nx = sqrt(Cs), ny = sqrt(Cs), spacing = 500/(sqrt(Cs) - 1))
rownames(traps) = seq(1,Cs)
camera_locations <- matrix(0, nrow = Cs, ncol = 2)
camera_locations[,1] <- 0.25 + 0.001*traps[,1]; camera_locations[,2] <- 0.25 + 0.001*traps[,2]
write.csv(camera_locations, paste("SESCRSims/Cameras_",N,".csv",sep=""))
write.csv(pars, paste("SESCRSims/Pars_",N,".csv",sep=""))



for(i in 1:nrow(pars)){
  identifier = paste(N,"_",i,sep="")
  bounds = c(pars$bounds1[i], pars$bounds2[i], pars$bounds3[i], pars$bounds4[i])
  
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


