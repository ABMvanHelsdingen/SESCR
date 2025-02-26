# Simulate from OU process and save output
# Last Updated: 31 October 2024
args <- commandArgs(trailingOnly = TRUE)
N <- as.numeric(args[1])
n_sims <- 8

## START OF SCRIPT ##
set.seed(2+10*N)
source("SimulateFunctions.R")
dir.create("OUSCRSims")

Ns <- 50; ts <- rep(seq(10,30,length.out=5),each=2)
Cs <- 25; taus <- 0.2; sigmas <- 0.15; epsilon = 0.004
bounds1 <- 0; bounds2 <- 1; bounds3 <- 0; bounds4 <- 1
  
pars <- data.frame(N = Ns, t = ts, C = Cs, tau = taus, sigma = sigmas,
                    epsilon = epsilon, bounds1 = bounds1, bounds2 = bounds2, 
                    bounds3 = bounds3, bounds4 = bounds4)
  
traps <- secr::make.grid(nx = sqrt(Cs), ny = sqrt(Cs), spacing = 500/(sqrt(Cs) - 1))
rownames(traps) = seq(1,Cs)
camera_locations <- matrix(0, nrow = Cs, ncol = 2)
camera_locations[,1] <- 0.25 + 0.001*traps[,1]; camera_locations[,2] <- 0.25 + 0.001*traps[,2]
write.csv(camera_locations, paste("OUSCRSims/Cameras_",N,".csv",sep=""))
write.csv(pars, paste("OUSCRSims/Pars_",N,".csv",sep=""))




for(i in 1:nrow(pars)){
  identifier = paste(N,"_",i,sep="")
  bounds = c(pars$bounds1[i], pars$bounds2[i], pars$bounds3[i], pars$bounds4[i])

  centers <- matrix(0, nrow = pars$N[i], ncol = 2)
  centers[,1] <- runif(pars$N[i], bounds[1], bounds[2])
  centers[,2] <- runif(pars$N[i], bounds[3], bounds[4])
  obs <- simOU(pars$N[i], pars$epsilon[i], pars$tau[i], pars$sigma[i], camera_locations, 
                  pars$t[i], s=centers, freq = 100)
  # Save simulation
  sim_data <- data.frame(times = obs$times, ids = obs$ids, cameras = obs$cameras)
  write.csv(sim_data, paste("OUSCRSims/Sim_Data_",identifier,".csv",sep=""))
}


