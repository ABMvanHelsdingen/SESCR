# Simulate from OU process and save output
# Last Updated: 28 May 2025

n_sims <- 1000
S <- 1 # study number

## START OF SCRIPT ##
set.seed(S)
source("SimulateFunctions.R")
dir.create("OUSCRSims")

N <- 50; C <- 25
t <- rep(seq(10,30,length.out=5),each=n_sims/5)
bounds <- c(0,1,0,1)
epsilon <- 0.004
tau <- 0.2
sigma <- 0.15
  
pars <- data.frame(N = N, t = t, C = C, tau = tau, sigma = sigma,
                    epsilon = epsilon, bounds1 = bounds[1], bounds2 = bounds[2], 
                    bounds3 = bounds[3], bounds4 = bounds[4])
  
locs <- seq(0.25,0.75,length.out=sqrt(C))
camera_locations <- matrix(unlist(expand.grid(locs, locs)), nrow = C, ncol = 2)
write.csv(camera_locations, paste("OUSCRSims/Cameras_",S,".csv",sep=""))
write.csv(pars, paste("OUSCRSims/Pars_",S,".csv",sep=""))


for(i in 1:n_sims){
  identifier = paste(S,"_",i,sep="")

  centers <- matrix(0, nrow = pars$N[i], ncol = 2)
  centers[,1] <- runif(pars$N[i], bounds[1], bounds[2])
  centers[,2] <- runif(pars$N[i], bounds[3], bounds[4])
  obs <- simOU(pars$N[i], pars$epsilon[i], pars$tau[i], pars$sigma[i], camera_locations, 
                  pars$t[i], s=centers, freq = 100)
  # Save simulation
  sim_data <- data.frame(times = obs$times, ids = obs$ids, cameras = obs$cameras)
  write.csv(sim_data, paste("OUSCRSims/Sim_Data_",identifier,".csv",sep=""))
}


