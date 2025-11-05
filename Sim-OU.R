# Simulate from OU process and save output
# Last Updated: 28 May 2025

n_sims <- 10
S <- 1 # study number

## START OF SCRIPT ##
set.seed(S)
source("SimulateFunctions.R")
dir.create("OUSCRSims")


N <- 50; C <- 25
t <- 10
bounds <- c(0,1,0,1)
epsilon <- 0.004
tau <- 1/rep(c(0.2,0.4,0.6,0.8,1),times=2)
sigma <- ifelse(S %% 2 == 0, 0.12, 0.18)
freq <- rep(c(135,270),each=5)
  
pars <- data.frame(N = N, t = t, C = C, beta = beta, sigma = sigma,
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
  obs <- simOU(pars$N[i], pars$epsilon[i], pars$beta[i], pars$sigma[i], camera_locations, 
                  pars$t[i], s=centers, freq = 100)
  # Save simulation
  sim_data <- data.frame(times = obs$times, ids = obs$ids, cameras = obs$cameras)
  #write.csv(sim_data, paste("OUSCRSims/Sim_Data_",identifier,".csv",sep=""))
}


