# Simulate from SESCR Model and save output
# Last Updated: 14 August 2025

n_sims <- 10
S <- 1 # set number

## START OF SCRIPT ##
set.seed(S)
source("SimulateFunctions.R")
dir.create("SESCRSims")

# Create grid of all parameter values and select 10 for this set
pars <- expand.grid(beta = c(0.2,0.4,0.6,0.8,1), sigma = c(0.12,0.18),
                    mu0 = c(0.005,0.01), r = c(0.25,0.5))

rem <- (S %% 4) + 1
pars <- pars[(1 + 10*(rem - 1)): (10*rem), ]

pars$N <- 50
pars$t <- 10
pars$C <- 25

bounds <- c(0,1,0,1)
pars$bounds1 <- bounds[1]
pars$bounds2 <- bounds[2]
pars$bounds3 <- bounds[3]
pars$bounds4 <- bounds[4]

# Camera locations
locs <- seq(0.25,0.75,length.out=sqrt(pars$C[1]))
camera_locations <- matrix(unlist(expand.grid(locs, locs)), nrow = pars$C[1], ncol = 2)

# Save parameter values and camera locations
write.csv(camera_locations, paste("SESCRSims/Cameras_",S,".csv",sep=""))
write.csv(pars, paste("SESCRSims/Pars_",S,".csv",sep=""))


for(i in 1:n_sims){
  identifier = paste(S,"_",i,sep="")
  
  centers <- matrix(0, nrow = pars$N[i], ncol = 2)
  centers[,1] <- runif(pars$N[i], bounds[1], bounds[2])
  centers[,2] <- runif(pars$N[i], bounds[3], bounds[4])
  Sigma = matrix(c(pars$sigma[i] ^ 2, 0,0,pars$sigma[i] ^2), nrow = 2)
  # Here d refers to the variance of the self-excitement in space
  d = (pars$r[i] ^2) * Sigma
  obs <- simSESCR(N = pars$N[i], mu0 = pars$mu0[i], beta = pars$beta[i], 
                        d = d, sigma = Sigma, camera_locations = camera_locations, t = pars$t[i],
                        s = centers)
  # Save simulation
  sim_data <- data.frame(times = obs$times, ids = obs$ids, cameras = obs$cameras)
  write.csv(sim_data, paste("SESCRSims/Sim_Data_",identifier,".csv",sep=""))
  print(length(obs$times))
}


