library(nimble)

# This function takes in raw data, checks it for inconsistencies
# Outputs several useful derived quantities

### DATA PROCESSING ###

prepare_SESCR_NIMBLE <- function(camera_locs, times, cameras, ids, params,
                                 bounds = c(0,1,0,1), nrand = 2000, 
                                 tmb_silent = TRUE, nlminb_silent = TRUE) {
  # cameras and ids must be 1-indexed (and this is used in the NIMBLE function)
  
  ## error checking
  for (i in 2:length(times)) {
    if ((times[i] - times[i - 1]) < 1.e-10)
      stop("times must be in ascending order with no simultaneous events")
  }
  
  m = max(ids) # number of observed animals
  N = nrow(camera_locs)
  
  # Calculate number of camera-animal pairs with detections
  streams = numeric(length=length(times))
  events_per_stream = numeric(length=m*N)
  
  for(i in 1:length(times)){
    sn = ((ids[i] - 1) * N) + cameras[i]
    streams[i] = sn
    events_per_stream[sn] = events_per_stream[sn] + 1
  }
  
  s = sum(events_per_stream > 0)
  
  # Mask
  # Generate Random points within rectangular study area
  rpts = matrix(0, nrow = nrand, ncol = 2)
  rpts[, 1] = runif(nrand, bounds[1], bounds[2])
  rpts[, 2] = runif(nrand, bounds[3], bounds[4])
  
  
  # Calculate area of Survey region
  area = (bounds[2] - bounds[1]) * (bounds[4] - bounds[3])
  
  # initial estimates for activity centers are the mean of the observation locations
  a = matrix(0, nrow = m, ncol= 3)
  for (i in 1:length(ids)){
    An = ids[i]
    C = cameras[i]
    a[An, 1:2] = a[An, 1:2] + camera_locs[C, ]
    a[An, 3] = a[An, 3] + 1
  }
  a[, 1:2] = a[, 1:2] / a[,3]
  a = a[, 1:2]
  
  output <- list(N = N, m = m, area = area, s = s, obsA = a, rpts = rpts, nrand = nrand)
  
  return(output)
}

### FIT SESCR IN NIMBLE ###

dSESCR_DA <- nimbleFunction(
  run = function(x = double(2),
                 N = integer(0), # number of cameras
                 m = integer(0), # number of detected animals
                 Nobs = integer(0), # total number of detections
                 area = double(0), # area of study area
                 camera_locations = double(2), # 2D camera locations
                 z = double(1), # binary indicators for unobserved animals
                 M = integer(0), # super-population size
                 psi = double(0), # probability that unobserved animals exist
                 lambda0 = double(0), # detection rate at ac for the baseline
                 beta = double(0), # self-excitement decay rate
                 Dratio = double(0), # d/sigma
                 a = double(2), # activity centers
                 sigma = double(0), # Home Range Covariance
                 log = integer(0, default = 0)) {
    returnType(double(0))
    times <- x[,1] # entry Nobs + 1 is the survey duration
    cameras <- x[,2] # cameras of each detection
    animals <- x[,3] # individuals of each detection
    last_cameras <- x[,4] # last detection camera
    
    
    #Indexing scheme:
    #i and k for cameras, j for times/events, ani for individual animals
    
    
    # Background detection rates (i.e. the home range) is a bivariate normal
    loci <- numeric(2)
    DF <- matrix(0, nrow = m, ncol = N) # Detectability at each animal-camera combination
    DFsums <- numeric(m) # Sum of DF for each animal
    
    for(ani in 1:m){
      for(i in 1:N){
        loci = camera_locations[i, ] - a[ani, ]
        dist2 = loci[1]^2 + loci[2]^2
        DF[ani, i] = exp(-dist2/(2*sigma^2))
        DFsums[ani] = DFsums[ani] + DF[ani, i]
      }
    }
    
    
    #Matrix of spatial self-excitement between cameras
    sse = matrix(0, nrow = N, ncol = N)
    d = Dratio * sigma
    for(i in 1:N){
      for(k in i:N){
        loci = camera_locations[i, ] - camera_locations[k, ]
        dist2 = loci[1]^2 + loci[2]^2
        sse[i,k] = exp(-dist2/(2*d^2)) / (Dratio^2)
        sse[k,i] = sse[i,k] # Matrix is symmetric
      }
    }
    
    
    duration = times[(Nobs+1)]
    
    # A[ani, j] determines the self-exciting effect of animal ani at event j
    # The other variables are used to calculate the sums of the integral for each stream
    A = matrix(0, nrow = m, ncol = Nobs + 1)
    decay_marks = numeric(m); spike_marks = numeric(m)
    decay_sums = numeric(m); spike_sums = numeric(m)
    
    for(j in 2:(Nobs + 1)){
      An = animals[(j-1)] # most recent detected animal
      for(ani in 1:m){
        if (ani == An){
          A[ani,j] = exp(-beta * (times[j] - times[(j - 1)])) # Markovian, no dependency on previous events
        } else {
          A[ani,j] = exp(-beta * (times[j] - times[(j - 1)])) * A[ani,(j-1)]
        }
      }
      
      # Areas of spikes and decays
      # Animal seen at observation j - 1 needs a new (truncated) exponential decay curve to begin
      # If the animal was detected previously, a curve needs to finish and the area be calculated
      # Integral = sum of marks - sum(mark_i * A_i+1) up to i = N -1
      if (decay_marks[An] > 0) { # Animal has been detected before
        decay_sums[An] = decay_sums[An] - decay_marks[An] * A[An,(j-1)]
        spike_sums[An] = spike_sums[An] - spike_marks[An] * A[An,(j-1)]
      }
      C = cameras[(j-1)]
      decay_marks[An] = DFsums[An]
      spike_marks[An] = sum(sse[,C])
      
      decay_sums[An] = decay_sums[An] + decay_marks[An]
      spike_sums[An] = spike_sums[An] + spike_marks[An]
    }
    
    # At end of survey, all exponential curves remaining end
    for(ani in 1:m){
      decay_sums[ani] = decay_sums[ani] - decay_marks[ani] * A[ani,(Nobs+1)]
      spike_sums[ani] = spike_sums[ani] - spike_marks[ani] * A[ani,(Nobs+1)]
    }
    
    
    #Third term of NLL: Sum of the log of detection rates when detections occur. 
    nll = 0
    for(j in 1:Nobs){
      An = animals[j]
      C = cameras[j]
      # Calculate lambda at the time of detection
      if (last_cameras[j] == 0){
        nll = nll - log(DF[An,C])
      } else{
        nll = nll - log(((1 - A[An,j])*DF[An,C]) + (A[An,j] * sse[last_cameras[j], C]))
      }
    }
    
    nll = nll - Nobs * log(lambda0)
    
    
    
    # First two terms of NLL: Sum of the Integral of Lambda for all streams
    nll = nll + (sum(spike_sums) - sum(decay_sums)) * (lambda0 / beta) # spikes and decay areas
    nll = nll + sum(DF) * (lambda0 * duration) # baseline
    
    #CONTRIBUTION TO LIKELIHOOD FROM UNDETECTED ANIMALS
    
    for(uani in 1:(M-m)){
      if (z[uani] > 0){
        DFR = 0
        for(i in 1:N){
          loci = camera_locations[i, ] - a[(m+uani), ]
          dist2 = loci[1]^2 + loci[2]^2
          DFR = DFR + exp(-dist2/(2*sigma^2))
        }
        nll = nll + lambda0*duration*DFR
      }
    }
    
    
    # AC point process
    nll = nll - m*log(psi)
    
    
    if(log) return(-1*nll)
    else return(exp(-1*nll))
  })


assign('dSESCR_DA', dSESCR_DA, .GlobalEnv)


mcSESCR_DA <- nimbleCode({
  # Priors
  for(i in 1:M){
    a[i, 1] ~ dunif(bounds[1], bounds[2])
    a[i, 2] ~ dunif(bounds[3], bounds[4])
  }
  
  
  psi ~ dbeta(1,1)
  lambda0 ~ dunif(1e-20, 1)
  beta ~ dunif(1e-10, 1000)
  Dratio ~ dunif(1e-10, 1)
  sigma ~ dunif(log(0.5*sqrt(area)), 10)
  
  for(i in 1:(M-m)){
    z[i] ~ dbern(psi)
  }
  
  
  # Likelihood
  Sigma <- sqrt(1/(2*exp(sigma)))
  Beta <- 1/beta
  events[,] ~ dSESCR_DA(N = N, m = m, Nobs = Nobs, area = area,
                         camera_locations = camera_locations[,],
                         z = z[], M = M, psi = psi, lambda0 = lambda0, beta = Beta,
                         Dratio = Dratio, a = a[,], sigma = Sigma)
  
  Nind <- sum(z[1:(M-m)]) + m
  
})
