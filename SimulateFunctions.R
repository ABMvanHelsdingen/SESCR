# Functions for simulating SCR, OU and SESCR detections
# Chen (2016a): https://www.math.fsu.edu/~ychen/research/Thinning%20algorithm.pdf
# Chen (2016b): https://www.math.fsu.edu/~ychen/research/multiHawkes.pdf

library(IHSEP)
library(mvtnorm)
library(spatstat.geom)

### SESCR ###

SESCROgata <- function(mu, alpha, beta, d, stream, mu0, N, t){
  # Custom Ogata algorithm for SESCR model
  # following Chen (2016a)
  # using the definition of A in Laub (2014) and stelfi
  
  # mu: baseline rate for each stream
  # alpha: strength of self-excitement, scalar between 0 and 1
  # beta: temporal decay of self-excitement
  # d: self-excitement over space (NxN)
  # mu0: volume of baseline rate across all streams
  # N: number of cameras and streams
  
  # check arguments are correct size
  if(length(mu) != N || length(d) != N*N){
    stop("mu and d must be length N and N^2 respectively")
  }
  if(length(alpha) != 1 || length(beta) != 1 || length(mu0) != 1){
    stop("alpha, beta, mu0 must be a scalar")
  }
  
  # following Chen (2016b)
  # using the definition of A in Laub (2014) and stelfi
  
  # calculate maximum intensity of each stream
  max_lambda = numeric(N)
  for(i in 1:N){
    max_lambda[i] = mu[i] + (alpha * mu0 * d[i,i])
  }
  lambda_bar = sum(max_lambda)
  times = numeric(0)
  streams = numeric(0)
  s = 0; n = 0; A = 0
  lambda_s = mu
  last_event = 0 # stream of last event
  while (s < t){
    u <- runif(1)
    w <- -1*log(u)/lambda_bar
    s <- s + w
    D <- runif(1)
    if (n > 0){
      A <- exp(-beta*(s-max(times))) # Markovian
      # lambda_s is intensity at candidate point s
      for(i in 1:N){
        lambda_s[i] = ((1 - alpha * A) * mu[i]) + (alpha * A * mu0 * d[i, last_event])
      }
    }
    if (D*lambda_bar <= sum(lambda_s)){
      k = 1
      while(D*lambda_bar > sum(lambda_s[1:k])){
        k = k + 1
      }
      n <- n + 1
      last_event <- k
      times <- append(times, s)
      streams <- append(streams, k)
    }
  }
  # Check if last event > t, and return output
  if (n == 0 || (n == 1 && max(times) > t)){
    return(list(times=numeric(0), streams=numeric(0)))
  } else if (max(times)<t){
    return(list(times=times, streams=streams))
  } else {
    return(list(times=times[1:(n-1)], streams=streams[1:(n-1)]))
  }
}

# Simulate SESCR model (with exponential kernels)

simSESCR <- function(N, mu0, beta, d, sigma, camera_locations, t,
                       s = NULL){
  
  # N: true population
  # mu0: baseline background detection rate for the MHP
  # beta: rate of temporal decay in self-excitement
  # d: covariance for SE between cameras (2x2)
  # sigma: covariance for activity ranges (2x2)
  # camera_locations: locations of cameras (Jx2)
  # t: latest possible detection time
  # s: activity centers (Nx2)
  
  # Simulate activity centers (unless provided) 
  # Study area is a unit square from (0,0) to (1,1)
  if (length(s) == 0){
    s = matrix(runif(N*2,0,1), nrow=N, ncol=2)
  }
  
  # check arguments
  if(ncol(camera_locations) != 2){
    stop("camera_locations must be a Jx2 matrix")
  }
  if(nrow(sigma) != 2 || ncol(sigma) != 2){
    stop("sigma must be a 2x2 matrix")
  }
  if(sigma[1,1] != sigma[2,2] || sigma[1,2] != sigma[2,1]){
    stop("sigma must have equal variance in x and y directions
         and zero covariance")
  }
  if(nrow(d) != 2 || ncol(d) != 2){
    stop("d must be a 2x2 matrix")
  }
  if(d[1,1] != d[2,2] || d[1,2] != d[2,1]){
    stop("d must have equal variance in x and y directions
         and zero covariance")
  }
  
  J <- nrow(camera_locations)
  n <- 0 # number of animals with observations
  
  # Calculate spatial self-excitement
  sse <- matrix(0, nrow = J, ncol = J)
  for(x in 1:J){
    for(y in 1:x){
      distance <- camera_locations[x, ] - camera_locations[y, ]
      sse[x,y] <- mnormt::dmnorm(x = distance, varcov = d)
      sse[y,x] <- sse[x,y] # sse is symmetric
    }
  }
  
  # empty vectors for storing simulations
  times <- numeric(0)
  cameras <- numeric(0)
  animals <- numeric(0)
  observed <- numeric(0)
  
  # for each animal:
  for(i in 1:N){
    # Calculate mu for MHP
    mus <- numeric(J)
    for(j in 1:J){
      distance <- camera_locations[j, ] - s[i, ]
      mus[j] <- mnormt::dmnorm(x = distance, varcov = sigma)
    }
    mus <- mu0 * mus
    
    # Generate MHP
    result <- SESCROgata(mu = mus, alpha = 1, beta = beta, 
                         d = sse, mu0 = mu0, N = J, t = t)
    if (length(result$times) > 0){ # If animal is detected
      n <- n + 1
      observed <- append(observed, n)
      times <- append(times, result$times)
      cameras <- append(cameras, result$streams)
      animals <- append(animals, rep(n, length(result$times)))
    }
  }
  times_order = order(times, decreasing = FALSE)
  times <- times[times_order]
  cameras <- cameras[times_order]
  animals <- animals[times_order]
  return(list(times = times, cameras = cameras, ids = animals, m = n, observed = observed))
}

### SCR ###
simSCR <- function(N, mu0, sigma, camera_locations, t,
                   s = NULL){
  # Simulate Continuous-Time SCR
  # Study area is a unit square from (0,0) to (1,1)
  
  # see simSESCR() for parameter definitions
  
  
  # Simulate activity centers (unless provided)
  if (length(s) == 0){
    s = matrix(runif(N*2,0,1), nrow=N, ncol=2)
  }
  
  # check arguments
  if(ncol(camera_locations) != 2){
    stop("camera_locations must be a Jx2 matrix")
  }
  if(nrow(sigma) != 2 || ncol(sigma) != 2){
    stop("sigma must be a 2x2 matrix")
  }
  if(sigma[1,1] != sigma[2,2] || sigma[1,2] != sigma[2,1]){
    stop("sigma must have equal variance in x and y directions
         and zero covariance")
  }
  
  J <- nrow(camera_locations)
  n <- 0 # number of animals with observations
  
  
  # empty vectors for storing simulations
  times <- numeric(0)
  cameras <- numeric(0)
  animals <- numeric(0)
  
  # for each animal:
  for(i in 1:N){
    # Calculate mu for MHP
    mus <- numeric(J)
    for(j in 1:J){
      distance <- camera_locations[j, ] - s[i, ]
      mus[j] <- mnormt::dmnorm(x = distance, varcov = sigma)
    }
    mus <- mu0 * mus
    
    # Generate MHP
    result <- IHSEP::simPois(int = function(x){sum(mus)}, cens = t, int.M = sum(mus))
    
    if (length(result) > 0){ # If animal is detected
      n <- n + 1
      times <- append(times, result)
      animals <- append(animals, rep(n, length(result)))
      cameras <- append(cameras, sample(c(1:J), size = length(result), replace = TRUE,
                                        prob = mus/sum(mus)))
    }
  }
  times_order = order(times, decreasing = FALSE)
  times <- times[times_order]
  cameras <- cameras[times_order]
  animals <- animals[times_order]
  return(list(times = times, cameras = cameras, ids = animals, m = n))
}

### OU Process ###
sim.ou <- function(mu, tau, sigma, n.steps, start = NULL){
  if (is.null(start)){
    start <- rmvnorm(1, mu, sigma^2*diag(2))
  }
  ## Changing Theo's parameterisation.
  b <- -1/tau
  v <- sigma^2
  out <- matrix(0, nrow = n.steps, ncol = 2)
  out[1, ] <- start
  for (i in 2:n.steps){
    out[i, ] <- rmvnorm(1, mu + exp(b)*(out[i - 1, ] - mu),
                        v*(1 - exp(2*b))*diag(2))
    
  }
  out
}

simOU <- function(N, epsilon=0.004, tau=100, sigma=0.2, camera_locations, t,
                  s = NULL, freq = 10, mode = 1){
  
  # Simulate activity centers (unless provided)
  if (length(s) == 0){
    s = matrix(runif(N*2,0,1), nrow=N, ncol=2)
  }
  # empty vectors for storing simulations
  times <- numeric(0)
  cameras <- numeric(0)
  animals <- numeric(0)
  n <- 0
  if (mode == 1){ # Every 1/freq units, animals are observed, jitter applied
    for(i in 1:N){
      locs <- sim.ou(mu = s[i,], tau = tau * freq, sigma = sigma, n.steps = freq * t)
      dists <- crossdist(locs[, 1], locs[, 2], camera_locations[, 1], camera_locations[, 2])
      dets <- which(dists < epsilon, arr.ind = TRUE)
      dets[,1] <- (dets[,1] + runif(nrow(dets),-1,0))/freq
      if (nrow(dets) > 0){
        n <- n + 1
        times <- append(times, dets[,1])
        cameras <- append(cameras, dets[,2])
        animals <- append(animals, rep(n, nrow(dets)))
      }
    }
  } else { # The animals are observed at a rate of freq
    for(i in 1:N){
      locs <- sim.ou2(mu = s[i,], tau = freq * tau, sigma = sigma, n.steps = freq * t)
      dists <- crossdist(locs[, 1], locs[, 2], camera_locations[, 1], camera_locations[, 2])
      dets <- which(dists < epsilon, arr.ind = TRUE)
      if (nrow(dets) > 0){
        n <- n + 1
        times <- append(times, locs[dets[,1], 3]/freq)
        cameras <- append(cameras, dets[,2])
        animals <- append(animals, rep(n, nrow(dets)))
      }
    }
  }
  times_order = order(times, decreasing = FALSE)
  times <- times[times_order]
  cameras <- cameras[times_order]
  animals <- animals[times_order]
  return(list(times = times, cameras = cameras, ids = animals, m = n))
}