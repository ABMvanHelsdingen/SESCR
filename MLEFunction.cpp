// Semi-complete likelihood for SESCR
// Activity Centers marginalized from likelihood
// Realised not expected population
// Last Updated 26 February 2025

#include <TMB.hpp>
#include <math.h> //for pi
#include <vector>
#include <iostream>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace Eigen;
  using namespace density; // MVNORM
  
  //Data
  DATA_INTEGER(J); // Number of cameras
  DATA_INTEGER(m); // Number of detected individuals
  DATA_INTEGER(K); // Number of observations
  DATA_INTEGER(nrand); // Number of points in mask
  DATA_SCALAR(area); // area of survey area
  DATA_MATRIX(camera_locations); // coordinates of cameras
  DATA_MATRIX(mask); // mask for integration
  
  // Detection Data
  DATA_VECTOR(times); // detection times
  DATA_IVECTOR(cameras); // camera number for each detection
  DATA_IVECTOR(animals); // animal number for each detection
  DATA_IVECTOR(last_cameras); // last camera for that individual
  
  
  
  // Parameters
  PARAMETER(log_N); // log(N-m)
  PARAMETER(log_sigma); // log of home range size
  PARAMETER(logit_d); // logit of d/sigma
  PARAMETER(log_beta); // log of beta
  PARAMETER(log_lambda0); // log of lambda0
  
  // Convert parameters to standard scales
  Type N = exp(log_N) + Type(m);
  Type sigma = exp(log_sigma);
  Type Dratio = exp(logit_d) / (Type(1.) + exp(logit_d));
  Type beta = exp(log_beta);
  Type lambda0 = exp(log_lambda0);
  
  
  // Indexing scheme:
  // i and k for cameras, j for times/events, ani for individual animals
  
  vector<Type> loci(2);
  Type dist2 = 0;
  
  // Matrix of spatial self-excitement between cameras
  matrix<Type> sse(J,J);
  Type d = Dratio * sigma;
  
  for(int i = 0; i < J; i++){
    for(int k = i; k < J; k++){
      loci = camera_locations.row(i) - camera_locations.row(k);
      dist2 = pow(loci[0],2) + pow(loci[1],2);
      sse(i,k) = exp(-dist2/(2*pow(d,2))) / pow(Dratio,2);
      sse(k,i) = sse(i,k); // Matrix is symmetric
    }
  }
  
  
  
  Type duration = times[K]; // survey duration
  
  
  // A[ani, j] determines the self-exciting effect of animal ani at event j
  // The other variables are used to calculate the sums of the integral for each stream
  
  matrix<Type> A(m,(K+1)); A.setZero();
  vector<Type> decay_marks(m); decay_marks.setZero();
  vector<Type> spike_marks(m); spike_marks.setZero();
  vector<Type> decay_sums(m); decay_sums.setZero();
  vector<Type> spike_sums(m); spike_sums.setZero();
  
  for(int j = 1; j < (K+1); j++){
    int An = animals[(j-1)];
    for(int ani = 0; ani < m; ani++){
      if (ani == An){
        A(ani,j) = exp(-beta * (times[j] - times[(j - 1)])); // Markovian, no dependency on previous events
      } else {
        A(ani,j) = exp(-beta * (times[j] - times[(j - 1)])) * A(ani,(j-1));
      }
    }
  
              
    // Areas of spikes and decays
    // Animal seen at observation j - 1 needs a new (truncated) exponential decay curve to begin
    // If the animal was detected previously, a curve needs to finish and the area be calculated
    // Integral = sum of marks - sum(mark_i * A_i+1) up to i = N -1

    spike_sums[An] -= spike_marks[An] * A(An,(j-1));
    int C = cameras[(j-1)];
    spike_marks[An] = sse.col(C).sum();
    spike_sums[An] += spike_marks[An];
  }
  
            
  // At end of survey, all exponential curves remaining end
  for(int ani = 0; ani < m; ani++){
    spike_sums[ani] = spike_sums[ani] - spike_marks[ani] * A(ani,K);
  }
  
  
  // DF and decay marks depend on the activity center
  vector<Type> ILL(m); // Individual log-likelihoods
  vector<Type> ELK(m); ELK.setZero();// Expected likelihoods
  vector<Type> DF(J); // baseline rate at each camera
  //Type ENDP = 0; //Expected probability of non-detection
  Type IDP = 0; // Integrated Detection probability
  Type DFsum = 0; // sum of all baseline rates
  for(int x = 0; x < nrand; x++){
    ILL.setZero();
    // Baseline rates
    for(int i = 0; i < J; ++i){
      loci = camera_locations.row(i) - mask.row(x);
      dist2 = pow(loci[0],2) + pow(loci[1],2);
      DF[i] = exp(-dist2/(2*pow(sigma,2)));
    }
    DFsum = DF.sum();
    
    // Decay mark
    decay_sums.setZero(); decay_marks.setZero();
    for(int j = 1; j < (K+1); j++){
      int An = animals[(j-1)]; // Indivdual detected last
      // Areas of spikes and decays (spikes done before loop)
      // Animal seen at observation j - 1 needs a new (truncated) exponential decay curve to begin
      // If the animal was detected previously, a curve needs to finish and the area be calculated
      // Integral = sum of marks - sum(mark_i * A_i+1) up to i = N -1
      // If animal has not been detected before, decay_marks = 0
      decay_sums[An] -= decay_marks[An] * A(An,(j-1));
      decay_marks[An] = DFsum;
      decay_sums[An] += decay_marks[An];
    }
    
    // At end of survey, all exponential curves remaining end
    // The very last event does not create a spike and decay
    for(int ani = 0; ani < m; ani++){
      decay_sums[ani] -= decay_marks[ani] * A(ani, K);
    }
    
    
    // log(lambda)
    for(int j = 0; j < K; j++){
      int An = animals[j];
      ILL[An] -= log(lambda0);
      int C = cameras[j];
      // Calculate lambda at the time of detection
      if (last_cameras[j] == (-1)){
        ILL[An] -= log(DF[C]);
      } else{
        ILL[An] -= log(((1 - A(An,j))*DF[C]) + (A(An,j) * sse(last_cameras[j], C)));
      }
    }
    
    //Integrals of Lambdas
    for(int ani = 0; ani < m; ani++){
      ILL[ani] += DFsum * lambda0 * duration; //baseline
      ILL[ani] += (spike_sums[ani] - decay_sums[ani]) * (lambda0 / beta); // spikes and decreases
      ELK[ani] += exp(-ILL[ani]);
    }
    
    //Computations for likelihood of non-observed animals
    IDP += (Type(1.) - exp(-duration * DFsum * lambda0));
  }
  
  ELK *= (area / nrand);
  Type nll = -1 * sum(log(ELK));
  
  // Contribution to likelihood from undetected animals
  IDP /= nrand; // Average probability of detection
  
  
  // Calculate N choose m
  Type Mp1 = N + 1;
  Type mp1 = m + 1;
  Type Mlmp1 = N - m + 1;
  Type lMCm = lgamma(Mp1) - lgamma(mp1) - lgamma(Mlmp1);
  nll -= lMCm;
  nll -= (N - m) * log(Type(1.) - IDP);
  
  ADREPORT(N);
  ADREPORT(beta);
  ADREPORT(lambda0);
  ADREPORT(sigma);
  ADREPORT(Dratio);
  
  return nll;
}
