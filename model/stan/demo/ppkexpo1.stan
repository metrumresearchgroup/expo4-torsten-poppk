/*
Base PPK model
  - Linear 2 compartment
  - Centered parameterization
  - lognormal residual variation
  
Based on NONMEM PopPK FOCE example project
*/

data{
  int<lower = 1> nId; // number of individuals
  int<lower = 1> nt; // number of events (rows in data set)
  int<lower = 1> nObs; // number of PK observations
  array[nObs] int<lower = 1> iObs; // event indices for PK observations
  int<lower = 1> nBlq; // number of BLQ observations
  array[nBlq] int<lower = 1> iBlq; // event indices for BLQ observations
  array[nt] real<lower = 0> amt, time;
  array[nt] int<lower = 1> cmt;
  array[nt] int<lower = 0> evid;
  array[nId] int<lower = 1> start, end; // indices of first & last observations
  vector<lower = 0>[nObs] cObs;
  real<lower = 0> LOQ;
}

transformed data{
  int<lower = 1> nRandom = 5, nCmt = 3;
  array[nt] real<lower = 0> rate = rep_array(0.0, nt), 
                            ii = rep_array(0.0, nt);
  array[nt] int<lower = 0> addl = rep_array(0, nt), 
                           ss = rep_array(0, nt);
  array[nId, nCmt] real F = rep_array(1.0, nId, nCmt), 
                        tLag = rep_array(0.0, nId, nCmt);
  }

parameters{
  real<lower = 0> CLHat, QHat, V2Hat, V3Hat;
// To constrain kaHat > lambda1 uncomment the following
  real<lower = (CLHat / V2Hat + QHat / V2Hat + QHat / V3Hat +
		sqrt((CLHat / V2Hat + QHat / V2Hat + QHat / V3Hat)^2 -
		     4 * CLHat / V2Hat * QHat / V3Hat)) / 2> kaHat; // ka > lambda_1
  
  real<lower = 0> sigma;
  vector<lower = 0>[nRandom] omega;
  corr_matrix[nRandom] rho;
  array[nId] vector[nRandom] theta;
}

transformed parameters{
  vector[nRandom] thetaHat = log([CLHat, QHat, V2Hat, V3Hat, kaHat]');
  // Individual parameters
  array[nId] real<lower = 0> CL = exp(theta[,1]),
                             Q = exp(theta[,2]),
                             V2 = exp(theta[,3]),
                             V3 = exp(theta[,4]),
                             ka = exp(theta[,5]);

  // Covariance matrix
  cov_matrix[nRandom] Omega = quad_form_diag(rho, omega);

  row_vector[nt] cHat; // predicted concentration
  matrix[nCmt, nt] x; // mass in all compartments
  
  for(j in 1:nId){
    x[, start[j]:end[j]] = pmx_solve_twocpt(time[start[j]:end[j]],
                                            amt[start[j]:end[j]],
                                            rate[start[j]:end[j]],
                                            ii[start[j]:end[j]],
                                            evid[start[j]:end[j]],
                                            cmt[start[j]:end[j]],
                                            addl[start[j]:end[j]],
                                            ss[start[j]:end[j]],
                                            {CL[j], Q[j], V2[j], V3[j], ka[j]});

    cHat[start[j]:end[j]] = x[2, start[j]:end[j]] / V2[j];
  }
}

model{
  // priors
  CLHat ~ normal(0, 10); 
  QHat ~ normal(0, 10);
  V2Hat ~ normal(0, 50);
  V3Hat ~ normal(0, 100);
  kaHat ~ normal(0, 3);
  sigma ~ normal(0, 0.5);
  omega ~ normal(0, 0.5); 
  rho ~ lkj_corr(2);
  
  // interindividual variability
  theta ~ multi_normal(thetaHat, Omega);

  // likelihood
  cObs ~ lognormal(log(cHat[iObs]), sigma); // observed data likelihood
  target += lognormal_lcdf(LOQ | log(cHat[iBlq]), sigma); // BLQ data likelihood

}

generated quantities{
  // Individual parameters for hypothetical new individuals
  array[nId] vector[nRandom] thetaNew = 
    multi_normal_rng(rep_array(thetaHat, nId), Omega);
  array[nId] real<lower = 0> CLNew = exp(thetaNew[,1]),
                             QNew = exp(thetaNew[,2]),
                             V2New = exp(thetaNew[,3]),
                             V3New = exp(thetaNew[,4]),
                             kaNew = exp(thetaNew[,5]);

  row_vector[nt] cHatNew; // predicted concentration in new individuals
  matrix[nCmt, nt] xNew; // mass in all compartments in new individuals
  array[nt] real cPred, cPredNew;

  vector[nObs + nBlq] log_lik;

  // Simulations for posterior predictive checks
  for(j in 1:nId){
    xNew[, start[j]:end[j]] = pmx_solve_twocpt(time[start[j]:end[j]],
                                               amt[start[j]:end[j]],
                                               rate[start[j]:end[j]],
                                               ii[start[j]:end[j]],
                                               evid[start[j]:end[j]],
                                               cmt[start[j]:end[j]],
                                               addl[start[j]:end[j]],
                                               ss[start[j]:end[j]],
                                               {CLNew[j], QNew[j], V2New[j], 
                                               V3New[j], kaNew[j]});

    cHatNew[start[j]:end[j]] = xNew[2, start[j]:end[j]] / V2New[j];
    
    cPred[start[j]] = 0.0;
    cPredNew[start[j]] = 0.0;
    // new observations in existing individuals
    cPred[(start[j] + 1):end[j]] = lognormal_rng(log(cHat[(start[j] + 1):end[j]]), 
      sigma);
    // new observations in new individuals
    cPredNew[(start[j] + 1):end[j]] = lognormal_rng(log(cHatNew[(start[j] + 1):end[j]]), 
      sigma); 
  }

  // Calculate log-likelihood values to use for LOO CV
  
  // Observation-level log-likelihood
  for(i in 1:nObs)
    log_lik[i] = lognormal_lpdf(cObs[i] | log(cHat[iObs[i]]), sigma);
  for(i in 1:nBlq)
    log_lik[i + nObs] = lognormal_lcdf(LOQ | log(cHat[iBlq[i]]), sigma);
    
}
