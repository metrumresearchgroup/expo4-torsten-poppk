/*
Final PPK model
  - Linear 2 compartment
  - Non-centered parameterization
  - lognormal residual variation
  - Allometric scaling
  - Additional covariates: EGFR, age, albumin
  
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
  vector<lower = 0>[nId] weight, EGFR, age, albumin;
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
	real EGFR_CL, age_CL, albumin_CL;
  
  real<lower = 0> sigma;
  vector<lower = 0>[nRandom] omega;
//  corr_matrix[nRandom] rho;
//  array[nId] vector[nRandom] theta;
  cholesky_factor_corr[nRandom] L_corr;
  matrix[nRandom, nId] z;
}

transformed parameters{
  vector[nRandom] thetaHat = log([CLHat, QHat, V2Hat, V3Hat, kaHat]');
  // Individual parameters
  matrix[nId, nRandom] theta = (rep_matrix(thetaHat, nId) + diag_pre_multiply(omega, L_corr * z))';
  vector<lower = 0>[nId] CL = exp(theta[,1] + 0.75 * log(weight / 70) +
                              EGFR_CL * log(EGFR / 90) +
                              age_CL * log(age / 35) +
                              albumin_CL * log(albumin / 4.5)),
                         Q = exp(theta[,2] + 0.75 * log(weight / 70)),
                         V2 = exp(theta[,3] + log(weight / 70)),
                         V3 = exp(theta[,4] + log(weight / 70)),
                         ka = exp(theta[,5]);
  corr_matrix[nRandom] rho =  L_corr * L_corr';

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
  EGFR_CL ~ normal(0, 2);
  age_CL ~ normal(0, 2);
  albumin_CL ~ normal(0, 2);
  
  sigma ~ normal(0, 0.5);
  omega ~ normal(0, 0.5); 
  L_corr ~ lkj_corr_cholesky(2);
  
  // interindividual variability
  to_vector(z) ~ std_normal();

  // likelihood
  cObs ~ lognormal(log(cHat[iObs]), sigma); // observed data likelihood
  target += lognormal_lcdf(LOQ | log(cHat[iBlq]), sigma); // BLQ data likelihood

}

generated quantities{
  // Individual parameters for hypothetical new individuals
  matrix[nRandom, nId] zNew;
  matrix[nId, nRandom] thetaNew;
  vector<lower = 0>[nId] CLNew, QNew, V2New, V3New, kaNew;

  row_vector[nt] cHatNew; // predicted concentration in new individuals
  matrix[nCmt, nt] xNew; // mass in all compartments in new individuals
  array[nt] real cPred, cPredNew;

  vector[nObs + nBlq] log_lik;

  // Simulations for posterior predictive checks
  for(j in 1:nId){
    zNew[,j] = to_vector(normal_rng(rep_vector(0, nRandom), 1));
  }
  thetaNew = (rep_matrix(thetaHat, nId) + diag_pre_multiply(omega, L_corr * zNew))';
  CLNew = exp(thetaNew[,1] + 0.75 * log(weight / 70) +
          EGFR_CL * log(EGFR / 90) +
          age_CL * log(age / 35) +
          albumin_CL * log(albumin / 4.5));
  QNew = exp(thetaNew[,2] + 0.75 * log(weight / 70));
  V2New = exp(thetaNew[,3] + log(weight / 70));
  V3New = exp(thetaNew[,4] + log(weight / 70));
  kaNew = exp(thetaNew[,5]);

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
