/*
Base PPK model
  - Linear 2 compartment
  - Non-centered parameterization
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
  array[nObs] int<lower = 1> idObs; // subject index associated with each cObs
  array[nBlq] int<lower = 1> idBlq; // subject index associated with each BLQ
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
//  corr_matrix[nRandom] rho;
//  array[nId] vector[nRandom] theta;
  cholesky_factor_corr[nRandom] L_corr;
  matrix[nRandom, nId] z;
}

transformed parameters{
  vector[nRandom] thetaHat = log([CLHat, QHat, V2Hat, V3Hat, kaHat]');
  // Individual parameters
  matrix[nId, nRandom] theta = (rep_matrix(thetaHat, nId) + diag_pre_multiply(omega, L_corr * z))';
  vector<lower = 0>[nId] CL = exp(theta[,1]),
                         Q = exp(theta[,2]),
                         V2 = exp(theta[,3]),
                         V3 = exp(theta[,4]),
                         ka = exp(theta[,5]);
  corr_matrix[nRandom] rho =  L_corr * L_corr';
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

}

generated quantities{
  // approximate individual log-likelihood
  vector[nId] log_lik_id;
  {
    int nSim = 1000;
    // simulated individual log-likelihood
    array[nId] vector[nSim] log_lik_id_sim;
    // Individual parameters for hypothetical new individuals
    array[nId] vector[nRandom] thetaNew;
    array[nId] real CLNew, QNew, V2New, V3New, kaNew;
    
    row_vector[nt] cHatNew; // predicted concentration in new individuals
    matrix[nCmt, nt] xNew; // mass in all compartments in new individuals

    for(m in 1:nSim){
      thetaNew = multi_normal_rng(rep_array(thetaHat, nId), Omega);
      CLNew = exp(thetaNew[,1]);
      QNew = exp(thetaNew[,2]);
      V2New = exp(thetaNew[,3]);
      V3New = exp(thetaNew[,4]);
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
      }
      
      // Calculate log-likelihood values to use for LOO CV
      
      // Individual-level log-likelihood
      for(j in 1:nId)
        log_lik_id_sim[j, m] = multi_normal_lpdf(thetaNew[j] | thetaHat, Omega);
      for(i in 1:nObs)
        log_lik_id_sim[idObs[i], m] += lognormal_lpdf(cObs[i] | log(cHatNew[iObs[i]]), sigma);
      for(i in 1:nBlq)
        log_lik_id_sim[idBlq[i], m] += lognormal_lcdf(LOQ | log(cHatNew[iBlq[i]]), sigma);
    }
    for(j in 1:nId){
      log_lik_id[j] = log_sum_exp(log_lik_id_sim[j,]) - log(nSim);
    }
  }
}
