/*
Final PPK model
  - Linear 2 compartment
  - Non-centered parameterization
  - lognormal residual variation
  - Allometric scaling
  - Additional covariates: EGFR, age, albumin
  - reduce_sum
  
Based on NONMEM PopPK FOCE example project
*/

functions{
  real partial_sum(array[] int id, 
  int id_start, int id_end,
  vector cObs,
  array[] int blq,
  real LOQ,
  array[] int start, array[] int end,
  int nCmt,
  array[] real time,
  array[] real amt,
  array[] real rate,
  array[] real ii,
  array[] int evid,
  array[] int cmt,
  array[] int addl,
  array[] int ss,
  vector CL,
  vector Q,
  vector V2,
  vector V3,
  vector ka,
  real sigma
  ){
    int nt = end[id_end] - start[id_start] + 1;
    row_vector[nt] cHat; // predicted concentration
    matrix[nCmt, nt] x; // mass in all compartments
    real result;

    for(j in id_start:id_end){
      x[, (start[j] - start[id_start] + 1):(end[j] - start[id_start] + 1)] = 
        pmx_solve_twocpt(time[start[j]:end[j]],
        amt[start[j]:end[j]],
        rate[start[j]:end[j]],
        ii[start[j]:end[j]],
        evid[start[j]:end[j]],
        cmt[start[j]:end[j]],
        addl[start[j]:end[j]],
        ss[start[j]:end[j]],
        {CL[j], Q[j], V2[j], V3[j], ka[j]});
      
      cHat[(start[j] - start[id_start] + 1):(end[j] - start[id_start] + 1)] = 
        x[2, (start[j] - start[id_start] + 1):(end[j] - start[id_start] + 1)] / V2[j];
    }
    
    // likelihood
    result = 0;
    for(i in start[id_start]:end[id_end]){
      if(cObs[i] > 0)
        result += lognormal_lpdf(cObs[i] | log(cHat[i - start[id_start] + 1]), sigma);
      else if(blq[i] == 1)
        result += lognormal_lcdf(LOQ | log(cHat[i - start[id_start] + 1]), sigma);
    }
    
    return result;
  }
}

data{
  int<lower = 1> nId; // number of individuals
  int<lower = 1> nt; // number of events (rows in data set)
  array[nt] real<lower = 0> amt, time;
  array[nt] int<lower = 1> cmt;
  array[nt] int<lower = 0> evid;
  array[nId] int<lower = 1> start, end; // indices of first & last events
  vector<lower = 0>[nId] weight, EGFR, age, albumin;
  vector<lower = -1>[nt] cObs;
  array[nt] int<lower = 0, upper = 1> blq;
  real<lower = 0> LOQ;
  int grainsize;
}

transformed data{
  int<lower = 1> nRandom = 5, nCmt = 3;
  array[nt] real<lower = 0> rate = rep_array(0.0, nt), 
                            ii = rep_array(0.0, nt);
  array[nt] int<lower = 0> addl = rep_array(0, nt), 
                           ss = rep_array(0, nt);
  array[nId, nCmt] real F = rep_array(1.0, nId, nCmt), 
                        tLag = rep_array(0.0, nId, nCmt);
  array[nId] int id;

  for(i in 1:nId) id[i] = i;
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
  to_vector(z) ~ normal(0, 1);

  // likelihood
  target += reduce_sum_static(partial_sum, id, grainsize, cObs, blq, LOQ, start, end,
            nCmt, time, amt, rate, ii, evid, cmt, addl, ss, CL, Q, V2, V3, ka, sigma);

}

generated quantities{

}
