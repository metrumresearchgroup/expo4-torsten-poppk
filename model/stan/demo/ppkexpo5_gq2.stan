/*
Generate simulations for patients with different EGFR levels
using final PPK model
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
  // Number of simulated individuals
  int nIdSim = 1000;
  
  // Data structure for each individual
  int ntSim = 193;
  int<lower = 1> nRandom = 5, nCmt = 3;
  array[ntSim] real<lower = 0> amtSim = rep_array(0.0, ntSim),
                            timeSim;
  array[ntSim] int<lower = 1> cmtSim = rep_array(1, ntSim);
  array[ntSim] int<lower = 0> evidSim = rep_array(0, ntSim);
  array[ntSim] real<lower = 0> rateSim = rep_array(0.0, ntSim), 
                            iiSim = rep_array(0.0, ntSim);
  array[ntSim] int<lower = 0> addlSim = rep_array(0, ntSim), 
                           ssSim = rep_array(0, ntSim);
  //array[nId, nCmt] real FSim = rep_array(1.0, nIdSim, nCmt), 
  //                      tLagSim = rep_array(0.0, nIdSim, nCmt);
  
  real weightSim = 70,
       ageSim = 35,
       albuminSim = 4.5;
  int nEGFR = 4;
  array[nEGFR] real EGFRSim = {20, 40, 60, 90};
  
  for(i in 1:ntSim){
    timeSim[i] = i - 1;
  }
  amtSim[{1, 25, 49, 73, 97, 121, 145}] = rep_array(25, 7);
  evidSim[{1, 25, 49, 73, 97, 121, 145}] = rep_array(1, 7);
  
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
  
}

model{

}

generated quantities{
  vector[nRandom] thetaHat = log([CLHat, QHat, V2Hat, V3Hat, kaHat]');
  // 5th, 50th & 95% percentiles of simulated concentrations
  array[nEGFR, ntSim] real cPred5, cPred50, cPred95;
  {
    // Code block used so following declared variables are not saved
    // to disk
    
    // Individual parameters for hypothetical new individuals
    matrix[nRandom, nIdSim] zSim;
    matrix[nIdSim, nRandom] thetaSim;
    vector[nIdSim] CLSim, QSim, V2Sim, V3Sim, kaSim;
    
    row_vector[ntSim] cHatSim; // predicted concentration in new individuals
    matrix[nCmt, ntSim] xSim; // mass in all compartments in new individuals
    matrix[nIdSim, ntSim] cPredSim;
    
    // Simulations for posterior predictive checks
    for(j in 1:nIdSim){
      zSim[,j] = to_vector(normal_rng(rep_vector(0, nRandom), 1));
    }
    thetaSim = (rep_matrix(thetaHat, nIdSim) + diag_pre_multiply(omega, L_corr * zSim))';
    QSim = exp(thetaSim[,2] + 0.75 * log(weightSim / 70));
    V2Sim = exp(thetaSim[,3] + log(weightSim / 70));
    V3Sim = exp(thetaSim[,4] + log(weightSim / 70));
    kaSim = exp(thetaSim[,5]);
    
    for(i in 1:nEGFR){
      CLSim = exp(thetaSim[,1] + 0.75 * log(weightSim / 70) +
      EGFR_CL * log(EGFRSim[i] / 90) +
      age_CL * log(ageSim / 35) +
      albumin_CL * log(albuminSim / 4.5));
      
      for(j in 1:nIdSim){
        xSim = pmx_solve_twocpt(timeSim,
                                amtSim,
                                rateSim,
                                iiSim,
                                evidSim,
                                cmtSim,
                                addlSim,
                                ssSim,
                                {CLSim[j], QSim[j], V2Sim[j], 
                                V3Sim[j], kaSim[j]});
        
        cHatSim = xSim[2,] / V2Sim[j];
        cPredSim[j, 1] = 0.0;
        cPredSim[j, 2:ntSim] = to_row_vector(lognormal_rng(log(cHatSim[2:ntSim]), 
                                          sigma)); 
      }
      for(k in 1:ntSim){
        cPred5[i, k] = quantile(cPredSim[,k], 0.05);
        cPred50[i, k] = quantile(cPredSim[,k], 0.5);
        cPred95[i, k] = quantile(cPredSim[,k], 0.95);
      }
    }
  }
}

