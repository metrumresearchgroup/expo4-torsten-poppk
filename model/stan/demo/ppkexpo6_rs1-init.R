# Create Stan initial values
#
# This function must return something that can be passed to the `init` argument
#   of `cmdstanr::sample()`. There are several options; see `?cmdstanr::sample`
#   for details.
#
# `.data` represents the list returned from `make_standata()` for this model.
#   This is provided in case any of your initial values are dependent on some
#   aspect of the data (e.g. the number of rows).
#
# `.args` represents the list of attached arguments that will be passed through to
#   cmdstanr::sample(). This is provided in case any of your initial values are
#   dependent on any of these arguments (e.g. the number of chains).
#
# Note: you _don't_ need to pass anything to either of these arguments, you only
#   use it within the function. `bbr` will pass in the correct objects when it calls
#   `make_init()` under the hood.
#
make_init <- function(.data, .args) {
  ### create initial estimates
  gen1init <- function(.data){
    nRandom <- 5
    nId = .data$nId
    list(
      CLHat = rlnorm(1, log(5), 0.5), 
      QHat = rlnorm(1, log(5), 0.5), 
      V2Hat = rlnorm(1, log(50), 0.5), 
      V3Hat = rlnorm(1, log(100), 0.5), 
      kaHat = rlnorm(1, log(2), 0.5),
      EGFR_CL = rnorm(1, 0, 1),
      age_CL = rnorm(1, 0, 1),
      albumin_CL = rnorm(1, 0, 1),
      sex_CL = rnorm(1, 0, 1),
      AAG_CL = rnorm(1, 0, 1),
      AST_CL = rnorm(1, 0, 1),
      ALT_CL = rnorm(1, 0, 1),
      lambda_mix = rbeta(7, 5, 5),
      sigma = rlnorm(1, log(0.25), 0.5), 
      omega = rlnorm(nRandom, log(0.25), 0.5),
      L_corr = diag(nRandom),
      z = matrix(rep(0, nRandom * nId), nrow = nRandom))
  }
  inits <- NULL
  for (n in 1:.args$chains)
  {
    inits[[n]] <- gen1init(.data)
  }
  return(inits)
}
