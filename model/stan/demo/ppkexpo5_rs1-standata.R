# Create Stan data
#
# This function must return the list that will be passed to `data` argument
#   of `cmdstanr::sample()`
#
# The `.dir` argument represents the absolute path to the directory containing
#   this file. This is useful for building file paths to the input files you will
#   load. Note: you _don't_ need to pass anything to this argument, you only use
#   it within the function. `bbr` will pass in the correct path when it calls
#   `make_standata()` under the hood.
make_standata <- function(.dir) {
  # read in any input data
  data1 <- readr::read_csv(file.path(.dir, '..', '..', '..', "data", "derived", 
                                     "pk.csv"), na = ".") %>%
    mutate(AMT = if_else(is.na(AMT), 0, AMT))

  ## Write data1 to file in model directory 
  saveRDS(data1, file = file.path(.dir, "data1.RData"))
  
  grainsize <- readRDS(file.path(.dir, "grainsize.RDS"))
  
  nt <- nrow(data1)
  start <- (1:nt)[!duplicated(data1$ID)]
  end <- c(start[-1] - 1, nt)
  nId <- length(unique(data1$ID))
  
  stan_data <- with(data1,
                    list(nId = nId,
                         nt = nt,
                         amt = 1000 * AMT,
                         cmt = CMT,
                         evid = EVID,
                         time = TIME,
                         start = start,
                         end = end,
                         weight = WT[start],
                         EGFR = EGFR[start],
                         age = AGE[start],
                         albumin = ALB[start],
                         ## cObs = -1 if missing or not an observation record
                         cObs = if_else(BLQ == 0 & EVID == 0, DV, -1),
                         blq = BLQ,
                         LOQ = 10,
                         grainsize = grainsize
                    ))
  return(stan_data)
}

