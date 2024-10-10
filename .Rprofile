local({
  
  r_version <- "4.3"
  
  correct_r <- grepl(
    paste0("R version ", r_version),
    R.version[["version.string"]], 
    fixed = TRUE
  )
  
  if (!correct_r) {
    stop(paste0("This project only works with R ", r_version), call. = FALSE)
  }
  
  options(
    
    # set up Box project folder(s)
    mrgboxr.working_dirs = list(
      "Put your project name here (abc-123)" = "Put your project's Box ID here (12 digit number)",
      "report" = "157358260454"
    ),
    
    # set some bbr opinionated defaults, these won't impact users who don't use bbr
    'bbr.bbi_exe_path' = file.path(getwd(), "bin", "bbi"),
    
    # do not check that the project library is synced with lockfile on load
    renv.config.synchronized.check = FALSE,
    
    # do not print renv startup messages - we have our own
    renv.config.startup.quiet = TRUE
  )
  
  source("renv/activate.R")
  
  libpath_1 <- .libPaths()[1]
  lib_pkgs <- list.files(libpath_1)
  
  # This is the case where only renv itself is installed
  if (identical(lib_pkgs, "renv")) {
    
    warning(
      "No packages installed, run 'pkgr install' in the terminal to install.",
      call. = FALSE
    )
    
  } else {
    
    # Attempt to set options(repos) from pkgr.yml "Repos" section.
    # This is primarily to ensure that creating an renv.lock file
    # with renv::snapshot() will point to the same repos as pkgr.
    pkgr_list <- tryCatch(
      suppressWarnings(yaml::read_yaml("pkgr.yml")),
      error = identity
    )
    
    if (inherits(pkgr_list, "error")) {
      
      warning(
        "Extracting repos from pkgr.yml failed: ", conditionMessage(pkgr_list),
        call. = FALSE
      )
      
    } else {
      
      # if Repos is empty, repos_to_set will be "NULL"
      repos_to_set <- unlist(pkgr_list$Repos)
      
      if (!is.null(repos_to_set)) {
        
        # a value must be set to CRAN or R will complain
        # point to MPN if it's there, else grab the first one (it's arbitrary)
        if (!("CRAN" %in% names(repos_to_set))) {
          repos_to_set[["CRAN"]] <- 
            if ("MPN" %in% names(repos_to_set)) {
              repos_to_set[["MPN"]]
            } else {
              repos_to_set[[1]]
            }
        }
        
        options(repos = repos_to_set)
        
      } else {
        
        warning(
          "No repos found in pkgr.yml",
          call. = FALSE
        )
        
      }
    }
  }
  
  if (interactive()) {
    
    # project info ------------------------------------------------------------
    message("Project:   ", basename(getwd()))                                                                         
    message("R version: ", r_version)                                                                                 
    message("Library:   ", libpath_1) 
  }
  
})
