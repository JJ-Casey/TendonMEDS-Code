library(tidyverse)

read_and_group <-
  function(group,
           type,
           filename,
           screen_data,
           fidels_loc,
           threshold = 0.1) {
    data <- read.csv(paste0(screen_data, filename))
    
    fidel.file.name <- paste0(fidels_loc, group, "_" , type, ".csv")
    fidels <- read.csv(fidel.file.name)
    fid.len <- nrow(fidels)
    
    fidels.mask <- fidels < threshold
    first_match <- match(TRUE, fidels.mask)[1]
    
    if (is.na(first_match)) {
      first_match <- nrow(fidels)
    }
    #first_match <- min(num_data_points, first_match)
    
    # if (first_match < 0.2 * nrow(fidels)) {
    #   return(data.frame())
    # }
    
    fidels <- fidels[1:first_match,]
    data <- data[1:length(fidels),]
    errs <- 0.15 * sqrt(1 / fidels)
    
    names(data) <- c('stretch', 'stress')
    
    stretchm1 <- data[, 1] ^ -1
    stretch2 <- data[, 1] ^ 2
    stretch2m1 <- stretch2 ^ -1
    stretchmstretch2m1 <- data[, 1] - stretch2 ^ -1
    stretchLogStretch <- data[, 1] * log(data[, 1])
    
    data %>% mutate(
      group = group,
      type = rep(type, nrow(data)),
      measureNo = 1:nrow(.),
      errs = errs,
      stretchm1 = stretchm1,
      stretch2 = stretch2,
      stretchmstretch2m1 = stretchmstretch2m1,
      stretchLogStretch = stretchLogStretch
    ) %>%
      select(group, measureNo, errs, everything()) #%>% filter(row_number() %% 3 == 1)
  }

load_tendon_data <-
  function(screen_data = "E:/Google Drive/Uni!/PhD/Stan/screen data/csv/",
           fidels_loc = "E:/Google Drive/Uni!/PhD/Stan/fidels/",
           threshold = 0.1) {
    ## Glob for the data files
    files.sdft <-
      list.files(screen_data, pattern = "/* sdft unfiltered\\.csv")
    files.cdet <-
      list.files(screen_data, pattern = "/* cdet unfiltered\\.csv")
    
    tendon.types <-
      c(rep("sdft", length(files.sdft)), rep("cdet", length(files.cdet)))
    
    ## Grab the horse names (HXX)
    files <- c(files.sdft, files.cdet)
    subjects <-
      files %>% substr(., start = 1, stop = 3) %>% tolower()
    
    ## Load all data
    files.paths <-
      data.frame(group = subjects,
                 type = tendon.types,
                 filename = files)
    
    ## Load and sort the data
    data <-
      files.paths %>% rowwise() %>% do(.,
                                       read_and_group(
                                         .$group,
                                         .$type,
                                         .$filename,
                                         screen_data,
                                         fidels_loc,
                                         threshold
                                       )) %>% ungroup()
    
    ## Return the data
    data
  }

gen_synthetic_data.me <-
  function(mu.pop = c(1, 1, 0, 0, 0),
           sigma.pop = diag(c(0.05, 0.05, 0.1, 0.1, 0.1), nrow = 5, ncol = 5)) {
    source("./modules/tendonModel_nophi_scaled.R")
    
    num_exps <- 2
    num_obs_per_exp <- 30
    
    sigma.chol <- chol(sigma.pop)
    
    gen_obs_data <- function(index, theta_val = NULL) {
      group <- paste0("group", index)
      
      if (is.null(theta_val)) {
        theta <- rnorm(5, mu.pop, diag(sigma.chol))
      } else {
        group <- "test"
        theta <- theta_val
      }
      stretch <-
        seq(1.0, 1.1 + runif(1, -0.05, 0.05), length.out = num_obs_per_exp)
      errs <- 0.15# * sqrt(1 / fidels)
      
      stretchm1 <- stretch ^ -1
      stretch2 <- stretch ^ 2
      stretch2m1 <- stretch ^ -1
      stretchLogStretch <- stretch * log(stretch)
      
      stress <-
        tendonModel(
          stretch,
          stretchm1,
          stretch2,
          stretch2m1,
          stretchLogStretch,
          theta[1],
          theta[2],
          theta[3],
          theta[4],
          theta[5]
        ) + rnorm(num_obs_per_exp, sd = errs)
      
      data <- data.frame(stretch = stretch, stress = stress)
      
      data %>% mutate(
        group = group,
        type = rep("sdft", nrow(data)),
        measureNo = 1:nrow(.),
        errs = errs,
        stretchm1 = stretchm1,
        stretch2 = stretch2,
        stretch2m1 = stretch2m1,
        stretchLogStretch = stretchLogStretch,
        nu = theta[1],
        eta = theta[2],
        tau = theta[3],
        kappa = theta[4],
        rho = theta[5],
        phimu = exp(2.50531765 + 0.38510136 * theta[1]),
        phiE = exp(6.10303632 + 0.64387023 * theta[2]),
        a = exp(-3.80045123 + 0.64387023 * theta[3]) + 1,
        c = exp(-3.59771868 + 0.7310165 * theta[4]) + exp(-3.80045123 + 0.64387023 * theta[3]) + 1,
        b = exp(-3.59771868 + 0.7310165 * theta[5]) + exp(-3.59771868 + 0.7310165 * theta[4]) + exp(-3.80045123 + 0.64387023 * theta[3]) + 1
      ) %>%
        select(group, measureNo, errs, everything())
    }
    dfs <- lapply(seq.int(num_exps), gen_obs_data)
    dfs <- list(dfs, gen_obs_data(0, mu.pop))
    bind_rows(dfs)
  }

gen_synthetic_data.me.st <-
  function(mu.pop = c(1, 1, 0, 0),
           sigma.pop = diag(c(0.05, 0.05, 0.1, 0.1), nrow = 4, ncol = 4)) {
    source("./modules/tendonModel_nophi_scaled_st.R")
    
    num_exps <- 10
    num_obs_per_exp <- 30
    
    sigma.chol <- chol(sigma.pop)
    
    gen_obs_data <- function(index, theta_val = NULL) {
      group <- paste0("group", index)
      
      if (is.null(theta_val)) {
        theta <- rnorm(4, mu.pop, diag(sigma.chol))
      } else {
        group <- "test"
        theta <- theta_val
      }
      stretch <-
        seq(1.0, 1.1 + runif(1, -0.05, 0.05), length.out = num_obs_per_exp)
      errs <- 0.15# * sqrt(1 / fidels)
      
      stretchm1 <- stretch ^ -1
      stretch2 <- stretch ^ 2
      stretch2m1 <- stretch ^ -1
      stretchLogStretch <- stretch * log(stretch)
      
      stress <-
        tendonModel(
          stretch,
          stretchm1,
          stretch2,
          stretch2m1,
          stretchLogStretch,
          theta[1],
          theta[2],
          theta[3],
          theta[4]
        ) + rnorm(num_obs_per_exp, sd = errs)
      
      data <- data.frame(stretch = stretch, stress = stress)
      
      data %>% mutate(
        group = group,
        type = rep("sdft", nrow(data)),
        measureNo = 1:nrow(.),
        errs = errs,
        stretchm1 = stretchm1,
        stretch2 = stretch2,
        stretch2m1 = stretch2m1,
        stretchLogStretch = stretchLogStretch,
        nu = theta[1],
        eta = theta[2],
        tau = theta[3],
        rho = theta[4],
        phimu = exp(2.50531765 + 0.38510136 * theta[1]),
        phiE = exp(6.10303632 + 0.64387023 * theta[2]),
        a = exp(-3.80045123 + 0.64387023 * theta[3]) + 1,
        b = exp(-3.59771868 + 0.7310165 * theta[4]) + exp(-3.80045123 + 0.64387023 * theta[3]) + 1
      ) %>%
        select(group, measureNo, errs, everything())
    }
    dfs <- lapply(seq.int(num_exps), gen_obs_data)
    dfs <- list(dfs, gen_obs_data(0, mu.pop))
    bind_rows(dfs)
  }