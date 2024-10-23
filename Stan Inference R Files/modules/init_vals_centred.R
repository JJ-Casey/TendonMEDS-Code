library(trialr)

get_init_vals.me.nophi.rescaled.st.previous_run <-
  function(num_chains, num_data) {
    inits <- list()
    
    prev_samples <- read.csv("./modules/inits.csv")
    
    nu_samps <- prev_samples$b_Nu_Intercept
    eta_samps <- prev_samples$b_eta_Intercept
    tau_samps <- prev_samples$b_tau_Intercept
    rho_samps <- prev_samples$b_rho_Intercept
    
    for (i in 1:num_chains) {
      inits[[i]] <- list(
        pop_mean = c(
          sample(nu_samps, size = 1),
          sample(eta_samps, size = 1),
          sample(tau_samps, size = 1),
          sample(rho_samps, size = 1)
        ),
        theta = t(sapply(1:num_data, function(i)
          c(
            sample(nu_samps, size = 1),
            sample(eta_samps, size = 1),
            sample(tau_samps, size = 1),
            sample(rho_samps, size = 1)
          ))),
        pop_sd = c(
          sample(prev_samples$sd_group__Nu_Intercept, size = 1),
          sample(prev_samples$sd_group__eta_Intercept, size = 1),
          sample(prev_samples$sd_group__tau_Intercept, size = 1),
          sample(prev_samples$sd_group__rho_Intercept, size = 1)
        ),
        L = chol(trialr::rlkjcorr(1, 4, 1))
      )
    }
    inits
  }