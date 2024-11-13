get_model_priors.me.nophi.rescaled.st <- function () {
  NuPrior <- prior(std_normal(), nlpar = "Nu")
  etaPrior <- prior(std_normal(), nlpar = "eta")
  tauPrior <- prior(std_normal(), nlpar = "tau")
  rhoPrior <- prior(std_normal(), nlpar = "rho")
  nuSDPrior <- prior(student_t(3, 0, 1),
                     class = "sd",
                     nlpar = "Nu",
                     lb = 0)
  etaSDPrior <- prior(student_t(3, 0, 1),
                      class = "sd",
                      nlpar = "eta",
                      lb = 0)
  tauSDPrior <- prior(student_t(3, 0, 1),
                      class = "sd",
                      nlpar = "tau",
                      lb = 0)
  rhoSDPrior <- prior(student_t(3, 0, 1),
                      class = "sd",
                      nlpar = "rho",
                      lb = 0)
  corr_prior <- prior(lkj(1),
                      nlpar= "L_1")
  
  NuPrior + etaPrior + tauPrior + rhoPrior + nuSDPrior + etaSDPrior + tauSDPrior + rhoSDPrior
}
