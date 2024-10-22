get_model_priors.me.phi <- function () {
  NuPrior <- prior(normal(-1.9560115 ,  2.12854828), nlpar = "Nu")
  etaPrior <- prior(normal(-0.33471533,  1.33570051), nlpar = "eta")
  betaPrior <- prior(normal(0, 0.5), nlpar = "beta")
  tauPrior <- prior(normal(-3.80045123,  0.64387023), nlpar = "tau")
  kappaPrior <-
    prior(normal(-3.59771868,  0.7310165), nlpar = "kappa")
  rhoPrior <- prior(normal(-3.59771868,  0.7310165), nlpar = "rho")
  betaSDPrior <- prior(student_t(3, 0, 0.1),
                       class = "sd",
                       nlpar = "beta",
                       lb = 0)
  tauSDPrior <- prior(student_t(3, 0, 0.1),
                      class = "sd",
                      nlpar = "tau",
                      lb = 0)
  kappaSDPrior <- prior(student_t(3, 0, 0.1),
                        class = "sd",
                        nlpar = "kappa",
                        lb = 0)
  rhoSDPrior <- prior(student_t(3, 0, 0.1),
                      class = "sd",
                      nlpar = "rho",
                      lb = 0)
  
  NuPrior + etaPrior + tauPrior + kappaPrior + rhoPrior + betaPrior + betaSDPrior + tauSDPrior + kappaSDPrior + rhoSDPrior
}

get_model_priors.me.nophi <- function () {
  NuPrior <- prior(normal(-3.00583363 ,  2.12854828), nlpar = "Nu")
  etaPrior <- prior(normal(-0.7654982,  1.33570051), nlpar = "eta")
  tauPrior <- prior(normal(-3.80045123,  0.64387023), nlpar = "tau")
  kappaPrior <-
    prior(normal(-3.59771868,  0.7310165), nlpar = "kappa")
  rhoPrior <- prior(normal(-3.59771868,  0.7310165), nlpar = "rho")
  nuSDPrior <- prior(student_t(3, 0, 0.1),
                     class = "sd",
                     nlpar = "Nu",
                     lb = 0)
  etaSDPrior <- prior(student_t(3, 0, 0.1),
                      class = "sd",
                      nlpar = "eta",
                      lb = 0)
  tauSDPrior <- prior(student_t(3, 0, 0.1),
                      class = "sd",
                      nlpar = "tau",
                      lb = 0)
  kappaSDPrior <- prior(student_t(3, 0, 0.1),
                        class = "sd",
                        nlpar = "kappa",
                        lb = 0)
  rhoSDPrior <- prior(student_t(3, 0, 0.1),
                      class = "sd",
                      nlpar = "rho",
                      lb = 0)
  
  ## Maybe set a prior for L?
  
  NuPrior + etaPrior + tauPrior + kappaPrior + rhoPrior + nuSDPrior + etaSDPrior + tauSDPrior + kappaSDPrior + rhoSDPrior
}


get_model_priors.me.nophi.rescaled <- function () {
  NuPrior <- prior(normal(2.50531765, 0.38510136), nlpar = "Nu")
  etaPrior <- prior(normal(6.10303632 , 0.64387023), nlpar = "eta")
  tauPrior <- prior(normal(-3.80045123 , 0.64387023), nlpar = "tau")
  kappaPrior <- prior(normal(-3.59771868 , 0.7310165), nlpar = "kappa")
  rhoPrior <- prior(normal(-3.59771868 , 0.7310165), nlpar = "rho")
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
  kappaSDPrior <- prior(student_t(3, 0, 1),
                        class = "sd",
                        nlpar = "kappa",
                        lb = 0)
  rhoSDPrior <- prior(student_t(3, 0, 1),
                      class = "sd",
                      nlpar = "rho",
                      lb = 0)
  
  NuPrior + etaPrior + tauPrior + kappaPrior + rhoPrior + nuSDPrior + etaSDPrior + tauSDPrior + kappaSDPrior + rhoSDPrior
}
get_model_priors.me.nophi.rescaled.old <- function () {
  NuPrior <- prior(std_normal(), nlpar = "Nu")
  etaPrior <- prior(std_normal(), nlpar = "eta")
  tauPrior <- prior(std_normal(), nlpar = "tau")
  kappaPrior <- prior(std_normal(), nlpar = "kappa")
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
  kappaSDPrior <- prior(student_t(3, 0, 1),
                        class = "sd",
                        nlpar = "kappa",
                        lb = 0)
  rhoSDPrior <- prior(student_t(3, 0, 1),
                      class = "sd",
                      nlpar = "rho",
                      lb = 0)
  
  NuPrior + etaPrior + tauPrior + kappaPrior + rhoPrior + nuSDPrior + etaSDPrior + tauSDPrior + kappaSDPrior + rhoSDPrior
}
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

get_model_priors.me.nophi.rescaled.st.norm <- function () {
  NuPrior <- prior(std_normal(), nlpar = "Nu")
  etaPrior <- prior(std_normal(), nlpar = "eta")
  tauPrior <- prior(std_normal(), nlpar = "tau")
  rhoPrior <- prior(std_normal(), nlpar = "rho")
  nuSDPrior <- prior(normal(0, 2),
                     class = "sd",
                     nlpar = "Nu",
                     lb = 0)
  etaSDPrior <- prior(normal(0, 2),
                      class = "sd",
                      nlpar = "eta",
                      lb = 0)
  tauSDPrior <- prior(normal(0, 2),
                      class = "sd",
                      nlpar = "tau",
                      lb = 0)
  rhoSDPrior <- prior(normal(0, 2),
                      class = "sd",
                      nlpar = "rho",
                      lb = 0)
  
  NuPrior + etaPrior + tauPrior + rhoPrior + nuSDPrior + etaSDPrior + tauSDPrior + rhoSDPrior
}

get_model_priors.me.nophi.rescaled.st.multinorm <- function () {
  thetaPrior <- prior(multi_normal(prior_mean, prior_covar), nlpar = "theta")
  thetaSDPrior <- prior(student_t(3, 0, 1),
                      class = "sd",
                      nlpar = "theta",
                      lb = 0)
  corr_prior <- prior(lkj(1),
                      class= "cor")
  
  thetaPrior + thetaSDPrior #+ corr_prior
}

get_model_priors.me.phi_constant <- function () {
  NuPrior <- prior(normal(-1.9560115 ,  2.12854828), nlpar = "Nu")
  etaPrior <- prior(normal(-0.33471533,  1.33570051), nlpar = "eta")
  betaPrior <- prior(normal(0, 0.5), nlpar = "beta")
  betaPrior <- prior(constant(0), nlpar = "beta")
  tauPrior <- prior(normal(-3.80045123,  0.64387023), nlpar = "tau")
  kappaPrior <-
    prior(normal(-3.59771868,  0.7310165), nlpar = "kappa")
  rhoPrior <- prior(normal(-3.59771868,  0.7310165), nlpar = "rho")
  tauSDPrior <- prior(student_t(3, 0, 0.1),
                      class = "sd",
                      nlpar = "tau",
                      lb = 0)
  kappaSDPrior <- prior(student_t(3, 0, 0.1),
                        class = "sd",
                        nlpar = "kappa",
                        lb = 0)
  rhoSDPrior <- prior(student_t(3, 0, 0.1),
                      class = "sd",
                      nlpar = "rho",
                      lb = 0)
  
  NuPrior + etaPrior + tauPrior + kappaPrior + rhoPrior + betaPrior + tauSDPrior + kappaSDPrior + rhoSDPrior
}

get_model_priors.nome.phi <- function () {
  NuPrior <- prior(normal(-1.9560115 ,  2.12854828), nlpar = "Nu")
  etaPrior <- prior(normal(-0.33471533,  1.33570051), nlpar = "eta")
  betaPrior <- prior(normal(0, 0.5), nlpar = "beta")
  tauPrior <- prior(normal(-3.80045123,  0.64387023), nlpar = "tau")
  kappaPrior <-
    prior(normal(-3.59771868,  0.7310165), nlpar = "kappa")
  rhoPrior <- prior(normal(-3.59771868,  0.7310165), nlpar = "rho")
  
  NuPrior + etaPrior + tauPrior + kappaPrior + rhoPrior + betaPrior
}

get_model_priors.nome.nophi <- function () {
  NuPrior <- prior(normal(-3.00583363 ,  2.12854828), nlpar = "Nu")
  etaPrior <- prior(normal(-0.7654982,  1.33570051), nlpar = "eta")
  tauPrior <- prior(normal(-3.80045123,  0.64387023), nlpar = "tau")
  kappaPrior <-
    prior(normal(-3.59771868,  0.7310165), nlpar = "kappa")
  rhoPrior <- prior(normal(-3.59771868,  0.7310165), nlpar = "rho")
  
  NuPrior + etaPrior + tauPrior + kappaPrior + rhoPrior
}

get_model_priors.nome.nophi.scaled <- function () {
  NuPrior <- prior(std_normal(), nlpar = "Nu")
  etaPrior <- prior(std_normal(), nlpar = "eta")
  tauPrior <- prior(std_normal(), nlpar = "tau")
  rhoPrior <- prior(std_normal(), nlpar = "rho")
  
  NuPrior + etaPrior + tauPrior + rhoPrior
}