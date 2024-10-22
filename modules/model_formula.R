library(brms)

mixed_effects_model.phi <- function() {
  bf(
    stress |
      se(errs) ~ tendonModel(
        stretch,
        stretchm1,
        stretch2,
        stretchmstretch2m1,
        stretchLogStretch,
        Nu,
        eta,
        beta,
        tau,
        kappa,
        rho
      ),
    Nu ~ 1,
    eta ~ 1,
    beta ~ (1 | p | group),
    tau ~ (1 | p | group),
    kappa ~ (1 | p | group),
    rho ~ (1 | p | group),
    nl = TRUE,
    decomp = "QR"
  )
}

mixed_effects_model.nophi <- function() {
  bf(
    stress |
      se(errs) ~ tendonModel(
        stretch,
        stretchm1,
        stretch2,
        stretchmstretch2m1,
        stretchLogStretch,
        Nu,
        eta,
        tau,
        kappa,
        rho
      ),
    Nu ~ (1 | p | group),
    eta ~ (1 | p | group),
    tau ~ (1 | p | group),
    kappa ~ (1 | p | group),
    rho ~ (1 | p | group),
    nl = TRUE,
    decomp = "QR"
  )
}

mixed_effects_model.nophi.st <- function() {
  bf(
    stress |
      se(errs) ~ tendonModel(
        stretch,
        stretchm1,
        stretch2,
        stretchmstretch2m1,
        stretchLogStretch,
        Nu,
        eta,
        tau,
        rho
      ),
    Nu ~ (1 | p | group),
    eta ~ (1 | p | group),
    tau ~ (1 | p | group),
    rho ~ (1 | p | group),
    nl = TRUE,
    decomp = "QR"
  )
}

mixed_effects_model.nophi.st.multinorm <- function() {
  bf(
    stress |
      se(errs) ~ tendonModel(
        stretch,
        stretchm1,
        stretch2,
        stretchmstretch2m1,
        stretchLogStretch,
        theta
      ),
    theta ~ (1 | p | group),
    nl = TRUE,
    decomp = "QR"
  )
}

mixed_effects_model.nophi.st.centred <- function() {
  bf(
    stress |
      se(errs) ~ tendonModel(
        stretch,
        stretchm1,
        stretch2,
        stretchmstretch2m1,
        stretchLogStretch,
        Nu,
        eta,
        tau,
        rho
      ),
    Nu ~ (1 | p | group),
    eta ~ (1 | p | group),
    tau ~ (1 | p | group),
    rho ~ (1 | p | group),
    nl = TRUE,
    decomp = "QR"
  )
}

mixed_effects_model.nophi.st.nocorr <- function() {
  bf(
    stress |
      se(errs) ~ tendonModel(
        stretch,
        stretchm1,
        stretch2,
        stretchmstretch2m1,
        stretchLogStretch,
        Nu,
        eta,
        tau,
        rho
      ),
    Nu ~ (1 | group),
    eta ~ (1 | group),
    tau ~ (1 | group),
    rho ~ (1 | group),
    nl = TRUE,
    decomp = "QR"
  )
}

single_data_model.phi <- function() {
  bf(
    stress |
      se(errs) ~ tendonModel(
        stretch,
        stretchm1,
        stretch2,
        stretchmstretch2m1,
        stretchLogStretch,
        Nu,
        eta,
        beta,
        tau,
        kappa,
        rho
      ),
    Nu ~ 1,
    eta ~ 1,
    beta ~ 1,
    tau ~ 1,
    kappa ~ 1,
    rho ~ 1,
    nl = TRUE,
    decomp = "QR"
  )
}

single_data_model.nophi <- function() {
  bf(
    stress |
      se(errs) ~ tendonModel(
        stretch,
        stretchm1,
        stretch2,
        stretchmstretch2m1,
        stretchLogStretch,
        Nu,
        eta,
        tau,
        kappa,
        rho
      ),
    Nu ~ 1,
    eta ~ 1,
    tau ~ 1,
    kappa ~ 1,
    rho ~ 1,
    nl = TRUE,
    decomp = "QR"
  )
}

single_data_model.st.nophi <- function() {
  bf(
    stress |
      se(errs) ~ tendonModel(
        stretch,
        stretchm1,
        stretch2,
        stretchmstretch2m1,
        stretchLogStretch,
        Nu,
        eta,
        tau,
        rho
      ),
    Nu ~ 1,
    eta ~ 1,
    tau ~ 1,
    rho ~ 1,
    nl = TRUE,
    decomp = "QR"
  )
}