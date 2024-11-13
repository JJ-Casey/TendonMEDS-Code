library(brms)

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

