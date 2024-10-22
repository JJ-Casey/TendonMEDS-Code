#Define the piecewise functions A, B, C, D
pieceA <- function (i4, a, b, c) {
  (i4 < a ^ 2) * (0) +
    ((a ^ 2 <= i4) &
       (i4 <= c ^ 2)) * (-a ^ 2 / ((b - a) * (c - a))) +
    ((c ^ 2 < i4) &
       (i4 <= b ^ 2)) * (c ^ 2 / ((c - a) * (b - c)) - a ^ 2 / ((b - a) * (c -
                                                                             a))) +
    (b ^ 2 < i4) * (-1)
}
pieceA.stanCode <- "
 real pieceA(real i4, real a, real b, real c) {
  return (i4 < a^2) * 0
  + ((a^2 <= i4) * (i4 < c^2)) * (-a^2 / ((b-a) * (c-a)))
  + ((c^2 <= i4) * (i4 < b^2)) * (c^2 / ((c-a) * (b-c)) - a^2 / ((b-a) * (c-a)))
  + (b^2 <= i4) * (-1);
}
"

pieceB <- function (i4, a, b, c) {
  (i4 < a ^ 2) * (0) +
    ((a ^ 2 <= i4) &
       (i4 <= c ^ 2)) * (2 * a * log(a) / ((b - a) * (c -
                                                        a))) +
    ((c ^ 2 < i4) &
       (i4 <= b ^ 2)) * (2 * a * log(a) / ((b - a) * (c - a)) - 2 * c * log(c) / ((b -
                                                                                     c) * (c - a))) +
    (b ^ 2 < i4) * (2 * a * log(a) / ((b - a) * (c - a)) + 2 * b * log(b) / ((b -
                                                                                c) * (b - a)) - 2 * c * log(c) / ((b - c) * (c - a)))
}
pieceB.stanCode <- "
real pieceB(real i4, real a, real b, real c) {
  return (i4 < a^2) * (0)
  + ((a^2 <= i4) * (i4 < c^2)) * (2*a*log(a) / ((b-a) * (c-a)))
  + ((c^2 <= i4) * (i4 <= b^2)) * (2*a*log(a) / ((b-a) * (c-a)) - 2*c*log(c) / ((b-c) * (c-a)))
  + (b^2 <= i4) * (2*a*log(a) / ((b-a) * (c-a)) + 2*b*log(b) / ((b-c) * (b-a)) - 2*c*log(c) / ((b-c) * (c-a)));
}
"

pieceC <- function (i4, a, b, c) {
  (i4 < a ^ 2) * (0) +
    ((a ^ 2 <= i4) & (i4 <= c ^ 2)) * (1 / ((b - a) * (c - a))) +
    ((c ^ 2 < i4) & (i4 <= b ^ 2)) * (-1 / ((b - a) * (b - c))) +
    (b ^ 2 < i4) * (0)
}
pieceC.stanCode <- "
real pieceC(real i4, real a, real b, real c) {
  return (i4 < a^2) * 0
  + ((a^2 <= i4) * (i4 < c^2)) * (1 / ((b-a) * (c-a)))
  + ((c^2 <= i4) * (i4 < b^2)) * (-1 / ((b-a) * (b-c)))
  + (b^2 <= i4) * 0;
}
"

pieceD <- function (i4, a, b, c) {
  (i4 < a ^ 2) * (0) +
    ((a ^ 2 <= i4) &
       (i4 <= c ^ 2)) * (-2 * a / ((b - a) * (c - a))) +
    ((c ^ 2 < i4) & (i4 <= b ^ 2)) * (2 * b / ((b - a) * (b - c))) +
    (b ^ 2 < i4) * (0)
}
pieceD.stanCode <- "
real pieceD(real i4, real a, real b, real c) {
  return (i4 < a^2) * 0
  + ((a^2 <= i4) * (i4 < c^2)) * (-2*a / ((b-a) * (c-a)))
  + ((c^2 <= i4) * (i4 < b^2)) * (2*b / ((b-a) * (b-c)))
  + (b^2 <= i4) * 0;
}
"

tendonModel <-
  function (stretch,
            stretchm1,
            stretch2,
            stretch2m1,
            stretchLogStretch,
            nu,
            eta,
            tau,
            rho) {
    phimu <- exp(1.05309738 + 1.30927056 * nu)
    phiE <- exp(6.83672018 + 0.47191773 * eta)
    a <- exp(-3.80045123 + 0.64387023 * tau) + 1
    b <- exp(-3.59771868 + 0.7310165 * rho) + a
    c <- (a + b) / 2
    
    phimu * (stretch - stretch2m1) +
      phiE * stretchm1 * (
        pieceA(stretch2, a, b, c) +
          pieceB(stretch2, a, b, c) * stretch +
          pieceC(stretch2, a, b, c) * stretch2 +
          pieceD(stretch2, a, b, c) * stretchLogStretch
      )
  }

region1 <- "

"

tendonModel.stanCode <- "
real tendonModel (real stretch, real stretchm1, real stretch2, real stretchmstretch2m1, real stretchLogStretch, real nu, real eta, real tau, real rho) {
  real phimu = exp(1.05309738 + 1.30927056 * nu);
  real phiE = exp(6.83672018 + 0.47191773 * eta);
  real a = exp(-3.80045123 + 0.64387023 * tau) + 1;
  real b = exp(-3.59771868 + 0.7310165 * rho) + a;
  real c = (a + b) / 2;

  real index1 = (stretch2 < a^2);
  real index2 = (stretch2 >= a^2) * (stretch2 < c^2);
  real index3 = (stretch2 >= c^2) * (stretch2 < b^2);
  real index4 = (stretch2 >= b^2);

  // A + Bl + Cl^2 + DlLogl
  real part1 = phimu * stretchmstretch2m1;
  real part2 = index2 * (-(a ^ 2) + 2 * a * log(a) * stretch + stretch2 + -2 * a * stretchLogStretch) / ((b - a) * (c - a));
  real part3 = index3 * ((c ^ 2 / ((c - a) * (b - c)) - (a ^ 2) / ((b - a) * (c - a))) +
                              (2 * a * log(a) / ((b - a) * (c - a)) - 2 * c * log(c) / ((b - c) * (c - a))) * stretch +
                              (-1 / ((b - a) * (b - c))) * stretch2 +
                              (2 * b / ((b - a) * (b - c))) * stretchLogStretch);
  real part4 = index4 * (-1 +
                              (2 * a * log(a) / ((b - a) * (c - a))
                                  + 2 * b * log(b) / ((b - c) * (b - a))
                                  - 2 * c * log(c) / ((b - c) * (c - a))) * stretch);

  return part1 + phiE * stretchm1 * (part2 + part3 + part4);
}
"
tendonModel.stanCode.multinorm <- "
real tendonModel (real stretch, real stretchm1, real stretch2, real stretchmstretch2m1, real stretchLogStretch, vector theta) {
  real phimu = exp(1.05309738 + 1.30927056 * theta[1]);
  real phiE = exp(6.83672018 + 0.47191773 * theta[2]);
  real a = exp(-3.80045123 + 0.64387023 * theta[3]) + 1;
  real b = exp(-3.59771868 + 0.7310165 * theta[4]) + a;
  real c = (a + b) / 2;

  real index1 = (stretch2 < a^2);
  real index2 = (stretch2 >= a^2) * (stretch2 < c^2);
  real index3 = (stretch2 >= c^2) * (stretch2 < b^2);
  real index4 = (stretch2 >= b^2);

  // A + Bl + Cl^2 + DlLogl
  real part1 = phimu * stretchmstretch2m1;
  real part2 = index2 * (-(a ^ 2) + 2 * a * log(a) * stretch + stretch2 + -2 * a * stretchLogStretch) / ((b - a) * (c - a));
  real part3 = index3 * ((c ^ 2 / ((c - a) * (b - c)) - (a ^ 2) / ((b - a) * (c - a))) +
                              (2 * a * log(a) / ((b - a) * (c - a)) - 2 * c * log(c) / ((b - c) * (c - a))) * stretch +
                              (-1 / ((b - a) * (b - c))) * stretch2 +
                              (2 * b / ((b - a) * (b - c))) * stretchLogStretch);
  real part4 = index4 * (-1 +
                              (2 * a * log(a) / ((b - a) * (c - a))
                                  + 2 * b * log(b) / ((b - c) * (b - a))
                                  - 2 * c * log(c) / ((b - c) * (c - a))) * stretch);

  return part1 + phiE * stretchm1 * (part2 + part3 + part4);
}
"

get_model_stancode <- function() {
  # stanvar(scode = pieceA.stanCode, block = "function") +
  #   stanvar(scode = pieceB.stanCode, block = "function") +
  #   stanvar(scode = pieceC.stanCode, block = "function") +
  #   stanvar(scode = pieceD.stanCode, block = "function") +
    stanvar(scode = tendonModel.stanCode, block = "function")
}

get_model_stancode.multinorm <- function() {
  # stanvar(scode = pieceA.stanCode, block = "function") +
  #   stanvar(scode = pieceB.stanCode, block = "function") +
  #   stanvar(scode = pieceC.stanCode, block = "function") +
  #   stanvar(scode = pieceD.stanCode, block = "function") +
  stanvar(scode = tendonModel.stanCode.multinorm, block = "function")
}