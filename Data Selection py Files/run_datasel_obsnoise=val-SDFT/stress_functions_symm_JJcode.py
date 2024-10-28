import numpy as np

# evaluation of modelled stress using code from JJ

def choose_stress_model(theta, lam_range):
    # theta is a vector of the true thetaeters
    ## theta[0] = (1-phi)mu
    ## theta[1] = phi E
    ## theta[2] = a - 1
    ## theta[3] = b - a
    #
    # lam is a vector of values for lamda

    a = theta[2] + 1
    b = theta[3] + a

    # Symmetric Model
    c = (a + b) / 2

    # Piecewise functions
    def A(i4):
        if i4 < a ** 2:
            return 0
        if a ** 2 <= i4 <= c ** 2:
            return -a ** 2 / ((b - a) * (c - a))
        if c ** 2 < i4 <= b ** 2:
            return c ** 2 / ((c - a) * (b - c)) - a ** 2 / ((b - a) * (c - a))
        if b ** 2 < i4:
            return -1

    def B(i4):
        if i4 < a ** 2:
            return 0
        if a ** 2 <= i4 <= c ** 2:
            return 2 * a * np.log(a) / ((b - a) * (c - a))
        if c ** 2 < i4 <= b ** 2:
            return 2 * a * np.log(a) / ((b - a) * (c - a)) - 2 * c * np.log(c) / ((b - c) * (c - a))
        if b ** 2 < i4:
            return 2 * a * np.log(a) / ((b - a) * (c - a)) + 2 * b * np.log(b) / ((b - c) * (b - a)) - 2 * c * np.log(
                c) / ((b - c) * (c - a))

    def C(i4):
        if i4 < a ** 2:
            return 0
        if a ** 2 <= i4 <= c ** 2:
            return 1 / ((b - a) * (c - a))
        if c ** 2 < i4 <= b ** 2:
            return -1 / ((b - a) * (b - c))
        if b ** 2 < i4:
            return 0

    def D(i4):
        if i4 < a ** 2:
            return 0
        if a ** 2 <= i4 <= c ** 2:
            return -2 * a / ((b - a) * (c - a))
        if c ** 2 < i4 <= b ** 2:
            return 2 * b / ((b - a) * (b - c))
        if b ** 2 < i4:
            return 0

    # Extract other theta
    phimu, phiE = theta[0], theta[1]

    # Return the model using list comprehension, using the equation in JH paper 3.2
    return [phimu * (lam - lam ** -2) + phiE * lam ** -1 *
            (A(lam ** 2) + B(lam ** 2) * lam + C(lam ** 2) * lam ** 2 + D(lam ** 2) * lam * np.log(lam))
            for lam in lam_range]






