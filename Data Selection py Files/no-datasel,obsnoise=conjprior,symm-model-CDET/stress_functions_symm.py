import numpy as np
import math


def choose_stress_model(theta, lamb):
    # this function chooses which model of the stress should be used depending on the value of the strain
    # ORDER: nu=ln((1-phi)*mu), eta=ln(E*phi), tau=ln(a-1), xi=ln(c-a),chi=ln(b-c),
    phi_mu = theta[0]  # this is actually phi_mu
    phi_E = theta[1]
    a = 1 + theta[2]
    b = a + theta[3]

    c = (a + b) / 2


    ncmterm = stressmodel_NCM(lamb, phi_mu)
    # when lamb**2>=a**2 && lamb**2<=c**2
    fibrilterm1 = np.where((lamb >= a) & (lamb <= c), stressmodel_region2(lamb, phi_E, a, b, c), lamb * 0)
    # when lamb**2<=b**2 && lamb**2>c**2
    fibrilterm2 = np.where((lamb <= b) & (lamb > c), stressmodel_region3(lamb, phi_E, a, b, c), lamb * 0)
    # when lamb**2>b**2
    fibrilterm3 = np.where(lamb > b, stressmodel_region4(lamb, phi_E, a, b, c), lamb * 0)

    return ncmterm + fibrilterm1 + fibrilterm2 + fibrilterm3


def stressmodel_NCM(lamb, phi_mu):
    # stress of non-collagen matrix (when lamb**2<a**2)
    # lamb=stretch, nu=ln((1-phi)*mu)
    val = phi_mu * (lamb - (1 / lamb ** 2))
    return val


def stressmodel_region2(lamb, phi_E, a, b, c):
    # region where lamb**2>a**2 && lamb**2<c**2
    term1 = lamb ** 2
    term2 = -2 * a * lamb * np.log(lamb)
    term2b = 2 * a * lamb * np.log(a)
    term3 = - (a) ** 2
    prelim = (term1 + term2 + term2b + term3) / ((b-a) * (c-a))

    val = (phi_E / lamb) * prelim
    return val


def stressmodel_region3(lamb, phi_E, a, b, c):
    # region where lamb**2<b**2 && lamb**2>c**2
    term1 = c ** 2 / ((c-a) * (b-c))
    term2 = -(a) ** 2 / ((b-a) * (c-a))
    term3 = (2 * a * lamb * np.log(a)) / ((b-a) * (c-a))
    term4 = -(2 * c * lamb * np.log(c)) / ((c-a) * (b-c))
    term5 = - lamb ** 2 / ((b-a) * (b-c))
    term6 = (2 * b * lamb * np.log(lamb)) / ((b-a) * (b-c))

    val = (phi_E / lamb) * (term1 + term2 + term3 + term4 + term5 + term6)
    return val


def stressmodel_region4(lamb, phi_E, a, b, c):
    # region where lamb**2>b**2
    term1 = -1
    term2 = (2 * a * lamb * np.log(a)) / ((b-a) * (c-a))
    term3 = (2 * b * lamb * np.log(b)) / ((b-a) * (b-c))
    term4 = -(2 * c * lamb * np.log(c)) / ((c-a) * (b-c))

    val = (phi_E / lamb) * (term1 + term2 + term3 + term4)
    return val
