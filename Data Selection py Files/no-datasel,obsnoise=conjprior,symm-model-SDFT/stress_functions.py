import numpy as np
import math


def choose_stress_model(trans_param, lamb):
    # this function chooses which model of the stress should be used depending on the value of the strain
    # ORDER: nu=ln((1-phi)*mu), eta=ln(E*phi), tau=ln(a-1), xi=ln(c-a),chi=ln(b-c),
    a = 1 + np.exp(trans_param[2])
    b = a + np.exp(trans_param[3])

    c = (a + b) / 2

    ncmterm = stressmodel_NCM(lamb, trans_param[0])
    # when lamb**2>=a**2 && lamb**2<=c**2
    fibrilterm1 = np.where((lamb >= a) & (lamb <= c), stressmodel_region2(lamb, trans_param[2], trans_param[3],
                                                                          trans_param[4], trans_param[1]), lamb * 0)
    # when lamb**2<=b**2 && lamb**2>c**2
    fibrilterm2 = np.where((lamb <= b) & (lamb > c), stressmodel_region3(lamb, trans_param[2], trans_param[3],
                                                                         trans_param[4], trans_param[1]), lamb * 0)
    # when lamb**2>b**2
    fibrilterm3 = np.where(lamb > b, stressmodel_region4(lamb, trans_param[2], trans_param[3],
                                                         trans_param[4], trans_param[1]), lamb * 0)

    return ncmterm + fibrilterm1 + fibrilterm2 + fibrilterm3


def stressmodel_NCM(lamb, nu):
    # stress of non-collagen matrix (when lamb**2<a**2)
    # lamb=stretch, nu=ln((1-phi)*mu)
    val = math.exp(nu) * (lamb - (1 / lamb ** 2))
    return val


def stressmodel_region2(lamb, tau, xi, chi, eta):
    # region where lamb**2>a**2 && lamb**2<c**2
    # lamb=stretch, tau=ln(a-1), xi=ln(c-a),chi=ln(b-c), eta=ln(E*phi)
    term1 = lamb ** 2
    term2 = -2 * (1 + math.exp(tau)) * lamb * np.log(lamb)
    term2b = 2 * (1 + math.exp(tau)) * lamb * np.log(math.exp(tau) + 1)
    term3 = - (1 + math.exp(tau)) ** 2
    prelim = (term1 + term2 + term2b + term3) / ((math.exp(xi) + math.exp(chi)) * math.exp(xi))

    val = (math.exp(eta) / lamb) * prelim
    return val


def stressmodel_region3(lamb, tau, xi, chi, eta):
    # region where lamb**2<b**2 && lamb**2>c**2
    # lamb=stretch, tau=ln(a-1), xi=ln(c-a),chi=ln(b-c), eta=ln(E*phi)
    term1 = (1 + math.exp(tau) + math.exp(xi)) ** 2 / (math.exp(xi) * math.exp(chi))
    term2 = -(1 + math.exp(tau)) ** 2 / ((math.exp(xi) + math.exp(chi)) * math.exp(xi))
    term3 = (2 * (1 + math.exp(tau)) * lamb * np.log(1 + math.exp(tau))) / \
            ((math.exp(xi) + math.exp(chi)) * math.exp(xi))
    term4 = -(2 * (1 + math.exp(tau) + math.exp(xi)) * lamb * np.log(1 + math.exp(tau)
                                                                     + math.exp(xi))) / (math.exp(xi) * math.exp(chi))
    term5 = - lamb ** 2 / ((math.exp(xi) + math.exp(chi)) * math.exp(chi))
    term6 = (2 * (1 + math.exp(tau) + math.exp(xi) + math.exp(chi)) * lamb *
             np.log(lamb)) / ((math.exp(xi) + math.exp(chi)) * math.exp(chi))

    val = (math.exp(eta) / lamb) * (term1 + term2 + term3 + term4 + term5 + term6)
    return val


def stressmodel_region4(lamb, tau, xi, chi, eta):
    # region where lamb**2>b**2
    # lamb=stretch, tau=ln(a-1), xi=ln(c-a),chi=ln(b-c), eta=ln(E*phi)
    term1 = -1
    term2 = (2 * (1 + math.exp(tau)) * lamb * np.log(1 + math.exp(tau))) / \
            ((math.exp(xi) + math.exp(chi)) * math.exp(xi))
    term3 = (2 * (1 + math.exp(tau) + math.exp(xi) + math.exp(chi)) * lamb
             * np.log(1 + math.exp(tau) + math.exp(xi)
                      + math.exp(chi))) / ((math.exp(xi) + math.exp(chi)) * math.exp(chi))
    term4 = -(2 * (1 + math.exp(tau) + math.exp(xi)) * lamb * np.log(1 + math.exp(tau)
                                                                     + math.exp(xi))) / (math.exp(xi) * math.exp(chi))

    val = (math.exp(eta) / lamb) * (term1 + term2 + term3 + term4)
    return val
