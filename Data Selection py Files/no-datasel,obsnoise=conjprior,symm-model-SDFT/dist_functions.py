import numpy as np
import scipy as sc


def prior_transtheta(prior_mu, prior_sig, x):
    # calculate the ratio of prior probabilities for each of the parameters, working in the transformed parameter space
    # (negative log of prior)
    term = sc.stats.norm.logpdf(x, loc=prior_mu, scale=prior_sig)
    return - np.sum(term)

def likelihood(prior_params, model_data, stress, gamma, num_data):
    # - log of student's t distribution post_pred (simplified and without constants)
    term1a = prior_params.alpha_sigma + (num_data/2)
    term1b = np.sum(np.multiply(gamma, np.power(stress-model_data, 2)))
    term1c = prior_params.beta_sigma + (term1b/2)
    term1 = term1a * np.log(term1c)
    term2 = - 0.5 * np.sum(np.log(gamma))  # this is the gamma dependent constant
    return term1 + term2

def ref_pullback_T(theta_param):
    # negative log of reference pullback of the forward transformation (model params)
    # input (1-phi)mu, phiE, a-1, c-a, b-c
    return np.sum(np.log(theta_param))

def get_nlog_posterior(current_vals):
    # checking the min negative log posterior ON THE UN-TRANSFORMED SPACE
    pullback_theta = ref_pullback_T(current_vals.theta)

    nlog_posterior = current_vals.prior_transtheta + current_vals.likelihood + pullback_theta
    return nlog_posterior

