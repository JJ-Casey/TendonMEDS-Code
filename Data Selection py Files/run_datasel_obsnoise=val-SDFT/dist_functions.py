import numpy as np
import scipy as sc


def prior_transtheta(prior_mu, prior_sig, x):
    # calculate the ratio of prior probabilities for each of the parameters, working in the transformed parameter space
    # (negative log of prior)
    term = sc.stats.norm.logpdf(x, loc=prior_mu, scale=prior_sig)
    return - np.sum(term)

def prior_transgamma(mu, invK, trans_gamma):
    # GP prior on transformed gamma terms
    # (negative log of prior) density logit-multivariate normal
    diff = np.asarray(trans_gamma - mu)
    val = 0.5 * (diff @ invK @ diff)
    return val

# def likelihood(prior_params, model_data, stress, gamma, num_data):
#     # - log of student's t distribution post_pred (simplified and without constants)
#     term1a = prior_params.alpha_sigma + (num_data/2)
#     term1b = np.sum(np.multiply(gamma, np.power(stress-model_data, 2)))
#     term1c = prior_params.beta_sigma + (term1b/2)
#     term1 = term1a * np.log(term1c)
#     term2 = - 0.5 * np.sum(np.log(gamma))  # this is the gamma dependent constant
#     return term1 + term2

def likelihood_gauss(stress, model_data, gamma, current_vals):
    # negative log of gaussian likelihood with chosen value of sigmasq
    diff = stress - model_data
    term1 = np.multiply(gamma / current_vals.obsnoise, np.power(diff, 2))
    term2 = np.log(gamma)
    return 0.5 * np.sum(term1 - term2)

def ref_pullback_T(theta_param):
    # negative log of reference pullback of the forward transformation (model params)
    # input (1-phi)mu, phiE, a-1, c-a, b-c
    return np.sum(np.log(theta_param))

def ref_pullback_T_gamma(gamma_param):
    # negative log of reference pullback of the forward transformation (model params)
    # input fidelity parameters (bounded on 0,1)
    term = np.log(1-gamma_param) + np.log(gamma_param)
    return np.sum(term)

def get_nlog_posterior(current_vals):
    # checking the min negative log posterior ON THE UN-TRANSFORMED SPACE
    pullback_theta = ref_pullback_T(current_vals.theta)
    pullback_gamma = ref_pullback_T_gamma(current_vals.gamma)

    nlog_posterior = current_vals.prior_transtheta + current_vals.prior_transgamma + current_vals.likelihood \
                     + pullback_theta + pullback_gamma
    return nlog_posterior

def prior_sigmasq(x, alpha, beta):
    # negative log of the inverse gamma prior on sigma^2 density
    return - sc.stats.invgamma.logpdf(x, alpha, scale=beta)



