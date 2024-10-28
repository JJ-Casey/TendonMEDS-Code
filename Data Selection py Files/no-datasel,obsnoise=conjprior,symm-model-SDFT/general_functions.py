import numpy as np
import math
import dist_functions as distf
import stress_functions_symm as stressf
import scipy as sc
import os

def sample_theta(it, prop_cov, strain,stress, num_data, prior_params, current_vals, arrays, MCMC_params):
    # this function performs sampling on the model parameters (theta)
    # make new proposal
    w = np.random.multivariate_normal(np.zeros(MCMC_params.num_param), prop_cov)
    v_trans_theta = current_vals.trans_theta + np.multiply(MCMC_params.beta_theta, w)
    v_theta = np.exp(v_trans_theta)

    # new model data
    v_model_data = stressf.choose_stress_model(v_theta, strain)
    # calculate prior density, transition ratio, and posterior predicitve
    v_prior_transtheta = distf.prior_transtheta(prior_params.prior_mu, prior_params.prior_sig, v_trans_theta)
    v_likelihood = distf.likelihood(prior_params, v_model_data, stress, np.ones(num_data), num_data)

    # calculate acceptance ratio
    alpha_term = min(0, current_vals.prior_transtheta - v_prior_transtheta + current_vals.likelihood - v_likelihood)
    alpha = min(1, math.exp(alpha_term))
    # accept or reject
    if alpha > np.random.random():
        # accept
        current_vals.trans_theta = v_trans_theta
        current_vals.theta = v_theta
        current_vals.model_data = v_model_data
        current_vals.likelihood = v_likelihood
        current_vals.prior_transtheta = v_prior_transtheta

    # append to arrays
    if it < MCMC_params.N_burnin:
        idx = it - ((it // (MCMC_params.num_samp_cov)) * MCMC_params.num_samp_cov) - 1
        arrays.trans_theta_proposal[int(idx), :] = current_vals.trans_theta

    # check acceptance rate, change beta if necessary
    [arrays.alpha_theta, arrays.acc_rate_theta, MCMC_params.beta_theta, arrays.beta_theta] = check_acc_rate(
        it, MCMC_params.check_freq, alpha, arrays.alpha_theta, arrays.acc_rate_theta, MCMC_params.beta_theta,
        arrays.beta_theta, MCMC_params.upper_tol, MCMC_params.lower_tol, MCMC_params.N_burnin)

    # adapt model proposal covariance matrix if within correct it range
    if MCMC_params.num_samp_cov < it < MCMC_params.N_burnin and (it % MCMC_params.num_samp_cov) == 0:
        # (every x' iterations update proposal covariance)
        prop_cov = learn_prop_cov(arrays.trans_theta_proposal, MCMC_params.num_param)
        arrays.prop_cov.append(np.reshape(prop_cov, MCMC_params.num_param ** 2).tolist())

    return [prop_cov, current_vals, MCMC_params, arrays]

def sample_obsnoise(stress, prior_params, current_vals, num_data):     #(stress, prior_params, current_vals):
    # sample directly from inv gamma distribution
    alpha_hat = prior_params.alpha_sigma + (num_data / 2)
    beta_hat = prior_params.beta_sigma + (0.5 * np.sum(np.power(stress-current_vals.model_data, 2)))
    current_vals.obs_noise = sc.stats.invgamma.rvs(alpha_hat, scale=beta_hat)

    return current_vals


def write_data(it, MCMC_params, current_vals, arrays, tosaveto):
    if it % MCMC_params.thin == 0 and it >= MCMC_params.N_burnin:
        idx = (it / MCMC_params.thin) - ((it // (MCMC_params.thin * MCMC_params.N_write)) * (MCMC_params.N_write))
        arrays.trans_theta[int(idx), :] = current_vals.trans_theta
        arrays.theta[int(idx), :] = current_vals.theta
        arrays.model_fit[int(idx), :] = current_vals.model_data
        arrays.obs_noise[int(idx)] = current_vals.obs_noise


        if idx == (MCMC_params.N_write - 1):
            # append data to text files
            thetatrans_txt = open(os.path.join(tosaveto, 'theta_trans.txt'), 'a')
            theta_txt = open(os.path.join(tosaveto, 'theta.txt'), 'a')
            modelfit_txt = open(os.path.join(tosaveto, 'modelfit.txt'), 'a')
            obsnoise_txt = open(os.path.join(tosaveto, 'obs_noise.txt'), 'a')

            np.savetxt(thetatrans_txt, arrays.trans_theta)
            np.savetxt(theta_txt, arrays.theta)
            np.savetxt(modelfit_txt, arrays.model_fit)
            np.savetxt(obsnoise_txt, arrays.obs_noise)

            thetatrans_txt.close()
            theta_txt.close()
            modelfit_txt.close()
            obsnoise_txt.close()

    return arrays


def check_acc_rate(it, check_freq, alpha, alpha_array, acc_rate_array, beta_step, beta_array,
                   upper_tol, lower_tol, N_burnin):

    # Check the acceptance rate of proposals and adjust step-size (beta) if neceesary.

    if (it % check_freq) != 0:
        alpha_array[int((it % check_freq) - 1)] = alpha
        # print(int((it % check_freq) - 1), alpha)
    else:
        alpha_array[check_freq - 1] = alpha
        acc_rate_array[int((it / check_freq) - 1)] = np.mean(alpha_array)
        if acc_rate_array[int((it / check_freq) - 1)] > upper_tol and it < N_burnin:
            beta_step = beta_step * 1.05
        elif acc_rate_array[int((it / check_freq) - 1)] < lower_tol and it < N_burnin:
            beta_step = beta_step * 0.95

        beta_array[int((it / check_freq) - 1)] = beta_step

    return [alpha_array, acc_rate_array, beta_step, beta_array]


def learn_prop_cov(trans_array, num_param):
    # learn the proposal covariance matrix (for the model parameters only), to improve proposals during sampling
    prop_cov = np.cov(trans_array, rowvar=False)
    prop_cov = prop_cov + (1e-5 * np.identity(num_param))  # add a small amount to make sure it is positive-semi-def
    return prop_cov


def thetadata_to_abc(data):
    data[:, 2] = data[:, 2] + 1  # these are now the values of the parameter a
    data[:, 3] = data[:, 3] + data[:, 2]  # these are now the values of parameter b
    return data


