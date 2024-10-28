import numpy as np
import dist_functions as distf
import stress_functions_symm as stressf

def import_data(fnametxt):
    # import stress-strain data from text file
    data = np.genfromtxt(fname=fnametxt)
    strain = data[:, 0]
    stress = data[:, 1]

    # select only data up to 5% strain
    id_lt5perc = np.where(strain < 1.05)
    strain = strain[id_lt5perc[0]]
    stress = stress[id_lt5perc[0]]
    print('Data trimmed to 0-5% strain only.')

    return [strain, stress]

def parameters(num_data, strain, stress, chain):
    class prior_params():  #structure to store all parameters relating to the priors (or hyper priors within model)
        pass
    # model parameter prior parameters
    [nu_sig, nu_mu] = [1.30927056, 1.05309738 ]    #(1-phi)mu
    [eta_sig, eta_mu] = [0.47191773, 6.83672018]    #phi E
    [tau_sig, tau_mu] = [0.64387023, -3.80045123]   # a - 1
    [rho_sig, rho_mu] = [0.7310165, -3.59771868]     # b - a
    prior_params.prior_mu = np.asarray([nu_mu, eta_mu, tau_mu, rho_mu])         # model parameter prior mu vals
    prior_params.prior_sig = np.asarray([nu_sig, eta_sig, tau_sig, rho_sig])  # model parameter prior sigma vals

    #hyper params
    prior_params.alpha_sigma = 3                    # hyper parameter for noise conjugate prior
    prior_params.beta_sigma = 0.3                   # hyper parameter for noise conjugate prior
    prior_params.df = 2 * prior_params.alpha_sigma  # parameter of the student's t-distribution

    # proposal covariance
    prop_cov = np.diag(np.power(prior_params.prior_sig,2))  # proposal covariance for model parameter RW

    class current_vals():  # structure to store all the current values of parameters, prior densities, likelihood etc.
        pass
    # Initialise start positions
    # model start
    current_vals.trans_theta = np.random.normal(prior_params.prior_mu, prior_params.prior_sig)
    current_vals.theta = np.exp(current_vals.trans_theta)

   # evaluate model data using these initial values
    current_vals.model_data = stressf.choose_stress_model(current_vals.theta, strain)
    current_vals.likelihood = distf.likelihood(prior_params, current_vals.model_data, stress, np.ones(num_data), num_data)
    #prior densities
    current_vals.prior_transtheta = distf.prior_transtheta(prior_params.prior_mu, prior_params.prior_sig, current_vals.trans_theta)

    return [prior_params, current_vals, prop_cov]

def MCMC_parameters():
    class MCMC_params():
        pass
    # MCMC parameters
    MCMC_params.N = int(2e6)
    MCMC_params.N_burnin = int(MCMC_params.N / 2)
    MCMC_params.N_write = int(5000)
    MCMC_params.num_samp_cov = 10000  # number of samples between proposal covariance updates
    MCMC_params.num_param = 4  # number of model parameters
    MCMC_params.beta_theta = 1  # RWMH step size for model module
    MCMC_params.check_freq = int(500)  # frequency to check acceptance rate
    MCMC_params.thin = 10  # thinning parameter
    MCMC_params.tol = 0.05  # tolerance of acceptance rate
    MCMC_params.upper_tol = 0.234 + MCMC_params.tol
    MCMC_params.lower_tol = 0.234 - MCMC_params.tol

    return [MCMC_params]

def arrays(MCMC_params, num_data):
    class arrays():
        pass
    arrays.alpha_theta = np.empty(MCMC_params.check_freq)
    arrays.acc_rate_theta = np.empty(int(MCMC_params.N / MCMC_params.check_freq))
    arrays.beta_theta = np.empty(int(MCMC_params.N / MCMC_params.check_freq))
    arrays.nlog_post = np.empty(int(MCMC_params.N / MCMC_params.check_freq))
    arrays.nlog_likelihood = np.empty(int(MCMC_params.N / MCMC_params.check_freq))

    arrays.trans_theta_proposal = np.empty([int(MCMC_params.num_samp_cov), MCMC_params.num_param], dtype=np.float32)
    arrays.trans_theta = np.empty([int(MCMC_params.N_write), MCMC_params.num_param], dtype=np.float32)
    arrays.theta = np.empty([int(MCMC_params.N_write), MCMC_params.num_param], dtype=np.float32)
    arrays.prop_cov = []
    arrays.model_fit = np.empty([int(MCMC_params.N_write), num_data], dtype=np.float32)
    arrays.obs_noise = np.empty([int(MCMC_params.N_write)], dtype=np.float32)

    return [arrays]

# def get_prior_params(upper_val, lower_val, upper_cpd, lower_cpd):
#     # calculate the prior means and standard deviaitions given the conditions in Table 1 in the Supplementary
#     # these are the means and std of the normal priors on the transformed parameters
#     numerator = np.log(upper_val) - np.log(lower_val)
#     denominator = math.sqrt(2) * (sc.special.erfinv((2 * upper_cpd) - 1) - sc.special.erfinv((2 * lower_cpd) - 1))
#     prior_sig = numerator / denominator
#     prior_mu = np.log(upper_val) - (math.sqrt(2) * prior_sig) * sc.special.erfinv((2 * upper_cpd) - 1)
#     return [prior_sig, prior_mu]

