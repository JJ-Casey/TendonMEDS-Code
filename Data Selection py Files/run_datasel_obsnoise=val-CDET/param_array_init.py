import numpy as np
import scipy as sc
import math
import dist_functions as distf
import stress_functions_symm as stressf
import general_functions as f

def import_data(fnametxt):
    # import stress-strain data from text file
    data = np.genfromtxt(fname=fnametxt)
    strain = data[:, 0]
    stress = data[:, 1]
    return [strain, stress]

def parameters(num_data, strain, stress, chain):
    class prior_params():  #structure to store all parameters relating to the priors (or hyper priors within model)
        pass
    # model parameter prior parameters
    [nu_sig, nu_mu] = [1.30927056, 1.05309738]  # (1-phi)mu
    [eta_sig, eta_mu] = [0.47191773, 6.83672018]  # phi E
    [tau_sig, tau_mu] = [0.64387023, -3.80045123]  # a - 1
    [rho_sig, rho_mu] = [0.7310165, -3.59771868]  # b - a
    prior_params.prior_mu = np.asarray([nu_mu, eta_mu, tau_mu, rho_mu])  # model parameter prior mu vals
    prior_params.prior_sig = np.asarray([nu_sig, eta_sig, tau_sig, rho_sig])  # model parameter prior sigma vals

    # fidelity prior parameters
    lengthscale = 0.05
    mu = sigmoid(strain, x0=1.1, rate=50, leftasymp=4, rightasymp=1)
    sigma = sigmoid(strain, x0=1.1, rate=50, leftasymp=0.25, rightasymp=1)

    [prior_params] = get_fidprior_params(lengthscale, mu, sigma, np.expand_dims(strain, 1), num_data, prior_params)

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

    # fidelity start
    current_vals.trans_gamma = f.quick_MVN_sample(prior_params.gamma_prior_mu, prior_params.CholDecomp_cov, num_data)
    current_vals.gamma = f.conv_to_gamma(current_vals.trans_gamma)
    current_vals.gamma = f.check_min_max_fid(current_vals.gamma, 1e-16)
    current_vals.trans_gamma = f.conv_to_transgamma(current_vals.gamma)

    # evaluate model data using these initial values
    current_vals.model_data = stressf.choose_stress_model(current_vals.theta, strain)

    # likelihood with chosen value of noise variance
    # choose observational noise variance
    current_vals.obsnoise = 0.15**2
    print('Observational noise variance = %.3f' % current_vals.obsnoise)
    current_vals.likelihood = distf.likelihood_gauss(stress, current_vals.model_data, current_vals.gamma, current_vals)


    #prior densities
    current_vals.prior_transtheta = distf.prior_transtheta(prior_params.prior_mu, prior_params.prior_sig, current_vals.trans_theta)
    current_vals.prior_transgamma = distf.prior_transgamma(prior_params.gamma_prior_mu, prior_params.inv_cov, current_vals.trans_gamma)

    return [prior_params, current_vals, prop_cov]

def MCMC_parameters():
    class MCMC_params():
        pass
    # MCMC parameters
    MCMC_params.N = int(5e6)
    MCMC_params.N_burnin = int(MCMC_params.N / 2) # int(1e6) # int(5e5)  # change to 5e5
    MCMC_params.N_write = int(5000)
    MCMC_params.num_samp_cov = 10000  # number of samples between proposal covariance updates
    MCMC_params.num_param = 4 # number of model parameters
    MCMC_params.beta_theta = 1  # RWMH step size for model module
    MCMC_params.beta_gamma = 0.5 # 0.1  # RWMH step size for fidelity module

    MCMC_params.beta_sigmasq = 0.001

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
    arrays.alpha_gamma = np.empty(MCMC_params.check_freq)
    arrays.acc_rate_theta = np.empty(int(MCMC_params.N / MCMC_params.check_freq))
    arrays.acc_rate_gamma = np.empty(int(MCMC_params.N / MCMC_params.check_freq))
    arrays.beta_theta = np.empty(int(MCMC_params.N / MCMC_params.check_freq))
    arrays.beta_gamma = np.empty(int(MCMC_params.N / MCMC_params.check_freq))
    arrays.nlog_post = np.empty(int(MCMC_params.N / MCMC_params.check_freq))
    arrays.nlog_likelihood = np.empty(int(MCMC_params.N / MCMC_params.check_freq))

    arrays.trans_theta_proposal = np.empty([int(MCMC_params.num_samp_cov), MCMC_params.num_param], dtype=np.float32)
    arrays.trans_theta = np.empty([int(MCMC_params.N_write), MCMC_params.num_param], dtype=np.float32)
    arrays.theta = np.empty([int(MCMC_params.N_write), MCMC_params.num_param], dtype=np.float32)
    arrays.prop_cov = []
    arrays.gamma = np.empty([int(MCMC_params.N_write), num_data], dtype=np.float32)
    arrays.model_fit = np.empty([int(MCMC_params.N_write), num_data], dtype=np.float32)

    return [arrays]

def get_fidprior_params(l, mu, sigma, X, num_data, prior_params):
    prior_params.gamma_prior_mu = mu
    cov = gen_gamma_cov(X, X, l, sigma)
    # check our kernel is positive semi-definite
    while is_pos_def(cov) == False:
        # add on minimum eigenvalue
        cov += 1e-10 * np.identity(num_data)  # added diagonal must be v. small due to small eigenvalue terms
    prior_params.gamma_prior_cov = cov
    prior_params.inv_cov = np.linalg.inv(cov)
    prior_params.CholDecomp_cov = np.linalg.cholesky(cov)
    return [prior_params]

def get_prior_params(upper_val, lower_val, upper_cpd, lower_cpd):
    # calculate the prior means and standard deviaitions given the conditions in Table 1 in the Supplementary
    # these are the means and std of the normal priors on the transformed parameters
    numerator = np.log(upper_val) - np.log(lower_val)
    denominator = math.sqrt(2) * (sc.special.erfinv((2 * upper_cpd) - 1) - sc.special.erfinv((2 * lower_cpd) - 1))
    prior_sig = numerator / denominator
    prior_mu = np.log(upper_val) - (math.sqrt(2) * prior_sig) * sc.special.erfinv((2 * upper_cpd) - 1)
    return [prior_sig, prior_mu]

def gen_gamma_cov(x1, x2, l, sigma):
    sq_norm = -0.5 * sc.spatial.distance.cdist(x1, x2, 'sqeuclidean')
    C = np.zeros([len(x1), len(x1)])
    for i in range(len(x1)):
        for j in range(len(x1)):
            C[i, j] = sigma[i] * sigma[j] * np.exp(sq_norm[i, j] / (l ** 2))

    return C

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

def sigmoid(x, x0, rate, leftasymp, rightasymp):
    # x0 is the x coord of the mid point of the range of y
    return leftasymp + np.divide(rightasymp-leftasymp, 1 + np.exp(-np.multiply(rate, (x-x0))))



