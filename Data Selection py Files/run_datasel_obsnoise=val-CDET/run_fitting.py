# This script will be the central script that calls all other elements of the code
# GT version of model
# JFORSYTH OCTOBER 2022

def run(chain, datafile, file):

    import numpy as np
    import general_functions as f
    import param_array_init as ini
    import matplotlib.pyplot as plt
    import output_processing as outp
    import dist_functions as distf
    import random
    import os
    import shutil
    import sys
    from datetime import datetime

    # close all figures
    plt.close('all')
    # raise all errors
    # np.seterr(all='raise')
    # set random seed (default is time - no good, instead base on chain number)
    np.random.seed(chain)

    # import data here : if testing multiple data sets, this can be done according to chain numbers here
    [strain, stress] = ini.import_data(datafile)

    num_data = len(strain)
    ###################################################################################################################
    # load parameters
    ###################################################################################################################
    [prior_params, current_vals, prop_cov] = ini.parameters(num_data, strain, stress, chain)

    [MCMC_params] = ini.MCMC_parameters()

    [arrays] = ini.arrays(MCMC_params, num_data)

    ###################################################################################################################
    # Create directories to save to
    ###################################################################################################################
    now = datetime.now()
    fname = now.strftime("_%d%m%Y_%H-%M_")
    fname = os.path.join('results', file[0:len(file) - 4] + fname + 'ch' + str(chain))
    try:
        os.mkdir(fname)
    except:
        fname = fname + '_2'
        try:
            os.mkdir(fname)
            print('New directory made with _2.')
        except:
            print('New directory not made, please remove existing directories.')
    tosaveto = os.path.join(os.getcwd(), fname)

    # make a copy of parameters file within tosaveto
    paramfile = os.path.join(tosaveto, 'parameters.txt')
    shutil.copyfile('param_array_init.py', paramfile)
    np.savetxt(os.path.join(tosaveto, 'strain.txt'), strain)
    np.savetxt(os.path.join(tosaveto, 'stress.txt'), stress)

    file = open(os.path.join(tosaveto, 'running_log.txt'), 'a')
    nowtime = now.strftime("_%d%m%Y_%H-%M-%S")
    file.write('Chain %d --- %s , %d data points\n' % (chain, datafile, num_data))
    file.write('Chain %d --- 0 %% @ %s\n' % (chain, nowtime))
    file.close()

    print('Chain %d --- loaded and MCMC initialised.' % chain)

    ###################################################################################################################
    # perform MCMC
    ###################################################################################################################
    for it in range(1, MCMC_params.N + 1):
        if it % (MCMC_params.N / 20) == 0:
            file = open(os.path.join(tosaveto, 'running_log.txt'), 'a')
            # update log file
            now = datetime.now()
            nowtime = now.strftime("_%d%m%Y_%H-%M-%S")
            file.write('Chain %d --- %d %% @ %s\n' % (chain, 100 * (it / MCMC_params.N), nowtime))
            file.close()

        # sample on model parameters here
        [prop_cov, current_vals, MCMC_params, arrays] = f.sample_theta(it, prop_cov, strain,stress, num_data,
                                                                       prior_params, current_vals, arrays, MCMC_params)
        #choose value of observational noise variance
        current_vals.sigmasq = 1.0

        # Now sample on the fidelity terms
        [current_vals, MCMC_params, arrays] = f.sample_gamma(it, num_data, stress, prior_params, MCMC_params,
                                                             current_vals, arrays)

        # Write data to arrays and text files
        arrays = f.write_data(it, MCMC_params, current_vals, arrays, tosaveto)

        # Append the negative log of the posterior and likelihood to array
        if (it % MCMC_params.check_freq) == 0:
            idx = int((it / MCMC_params.check_freq) - 1)
            nlog_post = distf.get_nlog_posterior(current_vals)
            arrays.nlog_post[idx] = nlog_post
            arrays.nlog_likelihood[idx] = current_vals.likelihood


    ###################################################################################################################
    # estimate parameter values to report + 95% CI limits
    ###################################################################################################################
    reported_vals = outp.get_parameter_vals_CI(tosaveto)
    ###################################################################################################################
    # make plots
    ###################################################################################################################
    ylabellisttrans = ['\u03BD', '\u03B7', '\u03C4', '\u03C1']  # 'nu','eta', 'tau', 'rho'
    ylabellisttrue = ['(1-\u03C6)\u03BC', '\u03C6 E', 'a', 'b']  # '(1-phi)mu', 'phiE', 'a', 'b'

    outp.plot_accr_beta(tosaveto, MCMC_params, arrays)
    outp.trace_plots(MCMC_params, ylabellisttrue, ylabellisttrans, tosaveto, reported_vals)
    outp.model_true_hist(ylabellisttrue, tosaveto, reported_vals, MCMC_params.num_param)
    outp.gamma_heatmaps(num_data, tosaveto, reported_vals)
    outp.plot_nlog_post_like(tosaveto, MCMC_params, arrays)
    outp.model_fit(strain, stress, reported_vals, tosaveto)

    if len(arrays.prop_cov) > 1:
        outp.proposal_variance_plot(arrays.prop_cov, MCMC_params.num_param, tosaveto)

    ###################################################################################################################
    # remove saved text files to save memory
    ###################################################################################################################
    os.remove(os.path.join(tosaveto, 'theta_trans.txt'))
    os.remove(os.path.join(tosaveto, 'gamma.txt'))
    os.remove(os.path.join(tosaveto, 'modelfit.txt'))
    os.remove(os.path.join(tosaveto, 'theta.txt'))

    print('Chain %d --- finished.' % chain)

    # plt.show()

    return
