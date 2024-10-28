import os
import numpy as np
import matplotlib.pyplot as plt
import general_functions as f
from scipy.stats import norm
from sklearn.metrics import r2_score


def get_parameter_vals_CI(tosaveto):
    class reported_vals():
        pass
    # model parameters
    data = np.loadtxt(os.path.join(tosaveto, 'theta.txt'), dtype=np.float32)
    data = f.thetadata_to_abc(data)

    reported_vals.theta_lowCI = np.round(np.percentile(data, 2.5, axis=0), 4)
    reported_vals.theta_highCI  = np.round(np.percentile(data, 97.5, axis=0), 4)

    CI95_intervals = [reported_vals.theta_lowCI, reported_vals.theta_highCI ]
    np.savetxt(os.path.join(tosaveto, '95CI_model_params.txt'), CI95_intervals)

    reported_vals.theta_mean = np.mean(data, axis=0)
    reported_vals.theta_median = np.median(data, axis=0)
    np.savetxt(os.path.join(tosaveto, 'mean_theta_params.txt'), reported_vals.theta_mean)
    np.savetxt(os.path.join(tosaveto, 'median_theta_params.txt'), reported_vals.theta_median)

    del data

    # model fit
    data = np.loadtxt(os.path.join(tosaveto, 'modelfit.txt'), dtype=np.float32)
    reported_vals.modelfit_lowCI = np.percentile(data, 2.5,axis=0)
    reported_vals.modelfit_highCI = np.percentile(data, 97.5, axis=0)

    reported_vals.modelfit_mean = np.mean(data, axis=0)
    reported_vals.modelfit_median = np.median(data, axis=0)

    del data

    return reported_vals

def plot_accr_beta(tosaveto, MCMC_params, arrays):
    # alpha and beta
    fig2, ax2 = plt.subplots(1, figsize=(8,3))
    itvals = np.linspace(1, int(MCMC_params.N / MCMC_params.check_freq), int(MCMC_params.N / MCMC_params.check_freq)) * MCMC_params.check_freq
    ax2.plot(itvals, arrays.acc_rate_theta, color='m', linewidth=0.75)
    ax2.axhline(y=MCMC_params.upper_tol, linestyle='dashed', color='k')
    ax2.axhline(y=MCMC_params.lower_tol, linestyle='dashed', color='k')
    ax2.set_xticks([])
    ax2.set_ylabel('Acceptance Rate', color='m')
    ax2.set_title('Model')
    ax2b = ax2.twinx()
    ax2b.plot(itvals, arrays.beta_theta, color='#666699', linewidth=0.75)
    ax2b.set_ylabel('\u03B2 (step-size)', color='#666699')  # beta
    ax2.set_xlabel('Iteration')
    fig2.tight_layout(pad=0.3)

    plt.savefig(os.path.join(tosaveto,'accr_beta.png'), dpi=600)  # windows

def trace_plots(MCMC_params, ylabellisttrue, ylabellisttrans, tosaveto, reported_vals):
    model_true_array_tmp = np.loadtxt(os.path.join(tosaveto, 'theta.txt'), dtype=np.float32)
    model_true_array = np.copy(model_true_array_tmp)
    model_true_array = f.thetadata_to_abc(model_true_array)

    model_trans_array = np.loadtxt(os.path.join(tosaveto,'theta_trans.txt'), dtype=np.float32)

    itall = np.linspace(MCMC_params.N_burnin, MCMC_params.N, int((MCMC_params.N-MCMC_params.N_burnin)/MCMC_params.thin))
    fig4, axs4 = plt.subplots(MCMC_params.num_param, figsize=(8,5))

    for i in range(MCMC_params.num_param):
        axs4[i].plot(itall, model_true_array[:, i], color= '#006699', linewidth=0.5)
        axs4[i].set_ylabel(ylabellisttrue[i],color='#006699')
        axs4[i].tick_params(axis='y', colors='#006699')

        ax4_2 = axs4[i].twinx()
        ax4_2.plot(itall, model_trans_array[:, i], color='#ff9933', linewidth=0.5)
        ax4_2.set_ylabel(ylabellisttrans[i], color='#ff9933')
        ax4_2.tick_params(axis='y', colors='#ff9933')

        axs4[i].axhline(y=reported_vals.theta_mean[i], linestyle='solid', color='k')
        axs4[i].axhline(y=reported_vals.theta_median[i], linestyle='dashed', color='#e06b0b')

        if i != MCMC_params.num_param - 1:
            axs4[i].tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=False,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelbottom=False)  # labels along the bottom edge are off
    axs4[i].set_xlabel('Iteration')

    fig4.align_ylabels(axs4)
    fig4.tight_layout(pad=0.3)

    plt.savefig(os.path.join(tosaveto,'traceplots.png'), dpi=600)  # windows

    del model_true_array
    del model_trans_array
    del model_true_array_tmp

def model_true_hist(ylabellisttrue, tosaveto, reported_vals, num_param):
    data = np.loadtxt(os.path.join(tosaveto, 'theta.txt'), dtype=np.float32)
    data = f.thetadata_to_abc(data)

    fig1, ax1 = plt.subplots(ncols=num_param, nrows=num_param, layout='tight', figsize=(10, 8))

    # plot histograms along diagonal
    for i in range(num_param):
        ax1[i, i].hist(data[:, i], bins=101, density=True, color='#6699ff')
        ax1[i, i].set_xlim(min(data[:, i]), max(data[:,i]))

        #plot mean and median model parameter values
        ax1[i, i].axvline(x=reported_vals.theta_mean[i], linestyle='solid', color='k')
        ax1[i, i].axvline(x=reported_vals.theta_median[i], linestyle='dashed', color='#e06b0b')

        if i == num_param - 1:
            ax1[i, i].set_xticks([reported_vals.theta_lowCI[i], reported_vals.theta_highCI[i]])
        else:
            ax1[i, i].set_xticks([reported_vals.theta_lowCI[i], reported_vals.theta_highCI[i]], labels=[' ', ' '])

    # Correlation heatmaps
    for irow in range(num_param):
        for icol in range(num_param):
            if icol == 0:
                ax1[irow, icol].set_ylabel(ylabellisttrue[irow], fontsize=14, fontweight='bold')
            if irow == 0:
                ax1[irow, icol].set_title(ylabellisttrue[icol], fontsize=14, fontweight='bold')

            if icol != irow:
                if icol < irow:
                    ax1[irow, icol].hist2d(data[:, icol], data[:, irow], bins=[101, 101], density=True)
                    ax1[irow, icol].set_xlim(min(data[:, icol]), max(data[:, icol]))
                    ax1[irow, icol].set_ylim(min(data[:, irow]), max(data[:, irow]))

                if icol == 0:
                    ax1[irow, icol].set_yticks([reported_vals.theta_lowCI[irow], reported_vals.theta_highCI[irow]])
                else:
                    ax1[irow, icol].set_yticks([reported_vals.theta_lowCI[irow], reported_vals.theta_highCI[irow]], labels=[' ', ' '])

                if irow == num_param - 1:
                    ax1[irow, icol].set_xticks([reported_vals.theta_lowCI[icol], reported_vals.theta_highCI[icol]])
                else:
                    ax1[irow, icol].set_xticks([reported_vals.theta_lowCI[icol], reported_vals.theta_highCI[icol]], labels=[' ', ' '])

                if icol > irow:
                    corr_array = np.corrcoef(data[:, icol], data[:, irow])
                    txt = 'Corr.: \n %.4f' % corr_array[0, 1]
                    ax1[irow, icol].text(0.5, 0.5, txt, horizontalalignment='center', verticalalignment='center',
                                         transform=ax1[irow, icol].transAxes, fontsize=11)
                    ax1[irow, icol].set_axis_off()

    fig1.align_ylabels(ax1)
    fig1.tight_layout(pad=0.3)

    plt.savefig(os.path.join(tosaveto, 'true_param_hist.png'), dpi=600)  # windows

    del data

def proposal_variance_plot(prop_cov_array, num_param, tosaveto):
    fig, axs = plt.subplots(num_param, num_param, figsize=(9, 9))
    prop_cov_array = np.array(prop_cov_array)
    count=0
    for i in range(num_param ):
        for j in range(num_param):
            axs[i, j].plot(prop_cov_array[:, count], linewidth=0.75)
            count += 1
    plt.suptitle('Proposal covariance')
    fig.tight_layout()
    plt.savefig(os.path.join(tosaveto,'covariance_entries.png'), dpi=600)  # windows

def model_fit(strain, stress, reported_vals, tosaveto):

    fig, axa = plt.subplots(1, figsize=(7, 5))
    # plot a, c, b as vertical dotted lines
    axa.axvline(x=reported_vals.theta_median[2], linestyle='dotted', color='k', linewidth=0.5, label='a, b')
    axa.axvline(x=reported_vals.theta_median[3], linestyle='dotted', color='k', linewidth=0.5)

    axa.plot(strain, stress, linestyle='-', label='Data', c='#ffcc00')
    axa.plot(strain, reported_vals.modelfit_mean, color='k', label='Model fit', linestyle='dashed')
    axa.fill_between(strain, reported_vals.modelfit_lowCI, reported_vals.modelfit_highCI, color='k', edgecolor='k', linewidth=0.2, label='95% CI', alpha=0.4)
    axa.set_xlabel('Strain')
    axa.set_ylabel('Stress')
    axa2 = axa.twinx()
    axa2.plot(strain, np.ones(len(strain)), linestyle='solid', color='#8b0be0', linewidth=0.75, label='\u03b3 mean')
    #axa2.plot(strain, np.ones(len(strain)), linestyle='dashed', color='#e00bbc', linewidth=0.75, label='\u03b3 median')
    #axa2.fill_between(strain, reported_vals.gamma_lowCI, reported_vals.gamma_highCI, color='#dbc0fa', label='95% CI', alpha=0.4)
    axa2.set_ylim([0,1])
    axa2.set_ylabel('Fidelity', color='k')

    axa.legend()
    axa2.legend()
    fig.tight_layout()

    plt.savefig(os.path.join(tosaveto,'modelfit.png'), dpi=600)  # windows

def plot_nlog_post_like(tosaveto, MCMC_params, arrays):
    # alpha and beta
    fig2, ax2 = plt.subplots(1, figsize=(8,3))
    itvals = np.linspace(1, int(MCMC_params.N / MCMC_params.check_freq), int(MCMC_params.N / MCMC_params.check_freq)) * MCMC_params.check_freq
    ax2.plot(itvals, arrays.nlog_post, color='m', linewidth=0.75)
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Negative log posterior', color='m')
    ax2b = ax2.twinx()
    ax2b.plot(itvals, arrays.nlog_likelihood, color='#666699', linewidth=0.75)
    ax2b.set_ylabel('Negative log likelihood', color='#666699')  # beta

    plt.savefig(os.path.join(tosaveto,'nlog_post_likelihood.png'), dpi=600)  # windows

def obs_noise(tosaveto):
    obs_noise = np.loadtxt(os.path.join(tosaveto, 'obs_noise.txt'), dtype=np.float32)
    mean_obs_noise = np.mean(obs_noise)
    p25 = np.percentile(obs_noise, 2.5)
    p975 = np.percentile(obs_noise, 97.5)
    fig, ax = plt.subplots()
    ax.hist(obs_noise, bins=101, density=True, color='#6699ff')
    lbl = 'mean' + str(mean_obs_noise)
    ax.axvline(x=mean_obs_noise, linestyle='dashed', color='k', label=lbl)
    ax.axvline(x=p25, linestyle='dotted', color='c', label='2.5th percentile')
    ax.axvline(x=p975, linestyle='dotted', color='b', label='97.5th percentile')
    ax.set_xlim(min(obs_noise), max(obs_noise))
    ax.set_xlabel('Observational noise variance')
    ax.set_ylabel('Density')
    ax.legend()

    plt.savefig(os.path.join(tosaveto, 'obs_noise.png'), dpi=600)  # windows

    # save mean and 95% CI to text file
    vals = [mean_obs_noise, p25, p975]
    np.savetxt(os.path.join(tosaveto, 'obsnoise_mean_95CI.txt'), vals)

def residual_analysis(tosaveto, stress, reported_vals, strain, strain_lowlim, strain_uplim, file_ext):
    to_include = (strain>strain_lowlim) & (strain<strain_uplim)
    num = sum(to_include)
    modelfits = np.loadtxt(os.path.join(tosaveto, 'modelfit.txt'), dtype=np.float32)
    modelfits = modelfits[:,to_include]
    experimental_data = np.tile(stress[to_include], (len(modelfits),1))
    data = experimental_data - modelfits

    fig, ax = plt.subplots((2))
    ax[0].hist(np.matrix.flatten(data), density=True, bins=101)
    # fit normal distribution to residuals
    obs_noise_mu, obs_noise_std = norm.fit(np.matrix.flatten(data), loc=0, scale=0.1)
    # calculate R-squared
    true_hist, bin_edges = np.histogram(np.matrix.flatten(data), density=True, bins=101)
    norm_model_hist = norm.pdf(bin_edges, loc=obs_noise_mu, scale=obs_noise_std)
    r2_value = r2_score(true_hist, norm_model_hist[1:])
    # strain dependent residuals
    mu_vals=np.zeros((num,1))
    std_vals=np.zeros((num,1))
    for i in range(num):
        mu_vals[i], std_vals[i] = norm.fit(data[:, i], loc=0,scale=0.1)
    ax[1].plot(strain[to_include], mu_vals,'.g', label='fitted mu vals')
    plt.legend()
    ax1b = ax[1].twinx()
    ax1b.plot(strain[to_include], std_vals,'.m', label='fitted std vals')
    plt.legend()
    vals = [obs_noise_mu, obs_noise_std, r2_value]
    fname = 'obsnoise_normal_fit_mu_std' + file_ext + '.txt'
    np.savetxt(os.path.join(tosaveto, fname), vals)

    x = np.linspace(np.min(data), np.max(data),101)
    y = norm.pdf(x,loc=obs_noise_mu, scale=obs_noise_std)
    ax[0].plot(x,y,'k-')
    ax[0].set_xlabel('Residual')
    ax[0].set_ylabel('Density')
    ax[0].set_title('Model versus data residuals')
    ax[1].set_title('Fitted mu and std w/r to strain')
    plt.tight_layout()
    fname = 'residuals' + file_ext + '.png'
    plt.savefig(os.path.join(tosaveto, fname), dpi=600)  # windows

    del data