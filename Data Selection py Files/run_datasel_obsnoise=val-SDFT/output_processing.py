import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import general_functions as f

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

    # gamma parameters
    data = np.loadtxt(os.path.join(tosaveto, 'gamma.txt'), dtype=np.float32)

    reported_vals.gamma_lowCI = np.round(np.percentile(data, 2.5, axis=0), 4)
    reported_vals.gamma_highCI = np.round(np.percentile(data, 97.5, axis=0), 4)

    CI95_intervals = [reported_vals.gamma_lowCI, reported_vals.gamma_highCI]
    np.savetxt(os.path.join(tosaveto, '95CI_gamma_params.txt'), CI95_intervals)

    reported_vals.gamma_mean = np.mean(data, axis=0)
    reported_vals.gamma_median = np.median(data, axis=0)
    np.savetxt(os.path.join(tosaveto, 'mean_gamma_params.txt'), reported_vals.gamma_mean)
    np.savetxt(os.path.join(tosaveto, 'median_gamma_params.txt'), reported_vals.gamma_median)

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
    fig2, ax2 = plt.subplots(2, figsize=(8,5))
    itvals = np.linspace(1, int(MCMC_params.N / MCMC_params.check_freq), int(MCMC_params.N / MCMC_params.check_freq)) * MCMC_params.check_freq
    ax2[0].plot(itvals, arrays.acc_rate_theta, color='m', linewidth=0.75)
    ax2[0].axhline(y=MCMC_params.upper_tol, linestyle='dashed', color='k')
    ax2[0].axhline(y=MCMC_params.lower_tol, linestyle='dashed', color='k')
    ax2[0].set_xticks([])
    ax2[0].set_ylabel('Acceptance Rate', color='m')
    ax2[0].set_title('Model')
    ax2b = ax2[0].twinx()
    ax2b.plot(itvals, arrays.beta_theta, color='#666699', linewidth=0.75)
    ax2b.set_ylabel('\u03B2 (step-size)', color='#666699')  # beta

    ax2[1].plot(itvals, arrays.acc_rate_gamma, color='m', linewidth=0.75)
    ax2[1].axhline(y=MCMC_params.upper_tol, linestyle='dashed', color='k')
    ax2[1].axhline(y=MCMC_params.lower_tol, linestyle='dashed', color='k')
    ax2[1].set_xlabel('Iteration')
    ax2[1].set_ylabel('Acceptance Rate', color='m')
    ax2[1].set_title('Data selection (gamma)')
    ax2b2 = ax2[1].twinx()
    ax2b2.plot(itvals, arrays.beta_gamma, color='#666699',  linewidth=0.75)
    ax2b2.set_ylabel('\u03B2 (step-size)', color='#666699')  # beta

    fig2.align_ylabels(ax2)
    fig2.align_ylabels(ax2b2)
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
    for i in range(num_param):
        for j in range(num_param):
            axs[i, j].plot(prop_cov_array[:, count], linewidth=0.75)
            count += 1
    plt.suptitle('Proposal covariance')
    fig.tight_layout()
    plt.savefig(os.path.join(tosaveto,'covariance_entries.png'), dpi=600)  # windows

def gamma_heatmaps(num_data, tosaveto, reported_vals):
    gamma_array = np.loadtxt(os.path.join(tosaveto, 'gamma.txt'), dtype=np.float32)

    fig11, ax11a = plt.subplots(1, figsize=(9, 5))
    fig11.set_figheight(4)
    fig11.set_figwidth(5)

    gamma_hist_array = np.empty([101, num_data])
    for i in range(num_data):
        gamma_hist_array[:, i], edges = np.histogram(gamma_array[:, i], bins=101,
                                                     density=True, range=(0, 1.0))
    # rescale fid_hist_array
    gamma_hist_array_max = np.amax(gamma_hist_array, axis=0)  # max value in each row
    gamma_hist_array_norm = gamma_hist_array / np.reshape(gamma_hist_array_max, [1, num_data])  # max val in each row =1
    dataframe_gamma = pd.DataFrame(np.flip(gamma_hist_array_norm, axis=0),
                                   index=np.round(np.linspace(1, 0, 101), 2),
                                   columns=np.round(np.linspace(1, num_data, num_data, dtype=int), 0))

    hm1 = sns.heatmap(dataframe_gamma, yticklabels=10, xticklabels=round(num_data / 5), cmap='plasma', ax=ax11a)
    hm1.set_ylabel('Gamma (\u03b3)')
    hm1.set_xlabel('Data Point')

    ax11a2 = ax11a.twinx()
    ax11a2.plot(np.linspace(1, num_data, num_data), reported_vals.gamma_mean, linestyle='solid', color='k', linewidth=0.75, label='\u03b3 mean')
    ax11a2.plot(np.linspace(1, num_data, num_data), reported_vals.gamma_median, linestyle='dashed', color='k', linewidth=0.75,
              label='\u03b3 median')
    ax11a2.fill_between(np.linspace(1, num_data, num_data), reported_vals.gamma_lowCI, reported_vals.gamma_highCI, color='k', label='95% CI', alpha=0.01, edgecolor='k', linewidth=1.0, linestyle='dotted')
    ax11a2.set_ylim([0, 1])
    ax11a2.set_yticks([])
    ax11a2.legend()

    fig11.tight_layout()

    plt.savefig(os.path.join(tosaveto, 'gamma_heatmap.png'), dpi=600)  # windows

    del gamma_array


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
    axa2.plot(strain, reported_vals.gamma_mean, linestyle='solid', color='#8b0be0', linewidth=0.75, label='\u03b3 mean')
    axa2.plot(strain, reported_vals.gamma_median, linestyle='dashed', color='#e00bbc', linewidth=0.75, label='\u03b3 median')
    axa2.fill_between(strain, reported_vals.gamma_lowCI, reported_vals.gamma_highCI, color='#dbc0fa', label='95% CI', alpha=0.4)
    axa2.set_ylim([0,1])
    axa2.set_ylabel('Mean/Median Fidelity', color='k')

    axa.legend()
    axa2.legend()
    fig.tight_layout()

    plt.savefig(os.path.join(tosaveto,'modelfit_withfid.png'), dpi=600)  # windows

def trace_plots_gamma(MCMC_params, tosaveto):
    gamma_array_tmp = np.loadtxt(os.path.join(tosaveto,'gamma.txt'), dtype=np.float32)
    itall = np.linspace(MCMC_params.N_burnin, MCMC_params.N, int((MCMC_params.N-MCMC_params.N_burnin)/MCMC_params.thin))
    fig4, axs4 = plt.subplots(5, figsize=(9,5))
    for i in range(5):
        axs4[i].plot(itall, gamma_array_tmp[:, i], linewidth=0.5)

    axs4[i].set_xlabel('Iteration')
    fig4.suptitle('Example fidelity trace plots')

    plt.savefig(os.path.join(tosaveto,'traceplots_fid.png'), dpi=600)  # windows

    del gamma_array_tmp

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
	
    avg_nlog_post = np.mean(arrays.nlog_post[len(arrays.nlog_post)//2 :]) 
    np.savetxt(os.path.join(tosaveto,'avg_nlog_posterior.txt'), [avg_nlog_post])
