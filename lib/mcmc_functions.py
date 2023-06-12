import numpy as np
import emcee
import corner
import pickle
import base64
import matplotlib.pyplot as plt
import scipy.interpolate as scp

import physical_functions as pf


def mcmc_full_run(freq, flux, fluxerr, fitparams, priors, nwalkers, niter, hold):

    # Define the initial theta array
    theta = np.array([ fitparams['mb_mass'], fitparams['mb_temp'], fitparams['mb_beta'] ])

    # Define stepping methodology
    ndim = len(theta)
    p0 = [theta+1e-7*np.random.randn(ndim) for i in range(nwalkers)]

    # Run emcee
    data = (freq, flux, fluxerr, fitparams['distance'], priors, hold)
    sampler, pos, prob, state = run_samples(p0, nwalkers, niter, ndim, lnprob, data)

    # Save relevant data
    ## Get sample data
    samples = sampler.get_chain(flat=True, discard=500)
    ## Find the most probable datum
    theta_max = samples[np.argmax(sampler.flatlnprobability)].copy()
    ## Draw some random samples and get statistics
    draw = np.floor(np.random.uniform(0,len(samples),size=200)).astype(int)
    theta_dist = samples[draw]
    mass_spread = np.std(theta_dist[:,0],axis=0)
    temp_spread = np.std(theta_dist[:,1],axis=0)
    beta_spread = np.std(theta_dist[:,2],axis=0)

    # Save the posterior distribution as a function of frequency
    basis = np.logspace(9, 14, num=500)
    models = []
    for itr in theta_dist:
        models.append(pf.mb_model(np.append(itr,fitparams['distance']), 1e-9*basis))
    model_spread = np.std(models,axis=0)
    model_median = np.median(models,axis=0)
    ## Encode posterior distribution to be saved as a string
    pos_dist_enc = base64.b64encode(pickle.dumps(np.array([basis,model_median,model_spread]))).decode('utf-8')

    # Show corner plot for reference
    cp_labels = ['Dust Mass ($M_{\odot}$)','Temp (K)','$\\beta$']
    ## Mask held variables for corner plot
    cp_labels = np.array([lab for lab,h in zip(cp_labels,hold) if h is np.inf])
    cp_samples = np.array([ [s for s,h in zip(samp,hold) if h is np.inf] for samp in samples])
    fig = corner.corner(
            cp_samples, 
            labels=cp_labels,
            quantiles=[0.16,0.5,0.84],
            show_titles=True, 
            plot_datapoints=True
            )
    plt.show()

    # Return
    return {
            "mass" : theta_max[0],
            "mass_spread" : mass_spread,
            "temp" : theta_max[1],
            "temp_spread" : temp_spread,
            "beta" : theta_max[2],
            "beta_spread" : beta_spread,
            "posterior_spread_obj" : pos_dist_enc
            }


def run_samples(p0, nwalkers, niter, ndim, lnprob, data):

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data)

    pos, prob, state = sampler.run_mcmc(p0, niter+500)

    return sampler, pos, prob, state


def lnlike(theta,x,y,yerr,distance):

    LnLike = -0.5*np.sum(((y-pf.mb_model(np.append(theta,distance),x))/yerr)**2)

    return LnLike


def lnprior(theta, priors):

    prior_test = [
            theta[0] < priors['mb_mass_uplim'],
            theta[0] > priors['mb_mass_lolim'],
            theta[1] < priors['mb_temp_uplim'],
            theta[1] > priors['mb_temp_lolim'],
            theta[2] < priors['mb_beta_uplim'],
            theta[2] > priors['mb_beta_lolim']
        ]

    if np.all(prior_test):
        return 0.0
    else:
        return -np.inf

def lnprob(theta,x,y,yerr,distance,priors,hold):

    # Hold parameters if necessary
    for itr in range(len(hold)):
        if hold[itr] is not np.inf:
            theta[itr] = hold[itr]

    lp = lnprior(theta, priors)

    if not lp == 0.0:
        return -np.inf
    else:
        return lp + lnlike(theta,x,y,yerr,distance)
