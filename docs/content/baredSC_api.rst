Use baredSC within python
=========================

.. contents::
    :local:

This part is not a proper documentation of baredSC as a package as we are missing a lot of documentation of our methods.
It just gives indications on how it can be used.

baredSC_1d
----------

Use the ``npz`` content
^^^^^^^^^^^^^^^^^^^^^^^

baredSC is mainly used with the command line interface.
By default, only the numpy compressed ``.npz`` file is output,
but if the ``--figure`` argument is used, it outputs much more.
All the outputs are described in the :doc:`outputs` page.

We provide two examples of plots using the text outputs in :doc:`tuto/tutorial_sim_compare_means` and :doc:`tuto/tutorial_sim_custom_simple_plot`.
However, sometimes the user may want to plot information which is not part of the text outputs.
We describe here a script which will use the baredSC_1d output to plot a custom error bars.

First we get the value of parameters at each step of MCMC:

.. ipython::

    In [1]: import numpy as np
       ...: import baredSC.oned
       ...: import matplotlib.pyplot as plt
       ...: 
       ...: output_baredSC = "../example/first_example_1d_2gauss.npz"
       ...: 
       ...: # First option, use the function extract_from_npz
       ...: mu, cov, ox, oxpdf, x, logprob_values, samples = \
       ...:   baredSC.oned.extract_from_npz(output_baredSC)
       ...: # mu is the mean of each parameter,
       ...: # cov is the covariance of parameters
       ...: # ox is the oversampled x to compute the poisson noise pdf
       ...: # oxpdf is the oversampled x to compute the pdf
       ...: # x is the x used to compute the likelihood
       ...: # logprob_values is the value of log likelihood at each step fo the MCMC
       ...: # samples is the value of each parameter at each step of the MCMC
       ...: 
       ...: # Second option, get the samples and x, oxpdf from the npz:
       ...: mcmc_res = np.load(output_baredSC, allow_pickle=True)
       ...: samples = mcmc_res['samples']
       ...: x = mcmc_res['x']
       ...: oxpdf = mcmc_res['oxpdf']
       ...: 
       ...: # From the sample size we deduce the number of Gaussians
       ...: nnorm = (samples.shape[1] + 1) // 3
       ...: # We display the parameter names:
       ...: p_names = [f'{pn}{i}' for i in range(nnorm) for pn in ['amp', 'mu', 'scale']][1:]
       ...: print(p_names)
       ...: 

Then we compute the pdf for each step of the MCMC and we plot the pdf with the custom error bars:

.. ipython::    

    @savefig plot_from_npz1d.png
    In [2]: # We assume x is equally spaced
       ...: dx = x[1] - x[0]
       ...: nx = x.size
       ...: noxpdf = oxpdf.size
       ...: # We assume oxpdf is equally spaced
       ...: odxpdf = oxpdf[1] - oxpdf[0]
       ...: 
       ...: # Compute the pdf for each sample
       ...: # This can be long
       ...: pdf = np.array([baredSC.oned.get_pdf(p, nx, noxpdf, oxpdf, odxpdf)
       ...:                 for p in samples])
       ...: 
       ...: xmin = x[0] - dx / 2
       ...: xmax = x[-1] + dx / 2
       ...: 
       ...: my_custom_quantiles = {'0.3': [0.25, 0.75], '0.1': [0.1, 0.9]}
       ...: 
       ...: plt.figure()
       ...: for alpha in my_custom_quantiles:
       ...:   pm = np.quantile(pdf, my_custom_quantiles[alpha][0], axis=0)
       ...:   pp = np.quantile(pdf, my_custom_quantiles[alpha][1], axis=0)
       ...:   plt.fill_between(x, pm, pp, color='g', alpha=float(alpha),
       ...:   rasterized=True)
       ...: # Mean
       ...: plt.plot(x, np.mean(pdf, axis=0), 'r', lw=2, rasterized=True)
       ...: 

Run baredSC_1d
^^^^^^^^^^^^^^

You can also run the MCMC from python directly.

However, it requires formating of the input:

.. ipython::    
    In [1]: import numpy as np
       ...: import pandas as pd
       ...: from scipy.stats import lognorm, truncnorm, poisson
       ...: from baredSC.baredSC_1d import gauss_mcmc
       ...: 
       ...: # I generate 200 cells with normal expression at 1.5 with scale of 0.2
       ...: # In the Seurat scale (log(1 + 10^4 X))
       ...: n_cells = 200
       ...: cur_loc = 1.5
       ...: cur_scale = 0.2
       ...: N = lognorm.rvs(s=0.3, scale=16000, size=n_cells, random_state=1)
       ...: expression = truncnorm.rvs(- cur_loc / cur_scale, np.inf,
       ...:                            loc=cur_loc, scale=cur_scale,
       ...:                            size=n_cells,
       ...:                            random_state=2)
       ...: 
       ...: ks = poisson.rvs(mu=N * 1e-4 * (np.exp(expression) - 1),
       ...:                  random_state=3)
       ...: 
       ...: # I need to put the ks and the N in a data frame:
       ...: # The column containing the total number of UMI per cell
       ...: # must be 'nCount_RNA'
       ...: data = pd.DataFrame({'my_gene': ks, 'nCount_RNA': N})
       ...: 

Then the actual MCMC can be run with:

.. ipython::    

In [2]: results = gauss_mcmc(data=data,
   ...:                      col_gene='my_gene', # Put here the colname you put in your data
   ...:                      nx=50, # Number of bins in x
   ...:                      osampx=10, # Oversampling factor of the Poisson distribution
   ...:                      osampxpdf=5, # Oversampling factor of the PDF
   ...:                      xmin=0,
   ...:                      xmax=3,
   ...:                      min_scale=0.1, # Minimal value of the scale
   ...:                      xscale="Seurat",
   ...:                      target_sum=10000,
   ...:                      nnorm=1, # We use models with a single Gaussian
   ...:                      nsamples_mcmc=100000, # Number of steps in the MCMC
   ...:                      nsamples_burn=25000, # Number of steps in the burning phase of MCMC (we recommand nsampMCMC / 4)
   ...:                      nsplit_burn=10, # The burning phase is splitted in multiple sub-phase where the temperature is decreasing
   ...:                      T0_burn=100.0,
   ...:                      output='temp', # Where the npz output should be stored
   ...:                      seed=1)
   ...: print(len(results))
   ...: # The results are:
   ...: # mu, cov, ox, oxpdf, x, logprob_values, samples
   ...: # mu is the mean of each parameter,
   ...: # cov is the covariance of parameters
   ...: # ox is the oversampled x to compute the poisson noise pdf
   ...: # oxpdf is the oversampled x to compute the pdf
   ...: # x is the x used to compute the likelihood
   ...: # logprob_values is the value of log likelihood at each step fo the MCMC
   ...: # samples is the value of each parameter at each step of the MCMC
   ...: 
