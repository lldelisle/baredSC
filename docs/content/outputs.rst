Outputs
=======

.. contents:: 
    :local:

We will describe here each of the output file.

MCMC output
-----------

The only output by defaut is a numpy compressed ``.npz`` file. This output contains the result of the MCMC:

* samples: the value of the parameter at each step of the MCMC

  * In the 1d case, the parameters are mu0, scale0, (amp1, mu1, scale1, ...). The first amplitude can be deduced as 1 minus the sum of all other amplitudes.
  * In the 2d case, the parameters are mux0, scalex0, muy0, scaley0, corr0, (amp1, mux1, scalex1, muy1, scaley1, corr1, ..). The first amplitude can be deduced as 1 minus the sum of all other amplitudes.

* diagnostics: a dictionary with the diagnostics at each step of the MCMC, among them:

  * logprob: the log probability at each step of the MCMC
  * mu: the final estimate of the mean of each parameter
  * cov: the final estimate of the covariance matrix of parameters

* some of the input values

When the tool is run while the output exists it will use it instead of rerunning it which is useful to get more plots.

Plots and txt outputs
---------------------

When ``--figure name.extension`` is given then some QC and results are given.

QC
^^

name_convergence.extension
""""""""""""""""""""""""""

This plot shows the autocorrelation between samples for all parameters (the solid line shows the median, the shaded area shows the min and max). If the MCMC converged, you should see a value of ACF close to 0 since a small value of T.

name_neff.txt
"""""""""""""

From the autocorrelation displayed above, we can evaluate the number of independent samples, also called effective sample size. The value is printed and stored in this text file.

name_p.extension
""""""""""""""""

This plot shows the value of each parameter and the log probability (y axis, one panel per parameter) for all samples (x-axis). When the MCMC did not converged, it can be helpful to see if it can be explained by the fact that the first samples considered where not around the final solution. In this case, it can be useful to rerun the plots using an increase value of ``--removeFirstSamples`` (by default it is 1/4 of the number of samples), or increase the number of samples.

name_corner.extension
"""""""""""""""""""""

This plot shows the distribution of value of each parameter in relationship the one with the other. It can help to see which parameters are correlated. Also, when the MCMC did not converged, it can help to identify if 2 or more solutions were explored.


Results
^^^^^^^

name.extension
""""""""""""""

This is the figure with the results. 

- When the 1d version is used, it displays the mean pdf in solid red line, the median in black dashed lines (/!\backslash the integral of the median is not equal to 1) with the confidence interval of 1 sigma (68%), 2 sigma (95%) and 3 sigma (99.7%) as well as in green, the kernel density estimate of the input values, the detected expression (``log(1 + targetSum * raw / total UMI)``).

- When the 2d version is used, it displays the pdf as a heatmap as well as a projection on the x and y axis. On the projection, the confidence interval 68% is indicated as a shaded area as well as the mean with a solid red line and the median with a dashed black line. On the top right corner, the correlation is indicated with the confidence interval 68% as well as a confidence interval on the one-sided p-value (the probability that the correlation is the opposite sign of the mean, one sigma confidence interval).

name_individuals.extension
""""""""""""""""""""""""""

- When the 1d version is used, it displays the pdf of 100 samples.

- When the 2d version is used, it displays the projection of the pdf of 100 samples.

name_p.txt
""""""""""

This is a tabulated delimited table with the 16 percentile (low), median, 84 percentile (high) value of each parameter.


name_pdf.txt (1d only)
""""""""""""""""""""""

For each value of x, the 16 percentile (low), mean, 84 percentile (high) and median, is given in a tabulated delimited file.

name_with_posterior.extension (1d only)
"""""""""""""""""""""""""""""""""""""""

Same as name.extension except that a new orange line is plotted showing the posterior density evaluated as the average of the posterior density of each cell.

name_posterior_individuals.extension (1d only)
""""""""""""""""""""""""""""""""""""""""""""""

Showing posterior density probability of 50 random cells.

name_posterior_per_cell.txt (1d only)
"""""""""""""""""""""""""""""""""""""

For each cell of the input, providing the posterior average and standard deviation of the density probability.

name_posterior_andco.extension (1d only)
""""""""""""""""""""""""""""""""""""""""

Showing the mean pdf, the median pdf, the density from raw counts normalized, the average of the posterior density from all cells,  the density and a histogram using only the average value of the posterior distribution of each cell and the posterior density approximating the pdf of each cell by a Gaussian using values in the "posterior_per_cell.txt" file.

name_means.txt.gz (1d only)
"""""""""""""""""""""""""""

Each line correspond to the value of the mean expression evaluated at each sample of the MCMC.

name_median.extension (2d only)
"""""""""""""""""""""""""""""""

Same as name.extension except that the median instead of the mean is used.

name_corr.txt (2d only)
"""""""""""""""""""""""

The mean, median, 16 percentile, 84 percentile, p-value and error on the p-value for the correlation (see above).


name_pdf2d.txt (2d only)
""""""""""""""""""""""""

The mean pdf and the x and y values stored in a tabulated delimited file in a matrix format. Different x values correspond to different columns while different y values correspond to different rows.

name_pdf2d_flat.txt (2d only)
"""""""""""""""""""""""""""""

The x, y, 16 percentile (low), mean, 84 percentile (high) and median of pdf in a tabulated delimited file.

Results when ``--splity`` is provided in 2d
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``--splity`` is provided the pdf above and below this threshold on the y axis are summed up, resulting in 2 pdf along the x axis.

name_splitX.extension
"""""""""""""""""""""
This plot shows the 2 pdfs. The ratio between the area represent the ratio of cells above and below the threshold of the gene y. The pdf for cells below the threshold is in red (with the shaded area for the 68% confidence interval) and the pdf for cells above the threshold is in green. In black is the pdf of all cells projected on the x axis (sum of the 2).

name_splitX_renorm.extension
""""""""""""""""""""""""""""
Same plot as above except that the pdf were renormalized so the area of each pdf is equal to 1. Also the median is added in dashed black lines.

name_splitX.txt
"""""""""""""""
This is a tabulated delimited table with the x values, the 16 percentile (low), mean, 84 percentile (high) values of each pdf (below and above the threshold) before normalization.


Evidence
--------
When ``--logevidence`` is set. The log evidence is calculated and stored in this file. This can be used to compare different models, here different number of gaussians.

