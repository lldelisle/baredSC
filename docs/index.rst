.. baredSC documentation master file, created by
   sphinx-quickstart on Wed Jan 27 07:35:21 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to baredSC's documentation!
====================================

BARED (Bayesian Approach to Retreive Expression Distribution of) Single Cell
----------------------------------------------------------------------------

baredSC is a tool that uses a Monte-Carlo Markov Chain to estimate a confidence interval on the probability density function (PDF) of expression of one or two genes from single-cell RNA-seq data. It uses the raw counts and the total number of UMI for each cell. The PDF is approximated by a number of 1d or 2d gaussians provided by the user. The likelihood is estimated using the asumption that the raw counts follow a Poisson distribution of parameter equal to the proportion of mRNA for the gene in the cell multiplied by the total number of UMI identified in this cell.

We encourage users which are not familiar with MCMC to start with the :doc:`content/installation` page to install the baredSC package. Then, follow the tutorial for a single gene (:doc:`content/tuto/tutorial_sim_default_1d`) which is a step by step example of application using an input provided on the github repository with default parameters. The other tutorials describe influences of some parameters or use a 2d example.

Users familiar with MCMC could directly go to the :doc:`content/usage` page and read the :doc:`content/outputs` page.

Advanced users who wants to use baredSC as a package within python should go to the :doc:`content/baredSC_api`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   content/installation.rst
   content/usage.rst
   content/outputs.rst
   content/tutorial_sim.rst
   content/baredSC_api.rst
   content/releases.rst
   content/citation.rst

