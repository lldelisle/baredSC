baredSC
========
[![PyPI Version](https://img.shields.io/pypi/v/baredsc.svg?style=plastic)](https://pypi.org/project/baredSC/) [![DOI](https://zenodo.org/badge/370966963.svg)](https://zenodo.org/badge/latestdoi/370966963)


baredSC (Bayesian Approach to Retreive Expression Distribution of Single Cell) is a tool that uses a Monte-Carlo Markov Chain to estimate a confidence interval on the probability density function (PDF) of expression of one or two genes from single-cell RNA-seq data. It uses the raw counts and the total number of UMI for each cell. The PDF is approximated by a number of 1d or 2d gaussians provided by the user. The likelihood is estimated using the asumption that the raw counts follow a Poisson distribution of parameter equal to the proportion of mRNA for the gene in the cell multiplied by the total number of UMI identified in this cell.

Documentation
-------------
Visit [our documentation](https://baredsc.readthedocs.io) to see the possible options and follow the tutorials.

Citation
--------
If you are using baredSC, please cite our [biorxiv paper](https://doi.org/10.1101/2021.05.26.445740).
