# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
"""
  Module with all functions common to 1d and 2d
"""
import argparse
import os
import time
from tempfile import NamedTemporaryFile
from shutil import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import poisson
from scipy.sparse import issparse
import anndata
from samsam import acf

from . _version import __version__

mpl.use('agg')


def get_data_txt(input_file, metadata1_cn, metadata1_v,
                 metadata2_cn, metadata2_v,
                 metadata3_cn, metadata3_v,
                 genes):
  """
  Create a dataframe from input file
  keep only rows matching metadata values

  Args:
    input_file (str): File path for input table
    metadata1_cn (str): Column name for the first metadata
    metadata1_v (str): Comma separated possible values for metadata1
    metadata2_cn (str): Column name for the second metadata
    metadata2_v (str): Comma separated possible values for metadata2
    metadata3_cn (str): Column name for the third metadata
    metadata3_v (str): Comma separated possible values for metadata3
    genes ([str]): List of column names with genes

  Returns:
    data (pandas.DataFrame): Data Frame with nCount_RNA and the genes of interest
                             for the cells matching metadata
  """
  data = pd.read_csv(input_file, sep="\t")
  filters = {metadata1_cn: metadata1_v,
             metadata2_cn: metadata2_v,
             metadata3_cn: metadata3_v}
  for col_name, values in filters.items():
    if col_name is None:
      continue
    data = data[data[col_name].astype(str).isin(values.split(','))]
  print(f"You have {data.shape[0]} cells.")
  if data.shape[0] == 0:
    raise ValueError("No cell with selected filters.")
  for col_name in genes + ['nCount_RNA']:
    if col_name not in data.columns:
      raise ValueError(f"{col_name} not in the data table.")
  return data[genes + ['nCount_RNA']]


def get_raw_mat_from_annData(adata, used_raw=False, round_threshold=0.01):
  """
  Retrieves the raw counts matrix from annData.
  We expect that if log and norm were applied
  it was applied first Norm and then Log
  """
  add_message = ""
  if used_raw:
    add_message = "in .raw"
  # Transform sparse to array
  if issparse(adata.X):
    X = adata.X.toarray()
  else:
    X = np.copy(adata.X)
  # Check if the data (X) are all positive (not scaled):
  if not np.all(X >= 0):
    raise ValueError(f"Negative values found {add_message}. Cannot go back to raw counts.")
  # Check if they are already raw_counts:
  if np.array_equal(X, X.astype(int)):
    return X
  # Check if a log was applied
  if 'log1p' in adata.uns:
    base = adata.uns['log1p']['base']
    if base is not None:
      X *= np.log(base)
    X = np.expm1(X)
    add_message += " delogged"
  # Check if normalization was done:
  # We assume that if all values
  # are really close to integer no norm was performed.
  if np.max(np.abs(X - np.round(X))) < round_threshold:
    return np.round(X)
  # We first check that all values are there:
  # Try to find the target_sum:
  N_t = np.median(X.sum(1))
  if np.max(np.abs(X.sum(1) - N_t)) > 0.1:
    # Maybe it was logged but the lop1p in uns has not been exported
    X = np.expm1(X)
    add_message += " delogged"
    if np.max(np.abs(X - np.round(X))) < round_threshold:
      return np.round(X)
    N_t = np.median(X.sum(1))
    if np.max(np.abs(X.sum(1) - N_t)) > 0.1:
      raise ValueError("The sum of expression for each cell"
                       f" {add_message} is not constant."
                       " Genes are probably missing.")

  # We could also check the number of genes detected with:
  # np.array_equal(np.sum(X > 0, axis=1), adata.obs['n_genes_by_counts'])
  if 'total_counts' not in adata.obs.keys():
    raise ValueError(f"The expression {add_message}seems normalized"
                     " but 'total_counts' is not in obs."
                     "Cannot get raw values.")
  # Denormalize:
  X = X * adata.obs['total_counts'][:, None] / N_t
  if np.max(np.abs(X - np.round(X))) >= round_threshold:
    raise ValueError(f"Could not extract raw values of expression {add_message}.")

  return np.round(X)


def get_data_annData(input_file, metadata1_cn, metadata1_v,
                     metadata2_cn, metadata2_v,
                     metadata3_cn, metadata3_v,
                     genes):
  """
  Create a dataframe from input file
  keep only rows matching metadata values
  """
  # First load the annData:
  # Deprecation warning:
  # FutureWarning: `anndata.read` is deprecated,
  # use `anndata.read_h5ad` instead. `ad.read` will be removed in mid 2024.
  adata = anndata.read(input_file)
  raw_mat = None
  used_raw = False
  # If a raw exists we will use it:
  if adata.raw is not None:
    adata = adata.raw.to_adata()
    used_raw = True
  raw_mat = get_raw_mat_from_annData(adata, used_raw)
  data = pd.DataFrame(data=raw_mat,
                      index=adata.obs_names,
                      columns=adata.var_names)
  data = pd.concat([data, adata.obs], axis=1)
  if 'total_counts' in data.keys():
    data['nCount_RNA'] = data['total_counts']
  else:
    data['nCount_RNA'] = raw_mat.sum(1)

  filters = {metadata1_cn: metadata1_v,
             metadata2_cn: metadata2_v,
             metadata3_cn: metadata3_v}
  for col_name, values in filters.items():
    if col_name is None:
      continue
    data = data[data[col_name].astype(str).isin(values.split(','))]
  print(f"You have {data.shape[0]} cells.")
  if data.shape[0] == 0:
    raise ValueError("No cell with selected filters.")
  for col_name in genes + ['nCount_RNA']:
    if col_name not in data.columns:
      raise ValueError(f"{col_name} not in the data table.")
  return data[genes + ['nCount_RNA']]


def permuted(p, pref, wdist, permut, n_param_per_gauss):
  """Permut p (parameters) to get the permutation of p the closest to pref

  Args:
    p (np.ndarray): Set of parameters to permute
    pref (np.ndarray): Set of parameters of reference
    wdist (np.ndarray): Weight to give to parameters in distance evaluation
    permut (np.ndarray): Array with all possible permutations
    n_param_per_gauss (int): Number of parameter per Gaussian

  Returns:
    min_perm (np.ndarray): A permuted version of p which is the closest to pref
  """
  amp = p[(n_param_per_gauss - 1)::n_param_per_gauss]
  # We complete the list of parameters:
  p_full = np.concatenate(([1 - np.sum(amp)], p))
  dp = p - pref
  min_dist = np.sum(wdist * dp * dp)
  min_perm = p
  for ks in permut[1:]:
    dpk = p_full[ks][1:] - pref
    distk = np.sum(wdist * dpk * dpk)
    if distk < min_dist:
      min_dist = distk
      min_perm = p_full[ks][1:]
  return min_perm


def get_Ax(data, col_gene, nx, dx, ox, xscale, target_sum):
  """Get p(ki|Ni,xj) * dx

  Args:
    data (pandas.DataFrame): Data Frame with gene values and nCount_RNA
    col_gene (str): Column label in `data` with the gene of interest
    nx (int): Number of values in x to check how your evaluated pdf is compatible with the model.
    dx (float): Distance between 2 bin centers
    ox (np.ndarray): The nx * osampx bin centers
    xscale (str): scale for the x-axis: Seurat (log(1+targetSum*X)) or log (log(X))
    target_sum (int): factor when Seurat scale is used: (log(1+targetSum*X))

  Raises:
      ValueError: If size of ox is not multiple of size of x
      NotImplementedError: If xscale is not implemented

  Returns:
      np.ndarray: Matrix with p(ki|Ni,xj) * dx
  """
  # Get all UMI nbs:
  N = data['nCount_RNA']
  # Get the raw counts:
  kx = data[col_gene]
  # To get a better approximation we will use
  # a finer grid for x (ox)
  if ox.size % nx != 0:
    raise ValueError("The size of ox is not multiple of size of x.")
  # The Poisson law is evaluated on the finest grid (ox)
  # Summed up on the grid (x)
  if xscale == 'Seurat':
    # We assume ki follow
    # a Poisson distribution of parameter (exp(x) - 1)*Ni / target_sum
    # as x is expressed in the axis used in Seurat:
    # log(1 + target_sum X)
    if target_sum == 0:
      target_sum = np.median(N)
    poissx = np.array([
        np.mean(poisson.pmf(k=kxi, mu=Ni / target_sum * (np.exp(ox) - 1)).reshape(nx, -1), axis=1)
        for kxi, Ni in zip(kx, N)])
  elif xscale == 'log':
    # or a Poisson distribution of parameter exp(x)*Ni
    # as x is expressed in the log axis:
    # log(X)
    poissx = np.array([
        np.mean(poisson.pmf(k=kxi, mu=Ni * np.exp(ox)).reshape(nx, -1), axis=1)
        for kxi, Ni in zip(kx, N)])
  else:
    raise NotImplementedError("Only Seurat and log xscale are implemented.")
  # To evaluate the integral more easily multiply by the delta x
  Ax = poissx * dx

  return Ax


def get_prefix_suffix(output):
  """Process the output to get prefix and suffix

  Args:
    output (str): File path from which prefix and suffix should be extracted

  Returns:
    file_prefix (str): Basename
    file_suffix (str): Suffix (extension)
  """
  name = output.split(".")
  file_suffix = name[-1]
  file_prefix = ".".join(name[:-1])
  return(file_prefix, file_suffix)


def get_Neff(samples):
  """Get the number of effective samples

  Args:
    samples (np.ndarray): matrix with parameters for each sample

  Returns:
    float: Evaluation of the number of effective samples
  """
  # Autocorrelation (to test the convergence of MCMC)
  R = np.array([acf.acf(pi) for pi in samples.T]).T
  # Estimate the number of independent samples:
  Neff = samples.shape[0] / np.max([acf.iat(R=r) for r in R.T])
  return Neff


def plot_QC(logprob_values, samples, title, output,
            removeFirstSamples, nsampInPlot,
            p_names):
  """Do some plots and exports in txt

  Args:
    logprob_values (np.ndarray): The logprob values
    samples (np.ndarray): matrix with parameters for each sample
    title (str): Title for plots
    output (str): Path for output plots
    removeFirstSamples (int): Number of samples to ignore before making the plots
    nsampInPlot (int): Approximate number of samples to use in plots.
    p_names ([str]): Array with parameter names
  """
  file_prefix, file_suffix = get_prefix_suffix(output)
  # Remove the first samples:
  samples = samples[removeFirstSamples:]
  print(f"Considering the last {samples.shape[0]} samples.")
  logprob_values = logprob_values[removeFirstSamples:]

  # Select nsampInPlot
  samples_selected = np.unique(np.linspace(0, samples.shape[0] - 1, nsampInPlot, dtype=int))
  samples = samples[samples_selected]
  print(f"Using {samples.shape[0]} samples for plots.")
  logprob_values = logprob_values[samples_selected]

  # Autocorrelation (to test the convergence of MCMC)
  R = np.array([acf.acf(pi) for pi in samples.T]).T
  Rmin = np.min(R, axis=1)
  Rmax = np.max(R, axis=1)
  Rmed = np.median(R, axis=1)
  x_values = samples_selected[1:]
  plt.figure()
  plt.fill_between(x_values, Rmin[1:], Rmax[1:], color='r', alpha=0.3, rasterized=True)
  plt.plot(x_values, Rmed[1:], 'r', lw=2, rasterized=True)
  plt.xscale('log')
  plt.xlim(x_values[0], x_values[-1])
  plt.xlabel(r'$\tau$')
  plt.ylabel('ACF')
  plt.title(title)
  plt.savefig(f'{file_prefix}_convergence.{file_suffix}')

  # Estimate the number of independent samples:
  Neff = samples.shape[0] / np.max([acf.iat(R=r) for r in R.T])

  print(f"Neff is {Neff}")
  with open(f'{file_prefix}_neff.txt', 'w', encoding="utf-8") as fo:
    fo.write(f'{Neff}\n')

  # Plot the value of each parameter and logprob_values along the samples:
  fig, axs = plt.subplots(samples.shape[1] + 1, 1,
                          sharex='col', figsize=(6, (samples.shape[1] + 1) * 1.5))
  for i, p_name in enumerate(p_names):
    ax = axs[i]
    ax.plot(samples_selected, samples[:, i])
    ax.set_ylabel(p_name)
  ax = axs[-1]
  ax.plot(samples_selected, logprob_values)
  ax.set_ylabel('logprob')
  plt.tight_layout()
  fig.savefig(f'{file_prefix}_p.{file_suffix}')

  # 1, 2, 3 sigma:
  qm = np.array([(1 - lvl) / 2 for lvl in [0.6827, 0.9545, 0.9973]])
  qp = 1 - qm

  # Export the parameters median and 1 sigma
  p_summary = pd.DataFrame(columns=['name', 'low', 'median', 'high'])
  p_summary['name'] = p_names
  p_summary['low'] = np.quantile(samples, qm[0], axis=0)
  p_summary['high'] = np.quantile(samples, qp[0], axis=0)
  p_summary['median'] = np.median(samples, axis=0)
  p_summary.to_csv(f'{file_prefix}_p.txt',
                   index=False,
                   header=True, sep='\t')

  # Try to plot corner
  try:
    import corner  # pylint: disable=C0415
    plt.figure()
    figure = corner.corner(samples,
                           labels=p_names,
                           quantiles=[0.15865, 0.5, 0.84135],
                           show_titles=True, title_kwargs={"fontsize": 12})
    figure.suptitle(title, ha='right', x=0.98)
    plt.savefig(f'{file_prefix}_corner.{file_suffix}', dpi=100)
  except ImportError as e:
    print("corner plot could not be loaded.")
    print(e)

  return(samples, Neff)


def get_bins_centers(xmin, xmax, nx):
  """Get the position of bin centers

  Args:
    xmin (int): Minimum value to consider.
    xmax (int): Maximum value to consider.
    nx (int): Number of values.

  Returns:
      np.ndarray: The nx bin centers
  """
  x_borders = np.linspace(xmin, xmax, nx + 1)
  return (x_borders[:-1] + x_borders[1:]) / 2


def parse_arguments(tool: str) -> argparse.ArgumentParser:
  """Argument parser for all 4 tools
  """
  docu_ref = "The full documentation is available at https://baredsc.readthedocs.io"
  if tool == "baredSC_1d":
    argp = argparse.ArgumentParser(
        description=("Run mcmc to get the pdf for a given gene using a"
                     f" fixed number normal distributions. {docu_ref}"))
  elif tool == "baredSC_2d":
    argp = argparse.ArgumentParser(
        description=("Run mcmc to get the pdf in 2D for 2 given genes using a"
                     f" fixed number normal distributions. {docu_ref}"))
  elif tool == "combineMultipleModels_1d":
    argp = argparse.ArgumentParser(
        description=("Combine mcmc 1D results from multiple models to get a"
                     " mixture using logevidence to infer weights."))
  elif tool == "combineMultipleModels_2d":
    argp = argparse.ArgumentParser(
        description=("Combine mcmc 2D results from multiple models to get a"
                     " mixture using logevidence to infer weights."))
  else:
    raise NotImplementedError("This tool has not been implemented")

  argprequired = argp.add_argument_group('Required arguments')
  argpopt_data = argp.add_argument_group('Optional arguments to select input data')
  argpopt_mcmc = argp.add_argument_group('Optional arguments to run MCMC')
  argpopt_plot = argp.add_argument_group('Optional arguments to customize plots and text outputs')
  argpopt_loge = argp.add_argument_group('Optional arguments to evaluate logevidence')

  # Get data:
  group = argprequired.add_mutually_exclusive_group(required=True)
  group.add_argument('--input', default=None,
                     help="Input table (tabular separated"
                     " with header) with one line per cell"
                     " columns with raw counts and one column"
                     " nCount_RNA with total number of UMI per cell"
                     " optionally other meta data to filter.")
  group.add_argument('--inputAnnData', default=None,
                     help="Input annData (for example from Scanpy).")
  if tool.endswith("1d"):
    argprequired.add_argument('--geneColName', default=None, required=True,
                              help="Name of the column with gene counts.")
  elif tool.endswith("2d"):
    argprequired.add_argument('--geneXColName', default=None, required=True,
                              help="Name of the column with gene counts for gene in x.")
    argprequired.add_argument('--geneYColName', default=None, required=True,
                              help="Name of the column with gene counts for gene in y.")
  argpopt_data.add_argument('--metadata1ColName', default=None,
                            help="Name of the column with metadata1 to filter.")
  argpopt_data.add_argument('--metadata1Values', default=None,
                            help="Comma separated values for metadata1 of cells to keep.")
  argpopt_data.add_argument('--metadata2ColName', default=None,
                            help="Name of the column with metadata2 to filter.")
  argpopt_data.add_argument('--metadata2Values', default=None,
                            help="Comma separated values for metadata2 of cells to keep.")
  argpopt_data.add_argument('--metadata3ColName', default=None,
                            help="Name of the column with metadata3 to filter.")
  argpopt_data.add_argument('--metadata3Values', default=None,
                            help="Comma separated values for metadata3 of cells to keep.")
  # MCMC
  if tool.startswith("combineMultipleModels"):
    argprequired.add_argument('--outputs', default=None, required=True, nargs='+',
                              help="Ouput files basename (will be npz)"
                              " with different results of mcmc to combine.")
  argpopt_mcmc.add_argument('--xmin', default=0, type=float,
                            help="Minimum value to consider in x axis.")
  argpopt_mcmc.add_argument('--xmax', default=2.5, type=float,
                            help="Maximum value to consider in x axis.")
  if tool.endswith("1d"):
    default_nx = 100
    default_osampxpdf = 5
  else:
    default_nx = 50
    default_osampxpdf = 4
  argpopt_mcmc.add_argument('--nx', default=default_nx, type=int,
                            help="Number of values in x to check how "
                            "your evaluated pdf is compatible with the model.")
  argpopt_mcmc.add_argument('--osampx', default=10, type=int,
                            help="Oversampling factor of x values when evaluating "
                            "pdf of Poisson distribution.")
  argpopt_mcmc.add_argument('--osampxpdf', default=default_osampxpdf, type=int,
                            help="Oversampling factor of x values when evaluating "
                            "pdf at each step of the MCMC.")

  if tool.endswith("1d"):
    argpopt_mcmc.add_argument('--minScale', default=0.1, type=float,
                              help="Minimal value of the scale of gaussians"
                              " (Default is 0.1 but cannot be smaller than "
                              "max of twice the bin size of pdf evaluation"
                              " and half the bin size).")
  else:
    argpopt_mcmc.add_argument('--minScalex', default=0.1, type=float,
                              help="Minimal value of the scale of gaussians on x"
                              " (Default is 0.1 but cannot be smaller than "
                              "max of twice the bin size of pdf evaluation"
                              " and half the bin size on x axis).")
    argpopt_mcmc.add_argument('--ymin', default=0, type=float,
                              help="Minimum value to consider in y axis.")
    argpopt_mcmc.add_argument('--ymax', default=2.5, type=float,
                              help="Maximum value to consider in y axis.")
    argpopt_mcmc.add_argument('--ny', default=50, type=int,
                              help="Number of values in y to check how "
                              "your evaluated pdf is compatible with the model.")
    argpopt_mcmc.add_argument('--osampy', default=10, type=int,
                              help="Oversampling factor of y values when evaluating "
                              "pdf of Poisson distribution.")
    argpopt_mcmc.add_argument('--osampypdf', default=4, type=int,
                              help="Oversampling factor of y values when evaluating "
                              "pdf at each step of the MCMC.")
    argpopt_mcmc.add_argument('--minScaley', default=0.1, type=float,
                              help="Minimal value of the scale of gaussians on yx"
                              " (Default is 0.1 but cannot be smaller than "
                              "max of twice the bin size of pdf evaluation"
                              " and half the bin size on y axis).")
    argpopt_mcmc.add_argument('--scalePrior', default=0.3, type=float,
                              help="Scale of the truncnorm used in the prior for "
                              "the correlation.")
  if tool.endswith("1d"):
    argpopt_mcmc.add_argument('--xscale', default="Seurat", choices=['Seurat', 'log'],
                              help="scale for the x-axis: Seurat (log(1+targetSum*X))"
                              " or log (log(X))")
  else:
    argpopt_mcmc.add_argument('--scale', default="Seurat", choices=['Seurat', 'log'],
                              help="scale for the x-axis and y-axis: Seurat"
                              " (log(1+targetSum*X)) or log (log(X))")
  argpopt_mcmc.add_argument('--targetSum', default=10**4, type=float,
                            help="factor when Seurat scale is used: (log(1+targetSum*X))"
                            " (default is 10^4, use 0 for the median of nRNA_Counts)")
  argpopt_mcmc.add_argument('--seed', default=1, type=int,
                            help="Change seed for another output.")

  if tool.startswith("baredSC_"):
    if tool == 'baredSC_1d':
      default_nnorm = 2
    else:
      default_nnorm = 1
    argpopt_mcmc.add_argument('--nnorm', default=default_nnorm, type=int,
                              help="Number of gaussians to fit.")
    argpopt_mcmc.add_argument('--nsampMCMC', default=100000, type=int,
                              help="Number of samplings (iteractions) of mcmc.")
    argpopt_mcmc.add_argument('--nsampBurnMCMC', default=None, type=int,
                              help="Number of samplings (iteractions) in the "
                              "burning phase of mcmc (Default is nsampMCMC / 4).")
    argpopt_mcmc.add_argument('--nsplitBurnMCMC', default=10, type=int,
                              help="Number of steps in the "
                              "burning phase of mcmc.")
    argpopt_mcmc.add_argument('--T0BurnMCMC', default=100.0, type=float,
                              help="Initial temperature in the "
                              "burning phase of mcmc.")
    argpopt_mcmc.add_argument('--minNeff', default=None, type=float,
                              help="Will redo the MCMC with 10 times more samples until "
                              "the number of effective samples that this value "
                              "(Default is not set so will not rerun MCMC).")
    # To save/get MCMC
    argpopt_mcmc.add_argument('--force', default=None, action='store_true',
                              help="Force to redo the mcmc even if output exists.")
    argprequired.add_argument('--output', default=None, required=True,
                              help="Ouput file basename (will be npz)"
                              " with results of mcmc.")

  # Plot
  argprequired.add_argument('--figure', default=None,
                            required=tool.startswith("combineMultipleModels"),
                            help="Ouput figure basename.")
  argpopt_plot.add_argument('--title', default=None,
                            help="Title in figures.")
  if tool.endswith("2d"):
    argpopt_plot.add_argument('--splity', default=None, nargs='+', type=float,
                              help="Threshold value to plot the density for genex"
                              " for 2 categories in geney values.")
  argpopt_plot.add_argument('--removeFirstSamples', default=None, type=int,
                            help="Number of samples to ignore before making the plots"
                            " (default is nsampMCMC / 4).")
  argpopt_plot.add_argument('--nsampInPlot', default=100000, type=int,
                            help="Approximate number of samples to use in plots.")
  if tool.endswith("1d"):
    argpopt_plot.add_argument('--prettyBins', default=None, type=int,
                              help="Number of bins to use in plots (Default is nx).")
  if tool.endswith("2d"):
    argpopt_plot.add_argument('--prettyBinsx', default=None, type=int,
                              help="Number of bins to use in x in plots (Default is nx).")
    argpopt_plot.add_argument('--prettyBinsy', default=None, type=int,
                              help="Number of bins to use in y in plots (Default is ny).")
    argpopt_plot.add_argument('--log1pColorScale', action='store_true',
                              help="Use log1p color scale instead of linear color scale.")
  if tool == "combineMultipleModels_2d":
    argpopt_plot.add_argument('--getPVal', action='store_true',
                              help="Use less samples to get an estimation of the p-value.")
  # Evidence
  if tool.startswith("combineMultipleModels"):
    argpopt_loge.add_argument('--logevidences', default=None, nargs='+',
                              help="Ouput files of precalculated log evidence values."
                              "(if not provided will be calculated).")
  if tool.startswith("baredSC_"):
    argpopt_loge.add_argument('--logevidence', default=None,
                              help="Ouput file to put logevidence value.")

  argpopt_loge.add_argument('--coviscale', default=1, type=float,
                            help="Scale factor to apply to covariance of parameters"
                            " to get random parameters in logevidence evaluation.")
  argpopt_loge.add_argument('--nis', default=1000, type=int,
                            help="Size of sampling of random parameters in logevidence evaluation.")
  # Version
  argp.add_argument('--version', action='version',
                    version=__version__)
  return argp


def checkNeff(args, original_args):
  """Check if the number of effective samples is compatible
  with args.minNeff and if it is not provides the new args.

  Args:
      args (namespace): arguments used to run MCMC
      original_args (namespace): arguments provided by user before checks

  Returns:
      namespace: New arguments to use to run or None if we don't need to rerun
  """

  # If minNeff is not set, no need to rerun
  if args.minNeff is None:
    return None

  # Process the output to get prefix and suffix and read Neff
  file_prefix, _ = get_prefix_suffix(args.figure)
  with open(f'{file_prefix}_neff.txt', 'r', encoding="utf-8") as f:
    neff = float(f.read())
  if neff < args.minNeff:
    print("Neff is below the minimum required.")
    with NamedTemporaryFile(delete=False, suffix=".npz") as temp:
      temp_file = temp.name
      print(f"Results are moved to {temp_file}")
      copy(f'{args.output}.npz', temp_file)
    os.remove(f'{args.output}.npz')
    # Change the args to multiply by 10 the nb of samples
    new_args = original_args
    if '--nsampMCMC' in new_args:
      i = [i for i, v in enumerate(new_args) if v == '--nsampMCMC'][0]
      new_args[i + 1] = f"{args.nsampMCMC}0"
    else:
      new_args.append('--nsampMCMC')
      new_args.append(f'{args.nsampMCMC}0')
    return new_args
  # neff is large enough
  return None


def args_check_baredSC(args):
  """_summary_

  Args:
      args (namespace): args from argparse

  Raises:
      ValueError: if incompatible options are set

  Returns:
      namespace: updated args
  """
  # Put default:
  if args.nsampBurnMCMC is None:
    args.nsampBurnMCMC = args.nsampMCMC // 4
  if args.removeFirstSamples is None:
    args.removeFirstSamples = args.nsampMCMC // 4
  # Remove potential suffix:
  args.output = args.output.removesuffix('.npz')
  # Check incompatibilities:
  if args.minNeff is not None and args.figure is None:
    raise ValueError("--minNeff requires --figure to be set.")
  return args


def get_data_from_args(args, genes):
  """Get data from input or inputAnnData

  Args:
    args (namespace): arguments
    genes ([str]): List of column names with genes

  Returns:
    pandas.DataFrame: Data Frame with selected rows and columns
  """
  print("Get raw data")
  start = time.time()
  if args.input is not None:
    input_file = args.input
    get_data = get_data_txt
  else:
    input_file = args.inputAnnData
    get_data = get_data_annData

  data = get_data(input_file,
                  args.metadata1ColName, args.metadata1Values,
                  args.metadata2ColName, args.metadata2Values,
                  args.metadata3ColName, args.metadata3Values,
                  genes)
  print(f"Got. It took {(time.time() - start):.2f} seconds.")
  return data
