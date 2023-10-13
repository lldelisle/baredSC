# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
"""
  Module with all functions common to baredSC_1d combineMultipleModels_1d
"""
from itertools import permutations
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import kde
from scipy.special import factorial
import pandas as pd
import samsam
import samsam.logprior

# Local imports
from . common import get_Ax, get_prefix_suffix, get_bins_centers

mpl.use('agg')


def logprob(p,
            pref, wdist, permutx,
            xmin, xmax, nx,
            noxpdf, oxpdf, odxpdf,
            min_scale,
            Ax, temp=1.0):
  """Define log probability function

  Args:
    p (np.ndarray): Set of parameters to get logprob
    pref (np.ndarray): Set of parameters of reference
    wdist (np.ndarray): Weight to give to parameters in distance evaluation
    permutx (np.ndarray): Array with all possible permutations on x
    xmin (int): Minimum value to consider in x axis.
    xmax (int): Maximum value to consider in x axis.
    nx (int): Number of values in x to check how your evaluated pdf is compatible with the model.
    noxpdf (int): Number of values in x to evaluate pdf
    oxpdf (np.ndarray): bin centers for pdf
    odxpdf (float): Distance between 2 bin centers in oxpdf
    min_scale (float): Minimal value of the scale of gaussians
    Ax (np.ndarray): Matrix with p(ki|Ni,xj) * dx
    temp (float, optional): Temperature. Defaults to 1.0.

  Returns:
    float: Log probability
  """
  # Log prior
  lp = logprior(p, pref, wdist, permutx, xmin, xmax, min_scale)
  if np.isfinite(lp):
    # If the parameters are possible
    # Log likelihood
    ll = loglike(p, nx, noxpdf, oxpdf, odxpdf, Ax)
    # Return sum of log prior and log likelihood
    return (lp + ll) / temp
  # If the parameters are impossible
  return -np.inf


def logprior(p, pref, wdist, permutx, xmin, xmax, min_scale):
  """Define the log prior

  Args:
    p (np.ndarray): Set of parameters to get logprob
    pref (np.ndarray): Set of parameters of reference
    wdist (np.ndarray): Weight to give to parameters in distance evaluation
    permutx (np.ndarray): Array with all possible permutations on x
    xmin (int): Minimum value to consider in x axis.
    xmax (int): Maximum value to consider in x axis.
    min_scale (float): Minimal value of the scale of gaussians

  Returns:
    float: Log prior value and -inf if it is not in prior
  """
  # Raise exception if one of the condition is not
  # satisfied
  try:
    amp = p[2::3]
    mu = p[::3]
    sig = p[1::3]
    # We complete the list of parameters:
    p_full = np.concatenate(([1 - np.sum(amp)], p))
    amp_full = np.concatenate(([1 - np.sum(amp)], amp))
    # Amplitude
    # Amplitudes should be positive
    if np.any(amp_full < 0):
      raise ValueError('Negative amplitude')
    # Because we impose that sum of amplitude is below 1
    # We add this term:
    lp = np.log(factorial(amp.size))
    # Check swapping (symmetry of the problem)
    # pref is a parameter set and we impose that pref
    # is closer to p than any permutation of p
    dp = p - pref
    dist = np.sum(wdist * dp * dp)
    for ks in permutx[1:]:
      dpk = p_full[ks][1:] - pref
      distk = np.sum(wdist * dpk * dpk)
      if distk <= dist:
        raise ValueError('Swapping')
    # Number of possible swappings
    lp += np.log(permutx.shape[0])
    # Finally we add the priors on other parameters
    # Location
    lp += sum(samsam.logprior.uniform(muk, xmin - 3 * sigk, xmax + 3 * sigk)
              for muk, sigk in zip(mu, sig))
    # Scale
    lp += sum(samsam.logprior.loguniform(sigk, np.log(min_scale), np.log(xmax - xmin))
              for sigk in sig)
  except Exception:  # pylint: disable=W0718
    return -np.inf
  return lp


def trunc_norm(x, dx, mu, sigma):
  """pdf of a trunc_norm as scipy.stats.truncnorm.pdf is very slow

  Args:
      x (np.ndarray): Regular grid where pdf should be evaluated
      dx (float): Space between 2 values in x
      mu (float): Mean of Gaussian
      sigma (float): Scale of Gaussian

  Returns:
      np.ndarray: pdf values normalized
  """
  u = (x - mu) / sigma
  lognorm = -0.5 * u * u
  lognorm -= np.max(lognorm)  # Avoid underflow
  norm = np.exp(lognorm)
  norm /= np.sum(norm * dx)
  return norm


def loglike(p, nx, noxpdf, oxpdf, odxpdf, Ax):
  """Log likelihood

  Args:
    p (np.ndarray): Set of parameters to get loglikelihood
    nx (int): Number of values in x to check how your evaluated pdf is compatible with the model.
    noxpdf (int): Number of values in x to evaluate pdf
    oxpdf (np.ndarray): bin centers for pdf
    odxpdf (float): Distance between 2 bin centers in oxpdf
    Ax (np.ndarray): Matrix with p(ki|Ni,xj) * dx
  
  Returns:
    float: Log likelihood
  """
  # Get the pdf from p:
  pdf = get_pdf(p, nx, noxpdf, oxpdf, odxpdf)
  # The likelihood for a cell i is
  # integral(p(ki|Ni,x)p(x|amp,mu,sig)dx)
  # We approximate by the sum
  # p(ki|Ni,xj)*dx is stored in Ax
  likecell = Ax @ pdf
  # The log likelihood for all cells is
  # sum of log likelihood for each cell
  ll = np.sum(np.log(likecell))
  return ll


def get_pdf(p, nx, noxpdf, oxpdf, odxpdf):
  """Get the pdf on x 1d using parameters

  Args:
    p (np.ndarray): Set of parameters to get pdf
    nx (int): Number of values in x.
    noxpdf (int): Number of values in x to evaluate pdf
    oxpdf (np.ndarray): bin centers for pdf
    odxpdf (float): Distance between 2 bin centers in oxpdf
    
  Returns:
    np.ndarray: The pdf evaluated on x
  """
  amp = p[2::3]
  mu = p[::3]
  sig = p[1::3]
  # We complete the list of parameters:
  amp_full = np.concatenate(([1 - np.sum(amp)], amp))
  # Build the oversampled pdf
  pdf_o = np.zeros(noxpdf)
  for ampk, muk, sigk in zip(amp_full, mu, sig):
    # Add the trunc_norm
    pdf_o += ampk * trunc_norm(oxpdf, odxpdf, muk, sigk)
  # I do the mean over the osampx
  pdf = np.mean(pdf_o.reshape(nx, -1), axis=1)
  return pdf


def extract_from_npz(output):
  """Return the same as gauss_mcmc but from a npz

  Args:
    output (str): file_path to numpy archive
  
  Returns:
    mu (np.ndarray): The average of parameters
    cov (np.ndarray): The covariance between parameters
    ox (np.ndarray): The nx * osampx bin centers
    oxpdf (np.ndarray): The osampxpdf * nx bin centers
    x (np.ndarray): The nx bin centers
    logprob (np.ndarray): The logprob values
    samples (np.ndarray): matrix with parameters for each sample
  """
  print("Reading")
  start = time.time()
  mcmc_res = np.load(output, allow_pickle=True)
  diags = mcmc_res['diagnostics']
  mu = diags.item().get('mu')
  cov = diags.item().get('cov')
  logprob_values = diags.item().get('logprob')
  samples = mcmc_res['samples']
  x = mcmc_res['x']
  ox = mcmc_res['ox']
  oxpdf = mcmc_res['oxpdf']
  print(f"Read. It took {(time.time() - start):.2f} seconds.")
  return(mu, cov, ox, oxpdf, x, logprob_values, samples)


def write_evidence(data, col_gene, mu, cov, ox, oxpdf, x, min_scale,
                   covis_scale, nis, logprob_function, output, seed,
                   xscale, target_sum):
  """Write evidence to file
  In order to compare multiple models we use an estimate of the log evidence

  Args:
    data (pandas.DataFrame): Data frame with 'nCount_RNA' and the gene of interest
    col_gene (str): Column label in `data` with the gene of interest
    mu (np.ndarray): The average of parameters
    cov (np.ndarray): The covariance between parameters
    ox (np.ndarray): The nx * osampx bin centers
    oxpdf (np.ndarray): The osampxpdf * nx bin centers
    x (np.ndarray): The nx bin centers
    min_scale (float): Minimal value of the scale of gaussians
    covis_scale (float): Scale factor to apply to covariance of parameters
                         to get random parameters in logevidence evaluation.
    nis (int): Size of sampling of random parameters in logevidence evaluation.
    logprob_function (function): Log probability of the distribution to sample.
    output (str): File path to write value of log evidence
    seed (int): Seed
    xscale (str): scale for the x-axis: Seurat (log(1+targetSum*X)) or log (log(X))
    target_sum (int): factor when Seurat scale is used: (log(1+targetSum*X))

  Raises:
      ValueError: If ox or oxpdf have incompatible size with x
  """
  nx = x.size
  dx = x[1] - x[0]
  xmin = x[0] - dx / 2
  xmax = x[-1] + dx / 2
  noxpdf = oxpdf.size
  odxpdf = oxpdf[1] - oxpdf[0]
  if ox.size % nx != 0:
    raise ValueError("The size of ox is not multiple of size of x.")
  if oxpdf.size % nx != 0:
    raise ValueError("The size of oxpdf is not multiple of size of x.")

  Ax = get_Ax(data, col_gene, nx, dx, ox, xscale, target_sum)

  nnorm = (mu.size + 1) // 3

  # We want to permut the gaussians (3 parameters for each)
  permutx = np.array([np.concatenate(([3 * i + np.arange(3) for i in v]))
                      for v in permutations(np.arange(nnorm))])

  # We calculate logevidence
  np.random.seed(seed)
  logevidence = samsam.covis(mu, covis_scale ** 2 * cov,
                             logprob_function, nsamples=nis,
                             pref=mu, wdist=1 / (np.diag(cov)), permutx=permutx,
                             xmin=xmin, xmax=xmax, nx=nx,
                             noxpdf=noxpdf, oxpdf=oxpdf, odxpdf=odxpdf,
                             min_scale=min_scale,
                             Ax=Ax)[2]['logevidence']

  # We write it into a file
  with open(output, 'w', encoding="utf-8") as fo:
    fo.write(f'{logevidence}\n')


def plots_from_pdf(x, pdf, title, output, data, col_gene, osampx, xscale, target_sum):
  """Compute plots for 1d

  Args:
    x (np.ndarray): The nx bin centers on which pdf has been evaluated
    pdf (np.ndarray): PDF values on x for each sample to plot
    title (str): Title for plots
    output (str): Path for output plots
    data (pd.DataFrame): Data frame with 'nCount_RNA' and the gene of interest
    col_gene (str): Column label in `data` with the gene of interest
    osampx (np.ndarray): bin centers for gaussian evaluation
    xscale (str): scale for the x-axis: Seurat (log(1+targetSum*X)) or log (log(X))
    target_sum (int): factor when Seurat scale is used: (log(1+targetSum*X))

  Raises:
      NotImplementedError: If xscale is not implemented
  """
  # We assume x is equally spaced
  dx = x[1] - x[0]
  xmin = x[0] - dx / 2
  xmax = x[-1] + dx / 2
  # Process the output to get prefix and suffix
  file_prefix, file_suffix = get_prefix_suffix(output)

  # 1, 2, 3 sigma:
  qm = np.array([(1 - lvl) / 2 for lvl in [0.6827, 0.9545, 0.9973]])
  qp = 1 - qm

  # Export the pdf mean and 1 sigma (68%) and median
  pdf_summary = pd.DataFrame(columns=['x', 'low', 'mean', 'high', 'median'])
  pdf_summary['x'] = x
  # Plot PDF with 1, 2, 3 sigma
  plt.figure()
  # 1,2,3 sigma
  for i, alpha in enumerate([0.2, 0.1, 0.05]):
    pm = np.quantile(pdf, qm[i], axis=0)
    pp = np.quantile(pdf, qp[i], axis=0)
    if i == 0:
      pdf_summary['low'] = pm
      pdf_summary['high'] = pp
    plt.fill_between(x, pm, pp, color='r', alpha=alpha, rasterized=True)
  # Mean
  plt.plot(x, np.mean(pdf, axis=0), 'r', lw=2, rasterized=True)
  pdf_summary['mean'] = np.mean(pdf, axis=0)
  # Median
  median = np.median(pdf, axis=0)
  pdf_summary['median'] = median
  plt.plot(x, median, 'k--', lw=2, rasterized=True)
  # KDE
  ks = data[col_gene]
  Ns = data['nCount_RNA']
  try:
    if xscale == 'Seurat':
      if target_sum == 0:
        target_sum = np.median(Ns)
      density = kde.gaussian_kde(np.log(1 + target_sum * ks / Ns))(x)
    elif xscale == 'log':
      Xs = ks / Ns
      # I artificially substitute all 0 by exp(xmin)
      Xs[Xs == 0] = np.exp(xmin)
      density = kde.gaussian_kde(np.log(Xs))(x)
    else:
      raise NotImplementedError("Only xscale Seurat and log are implemented.")
    plt.plot(x, density, color='green')
  except NotImplementedError:
    print("Could not compute density from input data.")
  plt.xlim(xmin, xmax)
  plt.ylim(0, )
  if xscale == 'Seurat':
    plt.xlabel(f'log(1 + {target_sum} * expression)')
  elif xscale == 'log':
    plt.xlabel('log(expression)')
  plt.ylabel('PDF')
  plt.title(title)
  plt.savefig(output, dpi=300)

  # Compute the posterior probability of each cell
  nx = x.size
  nox = nx * osampx
  ox = get_bins_centers(xmin, xmax, nox)
  Ax = get_Ax(data, col_gene, nx, dx, ox, xscale, target_sum)
  post_per_cell = Ax * np.mean(pdf, axis=0)
  # I renormalize each pdf:
  post_per_cell /= (np.sum(post_per_cell, axis=1)[:, None] * dx)
  # I average:
  post_all_cells = np.mean(post_per_cell, axis=0)
  # I plot:
  plt.plot(x, post_all_cells, color='orange')
  plt.savefig(f'{file_prefix}_with_posterior.{file_suffix}', dpi=300)

  # Export
  pdf_summary.to_csv(f'{file_prefix}_pdf.txt',
                     index=False,
                     header=True, sep='\t')

  # Plot some individual pdfs
  cmap = plt.cm.get_cmap('jet', 100)

  plt.figure()
  plt.title(f"{title}\n100 individual samples")
  for i, cur_pdf in enumerate(pdf[::(pdf.shape[0] // 100)][:100]):
    plt.plot(x, cur_pdf, color=cmap(i), alpha=0.3, rasterized=True)
  plt.savefig(f'{file_prefix}_individuals.{file_suffix}')

  # I plot posterior of first cells:
  ncells = min(50, post_per_cell.shape[0])
  cmap = plt.cm.get_cmap('jet', ncells)
  plt.figure()
  plt.title(f"{title}\nposterior of {ncells} cells")
  for i, cur_pdf in enumerate(post_per_cell[:ncells, :]):
    plt.plot(x, cur_pdf, color=cmap(i), alpha=0.3, rasterized=True)
  plt.savefig(f'{file_prefix}_posterior_individuals.{file_suffix}')

  average_per_cell = np.sum(post_per_cell * x, axis=1) * dx
  average_square_per_cell = np.sum(post_per_cell * x * x, axis=1) * dx
  var_per_cell = average_square_per_cell - average_per_cell * average_per_cell
  infered_pdf_per_cell = np.array([trunc_norm(x, dx, mu, np.sqrt(var))
                                   for mu, var in zip(average_per_cell, var_per_cell)])
  infered_pdf_global = np.mean(infered_pdf_per_cell, axis=0)
  # Plot
  plt.figure()
  # Mean
  plt.plot(x, np.mean(pdf, axis=0), 'r', lw=2, rasterized=True, label='baredSC mean pdf')
  # Median
  median = np.median(pdf, axis=0)
  plt.plot(x, median, 'k--', lw=2, rasterized=True, label='baredSC median pdf')
  # KDE
  ks = data[col_gene]
  Ns = data['nCount_RNA']
  try:
    if xscale == 'Seurat':
      density = kde.gaussian_kde(np.log(1 + target_sum * ks / Ns))(x)
    elif xscale == 'log':
      Xs = ks / Ns
      # I artificially substitute all 0 by exp(xmin)
      Xs[Xs == 0] = np.exp(xmin)
      density = kde.gaussian_kde(np.log(Xs))(x)
    else:
      raise NotImplementedError("Only xscale Seurat and log are implemented.")
    plt.plot(x, density, color='green', label='density from data')
  except NotImplementedError:
    print("Could not compute density from input data.")
  # Posterior using full pdf
  plt.plot(x, post_all_cells, color='orange', label='mean post pdf')
  # Only mean value per cell
  try:
    density = kde.gaussian_kde(average_per_cell)(x)
    plt.hist(average_per_cell, bins=40, color='magenta', label='histogram post mean', density=True)
    plt.plot(x, density, color='pink', label='density post mean')
  except Exception:  # pylint: disable=W0718
    print("Could not compute density from average per cell.")
  plt.plot(x, infered_pdf_global, color='purple', label='post gauss each cell')
  plt.legend()
  if xmax > 0:
    plt.xlim(xmin, 1.2 * xmax)
  else:
    plt.xlim(xmin, xmax + 1)
  plt.ylim(0, )
  if xscale == 'Seurat':
    plt.xlabel(f'log(1 + {target_sum} * expression)')
  elif xscale == 'log':
    plt.xlabel('log(expression)')
  plt.ylabel('PDF')
  plt.title(title)
  plt.savefig(f'{file_prefix}_posterior_andco.{file_suffix}', dpi=300)

  # Export the predicted average and sd for each cell
  per_cell = pd.DataFrame(columns=['mu', 'sd'])
  per_cell['mu'] = average_per_cell
  per_cell['sd'] = np.sqrt(var_per_cell)
  per_cell.to_csv(f'{file_prefix}_posterior_per_cell.txt',
                  index=False,
                  header=True, sep='\t')

  # Export values of mean found in each sample
  means = pdf @ (x * dx)
  np.savetxt(f'{file_prefix}_means.txt.gz', means)


def args_check(args):
  """Check args for 1d

  Args:
      args (Namespace): Namespace object to check

  Raises:
      ValueError: when arguments are not valid

  Returns:
      Namespace: Updated Namespace
  """
  if args.xscale == 'Seurat' and args.xmin < 0:
    raise ValueError("--xmin negative is not "
                     "compatible with --xscale Seurat "
                     "as it is log(1+targetSum*X).")
  if args.xscale == 'Seurat' and args.targetSum < 0:
    raise ValueError("--targetSum negative is not "
                     "compatible with --xscale Seurat "
                     "as it is log(1+targetSum*X).")
  if args.xscale == 'log' and args.xmax > 0:
    raise ValueError("--xmax positive is not "
                     "compatible with --xscale log "
                     "as it is log(X) and X < 1.")
  # Check the minScale
  dx = (args.xmax - args.xmin) / args.nx
  min_minScale = max(dx / args.osampxpdf * 2, dx / 2)
  if args.minScale < dx / 2:
    print(f"The value of minScale ({args.minScale})"
          f" must be above half the bin size ((xmax - xmin) / nx / 2 = {dx / 2}).\n"
          f"Minimum value will be used ({min_minScale}).")
    args.minScale = min_minScale
  if args.minScale < dx / args.osampxpdf * 2:
    print(f"The value of minScale ({args.minScale})"
          " must be above twice the bin size "
          "of pdf evaluation ((xmax - xmin) / (nx * osampxpdf) * 2 ="
          f" {dx / args.osampxpdf * 2}).\n"
          f"Minimum value will be used ({min_minScale}).")
    args.minScale = min_minScale
  return args
