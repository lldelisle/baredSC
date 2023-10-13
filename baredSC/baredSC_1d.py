# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
"""
Runner for baredSC_1d
"""

import sys
import os
from itertools import permutations
import time
import datetime
import numpy as np
from samsam import sam

# Local imports
from . common import permuted, get_Ax, \
  plot_QC, get_bins_centers, parse_arguments, checkNeff, args_check_baredSC, \
  get_data_from_args
from . oned import logprob, extract_from_npz, write_evidence, \
  get_pdf, plots_from_pdf, args_check


# Main function: running the mcmc:
def gauss_mcmc(data, col_gene,
               nx, osampx, osampxpdf,
               xmin, xmax, min_scale,
               xscale, target_sum, nnorm,
               nsamples_mcmc, nsamples_burn,
               nsplit_burn, T0_burn, output,
               seed):
  """Run MCMC with 1d Gaussians

  Args:
    data (pandas.DataFrame): Data frame with 'nCount_RNA' and the gene of interest
    col_gene (str): Column label in `data` with the gene of interest
    nx (int): Number of values in x to check how your evaluated pdf is compatible with the model.
    osampx (int): Oversampling factor of x values when evaluating pdf of Poisson distribution.
    osampxpdf (int): Oversampling factor of x values when evaluating pdf at each step of the MCMC.
    xmin (int): Minimum value to consider in x axis.
    xmax (int): Maximum value to consider in x axis.
    min_scale (float): Minimal value of the scale of gaussians
    xscale (str): scale for the x-axis: Seurat (log(1+targetSum*X)) or log (log(X))
    target_sum (int): factor when Seurat scale is used: (log(1+targetSum*X))
    nnorm (int): Number of gaussians to fit.
    nsamples_mcmc (int): Number of samplings (iteractions) of mcmc.
    nsamples_burn (int): Number of samplings (iteractions) in the burning phase of mcmc.
    nsplit_burn (int): Number of steps in the burning phase of mcmc.
    T0_burn (float): Initial temperature in the burning phase of mcmc.
    output (str): Ouput file basename (will be npz) with results of mcmc.
    seed (int): Seed
      
  Returns:
    mu (np.ndarray): The average of parameters
    cov (np.ndarray): The covariance between parameters
    ox (np.ndarray): The nx * osampx bin centers
    oxpdf (np.ndarray): The osampxpdf * nx bin centers
    x (np.ndarray): The nx bin centers
    logprob (np.ndarray): The logprob values
    samples (np.ndarray): matrix with parameters for each sample
  """
  # 1. Calculate the constants
  nox = nx * osampx
  noxpdf = nx * osampxpdf

  # Generate the grids using max and size
  x = get_bins_centers(xmin, xmax, nx)
  dx = x[1] - x[0]
  ox = get_bins_centers(xmin, xmax, nox)
  oxpdf = get_bins_centers(xmin, xmax, noxpdf)
  odxpdf = oxpdf[1] - oxpdf[0]

  # Evaluate p(ki|Ni,xj) * dx (independent of the pdf)
  Ax = get_Ax(data, col_gene, nx, dx, ox, xscale, target_sum)

  # We want to get all possible permutation of the gaussians (3 parameters for each)
  permutx = np.array([np.concatenate(([3 * i + np.arange(3) for i in v]))
                      for v in permutations(np.arange(nnorm))])

  np.random.seed(seed)
  # 2. Init params
  p0 = np.empty(3 * nnorm - 1)
  # We start with equal amplitudes for all gaussians:
  p0[2::3] = 1 / nnorm
  # We start loc (mean for gaussian) randomly between xmin and max
  p0[::3] = np.random.uniform(xmin, xmax, nnorm)
  # We start with all scales at half the range of x scale
  p0[1::3] = max(min_scale, (xmax - xmin) / 2)
  # MCMC
  # Burn phase to accurately evaluate the pref
  nsamples_split = nsamples_burn // nsplit_burn
  # Run a first mcmc from p0 with pref at p0
  # With a maximum temperature
  sk, dk = sam(p0, logprob, nsamples=nsamples_split,
               pref=p0, wdist=np.ones_like(p0), permutx=permutx,
               xmin=xmin, xmax=xmax, nx=nx,
               noxpdf=noxpdf, oxpdf=oxpdf, odxpdf=odxpdf,
               min_scale=min_scale,
               Ax=Ax, temp=T0_burn)
  # Use the last set of parameters as starting_point
  # but potentially permut it so it is the closest to
  # the new pref (the mean of the parameters in the last mcmc)
  starting_point = permuted(sk[-1], dk['mu'], 1 / (np.diag(dk['cov'])), permutx, 3)
  # Repeat the process to get a more accurate pref
  # Decreasing the temperature
  for ksplit in range(1, nsplit_burn):
    sk, dk = sam(starting_point, logprob, nsamples=nsamples_split,
                 mu0=dk['mu'], cov0=dk['cov'], scale0=dk['scale'],
                 pref=dk['mu'], wdist=1 / (np.diag(dk['cov'])), permutx=permutx,
                 xmin=xmin, xmax=xmax, nx=nx,
                 noxpdf=noxpdf, oxpdf=oxpdf, odxpdf=odxpdf,
                 min_scale=min_scale, Ax=Ax,
                 temp=T0_burn ** (1 - ksplit / nsplit_burn))
    starting_point = permuted(sk[-1], dk['mu'], 1 / (np.diag(dk['cov'])), permutx, 3)
  # Start actual MCMC with temperature at 1
  # Run MCMC
  samples, diags = sam(starting_point, logprob, nsamples=nsamples_mcmc,
                       mu0=dk['mu'], cov0=dk['cov'], scale0=dk['scale'],
                       pref=dk['mu'], wdist=1 / (np.diag(dk['cov'])), permutx=permutx,
                       xmin=xmin, xmax=xmax, nx=nx,
                       noxpdf=noxpdf, oxpdf=oxpdf, odxpdf=odxpdf,
                       min_scale=min_scale, Ax=Ax)
  # save data to output
  print("Saving")
  start = time.time()
  np.savez_compressed(output,
                      samples=samples,
                      diagnostics=diags,
                      x=x,
                      ox=ox,
                      oxpdf=oxpdf,
                      nsamples_burn=nsamples_burn,
                      nsplit_burn=nsplit_burn,
                      T0_burn=T0_burn)
  print(f"Saved. It took {(time.time() - start):.2f} seconds.")

  return(diags['mu'], diags['cov'], ox, oxpdf, x, diags['logprob'], samples)


# Do some plots and exports in txt
def plot(oxpdf, x, logprob_values, samples, title, output, data, col_gene,
         removeFirstSamples, nsampInPlot, pretty_bins, osampx, xscale, target_sum):
  """Plot baredSC_1d outputs

  Args:
    oxpdf (np.ndarray): bin centers for pdf
    x (np.ndarray): bin centers for likelihood
    logprob_values (np.ndarray): logprob values
    samples (np.ndarray): matrix with parameter values for each sample
    title (str): Title for plots
    output (str): Path for output plots
    data (pd.DataFrame): Data frame with 'nCount_RNA' and the gene of interest
    col_gene (str): Column label in `data` with the gene of interest
    removeFirstSamples (int): Number of samples to ignore before making the plots
    nsampInPlot (int): Approximate number of samples to use in plots.
    pretty_bins (int): Number of bins to use in plots
    osampx (np.ndarray): bin centers for gaussian evaluation
    xscale (str): scale for the x-axis: Seurat (log(1+targetSum*X)) or log (log(X))
    target_sum (int): factor when Seurat scale is used: (log(1+targetSum*X))
  """
  # We assume x is equally spaced
  dx = x[1] - x[0]
  if pretty_bins is None:
    nx = x.size
    noxpdf = oxpdf.size
    # We assume oxpdf is equally spaced
    odxpdf = oxpdf[1] - oxpdf[0]
  else:
    xmax = x[-1] + dx / 2
    xmin = x[0] - dx / 2
    nx = pretty_bins
    x = get_bins_centers(xmin, xmax, nx)
    dx = x[1] - x[0]
    noxpdf = nx
    oxpdf = x
    odxpdf = dx

  nnorm = (samples.shape[1] + 1) // 3
  p_names = [f'{pn}{i}' for i in range(nnorm) for pn in ['amp', 'mu', 'scale']][1:]

  samples, _ = plot_QC(logprob_values, samples, title, output,
                          removeFirstSamples, nsampInPlot,
                          p_names)

  print("Computing pdf.")
  start = time.time()
  pdf = np.array([get_pdf(p, nx, noxpdf, oxpdf, odxpdf) for p in samples])
  print(f"Done. It took {datetime.timedelta(seconds=int(time.time() - start))}.")

  plots_from_pdf(x, pdf, title, output, data, col_gene, osampx, xscale, target_sum)


def main(args=None):
  """Main function of baredSC_1d
  """
  original_args = sys.argv[1:]
  args = parse_arguments('baredSC_1d').parse_args(args)
  args = args_check_baredSC(args)
  # Update args and check
  args = args_check(args)
  # Load data
  data = get_data_from_args(args, [args.geneColName])

  # Check the directory of args.output is writtable:
  if os.path.dirname(args.output) != '' and not os.access(os.path.dirname(args.output), os.W_OK):
    raise OSError("The output is not writable,"
                  " check the directory exists.")
  if not os.path.exists(f'{args.output}.npz') or args.force:
    # Run mcmc:
    print("Run MCMC")
    start = time.time()
    results = gauss_mcmc(data, args.geneColName,
                         args.nx, args.osampx, args.osampxpdf,
                         args.xmin, args.xmax, args.minScale,
                         args.xscale, args.targetSum, args.nnorm,
                         args.nsampMCMC, args.nsampBurnMCMC,
                         args.nsplitBurnMCMC, args.T0BurnMCMC,
                         args.output, args.seed)
    print(f"MCMC + Save took {datetime.timedelta(seconds=int(time.time() - start))}.")
  else:
    results = extract_from_npz(f'{args.output}.npz')
    # return(mu, cov, ox, oxpdf, x, logprob_values, samples)
    x = results[4]
    dx = x[1] - x[0]
    xmin = x[0] - dx / 2
    xmax = x[-1] + dx / 2
    try:
      assert (results[0].size + 1) // 3 == args.nnorm, \
          "nnorm value do not match what is in output."
      assert results[2].size // results[4].size == args.osampx, \
          "osampx value do not match what is in output."
      assert results[3].size // results[4].size == args.osampxpdf, \
          "osampxpdf value do not match what is in output."
      assert x.size == args.nx, \
          "nx value do not match what is in output."
      assert abs(xmin - args.xmin) < dx, \
          "xmin value do not match what is in output."
      assert abs(xmax - args.xmax) < dx, \
          "xmax value do not match what is in output."
      assert (results[5].size == args.nsampMCMC + 1) or (args.minNeff is not None), \
          "nsampMCMC value do not match what is in output."
    except AssertionError as e:
      raise ValueError(f"Ouput file already exists and {e}"
                       " Use --force to rerun MCMC.") from e
    if args.minNeff is not None:
      args.nsampMCMC = results[5].size - 1

  if args.figure is not None:
    # Check the directory of args.figure is writtable:
    if os.path.dirname(args.figure) != '' and not os.access(os.path.dirname(args.figure), os.W_OK):
      print("The output figure is not writable no figure/export will be done,"
            " check the directory exists.")
    else:
      # Make the plots:
      plot(*results[3:], args.title, args.figure, data, args.geneColName,
           args.removeFirstSamples, args.nsampInPlot, args.prettyBins,
           args.osampx, args.xscale, args.targetSum)
      new_args = checkNeff(args, original_args)
      if new_args is not None:
        return main(new_args)
  if args.logevidence is not None:
    # Check the directory of args.logevidence is writtable:
    if os.path.dirname(args.logevidence) != '' and not os.access(os.path.dirname(args.logevidence),
                                                                 os.W_OK):
      print("The output logeveidence is not writable"
            " check the directory exists.")
    else:
      # Evaluate the logevid:
      write_evidence(data, args.geneColName, *results[:5], args.minScale, args.coviscale,
                     args.nis, logprob, args.logevidence, args.seed, args.xscale, args.targetSum)
  return None


if __name__ == "__main__":
    arguments = None
    if len(sys.argv) == 1:
      arguments = ["--help"]
    main(arguments)
