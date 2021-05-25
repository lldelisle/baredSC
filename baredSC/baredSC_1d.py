# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
import numpy as np
import sys
import argparse
import os
from itertools import permutations
import time
import datetime
from tempfile import NamedTemporaryFile
from shutil import copy
from samsam import sam

# Local imports
from . common import get_data, permuted, get_Ax, get_prefix_suffix, \
    plot_QC, get_bins_centers
from . oned import logprob, extract_from_npz, write_evidence, \
    get_pdf, plots_from_pdf
from . _version import __version__


# Main function: running the mcmc:
def gauss_mcmc(data, col_gene,
               nx, osampx, osampxpdf,
               xmin, xmax, min_scale,
               xscale, target_sum, nnorm,
               nsamples_mcmc, nsamples_burn,
               nsplit_burn, T0_burn, output,
               seed):

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
  Ax = get_Ax(data, col_gene, nx, x, dx, ox, xscale, target_sum)

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
               xmin=xmin, xmax=xmax, nx=nx, dx=dx,
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
                 xmin=xmin, xmax=xmax, nx=nx, dx=dx,
                 noxpdf=noxpdf, oxpdf=oxpdf, odxpdf=odxpdf,
                 min_scale=min_scale, Ax=Ax,
                 temp=T0_burn ** (1 - ksplit / nsplit_burn))
    starting_point = permuted(sk[-1], dk['mu'], 1 / (np.diag(dk['cov'])), permutx, 3)
  # Start actual MCMC with temperature at 1
  # Run MCMC
  samples, diags = sam(starting_point, logprob, nsamples=nsamples_mcmc,
                       mu0=dk['mu'], cov0=dk['cov'], scale0=dk['scale'],
                       pref=dk['mu'], wdist=1 / (np.diag(dk['cov'])), permutx=permutx,
                       xmin=xmin, xmax=xmax, nx=nx, dx=dx,
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

  samples, Neff = plot_QC(logprob_values, samples, title, output,
                          removeFirstSamples, nsampInPlot,
                          p_names)

  print("Computing pdf.")
  start = time.time()
  pdf = np.array([get_pdf(p, nx, noxpdf, oxpdf, odxpdf) for p in samples])
  print(f"Done. It took {datetime.timedelta(seconds=int(time.time() - start))}.")

  plots_from_pdf(x, pdf, title, output, data, col_gene, osampx, xscale, target_sum)


def parse_arguments(args=None):
  argp = argparse.ArgumentParser(
      description=("Run mcmc to get the pdf for a given gene using a"
                   " normal distributions."))
  argprequired = argp.add_argument_group('Required arguments')
  argpopt_data = argp.add_argument_group('Optional arguments to select input data')
  argpopt_mcmc = argp.add_argument_group('Optional arguments to run MCMC')
  argpopt_plot = argp.add_argument_group('Optional arguments to get plots and text outputs')
  argpopt_loge = argp.add_argument_group('Optional arguments to get logevidence')
  # Get data:
  argprequired.add_argument('--input', default=None, required=True,
                            help="Input table with one line per cell"
                            " columns with raw counts and one column"
                            " nCount_RNA with total number of UMI per cell"
                            " optionally other meta data to filter.")
  argprequired.add_argument('--geneColName', default=None, required=True,
                            help="Name of the column with gene counts.")
  argpopt_data.add_argument('--metadata1ColName', default=None,
                            help="Name of the column with metadata1 to filter.")
  argpopt_data.add_argument('--metadata1Values', default=None,
                            help="Comma separated values for metadata1.")
  argpopt_data.add_argument('--metadata2ColName', default=None,
                            help="Name of the column with metadata2 to filter.")
  argpopt_data.add_argument('--metadata2Values', default=None,
                            help="Comma separated values for metadata2.")
  argpopt_data.add_argument('--metadata3ColName', default=None,
                            help="Name of the column with metadata3 to filter.")
  argpopt_data.add_argument('--metadata3Values', default=None,
                            help="Comma separated values for metadata3.")
  # MCMC
  argpopt_mcmc.add_argument('--xmin', default=0, type=float,
                            help="Minimum value to consider in x axis.")
  argpopt_mcmc.add_argument('--xmax', default=2.5, type=float,
                            help="Maximum value to consider in x axis.")
  argpopt_mcmc.add_argument('--xscale', default="Seurat", choices=['Seurat', 'log'],
                            help="scale for the x-axis: Seurat (log(1+targetSum*X)) or log (log(X))")
  argpopt_mcmc.add_argument('--targetSum', default=10**4, type=float,
                            help="factor when Seurat scale is used: (log(1+targetSum*X)) (default is 10^4, use 0 for the median of nRNA_Counts)")
  argpopt_mcmc.add_argument('--nx', default=100, type=int,
                            help="Number of values in x to check how "
                            "your evaluated pdf is compatible with the model.")
  argpopt_mcmc.add_argument('--osampx', default=10, type=int,
                            help="Oversampling factor of x values when evaluating "
                            "pdf of Poisson distribution.")
  argpopt_mcmc.add_argument('--osampxpdf', default=5, type=int,
                            help="Oversampling factor of x values when evaluating "
                            "pdf at each step of the MCMC.")
  argpopt_mcmc.add_argument('--minScale', default=0.1, type=float,
                            help="Minimal value of the scale of gaussians"
                            " (Default is 0.1 but cannot be smaller than "
                            "max of twice the bin size of pdf evaluation"
                            " and half the bin size).")
  argpopt_mcmc.add_argument('--nnorm', default=2, type=int,
                            help="Number of gaussian to fit.")
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
                            "burning phase of mcmc (>1).")
  argpopt_mcmc.add_argument('--seed', default=1, type=int,
                            help="Change seed for another output.")
  argpopt_mcmc.add_argument('--minNeff', default=None, type=float,
                            help="Will redo the MCMC with 10 times more samples until "
                            "Neff is greater that this value (Default is not set so will not rerun MCMC).")
  # To save/get MCMC
  argpopt_mcmc.add_argument('--force', default=None, action='store_true',
                            help="Force to redo the mcmc even if output exists.")
  argprequired.add_argument('--output', default=None, required=True,
                            help="Ouput file basename (will be npz)"
                            " with results of mcmc.")
  # Plot
  argpopt_plot.add_argument('--figure', default=None,
                            help="Ouput figure filename.")
  argpopt_plot.add_argument('--title', default=None,
                            help="Title in figures.")
  argpopt_plot.add_argument('--removeFirstSamples', default=None, type=int,
                            help="Number of samples to ignore before making the plots"
                            " (default is nsampMCMC / 4).")
  argpopt_plot.add_argument('--nsampInPlot', default=100000, type=int,
                            help="Approximate number of samples to use in plots.")
  argpopt_plot.add_argument('--prettyBins', default=None, type=int,
                            help="Number of bins to use in plots (Default is nx).")
  # Calculate evidence
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
  return(argp)


def main(args=None):
  original_args = sys.argv[1:]
  args = parse_arguments().parse_args(args)
  # Put default:
  if args.nsampBurnMCMC is None:
    args.nsampBurnMCMC = args.nsampMCMC // 4
  if args.removeFirstSamples is None:
    args.removeFirstSamples = args.nsampMCMC // 4
  # Remove potential suffix:
  args.output = args.output.removesuffix('.npz')
  # Check incompatibilities:
  if args.minNeff is not None and args.figure is None:
    raise Exception("--minNeff requires --figure to be set.")
  if args.xscale == 'Seurat' and args.xmin < 0:
    raise Exception("--xmin negative is not "
                    "compatible with --xscale Seurat "
                    "as it is log(1+targetSum*X).")
  if args.xscale == 'Seurat' and args.targetSum < 0:
    raise Exception("--targetSum negative is not "
                    "compatible with --xscale Seurat "
                    "as it is log(1+targetSum*X).")
  if args.xscale == 'log' and args.xmax > 0:
    raise Exception("--xmax positive is not "
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

  # Load data
  print("Get raw data")
  start = time.time()
  data = get_data(args.input,
                  args.metadata1ColName, args.metadata1Values,
                  args.metadata2ColName, args.metadata2Values,
                  args.metadata3ColName, args.metadata3Values,
                  [args.geneColName])
  print(f"Got. It took {(time.time() - start):.2f} seconds.")

  # Check the directory of args.output is writtable:
  if os.path.dirname(args.output) != '' and not os.access(os.path.dirname(args.output), os.W_OK):
    raise Exception("The output is not writable,"
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
    except Exception as e:
      raise Exception(f"Ouput file already exists and {e}"
                      " Use --force to rerun MCMC.")
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
      if args.minNeff is not None:
        # Process the output to get prefix and suffix
        file_prefix, file_suffix = get_prefix_suffix(args.figure)
        with open(f'{file_prefix}_neff.txt', 'r') as f:
          neff = float(f.read())
        if neff < args.minNeff:
          print("Neff is below the minimum required.")
          temp_file = NamedTemporaryFile(delete=False, suffix=".npz").name
          print(f"Results are moved to {temp_file}")
          copy(f'{args.output}.npz', temp_file)
          os.remove(f'{args.output}.npz')
          # Change the args to multiply by 10 the nb of samples
          new_args = original_args
          if '--nsampMCMC' in new_args:
            i = [i for i, v in enumerate(new_args) if v == '--nsampMCMC'][0]
            new_args[i + 1] += '0'
          else:
            new_args.append('--nsampMCMC')
            new_args.append(f'{args.nsampMCMC}0')
          return main(new_args)
  if args.logevidence is not None:
    # Check the directory of args.logevidence is writtable:
    if os.path.dirname(args.logevidence) != '' and not os.access(os.path.dirname(args.logevidence), os.W_OK):
      print("The output logeveidence is not writable"
            " check the directory exists.")
    else:
      # Evaluate the logevid:
      write_evidence(data, args.geneColName, *results[:5], args.minScale, args.coviscale,
                     args.nis, logprob, args.logevidence, args.seed, args.xscale, args.targetSum)


if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
      args = ["--help"]
    main(args)
