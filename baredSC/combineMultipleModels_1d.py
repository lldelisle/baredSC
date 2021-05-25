# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
import numpy as np
import sys
import argparse
import os
import time
from tempfile import NamedTemporaryFile

# Local imports
from . common import get_data, get_bins_centers
from . oned import logprob, extract_from_npz, write_evidence, \
    get_pdf, plots_from_pdf
from baredSC._version import __version__


# Do some plots and exports in txt
def plot_combined(all_results, all_logevid, title, output, data, col_gene,
                  removeFirstSamples, nsampInPlot, pretty_bins, osampx,
                  xscale, target_sum):
  # Evaluate how many samples from each output you will get:
  all_logevid_delta = [le - max(all_logevid) for le in all_logevid]
  # Delog the log evid
  delog_all_logevid = np.array([np.exp(le) for le in all_logevid_delta])
  # Get the number of samples from each result
  evid_prop = delog_all_logevid / sum(delog_all_logevid)
  samples_per_results = np.array(evid_prop * nsampInPlot, dtype=int)
  # Compute the pdf array
  pdf = None
  for i, n_sample in enumerate(samples_per_results):
    if n_sample == 0:
      print(f"Using 0 sample from output {i}.")
    else:
      results = all_results[i]
      # (mu, cov, ox, oxpdf, x, logprob_values, samples)
      x = results[4]
      samples = results[6]
      # We assume x is equally spaced
      dx = x[1] - x[0]
      if pretty_bins is None:
        nx = x.size
        oxpdf = results[3]
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
      # Remove the first samples:
      if removeFirstSamples is None:
        samples = samples[(samples.shape[0] // 4):]
      else:
        samples = samples[removeFirstSamples:]
      print(f"Considering the last {samples.shape[0]} samples.")
      if samples.shape[0] < n_sample:
        new_nsampInPlot = nsampInPlot * samples.shape[0] // n_sample
        print(f"Could not sample {n_sample} samples from output {i} will plot {new_nsampInPlot} samples.")
        return(plot_combined(all_results, all_logevid, title, output, data, col_gene,
                             removeFirstSamples, new_nsampInPlot, pretty_bins, osampx,
                             xscale, target_sum))
      # Select n_sample
      samples_selected = np.unique(np.linspace(0, samples.shape[0] - 1, n_sample, dtype=int))
      samples = samples[samples_selected]
      print(f"Using {samples.shape[0]} samples from output {i}.")
      print("Computing pdf.")
      current_pdf = np.array([get_pdf(p, nx, noxpdf, oxpdf, odxpdf) for p in samples])
      if pdf is None:
        pdf = current_pdf
      else:
        pdf = np.concatenate((pdf, current_pdf))
  plots_from_pdf(x, pdf, title, output, data, col_gene, osampx, xscale, target_sum)


def parse_arguments(args=None):
  argp = argparse.ArgumentParser(
      description=("Combine mcmc results from multiple models to get a"
                   " mixture using logevidence to infer weights."))
  argprequired = argp.add_argument_group('Required arguments')
  argpopt_mcmc = argp.add_argument_group('Optional arguments used to run MCMC')
  argpopt_data = argp.add_argument_group('Optional arguments to select input data')
  argpopt_plot = argp.add_argument_group('Optional arguments to customize plots and text outputs')
  argpopt_loge = argp.add_argument_group('Optional arguments to evaluate logevidence')
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
  argprequired.add_argument('--outputs', default=None, required=True, nargs='+',
                            help="Ouput files basename (will be npz)"
                            " with different results of mcmc to combine.")
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
  argpopt_mcmc.add_argument('--seed', default=1, type=int,
                            help="Change seed for another output.")
  # Plot
  argprequired.add_argument('--figure', default=None, required=True,
                            help="Ouput figure basename.")
  argpopt_plot.add_argument('--title', default=None,
                            help="Title in figures.")
  argpopt_plot.add_argument('--removeFirstSamples', default=None, type=int,
                            help="Number of samples to ignore before making the plots"
                            " (default is nsampMCMC / 4).")
  argpopt_plot.add_argument('--nsampInPlot', default=100000, type=int,
                            help="Approximate number of samples to use in plots.")
  argpopt_plot.add_argument('--prettyBins', default=None, type=int,
                            help="Number of bins to use in plots (Default is nx).")
  # Evidence
  argpopt_loge.add_argument('--logevidences', default=None, nargs='+',
                            help="Ouput files of precalculated log evidence values."
                            "(if not provided will be calculated).")
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
  args = parse_arguments().parse_args(args)
  # Check incompatibilities:
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
  # Check the directory of args.figure is writtable:
  if os.path.dirname(args.figure) != '' and not os.access(os.path.dirname(args.figure), os.W_OK):
    raise Exception("The output figure is not writable no figure/export will be done,"
                    " check the directory exists.")

  # Load data
  print("Get raw data")
  start = time.time()
  data = get_data(args.input,
                  args.metadata1ColName, args.metadata1Values,
                  args.metadata2ColName, args.metadata2Values,
                  args.metadata3ColName, args.metadata3Values,
                  [args.geneColName])
  print(f"Got. It took {(time.time() - start):.2f} seconds.")

  # Check that output exists:
  for output in args.outputs:
    # Remove potential suffix:
    output = output.removesuffix('.npz')
    if not os.path.exists(f'{output}.npz'):
      raise Exception(f'{output}.npz does not exists.')
  if args.logevidences is not None:
    assert len(args.logevidences) == len(args.outputs), \
        "The number of logevid does not corresponds" \
        " to the number of outputs."
    for logevid_file in args.logevidences:
      if not os.path.exists(logevid_file):
        raise Exception(f'{logevid_file} does not exists.')

  # Get the results and the logevid:
  all_results = {}
  all_logevid_files = {}
  all_logevid = []
  for i, output in enumerate(args.outputs):
    # Remove potential suffix:
    output = output.removesuffix('.npz')
    all_results[i] = extract_from_npz(f'{output}.npz')
    # return(mu, cov, ox, oxpdf, x, logprob_values, samples)
    x = all_results[i][4]
    dx = x[1] - x[0]
    xmin = x[0] - dx / 2
    xmax = x[-1] + dx / 2
    try:
      assert all_results[i][2].size // all_results[i][4].size == args.osampx, \
          f"osampx value do not match what is in {output}.npz."
      assert all_results[i][3].size // all_results[i][4].size == args.osampxpdf, \
          f"osampxpdf value do not match what is in {output}.npz."
      assert x.size == args.nx, \
          "nx value do not match what is in output."
      assert abs(xmin - args.xmin) < dx, \
          "xmin value do not match what is in output."
      assert abs(xmax - args.xmax) < dx, \
          "xmax value do not match what is in output."
    except Exception as e:
      raise Exception(f"{e}\nChange args or rerun the MCMC.")
    # Get the logevid value if the file is provided:
    if args.logevidences is not None:
      all_logevid_files[i] = args.logevidences[i]
    else:
      # Else calculate it
      all_logevid_files[i] = NamedTemporaryFile(delete=False).name
      write_evidence(data, args.geneColName, *all_results[i][:5], args.minScale, args.coviscale,
                     args.nis, logprob, all_logevid_files[i], args.seed, args.xscale, args.targetSum)
    # Read the logevid value:
    try:
      with open(all_logevid_files[i], 'r') as f:
        logevid = f.read()
    except Exception as e:
      raise Exception(f"Could not read the file {all_logevid_files[i]}: {e}")
    try:
      logevid = float(logevid)
    except Exception:
      raise Exception(f"The value in the file {all_logevid_files[i]} is not a float.")
    all_logevid.append(logevid)
  # Make the plots:
  plot_combined(all_results, all_logevid, args.title, args.figure, data, args.geneColName,
                args.removeFirstSamples, args.nsampInPlot, args.prettyBins,
                args.osampx, args.xscale, args.targetSum)


if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
      args = ["--help"]
    main(args)
