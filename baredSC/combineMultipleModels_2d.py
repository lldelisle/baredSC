# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
import numpy as np
import sys
import argparse
import os
import time
from tempfile import NamedTemporaryFile

# Local imports
from . common import get_data, get_bins_centers, get_Neff
from . twod import logprob, extract_from_npz, write_evidence, \
    get_pdf, plots_from_pdf
from baredSC._version import __version__


# Do some plots and exports in txt
def plot_combined(all_results, all_logevid, title, output,
                  genex, geney, splitys,
                  removeFirstSamples, nsampInPlot,
                  log1pColorScale, pretty_bins_x, pretty_bins_y,
                  getPVal):
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
      # return(mu, cov, ox, oy, oxpdf, oypdf, x, y,
      #        logprob_values, samples)
      x = results[6]
      y = results[7]
      samples = results[9]
      # We assume x and y are equally spaced
      dx = x[1] - x[0]
      dy = y[1] - y[0]
      if pretty_bins_x is None:
        nx = x.size
        oxpdf = results[4]
        noxpdf = oxpdf.size
        # We assume oxpdf is equally spaced
        odxpdf = oxpdf[1] - oxpdf[0]
      else:
        xmin = x[0] - dx / 2
        xmax = x[-1] + dx / 2
        nx = pretty_bins_x
        x = get_bins_centers(xmin, xmax, nx)
        dx = x[1] - x[0]
        noxpdf = nx
        oxpdf = x
        odxpdf = dx
      if pretty_bins_y is None:
        ny = y.size
        oypdf = results[5]
        noypdf = oypdf.size
        # We assume oypdf is equally spaced
        odypdf = oypdf[1] - oypdf[0]
      else:
        ymin = y[0] - dy / 2
        ymax = y[-1] + dy / 2
        ny = pretty_bins_y
        y = get_bins_centers(ymin, ymax, ny)
        dy = y[1] - y[0]
        noypdf = ny
        oypdf = y
        odypdf = dy
      odxypdf = odxpdf * odypdf
      oxypdf = np.array(np.meshgrid(oxpdf, oypdf)).transpose(1, 2, 0)

      # Remove the first samples:
      if removeFirstSamples is None:
        samples = samples[(samples.shape[0] // 4):]
      else:
        samples = samples[removeFirstSamples:]
      print(f"Considering the last {samples.shape[0]} samples.")
      if getPVal:
        # We need to have independent samples:
        Neff = get_Neff(samples)
        print(f"Sampling {round(Neff)} independent samples.")
        # Select Neff samples
        samples_selected = np.unique(np.linspace(0, samples.shape[0] - 1, round(Neff), dtype=int))
        samples = samples[samples_selected]        
      if samples.shape[0] < n_sample:
        new_nsampInPlot = nsampInPlot * samples.shape[0] // n_sample
        print(f"Could not sample {n_sample} samples from output {i} will plot {new_nsampInPlot} samples.")
        return(plot_combined(all_results, all_logevid, title, output,
                             genex, geney, splitys,
                             removeFirstSamples, new_nsampInPlot,
                             log1pColorScale, pretty_bins_x, pretty_bins_y,
                             getPVal))
      # Select n_sample
      samples_selected = np.unique(np.linspace(0, samples.shape[0] - 1, n_sample, dtype=int))
      samples = samples[samples_selected]
      print(f"Using {samples.shape[0]} samples from output {i}.")
      print("Computing pdf.")
      current_pdf = np.array([get_pdf(p, nx, ny, noxpdf, noypdf, oxypdf, odxypdf) for p in samples])
      if pdf is None:
        pdf = current_pdf
      else:
        pdf = np.concatenate((pdf, current_pdf))
  if getPVal:
    Neff = pdf.shape[0]
  else:
    Neff = None
  plots_from_pdf(x, y, pdf, title, output,
                 genex, geney, log1pColorScale, splitys,
                 Neff)


def parse_arguments(args=None):
  argp = argparse.ArgumentParser(
      description=("Combine mcmc 2D results from multiple models to get a"
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
  argprequired.add_argument('--geneXColName', default=None, required=True,
                            help="Name of the column with gene counts for gene in x.")
  argprequired.add_argument('--geneYColName', default=None, required=True,
                            help="Name of the column with gene counts for gene in y.")
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
  argpopt_mcmc.add_argument('--nx', default=50, type=int,
                            help="Number of values in x to check how "
                            "your evaluated pdf is compatible with the model.")
  argpopt_mcmc.add_argument('--osampx', default=10, type=int,
                            help="Oversampling factor of x values when evaluating "
                            "pdf of Poisson distribution.")
  argpopt_mcmc.add_argument('--osampxpdf', default=4, type=int,
                            help="Oversampling factor of x values when evaluating "
                            "pdf at each step of the MCMC.")
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
  argpopt_mcmc.add_argument('--scale', default="Seurat", choices=['Seurat', 'log'],
                            help="scale for the x-axis and y-axis: Seurat (log(1+targetSum*X)) or log (log(X))")
  argpopt_mcmc.add_argument('--scalePrior', default=0.3, type=float,
                            help="Scale of the truncnorm used in the prior for "
                            "the correlation.")
  argpopt_mcmc.add_argument('--targetSum', default=10**4, type=float,
                            help="factor when Seurat scale is used: (log(1+targetSum*X)) (default is 10^4, use 0 for the median of nRNA_Counts)")
  argpopt_mcmc.add_argument('--seed', default=1, type=int,
                            help="Change seed for another output.")
  # Plot
  argprequired.add_argument('--figure', default=None, required=True,
                            help="Ouput figure basename.")
  argpopt_plot.add_argument('--title', default=None,
                            help="Title in figures.")
  argpopt_plot.add_argument('--splity', default=None, nargs='+', type=float,
                            help="Threshold value to plot the density for genex"
                            " for 2 categories in geney values.")
  argpopt_plot.add_argument('--removeFirstSamples', default=None, type=int,
                            help="Number of samples to ignore before making the plots"
                            " (default is nsampMCMC / 4).")
  argpopt_plot.add_argument('--nsampInPlot', default=100000, type=int,
                            help="Approximate number of samples to use in plots.")
  argpopt_plot.add_argument('--prettyBins', default=None, type=int,
                            help="Number of bins to use in plots (Default is nx).")
  argpopt_plot.add_argument('--prettyBinsx', default=None, type=int,
                            help="Number of bins to use in x in plots (Default is nx).")
  argpopt_plot.add_argument('--prettyBinsy', default=None, type=int,
                            help="Number of bins to use in y in plots (Default is ny).")
  argpopt_plot.add_argument('--log1pColorScale', action='store_true',
                            help="Use log1p color scale instead of linear color scale.")
  argpopt_plot.add_argument('--getPVal', action='store_true',
                            help="Use less samples to get an estimation of the p-value.")
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
  if args.scale == 'Seurat' and args.targetSum < 0:
    raise Exception("--targetSum negative is not "
                    "compatible with --xscale Seurat "
                    "as it is log(1+targetSum*X).")
  if args.scale == 'Seurat' and args.xmin < 0:
    raise Exception("--xmin negative is not "
                    "compatible with --scale Seurat "
                    "as it is log(1+targetSum*X).")
  if args.scale == 'log' and args.xmax > 0:
    raise Exception("--xmax positive is not "
                    "compatible with --scale log "
                    "as it is log(X) and X < 1.")
  if args.scale == 'Seurat' and args.ymin < 0:
    raise Exception("--ymin negative is not "
                    "compatible with --scale Seurat "
                    "as it is log(1+targetSum*X).")
  if args.scale == 'log' and args.ymax > 0:
    raise Exception("--ymax positive is not "
                    "compatible with --scale log "
                    "as it is log(X) and X < 1.")
  # Check the minScales
  dx = (args.xmax - args.xmin) / args.nx
  min_minScalex = max(dx / args.osampxpdf * 2, dx / 2)
  if args.minScalex < dx / 2:
    print(f"The value of minScale ({args.minScalex})"
          " must be above half the bin size in x"
          f"((xmax - xmin) / nx / 2 = {dx / 2}).\n"
          f"Minimum value will be used ({min_minScalex}).")
    args.minScalex = min_minScalex
  if args.minScalex < dx / args.osampxpdf * 2:
    print(f"The value of minScale ({args.minScalex})"
          " must be above twice the bin size "
          "of pdf evaluation in x((xmax - xmin) / (nx * osampxpdf) * 2 ="
          f" {dx / args.osampxpdf * 2}).\n"
          f"Minimum value will be used ({min_minScalex}).")
    args.minScalex = min_minScalex
  dy = (args.ymax - args.ymin) / args.ny
  min_minScaley = max(dy / args.osampypdf * 2, dy / 2)
  if args.minScaley < dy / 2:
    print(f"The value of minScale ({args.minScaley})"
          " must be above half the bin size in y"
          f"((ymax - ymin) / ny / 2 = {dy / 2}).\n"
          f"Minimum value will be used ({min_minScaley}).")
    args.minScaley = min_minScaley
  if args.minScaley < dy / args.osampypdf * 2:
    print(f"The value of minScale ({args.minScaley})"
          " must be above twice the bin size "
          "of pdf evaluation in y((ymax - ymin) / (ny * osampypdf) * 2 ="
          f" {dy / args.osampypdf * 2}).\n"
          f"Minimum value will be used ({min_minScaley}).")
    args.minScaley = min_minScaley

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
                  [args.geneXColName, args.geneYColName])
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
    # return(mu, cov, ox, oy, oxpdf, oypdf, x, y,
    #        logprob_values, samples)
    x = all_results[i][6]
    dx = x[1] - x[0]
    xmin = x[0] - dx / 2
    xmax = x[-1] + dx / 2
    y = all_results[i][7]
    dy = y[1] - y[0]
    ymin = y[0] - dy / 2
    ymax = y[-1] + dy / 2
    try:
      assert all_results[i][2].size // all_results[i][6].size == args.osampx, \
          "osampx value do not match what is in output."
      assert all_results[i][4].size // all_results[i][6].size == args.osampxpdf, \
          "osampxpdf value do not match what is in output."
      assert all_results[i][3].size // all_results[i][7].size == args.osampy, \
          "osampy value do not match what is in output."
      assert all_results[i][5].size // all_results[i][7].size == args.osampypdf, \
          "osampypdf value do not match what is in output."
      assert x.size == args.nx, \
          "nx value do not match what is in output."
      assert abs(xmin - args.xmin) < dx, \
          "xmin value do not match what is in output."
      assert abs(xmax - args.xmax) < dx, \
          "xmax value do not match what is in output."
      assert y.size == args.ny, \
          "ny value do not match what is in output."
      assert abs(ymin - args.ymin) < dy, \
          "ymin value do not match what is in output."
      assert abs(ymax - args.ymax) < dy, \
          "ymax value do not match what is in output."
    except Exception as e:
      raise Exception(f"{e}\nChange args or rerun the MCMC.")
    # Get the logevid value if the file is provided:
    if args.logevidences is not None:
      all_logevid_files[i] = args.logevidences[i]
    else:
      # Else calculate it
      all_logevid_files[i] = NamedTemporaryFile(delete=False).name
      write_evidence(data, args.geneXColName, args.geneYColName,
                     *all_results[i][:8], args.minScalex, args.minScaley,
                     args.scalePrior, args.coviscale,
                     args.nis, logprob, all_logevid_files[i], args.seed,
                     args.scale, args.targetSum)
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
  plot_combined(all_results, all_logevid, args.title, args.figure, args.geneXColName,
                args.geneYColName, args.splity, args.removeFirstSamples,
                args.nsampInPlot, args.log1pColorScale,
                args.prettyBinsx, args.prettyBinsy, args.getPVal)


if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
      args = ["--help"]
    main(args)
