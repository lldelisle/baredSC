# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
"""
Runner for combineMultipleModels_2d
"""
import sys
import os
from tempfile import NamedTemporaryFile

import numpy as np

# Local imports
from . common import get_bins_centers, \
  get_Neff, parse_arguments, get_data_from_args
from . twod import logprob, extract_from_npz, write_evidence, \
  get_pdf, plots_from_pdf, args_check


# Do some plots and exports in txt
def plot_combined(all_results, all_logevid, title, output,
                  genex, geney, splitys,
                  removeFirstSamples, nsampInPlot,
                  log1pColorScale, pretty_bins_x, pretty_bins_y,
                  getPVal):
  """Main function of combineMultipleModels_2d

  Args:
    all_results (dict of int: [np.ndarray]): Dictionary with results for each model
    all_logevid ([float]): Log evidence values for each model
    title (str): Title for plots
    output (str): Path for output plots
    genex (str): Column label in `data` with the gene of interest on x
    geney (str): Column label in `data` with the gene of interest on y
    splitys ([float]): Threshold value to plot the density for genex
                       for 2 categories in geney values.
    removeFirstSamples (int): Number of samples to ignore before making the plots
    nsampInPlot (int): Approximate number of samples to use in plots.
    log1pColorScale (bool): Use log1p color scale instead of linear color scale.
    pretty_bins_x (int): Number of bins to use in plots on x
    pretty_bins_y (int): Number of bins to use in plots on y
    getPVal (bool): Weither a p-value should be computed
  """
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
        print(f"Could not sample {n_sample} samples from output {i}"
              f" will plot {new_nsampInPlot} samples.")
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
  return None


def main(args=None):
  """Main function of combineMultipleModels_2d
  """
  args = parse_arguments('combineMultipleModels_2d').parse_args(args)
  # Update args and check
  args = args_check(args)

  # Check the directory of args.figure is writtable:
  if os.path.dirname(args.figure) != '' and not os.access(os.path.dirname(args.figure), os.W_OK):
    raise IOError("The output figure is not writable no figure/export will be done,"
                  " check the directory exists.")

  # Load data
  data = get_data_from_args(args, [args.geneXColName, args.geneYColName])

  # Check that output exists:
  for output in args.outputs:
    # Remove potential suffix:
    output = output.removesuffix('.npz')
    if not os.path.exists(f'{output}.npz'):
      raise ValueError(f'{output}.npz does not exists.')
  if args.logevidences is not None:
    assert len(args.logevidences) == len(args.outputs), \
        "The number of logevid does not corresponds" \
        " to the number of outputs."
    for logevid_file in args.logevidences:
      if not os.path.exists(logevid_file):
        raise ValueError(f'{logevid_file} does not exists.')

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
    except AssertionError as e:
      raise ValueError(f"{e}\nChange args or rerun the MCMC.") from e
    # Get the logevid value if the file is provided:
    if args.logevidences is not None:
      all_logevid_files[i] = args.logevidences[i]
    else:
      # Else calculate it
      with NamedTemporaryFile(delete=False) as temp:
        all_logevid_files[i] = temp.name
      write_evidence(data, args.geneXColName, args.geneYColName,
                     *all_results[i][:8], args.minScalex, args.minScaley,
                     args.scalePrior, args.coviscale,
                     args.nis, logprob, all_logevid_files[i], args.seed,
                     args.scale, args.targetSum)
    # Read the logevid value:
    try:
      with open(all_logevid_files[i], 'r', encoding="utf-8") as f:
        logevid = f.read()
    except Exception as e:
      raise ValueError(f"Could not read the file {all_logevid_files[i]}: {e}") from e
    try:
      logevid = float(logevid)
    except Exception as e:
      raise ValueError(f"The value in the file {all_logevid_files[i]} is not a float.") from e
    all_logevid.append(logevid)
  # Make the plots:
  plot_combined(all_results, all_logevid, args.title, args.figure, args.geneXColName,
                args.geneYColName, args.splity, args.removeFirstSamples,
                args.nsampInPlot, args.log1pColorScale,
                args.prettyBinsx, args.prettyBinsy, args.getPVal)


if __name__ == "__main__":
    arguments = None
    if len(sys.argv) == 1:
      arguments = ["--help"]
    main(arguments)
