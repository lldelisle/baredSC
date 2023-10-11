# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
"""
Runner for combineMultipleModels_1d
"""
import sys
import os
from tempfile import NamedTemporaryFile

import numpy as np

# Local imports
from . common import get_bins_centers, parse_arguments, \
  get_data_from_args
from . oned import logprob, extract_from_npz, write_evidence, \
    get_pdf, plots_from_pdf, args_check


# Do some plots and exports in txt
def plot_combined(all_results, all_logevid, title, output, data, col_gene,
                  removeFirstSamples, nsampInPlot, pretty_bins, osampx,
                  xscale, target_sum):
  """Main function of combineMultipleModels_1d

  Args:
    all_results (dict of int: [np.ndarray]): Dictionary with results for each model
    all_logevid ([float]): Log evidence values for each model
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
        print(f"Could not sample {n_sample} samples from output {i}"
              f" will plot {new_nsampInPlot} samples.")
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
  return None


def main(args=None):
  """Main function of combineMultipleModels_1d
  """
  args = parse_arguments('combineMultipleModels_1d').parse_args(args)
  # Update args and check
  args = args_check(args)
  # Check the directory of args.figure is writtable:
  if os.path.dirname(args.figure) != '' and not os.access(os.path.dirname(args.figure), os.W_OK):
    raise IOError("The output figure is not writable no figure/export will be done,"
                  " check the directory exists.")

  # Load data
  data = get_data_from_args(args, [args.geneColName])

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
    except AssertionError as e:
      raise ValueError(f"{e}\nChange args or rerun the MCMC.") from e
    # Get the logevid value if the file is provided:
    if args.logevidences is not None:
      all_logevid_files[i] = args.logevidences[i]
    else:
      # Else calculate it
      with NamedTemporaryFile(delete=False) as temp:
        all_logevid_files[i] = temp.name
      write_evidence(data, args.geneColName, *all_results[i][:5],
                     args.minScale, args.coviscale,
                     args.nis, logprob, all_logevid_files[i],
                     args.seed, args.xscale, args.targetSum)
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
  plot_combined(all_results, all_logevid, args.title, args.figure, data, args.geneColName,
                args.removeFirstSamples, args.nsampInPlot, args.prettyBins,
                args.osampx, args.xscale, args.targetSum)


if __name__ == "__main__":
    arguments = None
    if len(sys.argv) == 1:
      arguments = ["--help"]
    main(arguments)
