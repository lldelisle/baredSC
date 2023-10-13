# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
"""
Runner for baredSC_2d
"""

import sys
import os
from itertools import permutations
import time
import datetime
import numpy as np
import pandas as pd
from samsam import sam

# Local imports
from . common import permuted, get_Ax, get_prefix_suffix, \
  plot_QC, get_bins_centers, parse_arguments, checkNeff, args_check_baredSC, \
  get_data_from_args
from . twod import logprob, extract_from_npz, write_evidence, \
  get_pdf, plots_from_pdf, args_check


# Main function: running the mcmc:
def gauss_mcmc(data, genex, geney,
               nx, osampx, osampxpdf, xmin, xmax,
               ny, osampy, osampypdf, ymin, ymax,
               min_scale_x, min_scale_y,
               scale_prior, scale, target_sum, nnorm,
               nsamples_mcmc, nsamples_burn,
               nsplit_burn, T0_burn, output,
               seed):
  """Run MCMC with 2d Gaussians

  Args:
    data (pandas.DataFrame): Data frame with 'nCount_RNA' and the genes of interest
    genex (str): Column label in `data` with the gene of interest on x
    geney (str): Column label in `data` with the gene of interest on y
    nx (int): Number of values in x to check how your evaluated pdf is compatible with the model.
    osampx (int): Oversampling factor of x values when evaluating pdf of Poisson distribution.
    osampxpdf (int): Oversampling factor of x values when evaluating pdf at each step of the MCMC.
    xmin (int): Minimum value to consider in x axis.
    xmax (int): Maximum value to consider in x axis.
    ny (int): Number of values in y to check how your evaluated pdf is compatible with the model.
    osampy (int): Oversampling factor of y values when evaluating pdf of Poisson distribution.
    osampypdf (int): Oversampling factor of y values when evaluating pdf at each step of the MCMC.
    ymin (int): Minimum value to consider in y axis.
    ymax (int): Maximum value to consider in y axis.
    min_scale_x (float): Minimal value of the scale of gaussians in x
    min_scale_y (float): Minimal value of the scale of gaussians in y
    scale_prior (float): Scale of the truncnorm used in the prior for the correlation.
    scale (str): scale for the x-axis and y-axis: Seurat (log(1+targetSum*X)) or log (log(X))
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
    oy (np.ndarray): The ny * osampy bin centers
    oxpdf (np.ndarray): The osampxpdf * nx bin centers
    oypdf (np.ndarray): The osampypdf * ny bin centers
    x (np.ndarray): The nx bin centers
    y (np.ndarray): The ny bin centers
    logprob (np.ndarray): The logprob values
    samples (np.ndarray): The samples
  """
  # 1. Calculate the constants
  nox = nx * osampx
  noy = ny * osampy

  noxpdf = nx * osampxpdf
  noypdf = ny * osampypdf

  # Define nxy
  nxy = nx * ny

  # Generate the grids using max and size
  x = get_bins_centers(xmin, xmax, nx)
  dx = x[1] - x[0]
  ox = get_bins_centers(xmin, xmax, nox)
  oxpdf = get_bins_centers(xmin, xmax, noxpdf)
  odxpdf = oxpdf[1] - oxpdf[0]

  y = get_bins_centers(ymin, ymax, ny)
  dy = y[1] - y[0]
  oy = get_bins_centers(ymin, ymax, noy)
  oypdf = get_bins_centers(ymin, ymax, noypdf)
  odypdf = oypdf[1] - oypdf[0]

  odxypdf = odxpdf * odypdf
  oxypdf = np.array(np.meshgrid(oxpdf, oypdf)).transpose(1, 2, 0)

  # Evaluate p(ki|Ni,xj) * dx (independent of the pdf)
  Ax = get_Ax(data, genex, nx, dx, ox, scale, target_sum)
  Ay = get_Ax(data, geney, ny, dy, oy, scale, target_sum)

  # Evaluate p(ki,li|Ni,xj,yk)*dx*dy
  Axy = Ax[:, None, :] * Ay[:, :, None]
  Axy = Axy.reshape(-1, nxy)

  # We want to get all possible permutation of the gaussians (6 parameters for each)
  permutxy = np.array([np.concatenate(([6 * i + np.arange(6) for i in v]))
                       for v in permutations(np.arange(nnorm))])

  np.random.seed(seed)

  # 2. Init params
  p0 = np.empty(6 * nnorm - 1)
  # The first amplitude is 1 - the sum of others

  # We start with equal amplitudes for all gaussians:
  p0[5::6] = 1 / nnorm
  # We start loc (mean for gaussian) randomy between xmin and max
  p0[::6] = np.random.uniform(xmin, xmax, nnorm)
  p0[1::6] = np.random.uniform(ymin, ymax, nnorm)
  # We start with all scales at half the axis
  p0[2::6] = max(max(xmax - xmin, ymax - ymin) / 2, min_scale_x)
  p0[3::6] = max(max(xmax - xmin, ymax - ymin) / 2, min_scale_y)
  # And correlation at 0
  p0[4::6] = 0.0

  # MCMC
  # Burn phase to accuratey evaluate the pref
  nsamples_split = nsamples_burn // nsplit_burn
  # Run a first mcmc from p0 with pref at p0
  sk, dk = sam(p0, logprob, nsamples=nsamples_split,
               pref=p0, wdist=np.ones_like(p0),
               permutxy=permutxy,
               xmin=xmin, xmax=xmax, nx=nx,
               ymax=ymax, ymin=ymin, ny=ny,
               nxy=nxy, noxpdf=noxpdf, noypdf=noypdf,
               oxypdf=oxypdf, odxypdf=odxypdf,
               min_scale_x=min_scale_x, min_scale_y=min_scale_y,
               scale_prior=scale_prior,
               Axy=Axy, temp=T0_burn)
  # Use the last set of parameters as starting_point
  # but potentialy permut it so it is the closest to
  # the new pref (the mean of the parameters in the last mcmc)
  starting_point = permuted(sk[-1], dk['mu'], 1 / (np.diag(dk['cov'])),
                            permutxy, 6)
  # Repeat the process to get a more accurate pref
  for ksplit in range(1, nsplit_burn):
    sk, dk = sam(starting_point, logprob, nsamples=nsamples_split,
                 mu0=dk['mu'], cov0=dk['cov'], scale0=dk['scale'],
                 pref=dk['mu'], wdist=1 / (np.diag(dk['cov'])),
                 permutxy=permutxy,
                 xmin=xmin, xmax=xmax, nx=nx,
                 ymax=ymax, ymin=ymin, ny=ny,
                 nxy=nxy, noxpdf=noxpdf, noypdf=noypdf,
                 oxypdf=oxypdf, odxypdf=odxypdf,
                 min_scale_x=min_scale_x, min_scale_y=min_scale_y,
                 scale_prior=scale_prior,
                 Axy=Axy, temp=T0_burn ** (1 - ksplit / nsplit_burn))
    starting_point = permuted(sk[-1], dk['mu'], 1 / (np.diag(dk['cov'])),
                              permutxy, 6)
  # Start actual MCMC with temperature at 1
  # Run MCMC
  samples, diags = sam(starting_point, logprob, nsamples=nsamples_mcmc,
                       mu0=dk['mu'], cov0=dk['cov'], scale0=dk['scale'],
                       pref=dk['mu'], wdist=1 / (np.diag(dk['cov'])),
                       permutxy=permutxy,
                       xmin=xmin, xmax=xmax, nx=nx,
                       ymax=ymax, ymin=ymin, ny=ny,
                       nxy=nxy, noxpdf=noxpdf, noypdf=noypdf,
                       oxypdf=oxypdf, odxypdf=odxypdf,
                       min_scale_x=min_scale_x, min_scale_y=min_scale_y,
                       scale_prior=scale_prior,
                       Axy=Axy)
  # save data to output
  print("Saving")
  start = time.time()
  np.savez_compressed(output,
                      samples=samples,
                      diagnostics=diags,
                      x=x,
                      ox=ox,
                      oxpdf=oxpdf,
                      y=y,
                      oy=oy,
                      oypdf=oypdf,
                      nsamples_burn=nsamples_burn,
                      nsplit_burn=nsplit_burn,
                      T0_burn=T0_burn)
  print(f"Saved. It took {(time.time() - start):.2f} seconds.")

  return(diags['mu'], diags['cov'], ox, oy,
         oxpdf, oypdf, x, y,
         diags['logprob'], samples)


# Do some plots and exports in txt
def plot(oxpdf, oypdf, x, y, logprob_values, samples,
         title, output, genex, geney, splitys,
         removeFirstSamples, nsampInPlot,
         log1pColorScale, pretty_bins_x, pretty_bins_y):
  """Plot baredSC_2d outputs

  Args:
    oxpdf (np.ndarray): bin centers for pdf on x
    oypdf (np.ndarray): bin centers for pdf on y
    x (np.ndarray): bin centers for likelihood on x
    y (np.ndarray): bin centers for likelihood on y
    logprob_values (np.ndarray): logprob values
    samples (np.ndarray): matrix with parameter values for each sample
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
  """
  # We assume x and y are equally spaced
  dx = x[1] - x[0]
  dy = y[1] - y[0]
  if pretty_bins_x is None:
    nx = x.size
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

  p_names = [f'{pn}{i}' for i in range((samples.shape[1] + 1) // 6)
             for pn in ['xy_amp', 'xy_mux', 'xy_muy', 'xy_scalex', 'xy_scaley', 'xy_corr']][1:]

  samples, Neff = plot_QC(logprob_values, samples, title, output,
                          removeFirstSamples, nsampInPlot,
                          p_names)
  print("Computing pdf.")
  start = time.time()
  pdf = np.array([get_pdf(p, nx, ny, noxpdf, noypdf, oxypdf, odxypdf) for p in samples])
  print(f"Done. It took {datetime.timedelta(seconds=int(time.time() - start))}.")

  # Process the output to get prefix and suffix
  file_prefix, _ = get_prefix_suffix(output)

  # 1, 2, 3 sigma:
  qm = np.array([(1 - lvl) / 2 for lvl in [0.6827, 0.9545, 0.9973]])
  qp = 1 - qm

  # Export the parameters median and 1 sigma
  p_summary = pd.DataFrame(columns=['name', 'low', 'median', 'high'])
  p_summary['name'] = [f'{pn}{i}' for i in range(int((samples.shape[1] + 1) / 6))
                       for pn in
                       ['xy_amp', 'xy_mux', 'xy_muy', 'xy_scalex', 'xy_scaley', 'xy_corr']][1:]
  p_summary['low'] = np.quantile(samples, qm[0], axis=0)
  p_summary['high'] = np.quantile(samples, qp[0], axis=0)
  p_summary['median'] = np.median(samples, axis=0)
  p_summary.to_csv(f'{file_prefix}_p.txt',
                   index=False,
                   header=True, sep='\t')
  plots_from_pdf(x, y, pdf, title, output,
                 genex, geney, log1pColorScale, splitys,
                 Neff)


def main(args=None):
  """Main function of baredSC_2d
  """
  original_args = sys.argv[1:]
  args = parse_arguments('baredSC_2d').parse_args(args)
  args = args_check_baredSC(args)
  # Update args and check
  args = args_check(args)
  # Load data
  data = get_data_from_args(args, [args.geneXColName, args.geneYColName])

  # Check the directory of args.output is writtable:
  if os.path.dirname(args.output) != '' and not os.access(os.path.dirname(args.output), os.W_OK):
    raise OSError("The output is not writable,"
                  " check the directory exists.")
  if not os.path.exists(f'{args.output}.npz') or args.force:
    # Run mcmc:
    print("Run MCMC")
    start = time.time()
    results = gauss_mcmc(data, args.geneXColName, args.geneYColName,
                         args.nx, args.osampx, args.osampxpdf, args.xmin, args.xmax,
                         args.ny, args.osampy, args.osampypdf, args.ymin, args.ymax,
                         args.minScalex, args.minScaley,
                         args.scalePrior, args.scale, args.targetSum, args.nnorm,
                         args.nsampMCMC, args.nsampBurnMCMC,
                         args.nsplitBurnMCMC, args.T0BurnMCMC,
                         args.output, args.seed)
    print(f"MCMC + Save took {datetime.timedelta(seconds=int(time.time() - start))}.")
  else:
    results = extract_from_npz(f'{args.output}.npz')
    # return(mu, cov, ox, oy, oxpdf, oypdf, x, y,
    #        logprob_values, samples)
    x = results[6]
    dx = x[1] - x[0]
    xmin = x[0] - dx / 2
    xmax = x[-1] + dx / 2
    y = results[7]
    dy = y[1] - y[0]
    ymin = y[0] - dy / 2
    ymax = y[-1] + dy / 2
    try:
      assert (results[0].size + 1) // 6 == args.nnorm, \
          "nnorm value do not match what is in output."
      assert results[2].size // results[6].size == args.osampx, \
          "osampx value do not match what is in output."
      assert results[4].size // results[6].size == args.osampxpdf, \
          "osampxpdf value do not match what is in output."
      assert results[3].size // results[7].size == args.osampy, \
          "osampy value do not match what is in output."
      assert results[5].size // results[7].size == args.osampypdf, \
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
      assert (results[8].size == args.nsampMCMC + 1) or \
          (args.minNeff is not None), \
          "nsampMCMC value do not match what is in output."
    except AssertionError as e:
      raise ValueError(f"Ouput file already exists and {e}"
                       " Use --force to rerun MCMC.") from e
    if args.minNeff is not None:
      args.nsampMCMC = results[8].size - 1

  if args.figure is not None:
    # Check the directory of args.figure is writtable:
    if os.path.dirname(args.figure) != '' and not os.access(os.path.dirname(args.figure), os.W_OK):
      print("The output figure is not writable no figure/export will be done,"
            " check the directory exists.")
    else:
      # Make the plots:
      plot(*results[4:], args.title, args.figure, args.geneXColName,
           args.geneYColName, args.splity, args.removeFirstSamples,
           args.nsampInPlot, args.log1pColorScale,
           args.prettyBinsx, args.prettyBinsy)
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
      write_evidence(data, args.geneXColName, args.geneYColName,
                     *results[:8], args.minScalex, args.minScaley,
                     args.scalePrior, args.coviscale,
                     args.nis, logprob, args.logevidence, args.seed,
                     args.scale, args.targetSum)
  return None


if __name__ == "__main__":
    arguments = None
    if len(sys.argv) == 1:
        arguments = ["--help"]
    main(arguments)
