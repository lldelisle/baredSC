# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
"""
  Module with all functions common to baredSC_2d combineMultipleModels_2d
"""
from itertools import permutations
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.special import factorial
import pandas as pd
import samsam
import samsam.logprior

# Local imports
from . common import get_Ax, get_prefix_suffix

mpl.use('agg')


def logprob(p,
            pref, wdist, permutxy,
            xmin, xmax, nx,
            ymin, ymax, ny,
            nxy,
            noxpdf, noypdf,
            oxypdf, odxypdf,
            min_scale_x, min_scale_y,
            scale_prior,
            Axy,
            temp=1.0):
  """Define log probability function

  Args:
    p (np.ndarray): Set of parameters to get logprob
    pref (np.ndarray): Set of parameters of reference
    wdist (np.ndarray): Weight to give to parameters in distance evaluation
    permutxy (np.ndarray): Array with all possible permutations on x and y
    xmin (int): Minimum value to consider in x axis.
    xmax (int): Maximum value to consider in x axis.
    nx (int): Number of values in x to check how your evaluated pdf is compatible with the model.
    ymin (int): Minimum value to consider in y axis.
    ymax (int): Maximum value to consider in y axis.
    ny (int): Number of values in y to check how your evaluated pdf is compatible with the model.
    nxy (int): Number of values for a PDF in 2D
    noxpdf (int): Number of values in x to evaluate pdf
    noypdf (int): Number of values in y to evaluate pdf
    oxypdf (np.ndarray): Matrix with coordinates of centers of 2D bins
    odxypdf (float): Area of 2D bin
    min_scale_x (float): Minimal value of the scale of gaussians in x
    min_scale_y (float): Minimal value of the scale of gaussians in y
    scale_prior (float): Scale of the truncnorm used in the prior for the correlation.
    Axy (np.ndarray): Matrix with p(ki,li|Ni,xj,yk)*dx*dy
    temp (float, optional): Temperature. Defaults to 1.0.
  """
  # Log prior
  lp = logprior(p, pref, wdist,
                permutxy,
                xmin, xmax, ymin, ymax,
                min_scale_x, min_scale_y,
                scale_prior)
  if np.isfinite(lp):
    # If the parameters are possible
    # Log likelihood
    ll = loglike(p,
                 nx, ny, nxy, noxpdf, noypdf,
                 oxypdf, odxypdf, Axy)
    # Return sum of log prior and log likelihood
    return (lp + ll) / temp
  # If the parameters are impossible
  return -np.inf


def logprior(p,
             pref, wdist, permutxy,
             xmin, xmax, ymin, ymax,
             min_scale_x, min_scale_y,
             scale_prior):
  """Define the log prior

  Args:
    p (np.ndarray): Set of parameters to get logprob
    pref (np.ndarray): Set of parameters of reference
    wdist (np.ndarray): Weight to give to parameters in distance evaluation
    permutxy (np.ndarray): Array with all possible permutations on x and y
    xmin (int): Minimum value to consider in x axis.
    xmax (int): Maximum value to consider in x axis.
    ymin (int): Minimum value to consider in y axis.
    ymax (int): Maximum value to consider in y axis.
    min_scale_x (float): Minimal value of the scale of gaussians in x
    min_scale_y (float): Minimal value of the scale of gaussians in y
    scale_prior (float): Scale of the truncnorm used in the prior for the correlation.

  Returns:
    float: Log prior value and -inf if it is not in prior
  """
  # Raise exception if one of the condition is not
  # satisfied
  amp = p[5::6]
  mux = p[::6]
  muy = p[1::6]
  sigx = p[2::6]
  sigy = p[3::6]
  cor = p[4::6]
  # We complete the list of parameters:
  p_full = np.concatenate(([1 - np.sum(amp)], p))
  amp_full = np.concatenate(([1 - np.sum(amp)], amp))
  try:
    # Amplitudes
    # Amplitudes should be positive
    # Including the first one
    if np.any(amp_full < 0):
      raise ValueError('Negative amplitude')
    # Because we impose that sum of amplitude is equal to 1
    # We add this term:
    lp = np.log(factorial(amp.size))
    # Check swapping (symmetry of the problem)
    # pref is a parameter set and we impose that pref
    # is closer to p than any permutation of p
    dp = p - pref
    dist = np.sum(wdist * dp * dp)
    for ks in permutxy[1:]:
      dpk = p_full[ks][1:] - pref
      distk = np.sum(wdist * dpk * dpk)
      if distk <= dist:
        raise ValueError('Swapping')
    # Number of possible swappings
    lp += np.log(permutxy.shape[0])
    # Finaly we add the priors on other parameters
    # Location
    lp += sum(samsam.logprior.uniform(muk, xmin - 3 * sigk, xmax + 3 * sigk)
              for muk, sigk in zip(mux, sigx))
    lp += sum(samsam.logprior.uniform(muk, ymin - 3 * sigk, ymax + 3 * sigk)
              for muk, sigk in zip(muy, sigy))
    # Correlation
    # lp += sum(samsam.logprior.uniform(cor, -0.95, 0.95) for cor in cor)
    lp += sum(samsam.logprior.truncnormal(cork, mu=0, sig=scale_prior, a=-0.95, b=0.95)
              for cork in cor)
    # Scale
    lp += sum(samsam.logprior.loguniform(sigk, np.log(min_scale_x), np.log(xmax - xmin))
              for sigk in sigx)
    lp += sum(samsam.logprior.loguniform(sigk, np.log(min_scale_y), np.log(ymax - ymin))
              for sigk in sigy)
  except Exception:  # pylint: disable=W0718
    return -np.inf
  return lp


def trunc_norm2d(xy, dxy, mu, sigma, corr):
  """pdf of a trunc_norm in 2d

  Args:
    xy (np.ndarray): Matrix with coordinates of centers of 2D bins (regularly spaced)
    dxy (float): Area of 2D bin
    mu (np.ndarray): Means in x and y of the Gaussian
    sigma (np.ndarray): Scales in x and y of the Gaussian
    corr (float): Correlation value of the Gaussian
  """
  u = (xy - mu) / sigma
  lognorm = -0.5 * (u[..., 0] * u[..., 0] + u[..., 1] * u[..., 1] -
                    2 * corr * u[..., 0] * u[..., 1]) / (1 - corr * corr)
  lognorm -= np.max(lognorm)  # Avoid underflow
  norm = np.exp(lognorm)
  norm /= np.sum(norm * dxy)
  return norm


def loglike(p,
            nx, ny, nxy, noxpdf, noypdf,
            oxypdf, odxypdf, Axy):
  """Log likelihood

  Args:
    p (np.ndarray): Set of parameters to get loglikelihood
    nx (int): Number of values in x to check how your evaluated pdf is compatible with the model.
    ny (int): Number of values in y to check how your evaluated pdf is compatible with the model.
    nxy (int): Number of values for a PDF in 2D
    noxpdf (int): Number of values in x to evaluate pdf
    noypdf (int): Number of values in y to evaluate pdf
    oxypdf (np.ndarray): Matrix with coordinates of centers of 2D bins
    odxypdf (float): Area of 2D bin
    Axy (np.ndarray): Matrix with p(ki,li|Ni,xj,yk)*dx*dy
    
  Returns:
    float: Log likelihood
  """
  # Get the pdf from p:
  pdf = get_pdf(p, nx, ny, noxpdf, noypdf,
                oxypdf, odxypdf)
  # The likelihood for a cell i is
  # integral(p(ki,li,Ni|x, y)p(x, y|amp,mu,sig)dxdy)
  # We approximate the integral by the sum
  # p(ki,li,Ni|xj, yj)*dx*dy is stored in Axy
  likecell = Axy @ pdf.reshape(nxy)
  # The log likelihood for all cells is
  # sum of log likelihood for each cell
  ll = np.sum(np.log(likecell))
  return ll


def get_pdf(p,
            nx, ny,
            noxpdf, noypdf,
            oxypdf, odxypdf):
  """Get the pdf on xy 2d using parameters

  Args:
    p (np.ndarray): Set of parameters to get pdf
    nx (int): Number of values in x to check how your evaluated pdf is compatible with the model.
    ny (int): Number of values in y to check how your evaluated pdf is compatible with the model.
    nxy (int): Number of values for a PDF in 2D
    noxpdf (int): Number of values in x to evaluate pdf
    noypdf (int): Number of values in y to evaluate pdf
    oxypdf (np.ndarray): Matrix with coordinates of centers of 2D bins
    odxypdf (float): Area of 2D bin
    
  Returns:
    np.ndarray: The pdf evaluated on xy
  """
  amp = p[5::6]
  mux = p[::6]
  muy = p[1::6]
  sigx = p[2::6]
  sigy = p[3::6]
  cor = p[4::6]
  # We complete the list of parameters:
  amp_full = np.concatenate(([1 - np.sum(amp)], amp))
  # Build the oversampled pdf
  pdf_o = np.zeros((noypdf, noxpdf))
  for ampk, muxk, muyk, sigxk, sigyk, cork in zip(amp_full, mux, muy, sigx, sigy, cor):
    # Add the trunc_norm 2d
    muk = np.array([muxk, muyk])
    sigk = np.array([sigxk, sigyk])
    pdf_o += ampk * trunc_norm2d(oxypdf, odxypdf, muk, sigk, cork)
  # I do the mean in each bin:
  pdf = np.mean(pdf_o.reshape(ny, noypdf // ny, nx, noxpdf // nx), axis=(1, 3))
  return pdf


# Return the same as gauss_mcmc but from a npz
def extract_from_npz(output):
  """Return the same as gauss_mcmc but from a npz

  Args:
    output (str): file_path to numpy archive
      
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
  print("Reading")
  start = time.time()
  mcmc_res = np.load(output, allow_pickle=True)
  diags = mcmc_res['diagnostics']
  mu = diags.item().get('mu')
  cov = diags.item().get('cov')
  logprob_values = diags.item().get('logprob')
  samples = mcmc_res['samples']
  x = mcmc_res['x']
  y = mcmc_res['y']
  ox = mcmc_res['ox']
  oy = mcmc_res['oy']
  oxpdf = mcmc_res['oxpdf']
  oypdf = mcmc_res['oypdf']
  print(f"Read. It took {(time.time() - start):.2f} seconds.")
  return(mu, cov, ox, oy, oxpdf, oypdf, x, y,
         logprob_values, samples)


def write_evidence(data, genex, geney, mu, cov, ox, oy, oxpdf, oypdf, x, y,
                   min_scale_x, min_scale_y, scale_prior,
                   covis_scale, nis, logprob_function, output,
                   seed, scale, target_sum):
  """Write evidence to file
  In order to compare multiple models we use an estimate of the log evidence

  Args:
    data (pandas.DataFrame): Data frame with 'nCount_RNA' and the genes of interest
    genex (str): Column label in `data` with the gene of interest on x
    geney (str): Column label in `data` with the gene of interest on y
    mu (np.ndarray): The average of parameters
    cov (np.ndarray): The covariance between parameters
    ox (np.ndarray): The nx * osampx bin centers
    oy (np.ndarray): The ny * osampy bin centers
    oxpdf (np.ndarray): The osampxpdf * nx bin centers
    oypdf (np.ndarray): The osampypdf * ny bin centers
    x (np.ndarray): The nx bin centers
    y (np.ndarray): The ny bin centers
    min_scale_x (float): Minimal value of the scale of gaussians in x
    min_scale_y (float): Minimal value of the scale of gaussians in y
    scale_prior (float): Scale of the truncnorm used in the prior for the correlation.
    covis_scale (float): Scale factor to apply to covariance of parameters
                         to get random parameters in logevidence evaluation.
    nis (int): Size of sampling of random parameters in logevidence evaluation.
    logprob_function (function): Log probability of the distribution to sample.
    output (str): File path to write value of log evidence
    seed (int): Seed
    scale (str): scale for the x-axis and y-axis: Seurat (log(1+targetSum*X)) or log (log(X))
    target_sum (int): factor when Seurat scale is used: (log(1+targetSum*X))

  Raises:
      Exception: _description_
      Exception: _description_
      Exception: _description_
      Exception: _description_
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

  ny = y.size
  dy = y[1] - y[0]
  ymin = y[0] - dy / 2
  ymax = y[-1] + dy / 2
  noypdf = oypdf.size
  odypdf = oypdf[1] - oypdf[0]
  if oy.size % ny != 0:
    raise ValueError("The size of oy is not multiple of size of y.")
  if oypdf.size % ny != 0:
    raise ValueError("The size of oypdf is not multiple of size of y.")

  nxy = nx * ny
  odxypdf = odxpdf * odypdf
  oxypdf = np.array(np.meshgrid(oxpdf, oypdf)).transpose(1, 2, 0)

  # Poisson law
  Ax = get_Ax(data, genex, nx, dx, ox, scale, target_sum)
  Ay = get_Ax(data, geney, ny, dy, oy, scale, target_sum)

  Axy = Ax[:, None, :] * Ay[:, :, None]
  Axy = Axy.reshape(-1, nxy)

  nnorm = (mu.size + 1) // 6

  # We want to get all possible permutation of the gaussians (6 parameters for each)
  permutxy = np.array([np.concatenate(([6 * i + np.arange(6) for i in v]))
                       for v in permutations(np.arange(nnorm))])

  # We calculate logevidence
  np.random.seed(seed)
  logevidence = samsam.covis(mu, covis_scale ** 2 * cov,
                             logprob_function, nsamples=nis,
                             pref=mu, wdist=1 / np.diag(cov),
                             permutxy=permutxy,
                             xmin=xmin, xmax=xmax, nx=nx,
                             ymax=ymax, ymin=ymin, ny=ny,
                             nxy=nxy, noxpdf=noxpdf, noypdf=noypdf,
                             oxypdf=oxypdf, odxypdf=odxypdf,
                             min_scale_x=min_scale_x,
                             min_scale_y=min_scale_y,
                             scale_prior=scale_prior,
                             Axy=Axy)[2]['logevidence']

  with open(output, 'w', encoding="utf-8") as fo:
    fo.write(f'{logevidence}\n')


def plots_from_pdf(x, y, pdf, title, output,
                   gene_x, gene_y,
                   log1pColorScale, splitys, Neff=None):
  """Compute plots for 2d

  Args:
    x (np.ndarray): The nx bin centers on which pdf has been evaluated
    y (np.ndarray): The ny bin centers on which pdf has been evaluated
    pdf (np.ndarray): PDF values on x, y for each sample to plot
    title (str): Title for plots
    output (str): Path for output plots
    gene_x (str): Column label in `data` with the gene of interest on x
    gene_y (str): Column label in `data` with the gene of interest on y
    log1pColorScale (bool): Use log1p color scale instead of linear color scale.
    splitys ([float]): Threshold value to plot the density for genex
                       for 2 categories in geney values.
    Neff (float, optional): Value of Neff. Defaults to None.
  """
  # We assume x is equally spaced
  dx = x[1] - x[0]
  xmin = x[0] - dx / 2
  xmax = x[-1] + dx / 2
  dy = y[1] - y[0]
  ymin = y[0] - dy / 2
  ymax = y[-1] + dy / 2

  # Process the output to get prefix and suffix
  file_prefix, file_suffix = get_prefix_suffix(output)

  # 1, 2, 3 sigma:
  qm = np.array([(1 - lvl) / 2 for lvl in [0.6827, 0.9545, 0.9973]])
  qp = 1 - qm

  # Evaluate pdf1d:
  pdfx = np.sum(pdf, axis=1) * dy
  pdfy = np.sum(pdf, axis=2) * dx

  # Evaluate the correlation:
  to_export = []
  print("Evaluate mean")
  mux = pdfx @ x * dx
  muy = pdfy @ y * dy
  print("Done.")
  cx = x[None, :] - mux[:, None]
  cy = y[None, :] - muy[:, None]
  varx = np.sum(cx ** 2 * pdfx, axis=1) * dx
  vary = np.sum(cy ** 2 * pdfy, axis=1) * dy
  print("Evaluate covariance")
  covxy = np.sum(cx[:, None, :] * cy[:, :, None] * pdf, axis=(1, 2)) * dx * dy
  print("Done.")
  corrxy = covxy / np.sqrt(varx * vary)
  mucorrxy = np.mean(corrxy)
  print('correlation:')
  print('mean =', mucorrxy)
  to_export.append(mucorrxy)
  print('med =', np.median(corrxy))
  to_export.append(np.median(corrxy))
  for k in range(3):
    print(f'{k+1}sig =', np.quantile(corrxy, qm[k]), np.quantile(corrxy, qp[k]))
    if k == 0:
      to_export.append(np.quantile(corrxy, qm[k]))
      to_export.append(np.quantile(corrxy, qp[k]))

  # Error on correlation
  sigcorrxym = mucorrxy - np.quantile(corrxy, qm[0])
  sigcorrxyp = np.quantile(corrxy, qp[0]) - mucorrxy

  # Estimate the p-value if Neff is not None:
  if Neff is not None:
    # Estimate p-value for mucorrxy < 0
    # Get the fraction
    frac_cor_pos = np.sum(corrxy >= 0) / corrxy.size
    # The expected value of the fraction is:
    pval_pos = frac_cor_pos + 1 / Neff
    # Get the error
    pval_pos_error = np.sqrt(pval_pos / Neff)
    # Estimate p-value for mucorrxy > 0
    # Get the fraction
    frac_cor_neg = np.sum(corrxy <= 0) / corrxy.size
    # The expected value of the fraction is:
    pval_neg = frac_cor_neg + 1 / Neff
    # Get the error
    pval_neg_error = np.sqrt(pval_neg / Neff)

  # Plot 4 panels plot
  fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')
  # Heatmap:
  ax = axs[1, 0]
  x_borders = np.linspace(xmin, xmax, len(x) + 1)
  y_borders = np.linspace(ymin, ymax, len(y) + 1)
  if log1pColorScale:
    ax.pcolormesh(x_borders, y_borders, np.log1p(np.mean(pdf, axis=0)),
                  shading='flat', rasterized=True, cmap='Greys')
  else:
    ax.pcolormesh(x_borders, y_borders, np.mean(pdf, axis=0),
                  shading='flat', rasterized=True, cmap='Greys')
  ax.set_xlabel(gene_x)
  ax.set_ylabel(gene_y)

  # PDF 1d
  ax = axs[0, 0]
  ax.fill_between(x, np.quantile(pdfx, qm[0], axis=0), np.quantile(pdfx, qp[0], axis=0),
                  color='r', alpha=0.2)
  ax.plot(x, np.mean(pdfx, axis=0), 'r')
  ax.plot(x, np.median(pdfx, axis=0), 'k--')
  ax.set_ylim(0,)
  ax.set_ylabel('PDF')

  ax = axs[1, 1]
  ax.fill_betweenx(y, np.quantile(pdfy, qm[0], axis=0), np.quantile(pdfy, qp[0], axis=0),
                   color='r', alpha=0.2)
  ax.plot(np.mean(pdfy, axis=0), y, 'r')
  ax.plot(np.median(pdfy, axis=0), y, 'k--')
  ax.set_xlim(0,)
  ax.set_xlabel('PDF')

  # Write correlation
  ax = axs[0, 1]
  ax.axis('off')
  ax.text(0.5, 0.5, f'Corr = ${mucorrxy:.3f}_{{-{sigcorrxym:.3f}}}^{{+{sigcorrxyp:.3f}}}$',
          horizontalalignment='center', verticalalignment='center',
          transform=ax.transAxes)
  if Neff is not None:
    if mucorrxy > 0:
      text_pval = f'p-value = ${pval_neg:.2}^+_- {pval_neg_error:.2}$'
      to_export.append(pval_neg)
      to_export.append(pval_neg_error)
    else:
      text_pval = f'p-value = ${pval_pos:.2}^+_- {pval_pos_error:.2}$'
      to_export.append(pval_pos)
      to_export.append(pval_pos_error)
    ax.text(0.5, 0.25, text_pval,
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
  fig.suptitle(title)

  # Save
  plt.savefig(output)

  # Use median instead of mean:
  ax = axs[1, 0]
  if log1pColorScale:
    ax.pcolormesh(x_borders, y_borders, np.log1p(np.median(pdf, axis=0)),
                  shading='flat', rasterized=True, cmap='Greys')
  else:
    ax.pcolormesh(x_borders, y_borders, np.median(pdf, axis=0),
                  shading='flat', rasterized=True, cmap='Greys')

  # Save
  plt.savefig(f'{file_prefix}_median.{file_suffix}')

  # Plot the projection of individual pdfs:
  cmap = plt.cm.get_cmap('jet', 100)

  # Plot 4 panels plot
  fig, axs = plt.subplots(2, 2)
  fig.suptitle(f"{title}\n100 individual samples")
  axs[0, 1].axis('off')
  axs[1, 0].axis('off')
  # PDF 1d
  ax = axs[0, 0]
  for i, cur_pdf in enumerate(pdfx[::(pdfx.shape[0] // 100)][:100]):
    ax.plot(x, cur_pdf, color=cmap(i), alpha=0.3)
  ax.set_ylim(0,)
  ax.set_ylabel('PDF')
  ax.set_xlabel(gene_x)

  ax = axs[1, 1]
  for i, cur_pdf in enumerate(pdfy[::(pdfy.shape[0] // 100)][:100]):
    ax.plot(cur_pdf, y, color=cmap(i), alpha=0.3)
  ax.set_xlim(0,)
  ax.set_xlabel('PDF')
  ax.set_ylabel(gene_y)

  # Save
  plt.savefig(f'{file_prefix}_individuals.{file_suffix}')

  # Write info of correlation to file:
  with open(f'{file_prefix}_corr.txt', 'w', encoding="utf-8") as f:
    if Neff is None:
      f.write('mean\tmedian\tlow\thigh\n')
    else:
      f.write('mean\tmedian\tlow\thigh\tpval\terror\n')
    f.write('\t'.join([str(v) for v in to_export]) + '\n')

  # Export the mean pdf with colnames and rownames
  np.savetxt(f'{file_prefix}_pdf2d.txt',
             np.concatenate((np.array([y]).T, np.mean(pdf, axis=0)), axis=1),
             delimiter="\t", header="xy\t" + "\t".join([str(v) for v in x]),
             comments="")

  # Export the pdf mean and 1 sigma (68%) and median
  pdf_summary = pd.DataFrame(columns=['x', 'y', 'low', 'mean', 'high', 'median'])
  pdf_summary['x'] = np.tile(x, len(y))
  pdf_summary['y'] = np.repeat(y, len(x))
  pdf_summary['low'] = np.quantile(pdf, qm[0], axis=0).flatten()
  pdf_summary['mean'] = np.mean(pdf, axis=0).flatten()
  pdf_summary['high'] = np.quantile(pdf, qp[0], axis=0).flatten()
  pdf_summary['median'] = np.median(pdf, axis=0).flatten()
  # Export
  pdf_summary.to_csv(f'{file_prefix}_pdf2d_flat.txt',
                     index=False,
                     header=True, sep='\t')

  # Sum up pdf above and below threshold.
  if splitys is not None:
    for splity in splitys:
      # splity renorm
      pdfx1 = np.sum(pdf[:, y < splity], axis=1)
      pdfx1 /= np.sum(pdfx1, axis=1)[:, None] * dx
      pdfx2 = np.sum(pdf[:, y >= splity], axis=1)
      pdfx2 /= np.sum(pdfx2, axis=1)[:, None] * dx
      plt.figure()
      plt.fill_between(x, np.quantile(pdfx1, qm[0], axis=0), np.quantile(pdfx1, qp[0], axis=0),
                       color='r', alpha=0.2)
      plt.plot(x, np.mean(pdfx1, axis=0), 'r', label=f'{gene_y}$<${splity}')
      plt.plot(x, np.median(pdfx1, axis=0), 'k--')
      plt.fill_between(x, np.quantile(pdfx2, qm[0], axis=0), np.quantile(pdfx2, qp[0], axis=0),
                       color='g', alpha=0.2)
      plt.plot(x, np.mean(pdfx2, axis=0), 'g', label=f'{gene_y}$\\geq${splity}')
      plt.plot(x, np.median(pdfx2, axis=0), 'k--', label='median')
      plt.xlim(xmin, xmax)
      plt.ylim(0, 3)
      plt.xlabel(gene_x)
      plt.ylabel('PDF')
      plt.legend()
      plt.title(title)
      plt.savefig(f'{file_prefix}_split{splity}_renorm.{file_suffix}')

      # splity raw
      pdfx1 = np.sum(pdf[:, y < splity], axis=1) * dy
      pdfx2 = np.sum(pdf[:, y >= splity], axis=1) * dy
      plt.figure()
      plt.fill_between(x, np.quantile(pdfx, qm[0], axis=0), np.quantile(pdfx, qp[0], axis=0),
                       color='k', alpha=0.2)
      plt.plot(x, np.mean(pdfx, axis=0), 'k', label='all')
      plt.fill_between(x, np.quantile(pdfx1, qm[0], axis=0), np.quantile(pdfx1, qp[0], axis=0),
                       color='r', alpha=0.2)
      plt.plot(x, np.mean(pdfx1, axis=0), 'r', label=f'{gene_y}$<${splity}')
      plt.fill_between(x, np.quantile(pdfx2, qm[0], axis=0), np.quantile(pdfx2, qp[0], axis=0),
                       color='g', alpha=0.2)
      plt.plot(x, np.mean(pdfx2, axis=0), 'g', label=f'{gene_y}$\\geq${splity}')
      plt.xlim(xmin, xmax)
      plt.ylim(0, 3)
      plt.xlabel(gene_x)
      plt.ylabel('PDF')
      plt.legend()
      plt.title(title)
      plt.savefig(f'{file_prefix}_split{splity}.{file_suffix}')

      np.savetxt(f'{file_prefix}_split{splity}.txt',
                 np.array([x, np.quantile(pdfx1, qm[0], axis=0),
                           np.mean(pdfx1, axis=0),
                           np.quantile(pdfx1, qp[0], axis=0),
                           np.quantile(pdfx2, qm[0], axis=0),
                           np.mean(pdfx2, axis=0),
                           np.quantile(pdfx2, qp[0], axis=0)]).T,
                 delimiter="\t",
                 header="x\tlow1\tmean1\thigh1\tlow2\tmean2\thigh2", comments="")

def args_check(args):
  """Check args for 2d

  Args:
      args (Namespace): Namespace object to check

  Raises:
      ValueError: when arguments are not valid

  Returns:
      Namespace: Updated Namespace
  """

  # Check incompatibilities:
  if args.scale == 'Seurat' and args.targetSum < 0:
    raise ValueError("--targetSum negative is not "
                     "compatible with --xscale Seurat "
                     "as it is log(1+targetSum*X).")
  if args.scale == 'Seurat' and args.xmin < 0:
    raise ValueError("--xmin negative is not "
                     "compatible with --scale Seurat "
                     "as it is log(1+targetSum*X).")
  if args.scale == 'log' and args.xmax > 0:
    raise ValueError("--xmax positive is not "
                     "compatible with --scale log "
                     "as it is log(X) and X < 1.")
  if args.scale == 'Seurat' and args.ymin < 0:
    raise ValueError("--ymin negative is not "
                     "compatible with --scale Seurat "
                     "as it is log(1+targetSum*X).")
  if args.scale == 'log' and args.ymax > 0:
    raise ValueError("--ymax positive is not "
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
  return args
