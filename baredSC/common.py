# Copyright 2021 Jean-Baptiste Delisle and Lucille Delisle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import poisson
from scipy.sparse import issparse
import anndata
from samsam import acf

mpl.use('agg')


def get_data_txt(input_file, metadata1_cn, metadata1_v,
                 metadata2_cn, metadata2_v,
                 metadata3_cn, metadata3_v,
                 genes):
  """
  Create a dataframe from input file
  keep only rows matching metadata values
  """
  data = pd.read_csv(input_file, sep="\t")
  filters = {metadata1_cn: metadata1_v,
             metadata2_cn: metadata2_v,
             metadata3_cn: metadata3_v}
  for col_name in filters:
    if col_name is None:
      continue
    data = data[data[col_name].astype(str).isin(filters[col_name].split(','))]
  print(f"You have {data.shape[0]} cells.")
  if data.shape[0] == 0:
    raise Exception("No cell with selected filters.")
  for col_name in genes + ['nCount_RNA']:
    if col_name not in data.columns:
      raise Exception(f"{col_name} not in the data table.")
  return(data[genes + ['nCount_RNA']])


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
    raise Exception(f"Negative values found {add_message}. Cannot go back to raw counts.")
  # Check if they are already raw_counts:
  if np.array_equal(X, X.astype(int)):
    return(X)
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
    return(np.round(X))
  # We first check that all values are there:
  # Try to find the target_sum:
  N_t = np.median(X.sum(1))
  if np.max(np.abs(X.sum(1) - N_t)) > 0.1:
    # Maybe it was logged but the lop1p in uns has not been exported
    X = np.expm1(X)
    add_message += " delogged"
    if np.max(np.abs(X - np.round(X))) < round_threshold:
      return(np.round(X))
    N_t = np.median(X.sum(1))
    if np.max(np.abs(X.sum(1) - N_t)) > 0.1:
      raise Exception("The sum of expression for each cell"
                      f" {add_message} is not constant."
                      " Genes are probably missing.")

  # We could also check the number of genes detected with:
  # np.array_equal(np.sum(X > 0, axis=1), adata.obs['n_genes_by_counts'])
  if 'total_counts' not in adata.obs.keys():
    raise Exception(f"The expression {add_message}seems normalized"
                    " but 'total_counts' is not in obs."
                    "Cannot get raw values.")
  # Denormalize:
  X = X * adata.obs['total_counts'][:, None] / N_t
  if np.max(np.abs(X - np.round(X))) < round_threshold:
    return(np.round(X))
  else:
    raise Exception(f"Could not extract raw values of expression {add_message}.")


def get_data_annData(input_file, metadata1_cn, metadata1_v,
                     metadata2_cn, metadata2_v,
                     metadata3_cn, metadata3_v,
                     genes):
  """
  Create a dataframe from input file
  keep only rows matching metadata values
  """
  # First load the annData:
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
  for col_name in filters:
    if col_name is None:
      continue
    data = data[data[col_name].astype(str).isin(filters[col_name].split(','))]
  print(f"You have {data.shape[0]} cells.")
  if data.shape[0] == 0:
    raise Exception("No cell with selected filters.")
  for col_name in genes + ['nCount_RNA']:
    if col_name not in data.columns:
      raise Exception(f"{col_name} not in the data table.")
  return(data[genes + ['nCount_RNA']])


# Permut p to get the permutation of p the closest to pref
def permuted(p, pref, wdist, permut, n_param_per_gauss):
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
  return(min_perm)


# Get p(ki|Ni,xj) * dx
def get_Ax(data, col_gene, nx, x, dx, ox, xscale, target_sum):
  # Get all UMI nbs:
  N = data['nCount_RNA']
  # Get the raw counts:
  kx = data[col_gene]
  # To get a better approximation we will use
  # a finer grid for x (ox)
  if ox.size % nx != 0:
    raise Exception("The size of ox is not multiple of size of x.")
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
    raise Exception("Only Seurat and log xscale are implemented.")
  # To evaluate the integral more easily multiply by the delta x
  Ax = poissx * dx

  return(Ax)


# Process the output to get prefix and suffix
def get_prefix_suffix(output):
  name = output.split(".")
  file_suffix = name[-1]
  file_prefix = ".".join(name[:-1])
  return(file_prefix, file_suffix)


# Get the number of effective samples
def get_Neff(samples):
  # Autocorrelation (to test the convergence of MCMC)
  R = np.array([acf.acf(pi) for pi in samples.T]).T
  # Estimate the number of independent samples:
  Neff = samples.shape[0] / np.max([acf.iat(R=r) for r in R.T])
  return(Neff)


# Do some plots and exports in txt
def plot_QC(logprob_values, samples, title, output,
            removeFirstSamples, nsampInPlot,
            p_names):
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
  with open(f'{file_prefix}_neff.txt', 'w') as fo:
    fo.write(f'{Neff}\n')

  # Plot the value of each parameter and logprob_values along the samples:
  fig, axs = plt.subplots(samples.shape[1] + 1, 1, sharex='col', figsize=(6, (samples.shape[1] + 1) * 1.5))
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
    import corner
    plt.figure()
    figure = corner.corner(samples,
                           labels=p_names,
                           quantiles=[0.15865, 0.5, 0.84135],
                           show_titles=True, title_kwargs={"fontsize": 12})
    figure.suptitle(title, ha='right', x=0.98)
    plt.savefig(f'{file_prefix}_corner.{file_suffix}', dpi=100)
  except Exception as e:
    print("corner plot could not be plotted.")
    print(e)

  return(samples, Neff)


def get_bins_centers(xmin, xmax, nx):
  x_borders = np.linspace(xmin, xmax, nx + 1)
  return((x_borders[:-1] + x_borders[1:]) / 2)
