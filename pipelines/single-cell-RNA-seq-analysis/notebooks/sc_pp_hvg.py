# BSD 3-Clause License

# Copyright (c) 2022 Intel
# Copyright (c) 2017 F. Alexander Wolf, P. Angerer, Theis Lab
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# This code is based on scanpy/preprocessing/_highly_variable_genes.py
# Original Code Modified by Author Padmanabhan S. Pillai <padmanabhan.s.pillai@intel.com>

import warnings
from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from anndata import AnnData


#from .. import logging as logg
#from .._settings import settings, Verbosity
#from .._utils import sanitize_anndata, check_nonnegative_integers
#from .._compat import Literal
from scanpy.pp._utils import _get_mean_var
#from ._distributed import materialize_as_ndarray
#from ._simple import filter_genes
import time
import numba


#def _highly_variable_genes_single_batch(
def _get_hvg(
    adata: AnnData,
    nrows,          # real number of rows after initial filtering
    ncols,          # ignore column indices ncols and above
    colmask,        # prior column mask, if any
    n_top_genes: Optional[int] = None,
    n_bins: int = 20
) -> pd.DataFrame:
    """\
    See `highly_variable_genes`.

    Returns
    -------
    A DataFrame that contains the columns
    `highly_variable`, `means`, `dispersions`, and `dispersions_norm`.
    """

    @numba.njit(cache=True,parallel=True)
    def get_mean_var_disp(shp, indptr, indices, data, n, m, nthr):
        nr = len(indptr)-1
        s=np.zeros((nthr,m))
        ss=np.zeros((nthr,m))
        mean=np.zeros(m)
        var=np.zeros(m)
        for i in numba.prange(nthr):
            for r in range(i,nr,nthr):
                for j in range(indptr[r],indptr[r+1]):
                    c = indices[j]
                    if c>=m: continue
                    v = data[j]
                    s[i,c]+=v
                    ss[i,c]+=v*v
        for c in numba.prange(m):
            s0 = s[:,c].sum()
            mean[c] = s0/n
            var[c] = (ss[:,c].sum() - s0*s0/n)/(n-1)
            if mean[c]==0: mean[c] = 1e-12  # set entries equal to zero to small value
        # now actually compute the dispersion
        dispersion = var / mean
        return mean, var, dispersion


    t0 = time.time()
    X = adata.X

    #mean, var = materialize_as_ndarray(_get_mean_var(X))
    #mean, var = _get_mean_var(X)

    # now actually compute the dispersion
    #mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    #dispersion = var / mean
    nthr = numba.get_num_threads()
    mean, var, dispersion = get_mean_var_disp(adata.X.shape, adata.X.indptr, adata.X.indices, adata.X.data, nrows, ncols, nthr)
    
    print("got mean,var",time.time()-t0)

    # all of the following quantities are "per-gene" here
    df = pd.DataFrame()
    df['means'] = mean
    df['dispersions'] = dispersion

    from statsmodels import robust

    df['mean_bin'] = pd.cut(
        df['means'],
        np.r_[-np.inf, np.percentile(df['means'], np.arange(10, 105, 5)), np.inf],
    )
    disp_grouped = df.groupby('mean_bin')['dispersions']
    disp_median_bin = disp_grouped.median()
    # the next line raises the warning: "Mean of empty slice"
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        disp_mad_bin = disp_grouped.apply(robust.mad)
        df['dispersions_norm'] = (
            df['dispersions'].values - disp_median_bin[df['mean_bin'].values].values
        ) / disp_mad_bin[df['mean_bin'].values].values
        
    dispersion_norm = df['dispersions_norm'].values
    dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
    dispersion_norm[
        ::-1
    ].sort()  # interestingly, np.argpartition is slightly slower
    if n_top_genes > adata.n_vars:
        logg.info(f'`n_top_genes` > `adata.n_var`, returning all genes.')
        n_top_genes = adata.n_vars
    disp_cut_off = dispersion_norm[n_top_genes - 1]
    gene_subset = np.nan_to_num(df['dispersions_norm'].values) >= disp_cut_off
    #logg.debug(
    #    f'the {n_top_genes} top genes correspond to a '
    #    f'normalized dispersion cutoff of {disp_cut_off}'
    #)

    #df['highly_variable'] = gene_subset

    if ncols<adata.shape[1]:
        adata.var['highly_variable'] = np.zeros(adata.shape[1], dtype=bool)
        adata.var['means'] = np.zeros(adata.shape[1], dtype=np.float32)
        adata.var['dispersions'] = np.zeros(adata.shape[1], dtype=np.float32)
        adata.var['dispersions_norm'] = np.zeros(adata.shape[1], dtype=np.float32)
        adata.var['highly_variable'][colmask] = gene_subset
        adata.var['means'][colmask] = df['means'].values
        adata.var['dispersions'][colmask] = df['dispersions'].values
        adata.var['dispersions_norm'][colmask] = df['dispersions_norm'].values.astype('float32',copy=False)
    else:
        adata.var['highly_variable'] = gene_subset
        adata.var['means'] = df['means'].values
        adata.var['dispersions'] = df['dispersions'].values
        adata.var['dispersions_norm'] = df['dispersions_norm'].values.astype('float32',copy=False)


