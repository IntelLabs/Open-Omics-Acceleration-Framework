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


# THIS FILE IS BASED ON SCANPY SOURCE CODE
# Original Code Modified by Author Padmanabhan S. Pillai <padmanabhan.s.pillai@intel.com>

"""Simple Preprocessing Functions

Compositions of these functions are found in sc.preprocess.recipes.
"""
from functools import singledispatch
from numbers import Number
import warnings
from typing import Union, Optional, Tuple, Collection, Sequence, Iterable

import numba
import numpy as np
import scipy as sp
from scipy.sparse import issparse, isspmatrix_csr, csr_matrix, spmatrix
from sklearn.utils import sparsefuncs, check_array
from pandas.api.types import is_categorical_dtype
from anndata import AnnData
import os

regr_type = int(os.environ.get('FASTPP_REGRESS_STYLE', "2"))
assert(regr_type in [1,2])

#from .. import logging as logg
from scanpy import logging as logg
#from .._settings import settings as sett
from scanpy._settings import settings as sett
#from .._utils import sanitize_anndata, deprecated_arg_names, view_to_actual, AnyRandom, _check_array_function_arguments
from scanpy._utils import sanitize_anndata, deprecated_arg_names, view_to_actual, AnyRandom, _check_array_function_arguments
#from .._compat import Literal
#from ..get import _get_obs_rep, _set_obs_rep
#from ._distributed import materialize_as_ndarray
#from ._utils import _get_mean_var

import time
import psutil

def regress_out(
    adata: AnnData,
    keys: Union[str, Sequence[str]],
    n_jobs: Optional[int] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Regress out (mostly) unwanted sources of variation.

    Uses simple linear regression. This is inspired by Seurat's `regressOut`
    function in R [Satija15]. Note that this function tends to overcorrect
    in certain circumstances as described in :issue:`526`.

    Parameters
    ----------
    adata
        The annotated data matrix.
    keys
        Keys for observation annotation on which to regress on.
    n_jobs
        Number of jobs for parallel computation.
        `None` means using :attr:`scanpy._settings.ScanpyConfig.n_jobs`.
    copy
        Determines whether a copy of `adata` is returned.

    Returns
    -------
    Depending on `copy` returns or updates `adata` with the corrected data matrix.
    """
    start_time = time.time()
    print(psutil.Process().memory_info().rss // (1024*1024), "MB")
    start = logg.info(f'regressing out {keys}')
    if issparse(adata.X):
        logg.info(
            '    sparse input is densified and may '
            'lead to high memory use'
        )
    adata = adata.copy() if copy else adata

    sanitize_anndata(adata)

    # TODO: This should throw an implicit modification warning
    if adata.is_view:
        adata._init_as_actual(adata.copy())

    if isinstance(keys, str):
        keys = [keys]

    if issparse(adata.X):
        adata.X = adata.X.toarray()

    n_jobs = sett.n_jobs if n_jobs is None else n_jobs

    # regress on a single categorical variable
    variable_is_categorical = False
    if keys[0] in adata.obs_keys() and is_categorical_dtype(adata.obs[keys[0]]):
        if len(keys) > 1:
            raise ValueError(
                'If providing categorical variable, '
                'only a single one is allowed. For this one '
                'we regress on the mean for each category.'
            )
        logg.debug('... regressing on per-gene means within categories')
        regressors = np.zeros(adata.X.shape, dtype='float32')
        for category in adata.obs[keys[0]].cat.categories:
            mask = (category == adata.obs[keys[0]]).values
            for ix, x in enumerate(adata.X.T):
                regressors[mask, ix] = x[mask].mean()
        variable_is_categorical = True
    # regress on one or several ordinal variables
    else:
        # create data frame with selected keys (if given)
        if keys:
            regressors = adata.obs[keys]
        else:
            regressors = adata.obs.copy()

        # add column of ones at index 0 (first column)
        regressors.insert(0, 'ones', 1.0)

    len_chunk = np.ceil(min(1000, adata.X.shape[1]) / n_jobs).astype(int)
    n_chunks = np.ceil(adata.X.shape[1] / len_chunk).astype(int)

    print("REGRESSORS",regressors)

    tasks = []
    # split the adata.X matrix by columns in chunks of size n_chunk
    # (the last chunk could be of smaller size than the others)
    chunk_list = np.array_split(adata.X, n_chunks, axis=1)
    #chunk_list = [ adata.X[:,i*adata.X.shape[1]//n_chunks:(i+1)*adata.X.shape[1]//n_chunks].copy() for i in range(n_chunks)]
    print(psutil.Process().memory_info().rss // (1024*1024), "MB")
    if variable_is_categorical:
        regressors_chunk = np.array_split(regressors, n_chunks, axis=1)
    for idx, data_chunk in enumerate(chunk_list):
        # each task is a tuple of a data_chunk eg. (adata.X[:,0:100]) and
        # the regressors. This data will be passed to each of the jobs.
        if variable_is_categorical:
            regres = regressors_chunk[idx]
        else:
            regres = regressors
        tasks.append(tuple((data_chunk, regres, variable_is_categorical)))

    del(chunk_list)

    print(psutil.Process().memory_info().rss // (1024*1024), "MB")
    print ("prep done at:", time.time()-start_time)
    if n_jobs > 1 and n_chunks > 1:
        import multiprocessing
        pool = multiprocessing.Pool(n_jobs)
        res = pool.map_async(_regress_out_chunk, tasks).get(9999999)
        pool.close()

    else:
        res = list(map(_regress_out_chunk, tasks))

    print ("Main work completed at:", time.time()-start_time)
    # res is a list of vectors (each corresponding to a regressed gene column).
    # The transpose is needed to get the matrix in the shape needed
    print(psutil.Process().memory_info().rss // (1024*1024), "MB")
    #adata.X = np.vstack(res).T.astype(adata.X.dtype)
    i = 0
    for r in res:
        adata.X[:,i:i+r.shape[0]] = r.T
        i+=r.shape[0]
    print(psutil.Process().memory_info().rss // (1024*1024), "MB")
    logg.info('    finished', time=start)
    print ("Data collection completed at:", time.time()-start_time)
    return adata if copy else None


def _regress_out_chunk(data):
    # data is a tuple containing the selected columns from adata.X
    # and the regressors dataFrame
    data_chunk = data[0]
    regressors = data[1]
    variable_is_categorical = data[2]

    responses_chunk_list = []
    import statsmodels.api as sm
    from statsmodels.tools.sm_exceptions import PerfectSeparationError

    #output = np.zeros((data_chunk.shape[1],data_chunk.shape[0]))
    for col_index in range(data_chunk.shape[1]):

        # if all values are identical, the statsmodel.api.GLM throws an error;
        # but then no regression is necessary anyways...
        if not (data_chunk[:, col_index] != data_chunk[0, col_index]).any():
            responses_chunk_list.append(data_chunk[:, col_index])
            continue

        if variable_is_categorical:
            regres = np.c_[np.ones(regressors.shape[0]), regressors[:, col_index]]
        else:
            regres = regressors
        try:
            if regr_type==1:
                result = sm.GLM(data_chunk[:, col_index], regres, family=sm.families.Gaussian()).fit(maxiter=1,tol=1e-1)
            elif regr_type==2:
                result = sm.GLM(data_chunk[:, col_index], regres, family=sm.families.Gaussian()).fit(method='lbfgs',m=10,factr=1e24,maxfun=1)
            #result = sm.GLM(data_chunk[:, col_index], regres, family=sm.families.Gaussian()).fit(method='newton', tol=1e-1)
            new_column = result.resid_response
            #temp = np.zeros((data_chunk.shape[0],int(data_chunk.shape[0]**0.5))).T
            #temp += data_chunk[:,col_index]
            #new_column = temp[-1].copy()
            #del(temp)
            #new_column = result.resid_response.copy()
        except PerfectSeparationError:  # this emulates R's behavior
            logg.warning('Encountered PerfectSeparationError, setting to 0 as in R.')
            new_column = np.zeros(data_chunk.shape[0])

        responses_chunk_list.append(new_column)
        #output[col_index] = new_column
        #del(new_column)
        #del(result)

    return np.vstack(responses_chunk_list)
    #time.sleep(5)
    #return output

def dumb_regress_out(adata):
    
    @numba.njit(cache=True, parallel=True)
    def doit(X):
        for i in numba.prange(X.shape[1]):
            mean = np.mean(X[:,i])
            X[:,i] -= mean

    if issparse(adata.X):
        adata.X = adata.X.toarray()
    doit(adata.X)


@singledispatch
def scale(
    X: Union[AnnData, spmatrix, np.ndarray],
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
):
    """\
    Scale data to unit variance and zero mean.

    .. note::
        Variables (genes) that do not display any variation (are constant across
        all observations) are retained and (for zero_center==True) set to 0
        during this operation. In the future, they might be set to NaNs.

    Parameters
    ----------
    X
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    zero_center
        If `False`, omit zero-centering variables, which allows to handle sparse
        input efficiently.
    max_value
        Clip (truncate) to this value after scaling. If `None`, do not clip.
    copy
        Whether this function should be performed inplace. If an AnnData object
        is passed, this also determines if a copy is returned.
    layer
        If provided, which element of layers to scale.
    obsm
        If provided, which element of obsm to scale.

    Returns
    -------
    Depending on `copy` returns or updates `adata` with a scaled `adata.X`,
    annotated with `'mean'` and `'std'` in `adata.var`.
    """
    _check_array_function_arguments(layer=layer, obsm=obsm)
    return scale_array(data, zero_center=zero_center, max_value=max_value, copy=copy)


@scale.register(np.ndarray)
def scale_array(
    X,
    *,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
    return_mean_std: bool = False,
):
    if copy:
        X = X.copy()
    if not zero_center and max_value is not None:
        logg.info(  # Be careful of what? This should be more specific
            "... be careful when using `max_value` " "without `zero_center`."
        )

    if np.issubdtype(X.dtype, np.integer):
        logg.info(
            '... as scaling leads to float results, integer '
            'input is cast to float, returning copy.'
        )
        X = X.astype(float)

    mean, var = _get_mean_var(X)
    std = np.sqrt(var)
    std[std == 0] = 1
    if issparse(X):
        if zero_center:
            raise ValueError("Cannot zero-center sparse matrix.")
        sparsefuncs.inplace_column_scale(X, 1 / std)
    else:
        if zero_center:
            X -= mean
        X /= std

    # do the clipping
    if max_value is not None:
        logg.debug(f"... clipping at max_value {max_value}")
        X[X > max_value] = max_value

    if return_mean_std:
        return X, mean, std
    else:
        return X


@scale.register(spmatrix)
def scale_sparse(
    X,
    *,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
    return_mean_std: bool = False,
):
    # need to add the following here to make inplace logic work
    if zero_center:
        logg.info(
            "... as `zero_center=True`, sparse input is "
            "densified and may lead to large memory consumption"
        )
        X = X.toarray()
        copy = False  # Since the data has been copied
    return scale_array(
        X,
        zero_center=zero_center,
        copy=copy,
        max_value=max_value,
        return_mean_std=return_mean_std,
    )

@scale.register(AnnData)
def scale_anndata(
    adata: AnnData,
    *,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
) -> Optional[AnnData]:
    adata = adata.copy() if copy else adata
    view_to_actual(adata)
    X = _get_obs_rep(adata, layer=layer, obsm=obsm)
    X, adata.var["mean"], adata.var["std"] = scale(
        X,
        zero_center=zero_center,
        max_value=max_value,
        copy=False,  # because a copy has already been made, if it were to be made
        return_mean_std=True,
    )
    _set_obs_rep(adata, X, layer=layer, obsm=obsm)
    if copy:
        return adata


