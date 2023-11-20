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



#  Fast preprcoessing functions to use in place of scanpy.pp 
#  Original Code Modified by Author Padmanabhan S. Pillai <padmanabhan.s.pillai@intel.com>

# NOTE:  these are "fragile" in the sense they only work with sparse arrays i nthe SCR format;  there is no error checking though
# To quickly switch between fst and original functions, set the USE_FASTPP environment variable:
#        0 for scapy, non-zero for fast versions (default: fast)

import time
import numpy as np
import numba
import math
import os
import scanpy as sc
import anndata
import resource
import cProfile
import h5py
import pandas as pd
import scipy
import multiprocessing as mp
import threading

import mkl
import scanpy as sc
#import sc_pp_simple as myscpp
import sc_pp_hvg

use_fastpp = int(os.environ.get('USE_FASTPP', "1"))
if use_fastpp != 0:
    use_fastpp = True
else:
    use_fastpp = False
debug = int(os.environ.get('FASTPP_DEBUG', "0"))
strategy = int(os.environ.get('FASTPP_ADATA_COPY_STYLE', "11"))
assert(0<strategy and strategy<=11)

regr_type = int(os.environ.get('FASTPP_REGRESS_STYLE', "4"))
assert(regr_type in [0,1,2,3,4,5])

use_fastload = int(os.environ.get('FASTPP_FASTLOAD_STYLE', "5"))
assert(use_fastload in [0,1,2,3,4,5])

#  This seems to fail on some machines (perhaps those with newer software versions);  stick to int64
#if strategy<11:       # work around bug / expectation in scipy sparse array fancy indexing that assunes indptr type is same as indices type
#    indices_type = np.int64
#    indices_shm_type = "l"
#else:
#    indices_type = np.int32
#    indices_shm_type = "i"
indices_type = np.int64
indices_shm_type = "l"

def _load_helper(fname, s, e, dataArray, indicesArray):
    f = h5py.File(fname,'r')
    data = np.frombuffer(dataArray,dtype=np.float32)
    indices = np.frombuffer(indicesArray,dtype=indices_type)
    f['X']['data'].read_direct(data, np.s_[s:e], np.s_[s:e])
    f['X']['indices'].read_direct(indices, np.s_[s:e], np.s_[s:e])


def _load_helper2(fname, i, k, datalen, barrier, dataArray, indicesArray, startsArray, endsArray):
    f = h5py.File(fname,'r')
    dataA = np.frombuffer(dataArray,dtype=np.float32)
    indicesA = np.frombuffer(indicesArray,dtype=indices_type)
    startsA = np.frombuffer(startsArray,dtype=np.int64)
    endsA = np.frombuffer(endsArray,dtype=np.int64)
    for j in range(datalen//(k*1024*1024)+1):
        # compute start, end
        s = i*datalen//k + j*1024*1024
        e = min(s+1024*1024, (i+1)*datalen//k)
        length = e-s
        startsA[i]=s
        endsA[i]=e
        # read direct
        f['X']['data'].read_direct(dataA, np.s_[s:e], np.s_[i*1024*1024:i*1024*1024+length])
        f['X']['indices'].read_direct(indicesA, np.s_[s:e], np.s_[i*1024*1024:i*1024*1024+length])
        
        # coordinate workers
        barrier.wait()  # done data load
        barrier.wait()  # done data copy

semDataLoaded = None  # will be initialized later
semDataCopied = None  # will be initialized later

def _load_helper3(fname, i, k, datalen, dataArray, indicesArray, startsArray, endsArray):
    f = h5py.File(fname,'r')
    dataA = np.frombuffer(dataArray,dtype=np.float32)
    indicesA = np.frombuffer(indicesArray,dtype=indices_type)
    startsA = np.frombuffer(startsArray,dtype=np.int64)
    endsA = np.frombuffer(endsArray,dtype=np.int64)
    for j in range(datalen//(k*1024*1024)+1):
        # compute start, end
        s = i*datalen//k + j*1024*1024
        e = min(s+1024*1024, (i+1)*datalen//k)
        length = e-s
        startsA[i]=s
        endsA[i]=e
        # read direct
        f['X']['data'].read_direct(dataA, np.s_[s:e], np.s_[i*1024*1024:i*1024*1024+length])
        f['X']['indices'].read_direct(indicesA, np.s_[s:e], np.s_[i*1024*1024:i*1024*1024+length])
        
        # coordinate with copy threads
        semDataLoaded[i].release()  # done data load
        semDataCopied[i].acquire()  # wait until  data copied


@numba.njit(cache=True,parallel=True)
def _fast_copy(data,src,starts,ends,k):
    for i in numba.prange(k):
        length = ends[i]-starts[i]
        data[starts[i]:ends[i]] = src[i*1024*1024:i*1024*1024+length]

@numba.njit(cache=True,parallel=True,nogil=True)
def _fast_copy2(data,dataA,indices,indicesA,starts,ends,i):
    length = ends[i]-starts[i]
    data[starts[i]:ends[i]] = dataA[i*1024*1024:i*1024*1024+length]
    indices[starts[i]:ends[i]] = indicesA[i*1024*1024:i*1024*1024+length]

def _waitload(i):
    semDataLoaded[i].acquire()

def _signalcopy(i):
    semDataCopied[i].release()

@numba.njit(parallel=True)
def _fast_copy3(data,dataA,indices,indicesA,starts,ends,k,m):
    for i in numba.prange(k):
        for _ in range(m):
            with numba.objmode():
                _waitload(i)
            length = ends[i]-starts[i]
            data[starts[i]:ends[i]] = dataA[i*1024*1024:i*1024*1024+length]
            indices[starts[i]:ends[i]] = indicesA[i*1024*1024:i*1024*1024+length]
            with numba.objmode():
                _signalcopy(i)

if use_fastload==5:  # compile _fast_copy3
    d = np.zeros(10,dtype=np.float32)
    i = np.zeros(10,dtype=indices_type)
    s = np.zeros(10,dtype=np.int64)
    _fast_copy3(d,d,i,i,s,s,1,0)

def _copy_thread(data,dataA,indices,indicesA,starts,ends,i,m):
    for _ in range(m):
        semDataLoaded[i].acquire()  # wait for data load
        _fast_copy2(data,dataA,indices,indicesA,starts,ends,i)
        semDataCopied[i].release()  # signal data copied

def fastload(fname, firstn):
    t0 = time.time()
    f = h5py.File(fname,'r')
    assert ('X' in f.keys() and 'var' in f.keys() and 'obs' in f.keys())

    # get obs dataframe
    rows = f['obs'][ list(f['obs'].keys())[0] ].size
    print ("File has",rows,"rows")
    if firstn>1 and firstn<rows:
        rows = firstn
        print ("Loading",rows,"rows")
    if '_index' in f['obs'].keys():
        dfobsind = pd.Series(f['obs']['_index'].asstr()[0:rows])
        dfobs = pd.DataFrame(index=dfobsind)
    else:
        dfobs = pd.DataFrame()
    for k in f['obs'].keys():
        if k=='_index': continue
        dfobs[k] = f['obs'][k].asstr()[...]

    # get var dataframe
    if '_index' in f['var'].keys():
        dfvarind = pd.Series(f['var']['_index'].asstr()[...])
        dfvar = pd.DataFrame(index=dfvarind)
    else:
        dfvar = pd.DataFrame()
    for k in f['var'].keys():
        if k=='_index': continue
        dfvar[k] = f['var'][k].asstr()[...]
    print("Pandas load done",time.time()-t0)

    # load index pointers, prepare shared arrays
    indptr = f['X']['indptr'][0:rows+1]
    datalen = int(indptr[-1])

    if use_fastload==1:
        data = f['X']['data'][0:datalen]
        indices = f['X']['indices'][0:datalen]

    if use_fastload==2:
        f.close()
        dataArray = mp.Array('f',datalen,lock=False)     # should be in shared memory
        indicesArray = mp.Array(indices_shm_type,datalen,lock=False)  # should be in shared memory
        print("created shared memory buffers",time.time()-t0)

        #k = mkl.get_max_threads()
        k = numba.get_num_threads()
        procs = [mp.Process(target=_load_helper, args=(fname, i*datalen//k, (i+1)*datalen//k, dataArray, indicesArray))  for i in range(k)]
        for p in procs: p.start()
        for p in procs: p.join()
        print("multiproc load done",time.time()-t0)

        data = np.frombuffer(dataArray,dtype=np.float32)
        indices = np.frombuffer(indicesArray,dtype=indices_type)

    if use_fastload==3:
        f.close()
        #k = mkl.get_max_threads()
        k = numba.get_num_threads()
        dataArray = mp.Array('f',k*1024*1024,lock=False)     # should be in shared memory
        indicesArray = mp.Array(indices_shm_type,k*1024*1024,lock=False)  # should be in shared memory
        startsArray = mp.Array('l',k,lock=False)  # start index of data read
        endsArray = mp.Array('l',k,lock=False)    # end index (noninclusive) of data read
        barrier = mp.Barrier(k+1)            # k workers plus this process
        dataA = np.frombuffer(dataArray,dtype=np.float32)
        indicesA = np.frombuffer(indicesArray,dtype=indices_type)
        startsA = np.frombuffer(startsArray, dtype=np.int64)
        endsA = np.frombuffer(endsArray, dtype=np.int64)
        data = np.empty(datalen, dtype=np.float32)
        indices = np.empty(datalen, dtype=indices_type)
        print("created shared memory buffers",time.time()-t0)

        procs = [mp.Process(target=_load_helper2, args=(fname, i, k, datalen, barrier, dataArray, indicesArray, startsArray, endsArray))  for i in range(k)]
        for p in procs: p.start()

        for _ in range(datalen//(k*1024*1024)+1):
            barrier.wait()  # wait for data loads
            _fast_copy(data,dataA,startsA,endsA,k)
            _fast_copy(indices,indicesA,startsA,endsA,k)
            barrier.wait()  # finish data copy

        for p in procs: p.join()
        print("multiproc load done",time.time()-t0)
    
    if use_fastload==4 or use_fastload==5:
        f.close()
        #k = mkl.get_max_threads()
        k = numba.get_num_threads()
        dataArray = mp.Array('f',k*1024*1024,lock=False)     # should be in shared memory
        indicesArray = mp.Array(indices_shm_type,k*1024*1024,lock=False)  # should be in shared memory
        startsArray = mp.Array('l',k,lock=False)  # start index of data read
        endsArray = mp.Array('l',k,lock=False)    # end index (noninclusive) of data read
        global semDataLoaded
        global semDataCopied
        semDataLoaded = [mp.Semaphore(0) for _ in range(k)]
        semDataCopied = [mp.Semaphore(0) for _ in range(k)]
        dataA = np.frombuffer(dataArray,dtype=np.float32)
        indicesA = np.frombuffer(indicesArray,dtype=indices_type)
        startsA = np.frombuffer(startsArray, dtype=np.int64)
        endsA = np.frombuffer(endsArray, dtype=np.int64)
        data = np.empty(datalen, dtype=np.float32)
        indices = np.empty(datalen, dtype=indices_type)
        print("created shared memory buffers",time.time()-t0)

        procs = [mp.Process(target=_load_helper3, args=(fname, i, k, datalen, dataArray, indicesArray, startsArray, endsArray))  for i in range(k)]
        for p in procs: p.start()

        if use_fastload==4:
            threads = [threading.Thread(target=_copy_thread, args=(data,dataA,indices,indicesA,startsA,endsA,i,datalen//(k*1024*1024)+1)) for i in range(k)]
            for t in threads: t.start()
        else:
            _fast_copy3(data,dataA,indices,indicesA,startsA,endsA,k,datalen//(k*1024*1024)+1)

        for p in procs: p.join()
        if use_fastload==4:
            for t in threads: t.join()
        print("multiproc load done",time.time()-t0)
    
    #X = scipy.sparse.csr_matrix((data,indices,indptr))   #This is too slow
    X = scipy.sparse.csr_matrix((0,0))
    #X = scipy.sparse.csr_matrix((rows,dfvar.shape[0]))
    X.data = data
    X.indices = indices
    #X.data = np.empty(data.shape,dtype=data.dtype)
    #X.data[:] = data
    #X.indices = np.empty(indices.shape,dtype=indices.dtype)
    #X.indices[:] = indices
    X.indptr = indptr
    X._shape = ((rows, dfvar.shape[0]))
    print("Created scipy sparse array",time.time()-t0)

    # create AnnData
    adata = anndata.AnnData(X, dfobs, dfvar)
    print("anndata created",time.time()-t0)
    return adata    


def csr_subset(X, rows_to_keep=None, cols_to_keep=None):

    @numba.njit(cache=True)
    def _csr_subset(rmask, cmask, cols, indptr, data, indices):
        wrptr = 0
        rows = 0
        if cmask is not None:
            colmap = np.cumsum(cmask)
            cols = colmap[-1]
        for r in range(len(indptr)-1):
            if rmask is not None and not rmask[r]: continue
            rdstart = indptr[r]
            indptr[rows] = wrptr
            rows += 1
            for rdptr in range(rdstart,indptr[r+1]):
                c = indices[rdptr]
                if cmask is not None and not cmask[c]: continue
                indices[wrptr] = c if cmask is None else colmap[c]-1
                data[wrptr] = data[rdptr]
                wrptr += 1
        indptr[rows:] = wrptr
        indices[wrptr:]=cols
        return rows, cols, wrptr

    #assert (X.shape[0]==rows_to_keep.shape[0])
    rows, cols, datalen = _csr_subset(rows_to_keep, cols_to_keep, X.shape[1],X.indptr, X.data, X.indices)
    #X.indptr.resize(rows+1,refcheck=False)
    #X.data.resize(datalen,refcheck=False)
    #X.indices.resize(datalen,refcheck=False)
    X.resize((rows,cols))
    return X

#def csr_row_subset(X, rows_to_keep):
#
#    @numba.njit(cache=True)
#    def _csr_row_subset(rmask, indptr, data, indices):
#        wrptr = 0
#        rows = 0
#        for r in range(len(indptr)-1):
#            if not rmask[r]: continue
#            rdstart = indptr[r]
#            indptr[rows] = wrptr
#            rows += 1
#            for rdptr in range(rdstart,indptr[r+1]):
#                data[wrptr] = data[rdptr]
#                indices[wrptr] = indices[rdptr]
#                wrptr += 1
#        indptr[rows] = wrptr
#        indptr[rows+1:] = indptr[-1]
#        return rows, wrptr
#
#    #assert (X.shape[0]==rows_to_keep.shape[0])
#    rows, datalen = _csr_row_subset(rows_to_keep, X.indptr, X.data, X.indices)
#    #X.indptr.resize(rows+1,refcheck=False)
#    #X.data.resize(datalen,refcheck=False)
#    #X.indices.resize(datalen,refcheck=False)
#    X.resize((rows,X.shape[1]))
#    return X

#def csr_row_subset(X, rows_to_keep):
#
#    @numba.njit(cache=True)
#    def _csr_row_subset(rmask, indptr, data, indices):
#        wrptr = 0
#        rows = 0
#        for r in range(len(indptr)-1):
#            if not rmask[r]: continue
#            rdstart = indptr[r]
#            indptr[rows] = wrptr
#            rows += 1
#            for rdptr in range(rdstart,indptr[r+1]):
#                data[wrptr] = data[rdptr]
#                indices[wrptr] = indices[rdptr]
#                wrptr += 1
#        indptr[rows] = wrptr
#        indptr[rows+1:] = indptr[-1]
#        return rows, wrptr
#
#    #assert (X.shape[0]==rows_to_keep.shape[0])
#    rows, datalen = _csr_row_subset(rows_to_keep, X.indptr, X.data, X.indices)
#    X.resize((rows,X.shape[1]))
#    return X

def csr_row_subset(X, rows_to_keep):

    @numba.njit(cache=True)
    def _csr_row_subset(rmask, indptr, data, indices, junkcol):
        rows = 0
        for r in range(len(indptr)-1):
            if rmask[r]:
                rows += 1
                indptr[rows] = indptr[r+1]
            else:
                indices[indptr[r]:indptr[r+1]] = junkcol
        indptr[rows:] = indptr[-1]
        return rows

    #assert (X.shape[0]==rows_to_keep.shape[0])
    rows = _csr_row_subset(rows_to_keep, X.indptr, X.data, X.indices, X.shape[1])
    #rows = _csr_row_subset(rows_to_keep, X.indptr, X.data, X.indices, -1)
    # cleared out and shifted pointers for rows, but don't actually remove them yet
    return rows

def csr_row_subset2(X, rows_to_keep):

    @numba.njit(cache=True,parallel=True)
    def _csr_row_subset2(rmask, indptr, data, indices, junkcol):
        rows = 0
        #print(len(indptr),len(rmask),len(data),junkcol)
        for r in numba.prange(len(indptr)-1):
            if rmask[r]:
                rows += 1
            else:
                indices[indptr[r]:indptr[r+1]] = junkcol
                #data[indptr[r]:indptr[r+1]] = 0
        return rows

    rows = _csr_row_subset2(rows_to_keep, X.indptr, X.data, X.indices, X.shape[1])
    #rows = _csr_row_subset2(rows_to_keep, X.indptr, X.data, X.indices, 1000000000)
    # cleared out indices for "removed" rows, but don't actually remove them
    return rows

def csr_col_subset(X, cols_to_keep):

    @numba.njit(cache=True, parallel=True)
    def _csr_col_subset(cmask, indices, data):
        colmap = np.cumsum(cmask)
        cols = colmap[-1]
        junkcol = len(colmap)
        #junkcol = 1000000000
        #print (junkcol,len(colmap),cols,len(indices),len(data))
        for i in numba.prange(len(indices)):
            c = indices[i]
            if c<len(colmap) and cmask[c]:
                indices[i] = colmap[c]-1
            else:
                indices[i] = junkcol
                #data[i] = 0
            #indices[i] = colmap[c]-1 if cmask[c] else junkcol
        return cols

    #assert (X.shape[0]==rows_to_keep.shape[0])
    cols = _csr_col_subset(cols_to_keep, X.indices, X.data)
    # renumbered all columns, marking ones for deletion out of range;  don't actually remove yet
    return cols


def finalize_drop(adata,nrows,keeprows,ncols,keepcols):
    t0 = time.time()
    adata.X.resize((nrows,ncols))
    adata = anndata.AnnData(adata.X,adata.obs.iloc[keeprows,:],adata.var.iloc[keepcols,:])
    print("finalize drops: ",time.time()-t0)
    return adata

def finalize_drop_dense(adata, nrows, keeprows, keepcols0, ncols0, ncols, keepcols):
    @numba.njit(cache=True,parallel=True)
    def drop_to_dense(shp, indptr, indices, data, nrows, keeprows, ncols):
        X = np.zeros((nrows,ncols),dtype=data.dtype)
        rowmap = np.cumsum(keeprows)
        for r in numba.prange(shp[0]):
            #if not keeprows[r]: continue
            newr = rowmap[r]-1
            #newr = r
            sparserow = data[indptr[r]:indptr[r+1]]
            sparseind = indices[indptr[r]:indptr[r+1]]
            for j in range(len(sparserow)):
                c = sparseind[j]
                if c<ncols:
                    X[newr,c] = sparserow[j]
        return X

    t0 = time.time()
    newX = drop_to_dense(adata.X.shape, adata.X.indptr, adata.X.indices, adata.X.data, nrows, keeprows, ncols)
    adata = anndata.AnnData(newX, adata.obs.iloc[keeprows,:], adata.var.iloc[keepcols,:])
    #adata = anndata.AnnData(newX, adata.obs, adata.var.iloc[keepcols,:])
    #adata = anndata.AnnData(newX, adata.obs.iloc[keeprows,:], adata.var.iloc[keepcols0,:].iloc[keepcols[:ncols0],:])
    #adata = anndata.AnnData(newX, adata.obs, adata.var.iloc[keepcols0,:].iloc[keepcols[:ncols0],:])
    print("finalize drops, densify: ",time.time()-t0)
    return adata

def normtotal_log1p( adata, target, ncols ):
    t0 = time.time()
    if use_fastpp:
        #print (ncols)
        @numba.njit(cache=True,parallel=True)
        def _normtotal_log1p( data, ind, indices, tt, ncols ):
            for i in numba.prange(len(ind)-1):
                row = data[ ind[i]:ind[i+1] ]
                rowind = indices[ ind[i]:ind[i+1] ]
                #total = np.sum(row)
                total = 1e-12
                for j in range(len(row)):
                    if rowind[j]<ncols: total+=row[j]
                sf = tt/total
                for j in range(len(row)):
                    row[j] = math.log1p( row[j]*sf )

        _normtotal_log1p( adata.X.data, adata.X.indptr, adata.X.indices, target, ncols )

    else:
        sc.pp.normalize_total(adata, target_sum=target)
        print ("norm_log: normalize", time.time()-t0)
        sc.pp.log1p(adata)
    print ("norm_log: total", time.time()-t0)


def filter_cells(adata, ming, maxg):
    t0 = time.time()
    rows_to_keep = 0
    rows = 0

    if use_fastpp:
        @numba.njit(cache=True,parallel=True)
        def get_rows_to_keep( indptr, ming, maxg):
            lens = indptr[1:] - indptr[:-1]
            keep_rows = np.logical_and(lens>=ming, lens<=maxg)
            return lens, keep_rows

        nrows = adata.X.shape[0]
        nelems = adata.X.data.shape[0]
        nthr = numba.get_num_threads()
        print ("filter_cells:  prep ", time.time()-t0)
        row_lengths, rows_to_keep = get_rows_to_keep(adata.X.indptr, ming, maxg)
        print ("filter_cells:  compute ", time.time()-t0)
        adata.obs['n_genes'] = row_lengths
        print ("filter_cells:  set metadata ", time.time()-t0)
        if strategy==1: adata = adata[rows_to_keep]
        if strategy==2: adata._inplace_subset_obs(rows_to_keep)
        if strategy==3: adata = adata[rows_to_keep].copy()
        if strategy==4: adata = anndata.AnnData(adata.X[rows_to_keep],adata.obs.iloc[rows_to_keep,:],adata.var)
        if strategy==5: adata = anndata.AnnData(adata.X[rows_to_keep],adata.obs.drop(adata.obs.iloc[np.logical_not(rows_to_keep),:],inplace=True),adata.var)
        if strategy==6: adata._inplace_subset_obs(rows_to_keep)
        if strategy==7: adata = anndata.AnnData(csr_subset(adata.X,rows_to_keep),adata.obs.iloc[rows_to_keep,:],adata.var)
        if strategy==8: rows = csr_row_subset(adata.X,rows_to_keep)
        if strategy==9: adata._inplace_subset_obs(rows_to_keep)
        if strategy==10: adata = anndata.AnnData(adata.X[rows_to_keep],adata.obs.iloc[rows_to_keep,:],adata.var)
        if strategy==11: rows = csr_row_subset2(adata.X,rows_to_keep)
        #if strategy==11: adata = anndata.AnnData(adata.X[rows_to_keep],adata.obs.iloc[rows_to_keep,:],adata.var)
    else:
        sc.pp.filter_cells(adata, min_genes=ming)
        print ("filter_cells:  first call ", time.time()-t0)
        sc.pp.filter_cells(adata, max_genes=maxg)

    print ("filter_cells:  filter total", time.time()-t0)
    if rows==0: rows = adata.shape[0]
    return adata, rows_to_keep, rows

def filter_genes(adata, minc):
    t0 = time.time()
    cols_to_keep = 0
    cols = 0

    if use_fastpp:
        @numba.njit(cache=True,parallel=True)
        def get_cols_to_keep( indices, data, minc, colcount, nthr):
            counts = np.zeros((nthr,colcount),dtype=np.int32)
            for i in numba.prange(nthr):
                start = i*indices.shape[0]//nthr
                end = (i+1)*indices.shape[0]//nthr
                for j in range(start,end):
                    if data[j]!=0 and indices[j]<colcount:
                    #if indices[j]<colcount:
                        counts[i,indices[j]] += 1
            counts = np.sum(counts, axis=0)
            keep_cols = counts>=minc
            return counts, keep_cols

        ncols = adata.X.shape[1]
        nthr = numba.get_num_threads()
        print ("filter_genes:  prep ", time.time()-t0)
        counts, cols_to_keep = get_cols_to_keep(adata.X.indices, adata.X.data, minc, ncols, nthr)
        print ("filter_genes:  compute ", time.time()-t0)
        adata.var['n_cells'] = counts
        print ("filter_genes:  set metadata ", time.time()-t0)
        if strategy==1: adata = adata[:,cols_to_keep]
        if strategy==2: adata._inplace_subset_var(cols_to_keep)
        if strategy==3: adata = adata[:,cols_to_keep].copy()
        if strategy==4: adata = anndata.AnnData(adata.X[:,cols_to_keep],adata.obs,adata.var.iloc[cols_to_keep,:])
        if strategy==5: adata = anndata.AnnData(adata.X[:,cols_to_keep],adata.obs,adata.var.drop(adata.var.iloc[np.logical_not(cols_to_keep),:],inplace=True))
        if strategy==6: adata._inplace_subset_var(cols_to_keep)
        if strategy==7: adata = anndata.AnnData(csr_subset(adata.X,None,cols_to_keep),adata.obs,adata.var.iloc[cols_to_keep,:])
        if strategy==8: cols = csr_col_subset(adata.X, cols_to_keep)
        if strategy==9: adata._inplace_subset_var(cols_to_keep)
        if strategy==10: adata = anndata.AnnData(adata.X[:,cols_to_keep],adata.obs,adata.var.iloc[cols_to_keep,:])
        if strategy==11: cols = csr_col_subset(adata.X, cols_to_keep)
    else:
        sc.pp.filter_genes(adata, min_cells=minc)

    print ("filter_genes:  filter total", time.time()-t0)
    if cols==0: cols = adata.shape[1]
    return adata, cols_to_keep, cols


def scale(adata,maxv):
    #from scanpy.pp._utils import _get_mean_var
    #mean, var = _get_mean_var(adata.X)
    #print(time.time()-t0)
    #std = np.sqrt(var)
    #print(time.time()-t0)
    #adata.X -= mean
    #print(time.time()-t0)
    #adata.X /= std
    #print(time.time()-t0)
    #adata.X[adata.X > maxv] = maxv

    #@numba.njit(cache=True,parallel=True)
    #def do_scale(X, maxv):
    #  #nthr = numba.get_num_threads()
    #  #for i in numba.prange(nthr):
    #    for c in numba.prange(X.shape[1]):
    #    #for c in range(i,X.shape[1],nthr):
    #        s =0.0
    #        ss=0.0
    #        n = X.shape[0]
    #        for r in range(n):
    #            v=X[r,c]
    #            s+=v
    #            ss+=v*v
    #        mean = s/n
    #        var = (ss - s*s/n)/(n-1)
    #        invstd = 1.0/np.sqrt(var)
    #        for r in range(n):
    #            v = (X[r,c]-mean)*invstd
    #            X[r,c] = maxv if v>maxv else v

    #@numba.njit(cache=True,parallel=True)
    #def do_scale(X, maxv):
    #    s=np.zeros(X.shape[1])
    #    ss=np.zeros(X.shape[1])
    #    mean=np.zeros(X.shape[1])
    #    std=np.zeros(X.shape[1])
    #    nthr = numba.get_num_threads()
    #    n = X.shape[0]
    #    for i in numba.prange(nthr):
    #        for r in range(n):
    #            startc = i*X.shape[1]//nthr
    #            stopc = (i+1)*X.shape[1]//nthr
    #            for c in range(startc,stopc):
    #                v=X[r,c]
    #                s[c]+=v
    #                ss[c]+=v*v
    #        for c in range(startc,stopc):
    #            mean[c] = s[c]/n
    #            std[c] = np.sqrt((ss[c] - s[c]*s[c]/n)/(n-1))
    #        for r in range(n):
    #            for c in range(startc,stopc):
    #                v = (X[r,c]-mean[c])/std[c]
    #                X[r,c] = maxv if v>maxv else v

    @numba.njit(cache=True,parallel=True)
    def do_scale(X, maxv, nthr):
        with numba.objmode(t0='f8'):
            t0 = time.time()
        #nthr = numba.get_num_threads()
        s=np.zeros((nthr,X.shape[1]))
        ss=np.zeros((nthr,X.shape[1]))
        mean=np.zeros(X.shape[1])
        std=np.zeros(X.shape[1])
        n = X.shape[0]
        for i in numba.prange(nthr):
            for r in range(i,n,nthr):
                for c in range(X.shape[1]):
                    v=X[r,c]
                    s[i,c]+=v
                    ss[i,c]+=v*v
        for c in numba.prange(X.shape[1]):
            s0 = s[:,c].sum()
            mean[c] = s0/n
            std[c] = np.sqrt((ss[:,c].sum() - s0*s0/n)/(n-1))
        with numba.objmode():
            print ("finshed getting means, stddev", time.time()-t0)
        for (r,c) in numba.pndindex(X.shape):
            v = (X[r,c]-mean[c])/std[c]
            X[r,c] = maxv if v>maxv else v
        with numba.objmode():
            print ("finshed scaling values", time.time()-t0)

    #@numba.njit(cache=True,parallel=True)
    #def do_scale(X, maxv, nthr):
    #    with numba.objmode(t0='f8'):
    #        t0 = time.time()
    #    s=np.zeros(X.shape[1])
    #    ss=np.zeros(X.shape[1])
    #    mean=np.zeros(X.shape[1])
    #    std=np.zeros(X.shape[1])
    #    n = X.shape[0]
    #    for r,c in numba.pndindex(X.shape):
    #        v=X[r,c]
    #        s[c]+=v
    #        ss[c]+=v*v
    #    for c in numba.prange(X.shape[1]):
    #        mean[c] = s[c]/n
    #        std[c] = np.sqrt((ss[c] - s[c]*s[c]/n)/(n-1))
    #    with numba.objmode():
    #        print ("finshed getting means, stddev", time.time()-t0)
    #    for (r,c) in numba.pndindex(X.shape):
    #        v = (X[r,c]-mean[c])/std[c]
    #        X[r,c] = maxv if v>maxv else v
    #    with numba.objmode():
    #        print ("finshed scaling values", time.time()-t0)

    if use_fastpp:
        do_scale(adata.X, maxv, numba.get_num_threads())
    else:
        sc.pp.scale(adata, max_value=maxv)

    debug and print(adata.obs)
    debug and print(adata.var)
    debug and print(adata.X[0], sum(adata.X[0]))
    debug and print(adata.X[1], sum(adata.X[1]))
    debug and print(resource.getrusage(resource.RUSAGE_SELF))
    



def loadpreprocess(input_file,firstn,min_genes_per_cell,max_genes_per_cell,min_cells_per_gene,target_sum,markers,n_top_genes):
    debug and print(resource.getrusage(resource.RUSAGE_SELF))

    start = time.time()
    print ("Loading data")
    if use_fastpp and use_fastload:
        adata = fastload(input_file, firstn)
    else:
        adata = sc.read(input_file)
    #adata = anndata.read_h5ad(input_file)
    #cProfile.run("adata = sc.read(input_file)")
    #cProfile.runctx("adata = anndata.read_h5ad(input_file)",globals=globals(),locals=locals())
    print ("Load data time",time.time()-start)
    adata.var_names_make_unique()

    debug and print(resource.getrusage(resource.RUSAGE_SELF))
    if firstn is not None and 0<firstn and firstn<adata.shape[0]:
        t0 = time.time()
        if use_fastpp:
            if strategy==1: adata = adata[0:firstn]
            if strategy==2: adata._inplace_subset_obs(slice(0,firstn))
            if strategy==3: adata = adata[0:firstn].copy()
            if strategy==4: adata = anndata.AnnData(adata.X[0:firstn],adata.obs.iloc[0:firstn,:],adata.var)
            if strategy==5: adata = anndata.AnnData(adata.X[0:firstn],adata.obs.drop(adata.obs.iloc[firstn:,:],inplace=True),adata.var)
            if strategy==6: adata = adata[0:firstn]
            if strategy>=7: 
                adata.X.resize((firstn,adata.X.shape[1]))
                adata = anndata.AnnData(adata.X,adata.obs.iloc[0:firstn,:],adata.var)
        else:
            adata = adata[0:firstn]
        print ("downselect", time.time()-t0)
    print(adata.shape)

    if use_fastpp:
        debug and print(adata.obs)
        debug and print(adata.var)
        debug and print(resource.getrusage(resource.RUSAGE_SELF))
        filter_time = time.time()
        # We filter the count matrix to remove cells with an extreme number of genes expressed.
        #sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
        #sc.pp.filter_cells(adata, max_genes=max_genes_per_cell)
        adata, keeprows, nrows = filter_cells(adata, min_genes_per_cell, max_genes_per_cell)
        #cProfile.runctx("adata = filter_cells(adata, min_genes_per_cell, max_genes_per_cell)",globals(),locals())
        t0 = time.time()
        print ("filter cells time", t0-filter_time)
        debug and print (adata.shape)
        debug and print(resource.getrusage(resource.RUSAGE_SELF))
        # Some genes will now have zero expression in all cells. We filter out such genes.
        #sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
        adata, keepcols, ncols = filter_genes(adata, min_cells_per_gene)
        print ("filter genes time", time.time()-t0)
        if strategy==8:  adata = finalize_drop(adata,nrows,keeprows,ncols,keepcols)
        #adata.X.sort_indices()
        print(adata.shape)
        debug and print(adata.obs)
        debug and print(adata.var)
        debug and print(adata.X.data[0:1807], sum(adata.X.data[0:1807]))
        debug and print(adata.X.data[1807:3056], sum(adata.X.data[1807:3056]))
        debug and print(resource.getrusage(resource.RUSAGE_SELF))

        # Normalize
        normtotal_log1p(adata, target_sum, ncols)

        debug and print("after normalize, log1p")
        debug and print(adata.obs)
        debug and print(adata.var)
        debug and print(adata.X.data[0:1807], sum(adata.X.data[0:1807]))
        debug and print(adata.X.data[1807:3056], sum(adata.X.data[1807:3056]))
        debug and print(resource.getrusage(resource.RUSAGE_SELF))
        

        @numba.njit(cache=True,parallel=True)
        def get_dense_col(col, indptr, indices, data):
            rows = indptr.shape[0]-1
            X=np.zeros(rows,dtype=data.dtype)
            for r in numba.prange(rows):
                for i in range(indptr[r],indptr[r+1]):
                    if indices[i]==col:  X[r] = data[i]
            return X

        # Select Most Variable Genes
        t0=time.time()
        # Select highly variable genes
        #sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor = "cell_ranger")
        sc_pp_hvg._get_hvg(adata, nrows, ncols, keepcols, n_top_genes=n_top_genes)
        print (time.time()-t0)
        # Retain marker gene expression
        for marker in markers:
            #adata.obs[marker + "_raw"] = adata.X[:, adata.var.index == marker].todense()
            if strategy==11:
                col = np.argmax( adata.var.index[keepcols] == marker )
            else:
                col = np.argmax( adata.var.index == marker )
            adata.obs[marker + "_raw"] = get_dense_col(col, adata.X.indptr, adata.X.indices, adata.X.data)
        print (time.time()-t0)
        
        # Filter matrix to only variable genes
        if strategy==1: adata = adata[:,adata.var.highly_variable]
        if strategy==2: adata._inplace_subset_var(adata.var.highly_variable)
        if strategy==3: adata = adata[:,adata.var.highly_variable].copy()
        if strategy==4:
            adata = anndata.AnnData( adata.X[:,adata.var.highly_variable.values], adata.obs,adata.var.iloc[adata.var.highly_variable.values,:])
        if strategy==5:
            adata = anndata.AnnData( adata.X[:,adata.var.highly_variable.values], adata.obs,adata.var.drop(adata.var.iloc[np.logical_not(adata.var.highly_variable.values),:],inplace=True))
        if strategy==6: adata = adata[:,adata.var.highly_variable]
        if strategy==7: adata = anndata.AnnData(csr_subset(adata.X,None,adata.var.highly_variable.values),adata.obs,adata.var.iloc[adata.var.highly_variable.values,:])
        if strategy==8:
            adata = anndata.AnnData( adata.X[:,adata.var.highly_variable.values], adata.obs,adata.var.iloc[adata.var.highly_variable.values,:])
        if strategy==9: adata._inplace_subset_var(adata.var.highly_variable)
        if strategy==10:
            adata = anndata.AnnData( adata.X[:,adata.var.highly_variable.values], adata.obs,adata.var.iloc[adata.var.highly_variable.values,:])
        if strategy==11: 
            csr_col_subset(adata.X, adata.var.highly_variable.values[keepcols])
            adata = finalize_drop_dense(adata, nrows, keeprows, keepcols, ncols, n_top_genes, adata.var.highly_variable.values)
        print("downselect to top genes", time.time()-t0)
        print(adata.shape)

    else:  # Original Code
        debug and print(adata.obs)
        debug and print(adata.var)
        debug and print(resource.getrusage(resource.RUSAGE_SELF))
        filter_time = time.time()
        # We filter the count matrix to remove cells with an extreme number of genes expressed.
        sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
        sc.pp.filter_cells(adata, max_genes=max_genes_per_cell)
        t0 = time.time()
        print ("filter cells time", t0-filter_time)
        # Some genes will now have zero expression in all cells. We filter out such genes.
        sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
        print ("filter genes time", time.time()-t0)
        print(adata.shape)
        debug and print(adata.obs)
        debug and print(adata.var)
        debug and print(adata.X.data[0:1807], sum(adata.X.data[0:1807]))
        debug and print(adata.X.data[1807:3056], sum(adata.X.data[1807:3056]))
        debug and print(resource.getrusage(resource.RUSAGE_SELF))

        # Normalize
        t0=time.time()
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        print("norm, log time:",time.time()-t0)
        debug and print(adata.obs)
        debug and print(adata.var)
        debug and print(adata.X.data[0:10])
        debug and print(adata.X.data[1807:1817])
        debug and print(resource.getrusage(resource.RUSAGE_SELF))

        # Select Most Variable Genes
        t0=time.time()
        # Select highly variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor = "cell_ranger")
        print (time.time()-t0)
        # Retain marker gene expression
        for marker in markers:
            adata.obs[marker + "_raw"] = adata.X[:, adata.var.index == marker].todense()
        print (time.time()-t0)
        
        # Filter matrix to only variable genes
        adata = adata[:, adata.var.highly_variable]
        print("downselect to top genes", time.time()-t0)
        print(adata.shape)

    print("Total Filtering time : %s\n" % (time.time() - filter_time))
    return adata


def regress_out(adata, regr):
    debug and print(adata.obs)
    debug and print(adata.var)
    if use_fastpp and regr_type>0:
        #myscpp.regress_out(adata, regr)
        numpy_regress_out(adata, regr)
    else:
        kthr = mkl.get_max_threads()
        debug and print ("MKL threads was at", kthr,"setting to 1")
        mkl.set_num_threads(1)
        sc.pp.regress_out(adata, regr)
        mkl.set_num_threads(kthr)
        debug and print ("MKL threads restored to",kthr)

@numba.njit(cache=True,parallel=True)
def to_dense(shp, indptr, indices, data):
    X = np.empty(shp, dtype=data.dtype)
    for r in numba.prange(shp[0]):
        X[r] = 0
        for i in range(indptr[r],indptr[r+1]):
            X[r,indices[i]] = data[i]
    return X


from scipy.sparse import issparse
def numpy_regress_out(adata, regr):

    @numba.njit(cache=True,parallel=True)
    def get_resid(X, A, coeff):
        for i in numba.prange(X.shape[0]):
            X[i] -= A[i]@coeff
    
    #@numba.njit(cache=True,parallel=True)
    #def _ATxB(A, B, nthr):   # fast A.T@B, assumes A is large by tiny, B is large by meduium sized
    #    out = np.zeros((A.shape[1],B.shape[1]),dtype=A.dtype)
    #    for i in numba.prange(nthr):
    #        cstart = i*B.shape[1]//nthr
    #        cend = (i+1)*B.shape[1]//nthr
    #        for j in range(A.shape[0]):
    #            #for c in range(cstart,cend):
    #            #    bval = B[j,c]
    #            #    for r in range(A.shape[1]):
    #            #        out[r,c] += A[j,r]*bval
    #            bval = B[j,cstart:cend]
    #            for r in range(A.shape[1]):
    #                out[r,cstart:cend] += A[j,r]*bval
    #    return out

    @numba.njit(cache=True,parallel=True)
    def _ATxB(A, B, nthr):   # fast A.T@B, assumes A is large by tiny, B is large by meduium sized
        out = np.zeros((nthr,A.shape[1],B.shape[1]),dtype=A.dtype)
        for i in numba.prange(nthr):
            for j in range(i,A.shape[0],nthr):
                bval = B[j]
                for r in range(A.shape[1]):
                    out[i,r] += A[j,r]*bval
        return out.sum(axis=0)

    start = time.time()
    if issparse(adata.X):
        #adata.X = adata.X.toarray()
        adata.X = to_dense(adata.X.shape, adata.X.indptr, adata.X.indices, adata.X.data)
    print("convert to dense done at", time.time()-start)
    A = np.empty((adata.shape[0], len(regr)+1),dtype=adata.X.dtype)
    A[:,0]=1
    for i,c in enumerate(regr):
        A[:,i+1] = adata.obs[c]
    print("prep done at", time.time()-start)
    if regr_type in [1,2]:
        coeff = np.linalg.lstsq(A,adata.X)[0]
    else:
        #tmp = np.linalg.inv(A.T@A)
        if (regr_type>4):
            tmp0 = _ATxB(A,A, numba.get_num_threads())
        else:
            tmp0 = A.T@A
        print(time.time()-start)
        tmp = np.linalg.inv(tmp0)
        print(time.time()-start)
        #coeff = (tmp@A.T)@adata.X
        #coeff = tmp@(A.T@adata.X)  # this should produce a smaller intermediate array
        if (regr_type>4):
            tmp0 = _ATxB(A,adata.X, numba.get_num_threads())
        else:
            tmp0 = A.T@adata.X
        print(time.time()-start)
        coeff = tmp@tmp0
    #print (coeff)
    print("regression done at", time.time()-start)
    if regr_type in [1,3]:
        adata.X -= A@coeff
    else:
        get_resid(adata.X, A, coeff)
    print("residuals done at", time.time()-start)

def sum( X, axis=None ):
    @numba.njit(cache=True,parallel=True)
    def _sum(X):
        s = 0
        for i in numba.pndindex(X.shape):
            s+=X[i]
        return s

    @numba.njit(cache=True,parallel=True)
    def _sum0(X, nthr):
        s = np.empty((nthr,X.shape[1]), dtype=X.dtype)
        for i in numba.prange(nthr):
            for r in range(i,X.shape[0],nthr):
                s[i] = X[r]
        return s.sum(axis=0)

    @numba.njit(cache=True,parallel=True)
    def _sum1(X):
        s = np.empty(X.shape[0], dtype=X.dtype)
        for r in numba.prange(X.shape[0]):
            s[r] = X[r].sum()
        return s


    if issparse(X) or not use_fastpp:
        if axis is None:
            return np.array(X.sum())
        return np.array(X.sum(axis=axis))
    if axis is None:
        return _sum(X)
    if X.ndim==2:
        if axis==0:
            nthr = numba.get_num_threads()
            return _sum0(X, nthr)
        return _sum1(X)
    # if axis is None:
    #     return X.sum()
    return X.sum(axis=axis)

