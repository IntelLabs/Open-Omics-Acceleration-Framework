#===============================================================================
# Copyright 2020-2021 Intel Corporation
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#===============================================================================

# daal4py TSNE scikit-learn-compatible class

import warnings
from time import time
import numpy as np
from scipy.sparse import issparse
import daal4py
from daal4py.sklearn._utils import daal_check_version, sklearn_check_version

from sklearn.manifold import TSNE as BaseTSNE
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.utils.validation import check_non_negative
from sklearn.utils import check_random_state, check_array

from ..neighbors import NearestNeighbors
from .._device_offload import support_usm_ndarray

# if sklearn_check_version('0.22'):
#     from sklearn.manifold._t_sne import _joint_probabilities, _joint_probabilities_nn
# else:
#     from sklearn.manifold.t_sne import _joint_probabilities, _joint_probabilities_nn


from scipy.spatial.distance import squareform
from scipy.sparse import csr_matrix
import math
import numba
MACHINE_EPSILON = np.finfo(np.double).eps
# print(MACHINE_EPSILON)

@numba.njit(cache=True,parallel=True, nogil=True)
def _binary_search_perplexity(sqdistances, desired_perplexity,verbose):

    NPY_INFINITY = np.inf
    EPSILON_DBL = 1e-8
    PERPLEXITY_TOLERANCE = 1e-5
    """Binary search for sigmas of conditional Gaussians.
    This approximation reduces the computational complexity from O(N^2) to
    O(uN).
    Parameters
    ----------
    sqdistances : array-like, shape (n_samples, n_neighbors)
        Distances between training samples and their k nearest neighbors.
        When using the exact method, this is a square (n_samples, n_samples)
        distance matrix. The TSNE default metric is "euclidean" which is
        interpreted as squared euclidean distance.
    desired_perplexity : float
        Desired perplexity (2^entropy) of the conditional Gaussians.
    verbose : int
        Verbosity level.
    Returns
    -------
    P : array, shape (n_samples, n_samples)
        Probabilities of conditional Gaussian distributions p_i|j.
    """
    # Maximum number of binary search steps
    n_steps = 100

    n_samples = sqdistances.shape[0]
    n_neighbors = sqdistances.shape[1]
    using_neighbors = n_neighbors < n_samples
    # Precisions of conditional Gaussian distributions
    # double beta
    # double beta_min
    # double beta_max
    beta_sum = 0.0

    # Use log scale
    desired_entropy = math.log(desired_perplexity)
    # double entropy_diff

    # cdef double entropy
    # cdef double sum_Pi
    # cdef double sum_disti_Pi
    # cdef long i, j, k, l

    # This array is later used as a 32bit array. It has multiple intermediate
    # floating point additions that benefit from the extra precision
    P = np.zeros((n_samples, n_neighbors), dtype=np.float64)

    # for i in range(n_samples):
    for i in numba.prange(n_samples):
        beta_min = -NPY_INFINITY
        beta_max = NPY_INFINITY
        beta = 1.0

        # Binary search of precision for i-th conditional distribution
        for l in range(n_steps):
            # Compute current entropy and corresponding probabilities
            # computed just over the nearest neighbors or over all data
            # if we're not using neighbors
            sum_Pi = 0.0
            for j in range(n_neighbors):
                if j != i or using_neighbors:
                    P[i, j] = math.exp(-sqdistances[i, j] * beta)
                    sum_Pi += P[i, j]

            if sum_Pi == 0.0:
                sum_Pi = EPSILON_DBL
            sum_disti_Pi = 0.0

            for j in range(n_neighbors):
                P[i, j] /= sum_Pi
                sum_disti_Pi += sqdistances[i, j] * P[i, j]

            entropy = math.log(sum_Pi) + beta * sum_disti_Pi
            entropy_diff = entropy - desired_entropy

            if math.fabs(entropy_diff) <= PERPLEXITY_TOLERANCE:
                break

            if entropy_diff > 0.0:
                beta_min = beta
                if beta_max == NPY_INFINITY:
                    beta *= 2.0
                else:
                    beta = (beta + beta_max) / 2.0
            else:
                beta_max = beta
                if beta_min == -NPY_INFINITY:
                    beta /= 2.0
                else:
                    beta = (beta + beta_min) / 2.0

        beta_sum += beta

        # if verbose and ((i + 1) % 1000 == 0 or i + 1 == n_samples):
        #     print("[t-SNE] Computed conditional probabilities for sample "
        #           "%d / %d" % (i + 1, n_samples))

    # if verbose:
    #     print("[t-SNE] Mean sigma: %f"
    #           % np.mean(math.sqrt(n_samples / beta_sum)))
    
    # P = P.astype(sqdistances.dtype)
    # return P
    return P.astype(sqdistances.dtype)

def _joint_probabilities(distances, desired_perplexity, verbose):
    """Compute joint probabilities p_ij from distances.
    Parameters
    ----------
    distances : ndarray of shape (n_samples * (n_samples-1) / 2,)
        Distances of samples are stored as condensed matrices, i.e.
        we omit the diagonal and duplicate entries and store everything
        in a one-dimensional array.
    desired_perplexity : float
        Desired perplexity of the joint probability distributions.
    verbose : int
        Verbosity level.
    Returns
    -------
    P : ndarray of shape (n_samples * (n_samples-1) / 2,)
        Condensed joint probability matrix.
    """
    # Compute conditional probabilities such that they approximately match
    # the desired perplexity
    # distances = distances.astype(np.float32, copy=False)               # Changed by me
    conditional_P = _binary_search_perplexity(
        distances, desired_perplexity, verbose
    )
    P = conditional_P + conditional_P.T
    sum_P = np.maximum(np.sum(P), MACHINE_EPSILON)
    P = np.maximum(squareform(P) / sum_P, MACHINE_EPSILON)
    return P

def _joint_probabilities_nn(distances, desired_perplexity, verbose):
    """Compute joint probabilities p_ij from distances using just nearest
    neighbors.
    This method is approximately equal to _joint_probabilities. The latter
    is O(N), but limiting the joint probability to nearest neighbors improves
    this substantially to O(uN).
    Parameters
    ----------
    distances : sparse matrix of shape (n_samples, n_samples)
        Distances of samples to its n_neighbors nearest neighbors. All other
        distances are left to zero (and are not materialized in memory).
        Matrix should be of CSR format.
    desired_perplexity : float
        Desired perplexity of the joint probability distributions.
    verbose : int
        Verbosity level.
    Returns
    -------
    P : sparse matrix of shape (n_samples, n_samples)
        Condensed joint probability matrix with only nearest neighbors. Matrix
        will be of CSR format.
    """
    t0 = time()
    # Compute conditional probabilities such that they approximately match
    # the desired perplexity
    distances.sort_indices()
    n_samples = distances.shape[0]
    distances_data = distances.data.reshape(n_samples, -1)
    # distances_data = distances_data.astype(np.float32, copy=False)                # Changed by me
    conditional_P = _binary_search_perplexity(
        distances_data, desired_perplexity, verbose
    )
    assert np.all(np.isfinite(conditional_P)), "All probabilities should be finite"

    # Symmetrize the joint probability distribution using sparse operations
    P = csr_matrix(
        (conditional_P.ravel(), distances.indices, distances.indptr),
        shape=(n_samples, n_samples),
    )
    P = P + P.T
    
    # Normalize the joint probability distribution
    sum_P = np.maximum(P.sum(), MACHINE_EPSILON)
    P /= sum_P

    assert np.all(np.abs(P.data) <= 1.0)
    if verbose >= 2:
        duration = time() - t0
        print("[t-SNE] Computed conditional probabilities in {:.3f}s".format(duration))

    return P


class TSNE(BaseTSNE):
    @support_usm_ndarray()
    def fit_transform(self, X, y=None):
        return super().fit_transform(X, y)

    @support_usm_ndarray()
    def fit(self, X, y=None):
        return super().fit(X, y)

    def _daal_tsne(self, P, n_samples, X_embedded):
        """Runs t-SNE."""
        # t-SNE minimizes the Kullback-Leiber divergence of the Gaussians P
        # and the Student's t-distributions Q. The optimization algorithm that
        # we use is batch gradient descent with two stages:
        # * initial optimization with early exaggeration and momentum at 0.5
        # * final optimization with momentum at 0.8

        # N, nnz, n_iter_without_progress, n_iter
        size_iter = np.array([[n_samples], [P.nnz], [self.n_iter_without_progress],
                             [self.n_iter]], dtype=P.dtype)
        params = np.array([[self.early_exaggeration], [self._learning_rate],
                          [self.min_grad_norm], [self.angle]], dtype=P.dtype)
        results = np.zeros((3, 1), dtype=P.dtype)  # curIter, error, gradNorm

        if P.dtype == np.float64:
            # P.data = np.ascontiguousarray(P.data)
            # P.indices = np.ascontiguousarray(P.indices)
            # P.indptr = np.ascontiguousarray(P.indptr)
            print("running double precision")
            daal4py.daal_tsne_gradient_descent(
                X_embedded,
                P,
                size_iter,
                params,
                results,
                0)
        elif P.dtype == np.float32:
            print("running single-precision")
            daal4py.daal_tsne_gradient_descent(
                X_embedded,
                P,
                size_iter,
                params,
                results,
                1)
        else:
            raise ValueError("unsupported dtype of 'P' matrix")

        # Save the final number of iterations
        self.n_iter_ = int(results[0][0])

        # Save Kullback-Leiber divergence
        self.kl_divergence_ = results[1][0]

        return X_embedded

    def _fit(self, X, skip_num_points=0):
        """Private function to fit the model using X as training data."""
        if isinstance(self.init, str) and self.init == 'warn':
            warnings.warn("The default initialization in TSNE will change "
                          "from 'random' to 'pca' in 1.2.", FutureWarning)
            self._init = 'random'
        else:
            self._init = self.init

        if isinstance(self._init, str) and self._init == 'pca' and issparse(X):
            raise TypeError("PCA initialization is currently not suported "
                            "with the sparse input matrix. Use "
                            "init=\"random\" instead.")

        if self.method not in ['barnes_hut', 'exact']:
            raise ValueError("'method' must be 'barnes_hut' or 'exact'")
        if self.angle < 0.0 or self.angle > 1.0:
            raise ValueError("'angle' must be between 0.0 - 1.0")
        if self.learning_rate == 'warn':
            warnings.warn("The default learning rate in TSNE will change "
                          "from 200.0 to 'auto' in 1.2.", FutureWarning)
            self._learning_rate = 200.0
        else:
            self._learning_rate = self.learning_rate
        if self._learning_rate == 'auto':
            self._learning_rate = X.shape[0] / self.early_exaggeration / 4
            self._learning_rate = np.maximum(self._learning_rate, 50)
        else:
            if not (self._learning_rate > 0):
                raise ValueError("'learning_rate' must be a positive number "
                                 "or 'auto'.")

        if hasattr(self, 'square_distances'):
            if self.square_distances not in [True, 'legacy']:
                raise ValueError("'square_distances' must be True or 'legacy'.")
            if self.metric != "euclidean" and self.square_distances is not True:
                warnings.warn(("'square_distances' has been introduced in 0.24"
                               "to help phase out legacy squaring behavior. The "
                               "'legacy' setting will be removed in 0.26, and the "
                               "default setting will be changed to True. In 0.28, "
                               "'square_distances' will be removed altogether,"
                               "and distances will be squared by default. Set "
                               "'square_distances'=True to silence this warning."),
                              FutureWarning)

        if self.method == 'barnes_hut':
            if sklearn_check_version('0.23'):
                X = self._validate_data(X, accept_sparse=['csr'],
                                        ensure_min_samples=2,
                                        dtype=[np.float32, np.float64])
            else:
                X = check_array(X, accept_sparse=['csr'], ensure_min_samples=2,
                                dtype=[np.float32, np.float64])
        else:
            if sklearn_check_version('0.23'):
                X = self._validate_data(X, accept_sparse=['csr', 'csc', 'coo'],
                                        dtype=[np.float32, np.float64])
            else:
                X = check_array(X, accept_sparse=['csr', 'csc', 'coo'],
                                dtype=[np.float32, np.float64])

        if self.metric == "precomputed":
            if isinstance(self._init, str) and self._init == 'pca':
                raise ValueError("The parameter init=\"pca\" cannot be "
                                 "used with metric=\"precomputed\".")
            if X.shape[0] != X.shape[1]:
                raise ValueError("X should be a square distance matrix")

            check_non_negative(X, "TSNE.fit(). With metric='precomputed', X "
                                  "should contain positive distances.")

            if self.method == "exact" and issparse(X):
                raise TypeError(
                    'TSNE with method="exact" does not accept sparse '
                    'precomputed distance matrix. Use method="barnes_hut" '
                    'or provide the dense distance matrix.')

        if self.method == 'barnes_hut' and self.n_components > 3:
            raise ValueError("'n_components' should be inferior to 4 for the "
                             "barnes_hut algorithm as it relies on "
                             "quad-tree or oct-tree.")
        random_state = check_random_state(self.random_state)

        if self.early_exaggeration < 1.0:
            raise ValueError("early_exaggeration must be at least 1, but is {}"
                             .format(self.early_exaggeration))

        if self.n_iter < 250:
            raise ValueError("n_iter should be at least 250")

        n_samples = X.shape[0]

        neighbors_nn = None
        if self.method == "exact":
            # Retrieve the distance matrix, either using the precomputed one or
            # computing it.
            if self.metric == "precomputed":
                distances = X
            else:
                if self.verbose:
                    print("[t-SNE] Computing pairwise distances...")

                if self.metric == "euclidean":
                    # Euclidean is squared here, rather than using **= 2,
                    # because euclidean_distances already calculates
                    # squared distances, and returns np.sqrt(dist) for
                    # squared=False.
                    # Also, Euclidean is slower for n_jobs>1, so don't set here
                    distances = pairwise_distances(X, metric=self.metric,
                                                   squared=True)
                else:
                    distances = pairwise_distances(X, metric=self.metric,
                                                   n_jobs=self.n_jobs)

            if np.any(distances < 0):
                raise ValueError("All distances should be positive, the "
                                 "metric given is not correct")

            if self.metric != "euclidean" and \
                    getattr(self, 'square_distances', True) is True:
                distances **= 2

            # compute the joint probability distribution for the input space
            P = _joint_probabilities(distances, self.perplexity, self.verbose)
            assert np.all(np.isfinite(P)), "All probabilities should be finite"
            assert np.all(P >= 0), "All probabilities should be non-negative"
            assert np.all(P <= 1), ("All probabilities should be less "
                                    "or then equal to one")

        else:
            # Compute the number of nearest neighbors to find.
            # LvdM uses 3 * perplexity as the number of neighbors.
            # In the event that we have very small # of points
            # set the neighbors to n - 1.
            n_neighbors = min(n_samples - 1, int(3. * self.perplexity + 1))

            if self.verbose:
                print("[t-SNE] Computing {} nearest neighbors..."
                      .format(n_neighbors))

            # Find the nearest neighbors for every point
            knn = NearestNeighbors(
                algorithm='auto',
                n_jobs=self.n_jobs,
                n_neighbors=n_neighbors,
                metric=self.metric,
            )
            t0 = time()
            knn.fit(X)
            duration = time() - t0
            if self.verbose:
                print("[t-SNE] Indexed {} samples in {:.3f}s...".format(
                    n_samples, duration))

            t0 = time()
            distances_nn = knn.kneighbors_graph(mode='distance')
            duration = time() - t0
            if self.verbose:
                print("[t-SNE] Computed neighbors for {} samples "
                      "in {:.3f}s...".format(n_samples, duration))

            # Free the memory used by the ball_tree
            del knn

            if getattr(self, 'square_distances', True) is True or \
                    self.metric == "euclidean":
                # knn return the euclidean distance but we need it squared
                # to be consistent with the 'exact' method. Note that the
                # the method was derived using the euclidean method as in the
                # input space. Not sure of the implication of using a different
                # metric.
                distances_nn.data **= 2

            # compute the joint probability distribution for the input space
            P = _joint_probabilities_nn(distances_nn, self.perplexity,
                                        self.verbose)

        if isinstance(self._init, np.ndarray):
            X_embedded = self._init
        elif self._init == 'pca':
            pca = PCA(
                n_components=self.n_components,
                svd_solver='randomized',
                random_state=random_state,
            )
            X_embedded = pca.fit_transform(X).astype(np.float32, copy=False)
            warnings.warn("The PCA initialization in TSNE will change to "
                          "have the standard deviation of PC1 equal to 1e-4 "
                          "in 1.2. This will ensure better convergence.",
                          FutureWarning)
        elif self._init == 'random':
            # The embedding is initialized with iid samples from Gaussians with
            # standard deviation 1e-4.
            X_embedded = 1e-4 * random_state.randn(
                n_samples, self.n_components).astype(np.float32)
        else:
            raise ValueError("'init' must be 'pca', 'random', or "
                             "a numpy array")

        # Degrees of freedom of the Student's t-distribution. The suggestion
        # degrees_of_freedom = n_components - 1 comes from
        # "Learning a Parametric Embedding by Preserving Local Structure"
        # Laurens van der Maaten, 2009.
        degrees_of_freedom = max(self.n_components - 1, 1)

        daal_ready = self.method == 'barnes_hut' and self.n_components == 2 and \
            self.verbose == 0 and daal_check_version((2021, 'P', 500))

        if daal_ready:
            X_embedded = check_array(X_embedded, dtype=[np.float32, np.float64])
            return self._daal_tsne(
                P,
                n_samples,
                X_embedded=X_embedded
            )
        return self._tsne(
            P,
            degrees_of_freedom,
            n_samples,
            X_embedded=X_embedded,
            neighbors=neighbors_nn,
            skip_num_points=skip_num_points
        )
