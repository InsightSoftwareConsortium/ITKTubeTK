"""fsa.py

This modules implements fine-structure analysis of undirected
graphs with (numeric) vertex attributes. It further contains
functionality to estimate the feature distribution using Gaussian
mixture models, or to build a Bag-of-Words representation from
a collection of feature vectors.

The idea of fine-structure analysis was recently proposed in

[1] Macindoe, O. and W. Richards, "Graph Comparison Using Fine
    Structure Analysis". In: Social Computing '10

Note: We do not implement the LBG features of [1]; Our graph
features include a subset of the features proposed in [2]

[2] Li. G. et al., "Graph Classification via Topological and Label
    Attributes". In: MLG '11

as well as some additional generic features available in networkx.
"""


__license__ = "Apache License, Version 2.0 (see TubeTK)"
__author__  = "Roland Kwitt, Kitware Inc., 2013"
__email__   = "E-Mail: roland.kwitt@kitware.com"
__status__  = "Development"


# Graph handling
import networkx as nx
from networkx.algorithms import bipartite

# Machine learning
import sklearn.mixture.gmm as gm
from sklearn.cluster import KMeans

from collections import defaultdict

# Misc.
import logging
import numpy as np
import scipy.sparse
import time
import sys
import os


attr_list = [ #Average degree
             lambda g : np.mean([e for e in g.degree().values()]),
             # Average eccentricity
             lambda g : np.mean([i for i in nx.eccentricity(g).values()]),
             # Average closeness centrality
             lambda g : np.mean([e for e in nx.closeness_centrality(g).values()]),
             # Percentage of isolated points (i.e., degree(v) = 1)
             lambda g : float(len(np.where(np.array(nx.degree(g).values())==1)[0]))/g.order(),
             # Spectral radius (i.e., largest AM eigenvalue)
             lambda g : np.abs(nx.adjacency_spectrum(g))[0],
             # Spectral trace (i.e., sum of abs. eigenvalues)
             lambda g : np.sum(np.abs(nx.adjacency_spectrum(g))),
             # Label entropy, as defined in [2]
             lambda g : label_entropy([e[1]['type'] for e in g.nodes(data=True)]),
             # Mixing coefficient of attributes
             lambda g : np.linalg.det(nx.attribute_mixing_matrix(g,'type')),
             # Avg. #vertics with eccentricity == radius (i.e., central points)
             lambda g : np.mean(float(len(nx.center(g)))/g.order()),
             # Link impurity, as defined in [2]
             lambda g : link_impurity(g),
             # Diameter := max(eccentricity)
             lambda g : nx.diameter(g),
             # Radius := min(eccentricity)
             lambda g : nx.radius(g)]


def link_impurity(g):
    """Compute link impurity of vertex-labeled graph.

    Parameters
    ----------

    g : networkx Graph
        Input graph with vertex attribute stored as 'type'.

    Returns
    -------

    impurity : float
        Link impurity, see [2]
    """

    if len(g.nodes()) == 1:
        return 0
    edges = g.edges()
    u = np.array([g.node[a]['type'] for (a,b) in edges])
    v = np.array([g.node[b]['type'] for (a,b) in edges])
    return float(len(np.nonzero(u - v)[0]))/len(edges)


def label_entropy(labels):
    """Compute entropy of label vector.

    Parameters
    ----------

    labels : numpy array, shape (L,)
        The input labels.

    Returns
    -------

    entropy : float
        Entropy of the label vector, see [2]
    """

    H = np.bincount(labels)
    p = H[np.nonzero(H)].astype(float)/np.sum(H)
    return np.abs(-np.sum(p * np.log(p)))


def graph_from_file(graph_file, label_file=None, n_skip=0):
    """Load graph from an ASCII file containing adjacency information.

    Parameters
    ----------

    graph_file : string
        Filename of the file containing all the adjaceny information. Format of
        the adjaceny matrix file is as follows:

        [Header, optional]
        0 1 1
        1 0 0
        0 1 0

        Interpretation: 3x3 adjaceny matrix, e.g., with an edge between vertices
        (0,1) and (0,2), etc.

    label_file : string
        Filename of the label information file. Here is an example:

        [Header, optional]
        5
        2
        1

        Interpretation: 3 labels, v_0 label: 5, v_1 label: 2 and  v_2 label: 1.

    n_skip : int (default: 0)
        Skip n header lines.

    Returns
    -------
    G : networkx Graph
    """

    logger = logging.getLogger()
    if not os.path.exists(graph_file):
        raise Exception("Graph file %s not found!" % graph_file)

    # Load adjacency information and ensure (0,1) weights
    adj_info = np.genfromtxt(graph_file, skip_header=n_skip)
    adj_info[np.where(adj_info >= 1)] = 1
    G = nx.Graph(adj_info)

    if not label_file is None:
        if not os.path.exists(label_file):
            raise Exception("Label file %d not found!" % label_file)
        labels = np.genfromtxt(label_file, skip_header=n_skip)
        logger.debug("Loaded labelfile %s!" % label_file)

        if len(labels) != len(G):
            raise Exception("Size mismatch for labels!")

        for idx,l in enumerate(labels):
            G.node[idx]['type'] = int(l)

    logger.debug("Built graph from %s with %d vertices." %
                 (graph_file, len(G)))
    return G


def compute_graph_features(g, radius=2, sps=None, omit_degenerate=False):
    """Compute graph feature vector(s).

    Parameters
    ----------

    g : networkx input graph with N vertices
        The input graph on which we need to compute graph features.

    radius: int (default: 2)
        Compute graph features from local neighborhoods of vertices,
        where the notion of neighborhood is defined by the number of
        hops to the neighbor, i.e., the radius. This assumes that the
        initial edges weights when computing the shortest-paths are 1.

    sps: numpy matrix, shape (N, N) (default : None)
        Matrix of shortest-path information for the graph g.

    omit_degenerate : boolean (default: False)
        Currently, degenerate cases are subgraphs with just a single
        vertex. If 'omit_degenerate' is 'True', these subgraphs are
        not considered. Otherwise, the feature vector for such a sub-
        graph is just a vector of zeros.

    Returns
    -------

    v_mat : numpy matrix, shape (N, D)
        A D-dimensional feature matrix with one feature vector for
        each vertex. Features are computed for the given radius.
    """

    logger = logging.getLogger()

    # Recompute shortest paths if neccessary
    if sps is None:
        sps = nx.floyd_warshall_numpy(g)

    # Feature matrix representation of graph
    v_mat = np.zeros([len(g),len(attr_list)])

    # Iterate over all nodes
    degenerates = []
    for n in g.nodes():
        # Get n-th row of shortest path matrix
        nth_row = np.array(sps[n,:]).ravel()
        # Find elements within a certain radius
        within_radius = np.where(nth_row <= radius)
        # Build a subgraph from those nodes
        sg = g.subgraph(within_radius[0])
        # Single vertex sg is considered degenerate
        if len(sg.nodes()) == 1:
            # Keep track of degenerates
            degenerates.append(n)
            if omit_degenerate:
                continue
            # Feature vector is 0-vector
            v = np.zeros((len(attr_list),))
        else:
            v = [attr_fun(sg) for attr_fun in attr_list]
        v_mat[n,:] = np.asarray(v)

    logger.info("Found %d generate cases!" % len(degenerates))
    if len(degenerates):
        logger.info("Pruning %d degenerate cases ..." % len(degenerates))
        v_mat = np.delete(v_mat, degenerates, axis=0)
    logger.debug("Computed (%d x %d) feature matrix." %
                 (v_mat.shape[0], v_mat.shape[1]))
    return v_mat


def run_fsa(data, radii=None, recompute=True, out=None, skip=0,
        omit_degenerate=False):
    """Run (f)ine-(s)tructure (a)nalysis.

    Paramters
    ---------

    data : list of N 3-tuple of (graph files,label files, class
        indices). We iterate over this list and compute fine-structure
        topological features for each graph.

    radii : list of 'int'
        The desired neighborhood radii.

    recompute: bool (default : True)
        Recompote features, otherwise try to load them from disk.
        In case we try to load from disk, filenames are constructed
        based on the value of the 'out' parameter.

    out : string (default : None)
        Base file name for the generated data files, e.g.,
        '/tmp/data'. Two files will be written to disk:

        /tmp/data.mat
        /tmp/data.idx

        where 'data.mat' contains the feature matrix, i.e., one
        feature vector per vertex; 'data.idx' contains the indices
        that identify which graph each feature vector belongs to;

    skip : int (default : 0)
        Skip N header entries when loading graphs.

    omit_degenerate : boolean (default: False)
        Currently, degenerate cases are subgraphs with just a single
        vertex. If 'omit_degenerate' is 'True', these subgraphs are
        not considered. Otherwise, the feature vector for such a sub-
        graph is just a vector of zeros.


    Returns
    -------
        X : numpy matrix, shape (#vertices, len(radii)*D)
            Feature matrix, where D is the total number of
            features that are computed for one radius setting.

        L : numpy array, shape (#total vertices,)
            Identifies to which graph a feature vector belongs
            to.
    """

    logger = logging.getLogger()

    if radii is None:
        raise Exception("No radii given!")

    if not out is None:
        mat_file = "%s.mat" % out
        idx_file = "%s.idx" % out
        if not recompute:
            if (os.path.exists(mat_file) and
                os.path.exists(idx_file)):
                logger.info("Loading data from file(s).")
                data_mat = np.genfromtxt(mat_file)
                data_idx = np.genfromtxt(idx_file)
                return {'data_mat' : data_mat,
                        'data_idx' : data_idx}

    data_mat = []
    data_idx = []
    for idx, (cf, lf, lab) in enumerate(data):
        logger.info("Processing %d-th graph ..." % idx)

        T, x = graph_from_file(cf, lf, skip), []
        for r in radii:
            x.append(compute_graph_features(T, r, None, omit_degenerate))
        xs = np.hstack(tuple(x))
        data_mat.append(xs)
        data_idx.append(np.ones((xs.shape[0], 1))*idx)

    data_mat = np.vstack(tuple(data_mat))
    data_idx = np.vstack(tuple(data_idx))

    if not out is None:
        np.savetxt(mat_file, data_mat, delimiter=' ')
        np.savetxt(idx_file, data_idx, delimiter=' ',fmt="%d")

    return {'data_mat' : data_mat,
            'data_idx' : data_idx}


def estimate_gm(X,components=3,seed=None):
    """Estimate a Gaussian mixture model.

    Note: Uses diagonal covariance matrices.

    Parameters
    ----------

    X : numpy matrix, shape (N,D)
        Matrix of data samples (i-th row is i-th sample vector).

    c : int (default : 3)
        Number of desired mixture components.

    seed : int (default : None)
        Seed for the random number generator.

    Returns
    -------

    gm_obj : sklearn.mixture.gmm object
        Estimated GMM.
    """

    logger = logging.getLogger()

    n, d = X.shape
    logger.info("Estimating %d-comp. GMM from (%d x %d) ..." %
                (components, n, d))

    gm_obj = gm.GMM (n_components=components,
                     covariance_type='diag',
                     random_state=seed)

    gm_obj.fit(X)
    return gm_obj


def learn_codebook(X, codebook_size=200, seed=None):
    """Learn a codebook.

    Run K-Means clustering to compute a codebook. K-Means
    is initialized by K-Means++, uses a max. of 500 iter-
    ations and 10 times re-initialization.

    Paramters
    ---------

    X : numpy matrix, shape (N,D)
        Input data.

    codebook_size : int (default : 200)
        Desired number of codewords.

    seed : int (default : None)
        Seed for random number generator.

    Returns
    -------

    cb : sklearn.cluster.KMeans object
        KMeans object after fitting.
    """

    logger = logging.getLogger()
    logger.info("Learning codebook with %d words ..." % codebook_size)

    # Run vector-quantization
    cb = KMeans(codebook_size,
                init="k-means++",
                n_init=10,
                max_iter=500,
                random_state=seed)
    cb.fit(X)
    return cb


def bow(X, cb):
    """Compute a (normalized) BoW histogram.

    Parameters
    ----------

    X : numpy matrix, shape (N, D)
        Input data.

    cb : sklearn.cluster.KMeans
        Already estimated codebook with C codewords.

    Returns
    -------

    H : numpy array, shape (C,)
        Normalized (l2-norm) BoW histogram.
    """

    # Get nr. codewords
    n,d = cb.cluster_centers_.shape
    if d != X.shape[1]:
        raise Exception("Dimensionality mismatch!")
    # Compute closest cluster centers
    assignments = cb.predict(X)
    # Compute (normalized) BoW histogram
    B = range(0,n+1)
    return np.histogram(assignments,bins=B,density=True)[0]


def pp_gmm(X, models, argmax=True):
    """Compute the posterior probability of X under a set of GMM models.

    Parameters
    ----------
    X : numpy matrix, shape (N,D)
        Data samples.

    models : list of sklearn.mixture.gmm objects
        List of C estimated GMMs.

    argmax : boolean (default : True)
        If 'True', the index of the class (represented by
        it's model) with the highest a-posteriori probability
        is computed. If 'False', the a-posteriori probability
        if each class (represented by the model) is computed for
        each row in X. Note: We assume equal prior probabilities
        for each class.

    Returns
    -------

    maxp : numpy.int64, or np.array with shape (N, C)
        Depending on whether 'argmax' is 'True' or
        'False', the index of the class with the highest
        a-posteriori probability is returned, or the
        a-posteriori probabilities under each model (for
        each feature vector in X).
    """

    n,d = X.shape
    n_models = len(models)

    ll = np.zeros((n,n_models),dtype="float32")
    for i, model in enumerate(models):
        ll[:,i] = np.asarray(model.score(X)).ravel()

    if argmax:
        # Column-wise sum
        sump = np.sum(ll,axis=0)
        # LogSumExp to compute MAP
        t0 = np.max(sump)
        t1 = np.exp(sump - (np.log(np.sum(np.exp(sump - t0))) + t0))
        max_idx = np.argmax(t1)
        return max_idx
    else:
        # LogSumExp to compute row-wise MAP
        t0 = np.asmatrix(np.max(ll,axis=1)).transpose()
        t1 = np.log(np.sum(np.exp(ll - np.tile(t0,(1,n_models))),axis=1)) + t0
        prob = np.exp(np.asmatrix(ll) - t1)
        return prob
