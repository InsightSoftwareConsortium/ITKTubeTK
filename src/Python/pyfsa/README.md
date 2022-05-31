pyfsa
=====

Overview
--------

This package implements *Fine-structure analysis (FSA)* for (vertex-labeled)
undirected graphs (with unit edge weights).

Besides the core functionality of building general graph analysis tools upon
FSA, the package contains two graph-classification examples: 1) a **maximum
a-posteriori (MAP) classifier** that uses Gaussian mixtures to model the
class-conditional distribution(s) of FSA features and 2) a **SVM classifier**
that leverages a **Bag-of-Words (BoW)** representation. Both examples are set
up to *evaluate* the classifier(s) using cross-validation, but it should be
pretty straightforward to adapt those examples to build an actual application.

Requirements
------------

- **networkx**, available at
  [https://networkx.github.io](https://networkx.github.io)
- **scikit-learn**, available at
  [https://scikit-learn.org/stable/](https://scikit-learn.org/stable/)
- **numpy**, available at [https://www.numpy.org](https://www.numpy.org)

Fine-Structure Analysis
-----------------------

The idea of fine-structure analysis is rather simple. While many approaches in
graph analysis compute topological descriptors (e.g., avg. vertex degree,
centrality measures, etc.) from the *whole* graph and thus summarize the graph
in just 1 feature vector, FSA computes topological features for subgraphs
centered at each graph vertex.

A subgraph, centered at a specific vertex, is built from all neighbors of that
vertex reachable within a certain number of hops. This resembles the idea of
multi-scale analysis. The output of FSA for one graph is a collection of
feature vectors that describe local properties of the graph, very much like
SIFT descriptors in computer vision.

In the current implementation, varying the number of hops leads to an increase
in dimensionality of the feature vectors (i.e., column-wise concatenation), but
other strategies for combining multiple scales are obviously possible.

Once features are extracted, we have multiple options on how to build a robust,
yet discriminative representation of the graph (e.g., a BoW histogram).

Example - Classifying MUTAG graphs
----------------------------------

In this example, we show how to use the BoW discriminant classifier to classify
a popular mutagenity dataset (MUTAG, https://www.predictive-toxicology.org) with
188 graphs. The graphs represent chemical compounds that are labeled as either
having a mutagenic effect on the Gram-negative bacterium *Salmonella
typhimurium* or not (125 positive, 63 negative).

1.  **Input data:** Download the MUTAG graphs (already in ```pyfsa```
    compatible format) from MIDAS
    [here](http://midas3.kitware.com/midas/folder/1526).

       First, **unpack** the data using

       ```bash
       tar xvfz MUTAG.tar.gz
       ```

    The data for each graph is provided in two files: one file that contains
    the adjacency information (in the form of a {0,1} matrix), and one file
    containing one label for each vertex.

    Additionally, we have three auxiliary files: ```mutag.list```,
    ```mutag.labels``` and ```mutag.groups```. ```mutag.list``` contains the
    filenames of the files that contain the adjacency information.
    ```mutag.labels``` contains a numeric label for each graph that indicates
    its class membership (i.e., mutagen or not).  Finally, ```mutag.groups```
    contains, for each graph, grouping information, i.e., to which group the
    graph belongs to. In our case, this is unused and corresponds to the class
    labels. However, in cases where graphs have different attributes, this
    grouping can be used to establish group-conditional models for instance
    (example coming soon).

2.  **Execution:** For the sake of explanation, we assume that the extracted
    MUTAG data resides under ```/tmp/MUTAG```.

    The following command-line arguments to ```dcl.py``` specify that we want
    to compute FSA features for two different radii (i.e., the number of hops
    that define *reachable* neighbors) and that we want a codebook size of 32
    words for the computation of BoW histograms. Further, we perform 10-fold
    cross-validation (default testing size of 20%). The classifier is a SVM
    classifier with parameters cross-validated on the training portion of the
    current fold. The logging level is set to ```info```. We also specify that
    we want to recompute all features (for a detailed explanation of the CLI
    arguments, see ```python dcl.py --help```.

    ```bash
    python dcl.py \
      --graphList=/tmp/MUTAG/mutag.list \
      --baseDir /tmp/MUTAG \
      --labelPrefix vertexLabel \
      --classInfo /tmp/MUTAG/mutag.labels \
      --codewords 32 \
      --cvRuns 10 \
      --logLevel info \
      --writeAs /tmp/mutag \
      --groupInfo /tmp/MUTAG/mutag.groups \
      --recompute \
      --radii 1,2
    ```
