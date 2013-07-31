Python Robust PCA
=================

pyrpca implements two recent proposals for *robust PCA*:
```bibtex
@article{Candes11a,
    author = {E.J.~Cand\'es and X.~Li and Y.~Ma and J.~Wright},
    title = {Robust Principal Component Analysis?},
    year = 2011,
    volume = 58,
    number = 3,
    journal = {J. ACM},
    pages = {1-37}}
```
and
```bibtex
@article{Xu12a,
    author = {H.~Xu and C.~Caramanis and S.~Sanghavi},
    title = {Robust {PCA} via Outlier Pursuit},
    journal = {IEEE Trans. Inf. Theory},
    volume = 59,
    number = 5,
    pages = {3047-3064},
    year = 2012}
```
Please cite these articles in case you use this code. Note that the original
authors of those articles also provide MATLAB code. Note that the objectives
of the two works are different. Candes et al.'s approach assumes randomly
distributed corruptions throughout the dataset, while Xu et al.'s approach
assumes that full observations (i.e., column vectors of the data matrix) and
not just single entries are corrupted.

Requirements
------------

* [**numnpy**](http://www.numpy.org/)
* [**SimpleITK**](http://www.simpleitk.org) [Optional]

Problem Statement(s)
--------------------

See references (above) for the exact problem formulations of Candes
et al. and Xu et al.

Example
-------

An illustrative example for Candes et al.'s RPCA approach is to use a
checkerboard image (provided under the `examples` directory) which is,
by definition, low-rank and corrupt that image with randomly distributed
outliers. The task is then to recover the low-rank part and thus obtain
a *clean* version of the checkerboard image (as well as the sparsity
pattern).

The `examples` directory contains an example (`ex1.py`) that demonstrates
exactly this scenario.  (**Note:** The example requires
[SimpleITK](http://www.simpleitk.org)'s python wrapping for image loading and
image writing; it should be easy to replace these parts with your favorite
image handling library, though).

Run the code with
```bash
python ex1.py checkerboard.png 0.3 /tmp/outlierImage.png /tmp/lowRank.png
```
Two images will be written: `/tmp/outlierImage.png` (i.e., the image *with*
outliers) and `/tmp/lowRank.png` (i.e., the *low-rank* recovered part).

Using the IPhton Notebook
=========================

The following instructions were tested on an Ubuntu 12.04 machine (assuming you have ```virtualenv``` installed, e.g., using ```apt-get install python-virtualenv```):

```bash
cd ~
mkdir tutorial-env
virtualenv ~/tutorial-env --no-site-packages
~/tutorial-env/bin/pip install ipython
~/tutorial-env/bin/pip install ipython[zmq]
~/tutorial-env/bin/pip install tornado
~/tutorial-env/bin/pip install numpy
~/tutorial-env/bin/pip install matplotlib
~/tutorial-env/bin/easy_install SimpleITK
```
Next, launch the IPython notebook:
```bash
~/tutorial-env/bin/ipython notebook --pylab=inline
```











