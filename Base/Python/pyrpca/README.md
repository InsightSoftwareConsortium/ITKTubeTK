pyrpca
======

This module currently implements two recent proposals for *robust PCA*
(with different objectives):
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
authors of those articles also provide MATLAB code.

Problem Statement(s)
--------------------

tbd.

Example
-------
An illustrative example for Candes et al.'s RPCA approach is to use a
[checkerboard image](https://github.com/rkwitt/pyrpca/blob/master/examples/checkerboard.png) (which is inherently
low-rank), corrupted with random noise and try to recover the low-rank part (i.e. the checkerboard) as well as the
sparsity pattern (i.e., the noise).

The `examples` directory contains an example (`ex1.py`) that demonstrates exactly this scenario.
**Note:** The example requires [SimpleITK](http://www.simpleitk.org)'s python wrapping for image
loading and image writing (It should be easy to replace these parts with your favorite image handling
library, though). Run the code with
```bash
python ex1.py checkerboard.png 0.3 /tmp/outlierImage.png /tmp/lowRank.png
```
Two images will be written: `/tmp/outlierImage.png` (i.e., the image *with*
outliers) and `/tmp/lowRank.png` (i.e., the *low-rank* recovered part).