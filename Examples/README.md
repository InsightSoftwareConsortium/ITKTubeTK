TubeTK Examples
===============

This directory contains a collection of IPython Notebooks, designed for
demonstrations of TubeTK **applications** and **use-cases**. All notebooks use
ITK's Python wrapping (among other Python modules) in the `pylab`
environment for interactive processing.

## Papers

- **Kwitt13a-MICCAI** : Code to reproduce the experiments of the MICCAI '13 paper:

```
@inproceedings{Kwitt13a,
 author      = {R.~Kwitt and D.~Pace and M.~Niethammer and S.~Aylward},
 title       = {Studying Cerebral Vasculature Using Structure Proximity and Graph Kernels},
 booktitle   = {MICCAI},
 year        = {2013}}


## Quickstart

The easiest way to get started with TubeTK examples is to enable
`TubeTK_USE_PYTHON` and `TubeTK_USE_EXAMPLES_AS_TESTS` when configuring TubeTK.
To run the IPython (Jupyter) notebooks from the `Examples`
directory, you need to set the `ITK_BUILD_DIR` and the`TubeTK_BUILD_DIR` variables.
Run the notebook, e.g., for `MergeAdjacentImages` as

```bash
cd Examples/MergeAdjacentImages
export ITK_BUILD_DIR=<FullPathTo>/ITK-build
export TubeTK_BUILD_DIR=<FullPathTo>/TubeTK-build
jupyter notebook
```

In case you don't want to follow this way, you can obviously set up a Python
virtual environment yourself as it is described in the following sections.


## Setting Up Your Own Python Environment

Under the best of circumstances (tested on OSX 10.8 and 10.7.5, RH6, Ubuntu
12) a Python virtual environment to run the TubeTK examples can be set up as
follows:

```bash
mkdir ~/TubeTK-virtualenv
sudo pip install virtualenv
virtualenv ~/TubeTK-virtualenv --system-site-packages
~/TubeTK-virtualenv/bin/pip install jupyter
~/TubeTK-virtualenv/bin/pip install IPython
~/TubeTK-virtualenv/bin/pip install tornado
~/TubeTK-virtualenv/bin/pip install numpy
~/TubeTK-virtualenv/bin/pip install scipy
~/TubeTK-virtualenv/bin/pip install matplotlib
~/TubeTK-virtualenv/bin/pip install pyzmq
~/TubeTK-virtualenv/bin/pip install jinja2
~/TubeTK-virtualenv/bin/pip install tables
~/TubeTK-virtualenv/bin/pip install PyOpengl
~/TubeTK-virtualenv/bin/pip install PySide
~/TubeTK-virtualenv/bin/pip install pyqtgraph
```

**Note**: On Linux platforms you may be able to obtain many of these packages
as system packages which may suffice ( Ubuntu 12+). On Windows platforms some
of these packages should be obtained as binary downloads and installed.

In case an example needs any other Python package, installation instructions
will be given in the `README.md` file of the corresponding example directory.

### Install ITK

Please refer to the [ITK documentation](https://blog.kitware.com/itk-python-wrapping-now-available-for-the-latest-msvc-clang-and-gcc/).

### Run the environment

Launch a Jupyter notebook using
```bash
cd TubeTK/Examples/<ExampleDirectory>
~/TubeTK-virtualenv/bin/activate
jupyter notebook --pylab=inline
```
where `<ExampleDirectory>` is the example directory that contains the IPython (Jupyter) notebook that you would
like to run.
