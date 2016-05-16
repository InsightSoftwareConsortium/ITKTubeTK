TubeTK Examples
===============

This directory contains a collection of IPython Notebooks, designed for
demonstrations of TubeTK **applications** and **use-cases**. All notebooks use
SimpleITK's Python wrapping (among other Python modules) in the `pylab`
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
`TubeTK_USE_PYTHON` and `TubeTK_USE_PYTHON_EXAMPLES_AS_TESTS` when configuring TubeTK.
Only then, a Python virtual environment `PythonVirtualenv` will be automatically
set up in the `Temporary` folder of the TubeTK build directory (`TubeTK-build`).
This virtual environment **already** has all the required
packages installed. To run the IPython notebooks from the `Examples`
directory, you only need to set the `TubeTK_BINARY_DIR` variable and run the
notebook, e.g., for `MergeAdjacentImages` as

```bash
cd Examples/MergeAdjacentImages
export TubeTK_BINARY_DIR=<FullPathTo>/TubeTK-build
${TubeTK_BINARY_DIR}/Temporary/PythonVirtualenv/bin/ipython notebook --pylab inline
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
~/TubeTK-virtualenv/bin/pip install ipython[zmq]
~/TubeTK-virtualenv/bin/pip install tornado
~/TubeTK-virtualenv/bin/pip install numpy
~/TubeTK-virtualenv/bin/pip install matplotlib
~/TubeTK-virtualenv/bin/pip install sh
```

**Note**: On Linux platforms you may be able to obtain many of these packages
as system packages which may suffice ( Ubuntu 12+). On Windows platforms some
of these packages should be obtained as binary downloads and installed.

In case an example needs any other Python package, installation instructions
will be given in the `README.md` file of the corresponding example directory.

### Install SimpleITK

For many common platforms, a built distribution is available as an Python egg.
This can be downloaded and installed with the following command:

```bash
~/TubeTK-virtualenv/bin/easy_install SimpleITK
```

As of this writing, SimpleITK version >=0.6.1 is required to run these notebooks. This version currently needs to be
downloaded from [Source Forge](http://sourceforge.net/projects/simpleitk/files/SimpleITK/0.6.1/Python/). For general
information about installing SimpleITK please see the [SimpleITK wiki](http://www.itk.org/Wiki/ITK/Release_4/SimpleITK/GettingStarted).


### Run the environment

Lauch an IPython notebook using
```bash
cd TubeTK/Examples/<ExampleDirectory>
~/TubeTK-virtualenv/bin/ipython notebook --pylab=inline
```
where `<ExampleDirectory>` is the example directory that contains the IPython notebook that you would
like to run.
