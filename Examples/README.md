TubeTK Examples
===============

This directory contains a collection of IPython Notebooks, designed for demonstrations of TubeTK 
**applications** and **use-cases**. All notebooks use SimpleITK's Python wrapping (among other 
Python modules) in the `pylab` environment for interactive processing.

[SimpleITK](http://www.simpleitk.org) is an abstraction layer and wrapper around the 
[Insight Toolkit](http://www.itk.org) for many languages including Python.

# Getting Started

For general information about installing SimpleITK please see the 
[SimpleITK wiki](http://www.itk.org/Wiki/ITK/Release_4/SimpleITK/GettingStarted).

## Setting Up a Python Environment

It is recommended to setup a separate Python virtual environment to run through these notebooks.

Under the best of circumstances (tested on OSX 10.8 and 10.7.5, RH6, Ubuntu 12) this environment can be 
setup with the following:

```bash
mkdir ~/TubeTK-virtualenv
sudo pip install virtualenv
virtualenv ~/sitkpy --no-site-packages
~/TubeTK-virtualenv/bin/pip install ipython
~/TubeTK-virtualenv/bin/pip install ipython[zmq]
~/TubeTK-virtualenv/bin/pip install tornado
~/TubeTK-virtualenv/bin/pip install numpy
~/TubeTK-virtualenv/bin/pip install matplotlib
```

**Note**: On Linux platforms you may be able to obtain many of these packages as system packages which may 
suffice ( Ubuntu 12+). On Windows platforms some of these packages should be obtained as binary downloads and 
installed.

In case an example needs any other Python package, installation instructions will be given in the `README.md` 
file of the corresponding example directory.

### Install SimpleITK

For many common platforms, a built distribution is available as an Python egg. This can be downloaded 
and installed with the following command:

```bash
~/TubeTK-virtualenv/bin/easy_install SimpleITK
```

As of this writing, SimpleITK version >=0.6.1 is required to run these notebooks. This version currently needs to be 
downloaded from [Source Forge](http://sourceforge.net/projects/simpleitk/files/SimpleITK/0.6.1/Python/)


### Run the environment
 
Lauch an IPython notebook using
```bash
cd TubeTK/Examples/<ExampleDirectory>
~/TubeTK-virtualenv/bin/ipython notebook --pylab=inline
```
where `<ExampleDirectory>` is the example directory that contains the IPython notebook that you would
like to run.
