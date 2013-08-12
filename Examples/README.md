TubeTK Examples
===============

# SimpleITK Notebooks

This directory contains a collection of IPython Notebooks, designed for demonstrations of TubeTK examples and use-cases. All notebooks use SimpleITK 
in the `pylab` environment for interactive processing.

[SimpleITK](http://www.simpleitk.org) is an abstraction layer and wrapper around the [Insight Toolkit](http://www.itk.org) for many languages including Python.

# Getting Started

For general information about installing SimpleITK please see the [SimpleITK wiki](http://www.itk.org/Wiki/ITK/Release_4/SimpleITK/GettingStarted).


## Setting Up a Python Environment

It is recommended to setup a separate Python virtual environment to run through these notebooks.

Under the best of circumstances (tested on OSX 10.8 and 10.7.5, RH6, Ubuntu 12) this environment can be setup with the following:

    sudo pip install virtualenv
    virtualenv ~/sitkpy --no-site-packages
    ~/sitkpy/bin/pip install ipython
    ~/sitkpy/bin/pip install ipython[zmq]
    ~/sitkpy/bin/pip install tornado
    ~/sitkpy/bin/pip install numpy
    ~/sitkpy/bin/pip install matplotlib

**Note**: On Linux platforms you may be able to obtain many of these packages as system packages which may suffice ( Ubuntu 12+). On Windows platforms some of these packages should be obtained as binary downloads and installed.

### Install SimpleITK

For many common platforms, a built distribution is available as an Python egg. This can be downloading and installed with the following command:

    ~/sitkpy/bin/easy_install SimpleITK
 

As of this writing, SimpleITK version >=0.6.1 is required to run these notebooks. This version currently needs to be downloaded from [Source Forge](http://sourceforge.net/projects/simpleitk/files/SimpleITK/0.6.1/Python/)


### Run the environment
 
To launch:

```bash
cd SimpleITK-Notebooks
~/sitkpy/bin/ipython notebook --pylab=inline
```

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](http://www.kitware.com) and licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*
