ITKTubeTK: Tubular Object Extraction, Registration, and Analysis
================================================================

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/KitwareMedical/ITKTubeTK/blob/master/LICENSE.md)

![](https://github.com/InsightSoftwareConsortium/ITKTubeTK/workflows/Build,%20test,%20package/badge.svg)

Linux Debug: [![Build Status](https://dev.azure.com/KitwareMedical/ITKTubeTK/_apis/build/status/KitwareMedical.ITKTubeTK.LinuxGCCDebug?branchName=master)](https://dev.azure.com/KitwareMedical/ITKTubeTK/_build/latest?definitionId=4&branchName=master)

Linux Release: [![Build Status](https://dev.azure.com/KitwareMedical/ITKTubeTK/_apis/build/status/dashboards/AzureLinuxGCC/KitwareMedical.ITKTubeTK.Linux.GCC?branchName=master)](https://dev.azure.com/KitwareMedical/ITKTubeTK/_build/latest?definitionId=1&branchName=master)

Mac Release: [![Build Status](https://dev.azure.com/KitwareMedical/ITKTubeTK/_apis/build/status/dashboards/AzureMacOSGCC/KitwareMedical.ITKTubeTK.MacOS.GCC?branchName=master)](https://dev.azure.com/KitwareMedical/ITKTubeTK/_build/latest?definitionId=2&branchName=master)

Windows Release: [![Build Status](https://dev.azure.com/KitwareMedical/ITKTubeTK/_apis/build/status/dashboards/AxureWinVS/KitwareMedical.ITKTubeTK.WinVS2017?branchName=master)](https://dev.azure.com/KitwareMedical/ITKTubeTK/_build/latest?definitionId=3&branchName=master)

Overview
--------

[TubeTK](http://www.tubetk.org) is an open-source toolkit for the segmentation, registration, and analysis of tubes and surfaces in images, developed by [Kitware, Inc.](http://www.kitware.com)

Tubes and surfaces, as generalized 1D and 2D manifolds in N-dimensional images, are essential components in a variety of image analysis tasks. Instances of tubular structures in images include blood vessels in magnetic resonance angiograms and b-mode ultrasound images, wires in microscopy images of integrated circuits, roads in aerial photographs, and nerves in confocal microscopy.

A guiding premise of TubeTK is that by focusing on 1D and 2D manifolds we can devise methods that are insensitive to the modality, noise, contrast, and scale of the images being analyzed and to the arrangement and deformations of the objects in them. In particular, we propose that TubeTK's manifold methods offer improved performance for many applications, compared to methods involving the analysis of independent geometric measures (e.g., edges and corners) or requiring complete shape models.

TubeTK offers various interface layers:

* [TubeTK/src](TubeTK/src): This is the algorithms library.   It is the lowest level of access to the methods of TubeTK.  It is only available via C++, and it requires considerable expertise to effectively combine and call its methods to do anything useful.   Interfacing directly with these algorithms is not recommended and is not well supported. Unit-level testing is performed continuously on these methods.

* [TubeTK/include](TubeTK/include): This is the ITK interface to select methods in `TubeTK/src`.  This level of interface is intended for ITK users and Python scripts writers.  The methods exposed represent a level of modularization that invites experimentation, integration with other toolkits (e.g., Scikit-Learn), and development of processing pipelines that accomplish significant image analysis goals.  The interface is available as an ITK Extension and thereby available via Python using Wrapped ITK.

* [TubeTK/src/Applications](TubeTK/src/Applications): These are optional (TO BE DEPRICATED) command-line interface applications.  These applications are NOT distributed with the pip-installation of itk-tubetk, and they will be depicated in the near future - the python interface to TubeTK should be used instead.  Furthermore, building these applications requires a sequence of compilations, complex configuration changes, and third-party tools (e.g., SlicerExecutionModel). Continuous, unit-level testingis provided for these applications.   

Installing TubeTK via ITK
-------------------------
    > pip install itk-tubetk

A version of ITK with TubeTK already compiled into it is available via pip.    This provides the python interface to TubeTK via ITK.

Compiling TubeTK via ITK
-------------------------

If you wish to compile TubeTK from scratch (because you wish to modify it), then use the version of TubeTK that is bundled with ITK.   ITKTubeTK is available as a official ITK Remote Module, starting with [ITKv5.1rc02](https://github.com/InsightSoftwareConsortium/ITK/releases/tag/v5.1rc02).   `When you configure ITK using CMake`, set the options
* CMAKE_BUILD_TYPE = Release
* ITK_WRAP_PYTHON = On
* Module_TubeTK = On
and then when you build ITK, TubeTK will be automatically built as well.  Additionally, if you enable Python wrapping for ITK, that wrapping will include TubeTK.

Roadmap
-------

Our roadmap includes:
* Adding more Jupyter Notebook examples in ITKTubeTK/examples:
+ Sliding organ registration
+ Vessel-based registration
+ Tomosynthesis simulation
+ Additional vessel extraction demonstrations involving lungs, livers, and brains imaged via MRA, CT, and ultrasound.
+ Updating example/data directory to a wider variety of cases and enable its synchronization/download during cmake configuration.

For advanced developers
-----------------------

Select methods in TubeTK require VTK, and those methods are currently not distributed as pre-built binarys.   You must compile ITK with Module_TubeTK enabled and enable VTK, as follows:

*1) SlicerExecutionModel: needed for TubeTK's applications*

After compiling ITK with TubeTK enabled (see above section "Compiling TubeTK via ITK"), you can build SlicerExecutionModel and then re-configure ITK's TubeTK to build applications using SlicerExecutionModel.

*2) VTK: needed for Sliding organ registration and select other methods*

Optionally, you may also want to build VTK.   This is used by the Sliding
Organ Registration algorithm (anisotropic diffusion regularization and
registration).   Note that VTK must be build BEFORE building ITK's TubeTK.

Using CMake, you should configure VTK as follows:
* CMAKE_BUILD_TYPE = Release
* VTK_PYTHON_VERSION = 3
* VTK_WRAP_PYTHON = True
* You may also want to turn off building Tests and Examples

After configuring CMake and generating the build files, you should compile
VTK, and then re-configure ITK and enable TubeTK_USE_VTK, and set the VTK_DIR as appropriate.  You can then recompile ITK and it will include the VTK-based TubeTK methods.

Using a ITKTubeTK from Python
-----------------------------

First, to be able to run all the python tests and examples, the following packages are required:
* numpy
* scipy
* jupyter
* IPython
* tornado
* pyzmq
* jinja2
* tables
* matplotlib

Second, you will need to add the modules of python-wrapped ITK to your python environment.
This is accomplished by copying the files that specify the paths to their python modules into your
python site-packages directory.

To find the site-packages directory on your system, follow the directions on this link:
https://stackoverflow.com/questions/122327/how-do-i-find-the-location-of-my-python-site-packages-directory

Let's assume that path is /Python/Python36/site-packages.

First you will want to copy ITK's python paths file into that directory

    $ cp ~/src/ITK-Release/Wrapping/Generators/Python/WrapITK.pth /Python/Python36/site-packages

Then you can test your configuration:

    $ python -c "import itk"

and

    $ python -c "from itk import TubeTK"

Both of the above commands should execute and return without errors.


Acknowledgements
----------------

The development of TubeTK is supported in part by the

* [NCI](http://www.cancer.gov/) under award numbers R01CA138419, R01CA170665, R43CA165621, and R44CA143234;
* [NIBIB](http://www.nibib.nih.gov) (NBIB) of the National Institutes of Health (NIH) under award numbers R01EB014955, R41EB015775, R43EB016621, and U54EB005149;
* [NIBIB](http://www.nibib.nih.gov) and [NIGMS](http://www.nigms.nih.gov) R01EB021396;
* [NINDS](http://www.ninds.nih.gov) R42NS086295 and R41NS081792;
* [Defense Advanced Research Projects Agency](http://www.darpa.mil) (DARPA) under the TRUST program.

License
-------

This software is distributed under the Apache 2.0 license. Please see
the *LICENSE* file for details.

References
----------

*Publications*

( See also [http://www.aylward.org/biosketch/publications-1] )
* D.F. Pace, S.R. Aylward, M. Niethammer, "A Locally Adaptive Regularization Based on Anisotropic Diffusion for Deformable Image Registration of Sliding Organs," Medical Imaging, IEEE Transactions on , vol.32, no.11, pp.2114,2126, Nov. 2013 doi: 10.1109/TMI.2013.2274777
* E. Bullitt, D. Zeng, B. Mortamet, A. Ghosh, S. R. Aylward, W. Lin, B. L. Marks, and K. Smith, "The effects of healthy aging on intracerebral blood vessels visualized by magnetic resonance angiography," NEUROBIOLOGY OF AGING, vol. 31, no. 2, pp. 290-300, Feb. 2010.
* E. Bullitt, M. Ewend, J. Vredenburgh, A. Friedman, W. Lin, K. Wilber, D. Zeng, S. R. Aylward, and D. Reardon, "Computerized assessment of vessel morphological changes during treatment of glioblastoma multiforme: Report of a case imaged serially by MRA over four years," NEUROIMAGE, vol. 47, pp. T143-T151, Aug. 2009.
* E. Bullitt, K. Muller, I. Jung, W. Lin, and S. Aylward, "Analyzing attributes of vessel populations," MEDICAL IMAGE ANALYSIS, vol. 9, no. 1, pp. 39-49, Feb. 2005.
* S. Aylward, J. Jomier, S. Weeks, and E. Bullitt, "Registration and analysis of vascular images," INTERNATIONAL JOURNAL OF COMPUTER VISION, vol. 55, no. 2-3, pp. 123-138, Dec. 2003.
* E. Bullitt, G. Gerig, S. Pizer, W. Lin, and S. Aylward, "Measuring tortuosity of the intracerebral vasculature from MRA images," IEEE TRANSACTIONS ON MEDICAL IMAGING, vol. 22, no. 9, pp. 1163-1171, Sep. 2003.
* S. R. Aylward and E. Bullitt, "Initialization, noise, singularities, and scale in height ridge traversal for tubular object centerline extraction," Medical Imaging, IEEE Transactions on, vol. 21, no. 2, pp. 61-75, 2002.
* S. Aylward, S. Pizer, D. Eberly, and E. Bullitt, "Intensity Ridge and Widths for Tubular Object Segmentation and Description," in MMBIA '96: Proceedings of the 1996 Workshop on Mathematical Methods in Biomedical Image Analysis (MMBIA '96), Washington, DC, USA, 1996, p. 131.
