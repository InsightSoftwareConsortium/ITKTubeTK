ITKTubeTK: Tubular Object Extraction, Registration, and Analysis
================================================================

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/InsightSoftwareConsortium/ITKTubeTK/blob/master/LICENSE)

[![Build, test, package](https://github.com/InsightSoftwareConsortium/ITKTubeTK/actions/workflows/build-test-package.yml/badge.svg)](https://github.com/InsightSoftwareConsortium/ITKTubeTK/actions/workflows/build-test-package.yml)

[![Documentation Status](https://readthedocs.org/projects/tubetk/badge/?version=latest)](https://tubetk.readthedocs.io/en/latest/?badge=latest)

Available in C++ and Python for Linux, Windows, and MacOS.

Overview
--------

TubeTK is an open-source toolkit for the segmentation, registration, and analysis of tubes and surfaces in images, developed by [Kitware, Inc.](http://www.kitware.com)

Tubes and surfaces, as generalized 1D and 2D manifolds in N-dimensional images, are essential components in a variety of image analysis tasks. Instances of tubular structures in images include blood vessels in magnetic resonance angiograms and b-mode ultrasound images, wires in microscopy images of integrated circuits, roads in aerial photographs, and nerves in confocal microscopy.

A guiding premise of TubeTK is that by focusing on 1D and 2D manifolds we can devise methods that are insensitive to the modality, noise, contrast, and scale of the images being analyzed and to the arrangement and deformations of the objects in them. In particular, we propose that TubeTK's manifold methods offer improved performance for many applications, compared to methods involving the analysis of independent geometric measures (e.g., edges and corners) or requiring complete shape models.

TubeTK offers various interface layers:

* [TubeTK/src](src): This is the algorithms library.   It is the lowest level of access to the methods of TubeTK.  It is only available via C++, and it requires considerable expertise to effectively combine and call its methods to do anything useful.   Interfacing directly with these algorithms is not recommended and is not well supported. Unit-level testing is performed continuously on these methods.

* [TubeTK/include](include): This is the ITK interface to select methods in `TubeTK/src`.  This level of interface is intended for ITK users and Python scripts writers.  The methods exposed represent a level of modularization that invites experimentation, integration with other toolkits (e.g., Scikit-Learn), and development of processing pipelines that accomplish significant image analysis goals.  The interface is available as an ITK Extension and thereby available via Python using Wrapped ITK.

* [TubeTK/examples/Applications](examples/Applications): These are optional command-line interface applications.  These applications are mostly also available via the TubeTK/include interface, and thereby are available via python.  Expansion of ITK will focus on the TubeTK/include directory, and new applications will only rarely be added.  These applications are built when the cmake options BUILD_EXAMPLES is enabled.  These applications also require SlicerExecutionModel, see https://github.com/Slicer/SlicerExecutionModel.   

Installing TubeTK
=================

We recommend using TubeTK via Python.  To do so, the installation command is

    > pip install itk-tubetk

There may also be newer, experimental versions of TubeTK available via

    > pip install --pre itk-tubetk
    
For a list of present and past releases and pre-releases, see https://pypi.org/project/itk-tubetk/

Compiling TubeTK
================

We stronly reocmmend that you use the Python version of TubeTK, as described above.  However, if you wish to compile TubeTK from scratch (e.g., because you wish to modify it or use its C++ interface), then use the version of TubeTK that is bundled with ITK.   ITKTubeTK is available as a official ITK Remote Module, starting with [ITKv5.1.2](https://github.com/InsightSoftwareConsortium/ITK/releases/tag/v5.1.2).   

Details on compiling ITK (and optionally compiling it with VTK, compiling its example applications, and wrapping it for python) are described next.

Within ITK
----------

If you decide to compile TubeTK instead of using its convenient Python interface (see above), then when you configure ITK (https://github.com/InsightSoftwareConsortium/ITK) using CMake (https://cmake.org/), you must set the following options

* CMAKE_BUILD_TYPE = Release
* ITK_WRAP_PYTHON = On
* Module_TubeTK = On

and then, when you build ITK, TubeTK will be automatically built as well.  Additionally, if you enable Python wrapping for ITK, that wrapping will include TubeTK.

With VTK
--------

TubeTK will also default to requiring VTK (https://github.com/Kitware/VTK).  If you plan on wrapping for Python, you must build VTK with static libraries.  Typically, CMake will automatically find your VTK build directory, otherwise, you must set the VTK_DIR variable during ITK configuration:

* VTK_DIR = \<Path to VTK build directory\>


If you do not wish to have TubeTK use VTK, some functionality will be lost, but this can be accomplished by setting the ITK CMake variable

* TubeTK_USE_VTK = Off

Example Applications
--------------------

To build TubeTK's example applications, we recommend compiling with VTK (see above).   Additionally, you must do the following: 

1) Build Slicer Execution Model: https://github.com/Slicer/SlicerExecutionModel
2) Set the following configuration options in CMake for ITK:
    * BUILD_EXAMPLES = On
    * SlicerExecutionModel_DIR = \<Path to your build of Slicer Execution Model\>

We then recommend adding the following paths to your user environment:

For Python
----------

Again, our recommendation is to use the freely avaible and easy-to install Python wrapping of TubeTK that is available simply by issuing the following command:

    pip install itk-tubetk
    
However, if you are compiling your own version of ITK/TubeTK, and you have set ITK_WRAP_PYTHON = On, then when you compile ITK, you will generate the Python interface for ITK and TubeTK.

To use TubeTK from Python, you will also need the following packages on your build machine:
* numpy
* scipy
* jupyter
* matplotlib

Tou will also need to add the modules of python-wrapped ITK to your python environment. This is accomplished by copying the files that specify the paths to their python modules into your python site-packages directory.  To find the site-packages directory on your system, follow the directions on this link:
https://stackoverflow.com/questions/122327/how-do-i-find-the-location-of-my-python-site-packages-directory

If that reveals that your site-packages directory is /Python/Python36/site-packages. then copy ITK's python paths file into that directory, e.g.,

    $ cp ~/src/ITK-Release/Wrapping/Generators/Python/WrapITK.pth /Python/Python36/site-packages

Then you can test your configuration:

    $ python -c "import itk"

and

    $ python -c "from itk import TubeTK"

Both of the above commands should execute and return without errors.   Otherwise, please post a detailed description (of what you've done and what error you received) on the TubeTK issue tracker: https://github.com/InsightSoftwareConsortium/ITKTubeTK/issues

Roadmap
=======

Our roadmap includes:
* Adding more Jupyter Notebook examples in ITKTubeTK/examples:
    * Sliding organ registration
    * Vessel-based registration
    * Tomosynthesis simulation
    * Additional vessel extraction demonstrations involving lungs, livers, and brains imaged via MRA, CT, and ultrasound.

Acknowledgements
----------------

If you find TubeTK to be useful for your work, please cite the following publication when publishing your work:
* S. R. Aylward and E. Bullitt, "Initialization, noise, singularities, and scale in height ridge traversal for tubular object centerline extraction," Medical Imaging, IEEE Transactions on, vol. 21, no. 2, pp. 61-75, 2002.

The development of TubeTK has been supported, in part, by the following grants:

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

( See also [Stephen R. Aylward @ Google Scholar](https://scholar.google.com/citations?user=jHEgTSwAAAAJ&hl=en) )

* D.F. Pace, S.R. Aylward, M. Niethammer, "A Locally Adaptive Regularization Based on Anisotropic Diffusion for Deformable Image Registration of Sliding Organs," Medical Imaging, IEEE Transactions on , vol.32, no.11, pp.2114,2126, Nov. 2013 doi: 10.1109/TMI.2013.2274777
* E. Bullitt, D. Zeng, B. Mortamet, A. Ghosh, S. R. Aylward, W. Lin, B. L. Marks, and K. Smith, "The effects of healthy aging on intracerebral blood vessels visualized by magnetic resonance angiography," NEUROBIOLOGY OF AGING, vol. 31, no. 2, pp. 290-300, Feb. 2010.
* E. Bullitt, M. Ewend, J. Vredenburgh, A. Friedman, W. Lin, K. Wilber, D. Zeng, S. R. Aylward, and D. Reardon, "Computerized assessment of vessel morphological changes during treatment of glioblastoma multiforme: Report of a case imaged serially by MRA over four years," NEUROIMAGE, vol. 47, pp. T143-T151, Aug. 2009.
* E. Bullitt, K. Muller, I. Jung, W. Lin, and S. Aylward, "Analyzing attributes of vessel populations," MEDICAL IMAGE ANALYSIS, vol. 9, no. 1, pp. 39-49, Feb. 2005.
* S. Aylward, J. Jomier, S. Weeks, and E. Bullitt, "Registration and analysis of vascular images," INTERNATIONAL JOURNAL OF COMPUTER VISION, vol. 55, no. 2-3, pp. 123-138, Dec. 2003.
* E. Bullitt, G. Gerig, S. Pizer, W. Lin, and S. Aylward, "Measuring tortuosity of the intracerebral vasculature from MRA images," IEEE TRANSACTIONS ON MEDICAL IMAGING, vol. 22, no. 9, pp. 1163-1171, Sep. 2003.
* S. R. Aylward and E. Bullitt, "Initialization, noise, singularities, and scale in height ridge traversal for tubular object centerline extraction," Medical Imaging, IEEE Transactions on, vol. 21, no. 2, pp. 61-75, 2002.
* S. Aylward, S. Pizer, D. Eberly, and E. Bullitt, "Intensity Ridge and Widths for Tubular Object Segmentation and Description," in MMBIA '96: Proceedings of the 1996 Workshop on Mathematical Methods in Biomedical Image Analysis (MMBIA '96), Washington, DC, USA, 1996, p. 131.
