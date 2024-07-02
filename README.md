ITKTubeTK: Tubular Object Extraction, Registration, and Analysis
================================================================

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/InsightSoftwareConsortium/ITKTubeTK/blob/master/LICENSE)

[![Build, test, package](https://github.com/InsightSoftwareConsortium/ITKTubeTK/actions/workflows/build-test-package.yml/badge.svg)](https://github.com/InsightSoftwareConsortium/ITKTubeTK/actions/workflows/build-test-package.yml)

[![Documentation Status](https://readthedocs.org/projects/tubetk/badge/?version=latest)](https://tubetk.readthedocs.io/en/latest/?badge=latest)

Available in C++ and Python for Linux, Windows, and MacOS.

Overview
--------

TubeTK is an open-source toolkit for the segmentation, registration, and analysis of tubes and surfaces in images, developed by [Kitware, Inc.](https://www.kitware.com)

Tubes and surfaces, as generalized 1D and 2D manifolds in N-dimensional images, are essential components in a variety of image analysis tasks. Instances of tubular structures in images include blood vessels in magnetic resonance angiograms and b-mode ultrasound images, wires in microscopy images of integrated circuits, roads in aerial photographs, and nerves in confocal microscopy.

A guiding premise of TubeTK is that by focusing on 1D and 2D manifolds we can devise methods that are insensitive to the modality, noise, contrast, and scale of the images being analyzed and to the arrangement and deformations of the objects in them. In particular, we propose that TubeTK's manifold methods offer improved performance for many applications, compared to methods involving the analysis of independent geometric measures (e.g., edges and corners) or requiring complete shape models.

TubeTK offers various interface layers:

* [TubeTK/src](src): This is the algorithms library.   It is the lowest level of access to the methods of TubeTK.  It is only available via C++, and it requires considerable expertise to effectively combine and call its methods to do anything useful.   Interfacing directly with these algorithms is not recommended and is not well supported. Unit-level testing is performed continuously on these methods.

* [TubeTK/include](include): This is the ITK interface to select methods in `TubeTK/src`.  This level of interface is intended for ITK users and Python scripts writers.  The methods exposed represent a level of modularization that invites experimentation, integration with other toolkits (e.g., Scikit-Learn), and development of processing pipelines that accomplish significant image analysis goals.  The interface is available as an ITK Extension and thereby available via Python using Wrapped ITK.

* [TubeTK/examples/Applications](examples/Applications): These are optional command-line interface applications.  These applications are mostly also available via the TubeTK/include interface, and thereby are available via python.  Expansion of ITK will focus on the TubeTK/include directory, and new applications will only rarely be added.  These applications are built when the cmake options BUILD_EXAMPLES is enabled.  These applications also require SlicerExecutionModel, see https://github.com/Slicer/SlicerExecutionModel.   

Installing TubeTK
=================

There are now a full set of wheels for itk-TubeTK on Pypi, and they are
SOMEWHAT compatible with ITK v5.4. Regretfully getting them to work
required changes to ITK, and those changes aren't incorporated into
ITK v5.4. We will have to wait for the ITK v6.0 release, which will
likely be a few months from now.

So, for 90+% of the functionality working (certain Spatial Object I/O
methods will fail), you can use itk-tubetk 1.4.0 from pypi

    > pip install itk-tubetk

HOWEVER, for 100% functionality working, IF you are running Windows
with Python 3.11, I have created a full set of wheels for ITK,
ITKMinimumPathExtraction, and ITKTubeTK and uploaded them to the
TEST pypi server.  You can install them into a new Python venv (one
without ITK already installed), via

   > python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ itk itk-minimalpathextraction itk-tubetk

This should then install the following packages/wheels, ON WINDOWS machines
with Python 3.11 ONLY:

```
itk                       5.4.0.20240629
itk-core                  5.4.0.20240629
itk-filtering             5.4.0.20240629
itk-io                    5.4.0.20240629
itk-minimalpathextraction 2.0.0.20240629
itk-numerics              5.4.0.20240629
itk-registration          5.4.0.20240629
itk-segmentation          5.4.0.20240629
itk-tubetk                1.4.0.20240629
```

I will also upload Ubuntu 22.04 wheels to the PyPiTEST server. If you
need other wheels, please email me at stephen@aylward.org. However, my
time (like everyone's time) is quite booked so there will be a delay.
I will also get ITK master updated asap to build with the latest TubeTK
and MinimalPathExtraction...

Thanks for your patience!

Installing Additional Examples
==============================

For additional examples of how to use TubeTK, see

*** pytubetk 

This is a new project that will evolve to be the primary
source of examples for python tubetk.  It is very much a work-in-progress.
The scripts will be helpful, but they may require updates to work with
current itk-tubetk wheels:
    https://github.com/aylward/pytubetk/tree/main/examples

Compiling TubeTK
================

We stronly reocmmend that you use the Python version of TubeTK, as
described above.  However, if you wish to compile TubeTK from scratch
(e.g., because you wish to modify it or use its C++ interface), then
use the version of TubeTK that is bundled with ITK. ITKTubeTK is
available as a official ITK Remote Module, see instructions given
earlier in this document.

Details on compiling ITK (and optionally compiling its example applications and wrapping it for python) are described next.

Within ITK
----------

If you decide to compile TubeTK instead of using its convenient Python interface (see above), then when you configure ITK (https://github.com/InsightSoftwareConsortium/ITK) using CMake (https://cmake.org/), you must set the following options

* CMAKE_BUILD_TYPE = Release
* ITK_WRAP_PYTHON = On
* Module_TubeTK = On

and then, when you build ITK, TubeTK will be automatically built as well.  Additionally, if you enable Python wrapping for ITK, that wrapping will include TubeTK.

Example Applications
--------------------

To build TubeTK's example applications, you must do the following: 

1) Build Slicer Execution Model: https://github.com/Slicer/SlicerExecutionModel
2) Set the following configuration options in CMake for ITK:
    * BUILD_EXAMPLES = On
    * SlicerExecutionModel_DIR = \<Path to your build of Slicer Execution Model\>

We then recommend adding the following paths to your user environment:

For Python
----------

Again, our recommendation is to use the freely avaible and easy-to install
Python wrapping of TubeTK that is available as described earlier in this
document.
    
However, if you are compiling your own version of ITK/TubeTK, and you
have set ITK_WRAP_PYTHON = On, then when you compile ITK, you will
generate the Python interface for ITK and TubeTK. To use itk and 
itk-tubetk within you python environment, copy the file

    ITK-build/Wrapping/Generators/Python/WrapITK.pth 

into your python's site-packages directory, e.g.,

    venv/Lib/site-packages

For more info, see https://stackoverflow.com/questions/122327/how-do-i-find-the-location-of-my-python-site-packages-directory

Then you can test your configuration:

    $ python -c "import itk"

and

    $ python -c "from itk import TubeTK"

Both of the above commands should execute and return without errors.   Otherwise, please post a detailed description (of what you've done and what error you received) on the TubeTK issue tracker: https://github.com/InsightSoftwareConsortium/ITKTubeTK/issues

Roadmap
=======

Our roadmap includes:
* Continue to add examples to the pytubetk repo
* Use TubeTK segmentations to create MONAI models for vessel segmentation.

Acknowledgements
----------------

If you find TubeTK to be useful for your work, please cite the following publication when publishing your work:
* S. R. Aylward and E. Bullitt, "Initialization, noise, singularities, and scale in height ridge traversal for tubular object centerline extraction," Medical Imaging, IEEE Transactions on, vol. 21, no. 2, pp. 61-75, 2002.

The development of TubeTK has been supported, in part, by the following grants:

* [NCI](https://www.cancer.gov/) under award numbers R01CA138419, R01CA170665, R43CA165621, and R44CA143234;
* [NIBIB](https://www.nibib.nih.gov) (NBIB) of the National Institutes of Health (NIH) under award numbers R01EB014955, R41EB015775, R43EB016621, and U54EB005149;
* [NIBIB](https://www.nibib.nih.gov) and [NIGMS](https://www.nigms.nih.gov) R01EB021396;
* [NINDS](https://www.ninds.nih.gov) R42NS086295 and R41NS081792;
* [Defense Advanced Research Projects Agency](https://www.darpa.mil) (DARPA) under the TRUST program.

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
