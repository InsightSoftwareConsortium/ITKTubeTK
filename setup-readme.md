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

Using TubeTK
------------

Minimal
    $ python -c "from itk import TubeTK"

Recommended
    $ python -c "from itk import TubeTK as ttk"


Acknowledgements
----------------

If you find TubeTK to be useful for your work, please cite the following publication when publishing your work:
* S. R. Aylward and E. Bullitt, "Initialization, noise, singularities, and scale in height ridge traversal for tubular object centerline extraction," Medical Imaging, IEEE Transactions on, vol. 21, no. 2, pp. 61-75, 2002.
