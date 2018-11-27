TubeTK: Tubular Object Extraction, Registration, and Analysis
=============================================================

[![Circle CI](https://circleci.com/gh/KitwareMedical/TubeTK.svg?style=svg)](https://open.cdash.org/index.php?project=TubeTK)

[TubeTK](http://www.tubetk.org) is an open-source toolkit for the segmentation, registration, and analysis of tubes and surfaces in images, developed by [Kitware, Inc.](http://www.kitware.com)

Tubes and surfaces, as generalized 1D and 2D manifolds in N-dimensional images, are essential components in a variety of image analysis tasks. Instances of tubular structures in images include blood vessels in magnetic resonance angiograms and b-mode ultrasound images, wires in microscopy images of integrated circuits, roads in aerial photographs, and nerves in confocal microscopy.

A guiding premise of TubeTK is that by focusing on 1D and 2D manifolds we can devise methods that are insensitive to the modality, noise, contrast, and scale of the images being analyzed and to the arrangement and deformations of the objects in them. In particular, we propose that TubeTK's manifold methods offer improved performance for many applications, compared to methods involving the analysis of independent geometric measures (e.g., edges and corners) or requiring complete shape models.

TubeTK offers various interface layers:

* [TubeTK/Base][TubeTK/Base]: This is the algorithms library.   It is the lowest level of access to the methods of TubeTK.  It is only available via C++, and it requires considerable expertise to effectively combine and call its methods to do anything useful.   Interfacing directly with these algorithms is not recommended and is not well supported. Unit-level testing is performed continuously on these methods.

* [TubeTK/ITKModules][TubeTK/ITKModules]: This is the ITK interface to select methods in `TubeTK/Base`.  This level of interface is intended for ITK users and Python scripts writers.  The methods exposed represent a level of modularization that invites experimentation, integration with other toolkits (e.g., Scikit-Learn), and development of processing pipelines that accomplish significant image analysis goals.  The interface is available as an ITK Extension and thereby available via Python using Wrapped ITK.

* [TubeTK/Applications][TubeTK/Applications]: These are the command-line interface (CLI) equivalents to the methods available via `TubeTK/ITKModules`.  This is intended for bash, bat, and other system-call scripts.  The level of modularization and intended users are similar to those of `TubeTK/ITKModules`.  C++ and python-based CLIs are provided.  Continuous, unit-level testing of `TubeTK/ITKModules` is provided via these applications.

* [TubeTK/Experiments][TubeTK/Experiments]: These are Python Jupyter notebooks that combine many `TubeTK/ITKModules` into Python scripts that show how to accomplish high-level image analysis goals with TubeTK.  They are intended to be an interactive basis for exploring TubeTK.  Python and Jupyter notebooks packages must be installed on your computer to run these. These can also be (and are) run as tests to check performance (these test performance, whereas the unit-level tests focus on regression).


Python requirements
-------------------

To be able to run all the python tests and examples, the following packages are required:
* numpy
* scipy
* jupyter
* IPython
* tornado
* pyzmq
* jinja2
* tables
* matplotlib
* pyqtgraph
* PyOpengl
* PySide
* ITK

Installing most required packages can be done with the following command line:

```
pip install requirements.txt
```

To install ITK, please refer to the [ITK documentation](https://blog.kitware.com/itk-python-wrapping-now-available-for-the-latest-msvc-clang-and-gcc/) or define the environment variables `ITK_BUILD_DIR` and `TubeTK_BUILD_DIR` to point respectivally to the folder in which ITK and TubeTK are built.

License
-------

TubeTK is licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0) (the "License"). you may not use TubeTK except in compliance with the License. Unless required by applicable law or agreed to in writing, TubeTK is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

Acknowledgements
----------------

The development of TubeTK is supported in part by the

* [National Cancer Institute](http://www.cancer.gov‎) (NCI) of the [National Institutes of Health](http://www.nih.gov) (NIH) under award numbers R01CA138419, R01CA170665, R43CA165621, and R44CA143234;
* [National Institute of Biomedical Imaging and Bioengineering](http://www.nibib.nih.gov) (NBIB) of the National Institutes of Health (NIH) under award numbers R41EB015775, R43EB016621, and U54EB005149;
* [National Institute of Neurological Disorders and Stroke](http://www.ninds.nih.gov) (NINDS) of the National Institutes of Health (NIH) under award number R41NS081792;
* [Defense Advanced Research Projects Agency](http://www.darpa.mil) (DARPA) under the TRUST program.

[TubeTK/Base]: https://github.com/KitwareMedical/TubeTK/tree/master/Base
[TubeTK/ITKModules]: https://github.com/KitwareMedical/TubeTK/tree/master/ITKModules
[TubeTK/Applications]: https://github.com/KitwareMedical/TubeTK/tree/master/Applications
[TubeTK/Experiments]: https://github.com/KitwareMedical/TubeTK/tree/master/Experiments
