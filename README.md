TubeTK: Tubular Object Extraction, Registration, and Analysis
=============================================================

[![Build Status](https://travis-ci.org/TubeTK/TubeTK.png?branch=master)](https://travis-ci.org/TubeTK/TubeTK)

[TubeTK](http://www.tubetk.org) is an open-source toolkit for the segmentation, registration, and analysis of tubes and surfaces in images, developed by [Kitware, Inc.](http://www.kitware.com)

Tubes and surfaces, as generalized 1D and 2D manifolds in N-dimensional images, are essential components in a variety of image analysis tasks. Instances of tubular structures in images include blood vessels in magnetic resonance angiograms and b-mode ultrasound images, wires in microscopy images of integrated circuits, roads in aerial photographs, and nerves in confocal microscopy.

A guiding premise of TubeTK is that by focusing on 1D and 2D manifolds we can devise methods that are insensitive to the modality, noise, contrast, and scale of the images being analyzed and to the arrangement and deformations of the objects in them. In particular, we propose that TubeTK's manifold methods offer improved performance for many applications, compared to methods involving the analysis of independent geometric measures (e.g., edges and corners) or requiring complete shape models.

TubeTK makes extensive use of the [Insight Segmentation and Registration Toolkit](http://www.itk.org) (ITK) and the [Visualization Toolkit](http://www.vtk.org) (VTK). Select methods of TubeTK are provided as command-line applications and as extensions in [3D Slicer](http://www.slicer.org), an open-source medical imaging application.

License
-------

TubeTK is licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0) (the "License"); you may not use TubeTK except in compliance with the License. Unless required by applicable law or agreed to in writing, TubeTK is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

Acknowledgements
----------------

The development of TubeTK is supported in part by the

* [National Cancer Institute](http://www.cancer.gov‎) (NCI) of the [National Institutes of Health](http://www.nih.gov) (NIH) under award numbers R01CA138419, R01CA170665, R43CA165621, and R44CA143234;
* [National Institute of Biomedical Imaging and Bioengineering](http://www.nibib.nih.gov) (NBIB) of the National Institutes of Health (NIH) under award numbers R41EB015775, R43EB016621, and U54EB005149;
* [National Institute of Neurological Disorders and Stroke](http://www.ninds.nih.gov) (NINDS) of the National Institutes of Health (NIH) under award number R41NS081792;
* [Defense Advanced Research Projects Agency](http://www.darpa.mil) (DARPA) under the TRUST program.
