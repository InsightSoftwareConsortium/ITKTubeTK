TubeTK Enhance Edges Using Diffusion Application
================================================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](http://www.kitware.com) and licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*

#### Overview:

Performs edge enhancing anisotropic diffusion on an input image volume.
Smooths strongly in the direction parallel to an edge, while smoothing normal to the edge is inversely proportional
to the gradient. Implemented as described in Mendrik et al., Noise reduction in computed tomography scans using 3-D
anisotropic diffusion with continuous switch, IEEE Transactions on Medical Imaging, 28(10), pp. 1585-1594, 2009.

#### Command line usage:

```
USAGE:

 EnhanceEdgesUsingDiffusion  [--returnparameterfile <std::string>]
                             [--processinformationaddress <std::string>]
                             [--xml] [--echo] [-n <int>] [-t <double>]
                             [-e <double>] [-s <double>] [--] [--version]
                             [-h] <std::string> <std::string>

Where:

   --returnparameterfile <std::string>
     Filename in which to write simple return parameters (int, float,
     int-vector, etc.) as opposed to bulk return parameters (image,
     geometry, transform, measurement, table).

   --processinformationaddress <std::string>
     Address of a structure to store process information (progress, abort,
     etc.). (default: 0)

   --xml
     Produce xml description of command line arguments (default: 0)

   --echo
     Echo the command line arguments (default: 0)

   -n
     Number of iterations used in optimization. (default: 1)

   -t
     Time step used in optimization. (default: 0.11)

   -e
     Contrast of edges, as opposed to noise. If set too low, noise will not be filtered.
		 If set too high, plate-like structures will not be preserved. (default: 30.0)

   -s
     Scale at which first derivatives are calculated when determining the structure tensor and gradient magnitude.
		 If set too low, filter will be overly sensitive to noise. If set too high, small structures will not be
		 well preserved. (default: 1.0)

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input image for determining edge enhancement.

   <std::string>
     (required)  Output edge enhanced image.


   Description: Enhances edges in an image using anisotropic diffusion method.

   Author(s): Danielle F. Pace, Andinet Enquobahrie, Hua Yang, Stephen R. Aylward (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware
```
