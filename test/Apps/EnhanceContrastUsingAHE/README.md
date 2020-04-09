TubeTK Enhance Contrast Using AHE (Adaptive Histogram Equalization) Application
===============================================================================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](http://www.kitware.com) and licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*


#### Overview:

Power Law Adaptive Histogram Equalization.

Histogram equalization modifies the contrast in an image. The AdaptiveHistogramEqualizationImageFilter is a superset of many contrast enhancing filters. By modifying its parameters (alpha, beta, and window), the AdaptiveHistogramEqualizationImageFilter can produce an adaptively equalized histogram or a version of unsharp mask (local mean subtraction). Instead of applying a strict histogram equalization in a window about a pixel, this filter prescribes a mapping function (power law) controlled by the parameters alpha and beta.

The parameter alpha controls how much the filter acts like the classical histogram equalization method (alpha=0) to how much the filter acts like an unsharp mask (alpha=1).

The parameter beta controls how much the filter acts like an unsharp mask (beta=0) to much the filter acts like pass through (beta=1, with alpha=1).

The parameter window controls the size of the region over which local statistics are calculated.

By altering alpha, beta and window, a host of equalization and unsharp masking filters is available.

The boundary condition ignores the part of the neighborhood outside the image, and over-weights the valid part of the neighborhood.

For detail description, reference "Adaptive Image Contrast Enhancement using Generalizations of Histogram Equalization." J.Alex Stark. IEEE Transactions on Image Processing, May 2000.

#### Command line usage:

```
USAGE:

   EnhanceContrastUsingAHE  [--returnparameterfile <std::string>]
                              [--processinformationaddress <std::string>]
                              [--xml] [--echo] [-a <float>] [-b <float>]
                              [-w <int>]
                              [--] [--version] [-h] <std::string>
                              <std::string>


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

   -S <int>,  --seed <int>
     Random number seed (-1 = none). (default: -1)

   -i <int>,  --iterations <int>
     Number of iterations of optimization. (default: 100)

   -n <int>,  --maskBackgroundValue <int>
     Background label value in the mask. (default: 127)

   -o <int>,  --maskObjectValue <int>
     Object label value in the mask. (default: 255)

   -b <float>,  --backgroundScale <float>
     Scale of the background measure. (default: 20)

   -s <float>,  --objectScale <float>
     Scale of the image measure. (default: 1)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input volume.

   <std::string>
     (required)  Input mask.

   <std::string>
     (required)  Output, contrast-enhanced, image.


   Description: Given an image and a mask, process the image to maximize
   the SNR between the classes indicated by the mask.

   Author(s): Stephen R. Aylward (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware. It
   was funded in part by USC:EXPOSE.
```
