TubeTK Enhance Contrast Using Prior Application
===============================================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](http://www.kitware.com) and licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*


#### Overview:

Given an input image and a mask image, this application calculates an enhanced image.

#### Command line usage:

```
USAGE:

   EnhanceContrastUsingPrior  [--returnparameterfile <std::string>]
                              [--processinformationaddress <std::string>]
                              [--xml] [--echo] [-S <int>] [-i <int>] [-n
                              <int>] [-o <int>] [-b <float>] [-s <float>]
                              [--] [--version] [-h] <std::string>
                              <std::string> <std::string>


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