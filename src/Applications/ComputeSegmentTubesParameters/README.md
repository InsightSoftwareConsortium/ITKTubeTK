Compute Segment Tubes Parameters
================================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](http://www.kitware.com) and licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*

Overview:

Given input images this application create a parameter file which will be an input to segment tubes application.

Command line usage:
```
   ComputeSegmentTubesParameters  [--returnparameterfile <std::string>]
                                  [--processinformationaddress
                                  <std::string>] [--xml] [--echo] [-v
                                  <int>] [-b <int>] [--] [--version] [-h]
                                  <std::string> <std::string> <std::string
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

   -v <int>,  --maskTubeId <int>
     Value used to indicate vessel pixels in the mask. (default: 255)

   -b <int>,  --maskBackgroundId <int>
     Value used to indicate background pixels in the mask. (default: 127)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Mask indicates vessel and background labels

   <std::string>
     (required)  Scale of potential vessel at every point.

   <std::string>
     (required)  Input image.

   <std::string>
     (required)  Parameters for SegmentTubes are written to this file.


   Description: Demonstration of how to write a CLI application. Performs
   blurring.

   Author(s): Stephen R. Aylward (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.
```