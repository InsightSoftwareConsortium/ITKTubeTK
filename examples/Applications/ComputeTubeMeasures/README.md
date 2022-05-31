TubeTK Sample CLI Application
=============================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](https://www.kitware.com) and licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).*

#### Overview:

This application computes the ridgeness, roundness, curvature and levelness image for a given input image.

#### Command line usage:

```
USAGE:

   ComputeTubeMeasures  [--returnparameterfile <std::string>]
                        [--processinformationaddress <std::string>] [--xml]
                        [--echo] [-g <double>] [--] [--version] [-h]
                        <std::string> <std::string> <std::string>
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

   -g <double>,  --scale <double>
     Standard deviation of the ridge. (default: 4)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input image.

   <std::string>
     (required)  Output image.

   <std::string>
     (required)  Output image.

   <std::string>
     (required)  Output image.

   <std::string>
     (required)  Output image.


   Description: Returns the ridgeness, roundness, curvature, and levelness
   measures.

   Author(s): Stephen R. Aylward (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.

```