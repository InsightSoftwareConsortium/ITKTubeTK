TubeTK Tube Graph to Image Application
======================================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](https://www.kitware.com) and licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).*

#### Overview:

Given a spatial graph, this application computes an image from a graph representation and the CVT centroids.


#### Command line usage:

```
USAGE:

   ConvertSpatialGraphToImage  [--returnparameterfile <std::string>]
                               [--processinformationaddress <std::string>]
                               [--xml] [--echo] [--] [--version] [-h]
                               <std::string> <std::string> <std::string>


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

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  File with CVT centroids.

   <std::string>
     (required)  Input graph.

   <std::string>
     (required)  Output image filename (will be used as base name for more
     image output).


   Description: Computes an image from a graph representation and the CVT
   centroids.

   Author(s): Stephen R. Aylward, Roland Kwitt (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.