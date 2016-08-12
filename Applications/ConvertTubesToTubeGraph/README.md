TubeTK Tubes to Tube Graph Application
=====================================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](http://www.kitware.com) and licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*

#### Overview:

Given a TRE file containing a set of tubes/vessels, this module generates a
tube graph in terms of adjacency matrix using the information from Central Voronoi Tesselation input image.

#### Command line usage:

```
USAGE:

   ConvertTubesToTubeGraph  [--returnparameterfile <std::string>]
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
     (required)  Input tubes that are supposed to be converted.

   <std::string>
     (required)  Central Voronoi tesselation upon which the mapping is
     based.

   <std::string>
     (required)  Graph file that is about to be written.


   Description: Maps tubes into a graph structure using an central Voronoi
   tessellation.

   Author(s): Stephen R. Aylward, Roland Kwitt (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.
```