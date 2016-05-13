TubeTK ComputeTubeFlyThroughImage Application
=============================================

#### Overview:

Given an image, a TRE file, and the ID of a tube in the TRE file this
command-line modile outputs the Fly Through Image of the tube and a
binary mask indicating the tube voxels in the fly through image

#### Command line usage:

```
USAGE:

   ComputeTubeFlyThroughImage  [--returnparameterfile <std::string>]
                               [--processinformationaddress <std::string>]
                               [--xml] [--echo] [--] [--version] [-h]
                               <std::string> <std::string> <int>
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

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input Image

   <std::string>
     (required)  Input TRE file containing the tubes

   <int>
     (required)  ID of the tube in the input TRE file for which the fly
     through image should be generated

   <std::string>
     (required)  Output fly through image that will be generated

   <std::string>
     (required)  Output tube mask indicating the tube pixels in the
     generated fly through image


   Description: Computes a Fly Through Image for a Specified Tube in TRE
   File

   Author(s): Deepak Roy Chittajallu (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.
```
----
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by
[Kitware, Inc.](http://www.kitware.com) and licensed under the
[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*
