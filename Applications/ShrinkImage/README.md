TubeTK Shrink Image Application
=============================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](http://www.kitware.com) and licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*

#### Overview:

Project an image into fewer slices. Integration uses the 'max' function.

Author(s): Stephen R. Aylward (Kitware)

Acknowledgements: This work is part of the TubeTK project at Kitware.

#### Commandline Usage:

   ShrinkImage  [--returnparameterfile <std::string>]
                [--processinformationaddress <std::string>] [--xml]
                [--echo] [-p <std::string>] [-m] [-l] [-g] [-o
                <std::vector<int>>] [-n <std::vector<int>>] [-d
                <std::vector<int>>] [--] [--version] [-h] <std::string>
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

   -p <std::string>,  --mipPointImage <std::string>
     Save a vector image containing physical point of max pixel's original
     location

   -m,  --mean
     Generate new voxel using mean (instead of max) (default: 0)

   -l,  --log
     Compute the values in the log domain (default: 0)

   -g,  --gaussian
     Generate new voxel using a centered Gaussian weight (instead of max)
     (default: 0)

   -o <std::vector<int>>,  --overlap <std::vector<int>>
     Amount each sliding window region should overlap (in pixels)

   -n <std::vector<int>>,  --newSize <std::vector<int>>
     Target size of each dimension (size/newSize must be an integer)

   -d <std::vector<int>>,  --divideBy <std::vector<int>>
     Size of each dimension is divide by this much (must be integers)

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
