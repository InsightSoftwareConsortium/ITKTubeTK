TubeTK Compute Binary Image Similarity Metrics Application
==========================================================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](https://www.kitware.com) and
licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).*

#### Overview:

Given two binary volumes, this application calculates the measures to quantify the similarity between the two volumes.

#### Command line usage:

```
USAGE:

   ComputeBinaryImageSimilarityMetrics  [--returnparameterfile
                                        <std::string>]
                                        [--processinformationaddress
                                        <std::string>] [--xml] [--echo] [-o
                                        <std::string>] [--] [--version]
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

   -o <std::string>,  -- <std::string>
     Write results to a file instead of the standard output.

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input volume 1.

   <std::string>
     (required)  Input volume 2.


   Description: Return scores of how well two binary images
   match.

   Author(s): Stephen R. Aylward (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware. It
   was funded in part by USC:EXPOSE.
```