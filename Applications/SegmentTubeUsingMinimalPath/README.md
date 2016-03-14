TubeTK TubeMinimalPathExtraction Application
=============================================

#### Overview:

Given a speed function, path information and an optimizer, this module
generates an arrival function from which we extract the minimal path.

#### Usage:

```
   TubeMinimalPathExtraction  [--processinformationaddress <std::string>]
                              [--xml] [--echo] [--optimizer <std::string>]
                              [--pathPoints <point multiple="true">]
                              [--targetPath <string>]
                              [--terminationValue <double>]
                              [--numberOfIterations <int>] [--]
                              [--stepLengthFactor <double>]
                              [--stepLengthRelax <double>]
                              [--extractRadiusFromDistanceMap <bool>]
                              [--radiusImage <std::string>]
                              [--stepRadius <double>]
                              [--maxRadius <double>]
                              [--startRadius <double>]
                              [--endPointsFile <string>]
                              [--distance <double>]
                              [--hardBoundary <bool>]
                              [--version]
                              [-h] <std::string> <std::string>


Where:

   --processinformationaddress <std::string>
     Address of a structure to store process information (progress, abort,
     etc.). (default: 0)

   --xml
     Produce xml description of command line arguments. (default: 0)

   --echo
     Echo the command line arguments. (default: 0)

   --optimizer <std::string>
     Minimize arrival function at each point of the path.
     (default: Regular_Step_Gradient_Descent)

   --pathPoints <point multiple="true">
     List of path points from strart to end point.

   --targetPath <string>
     Target path file, starting from the path points, a new path is extracted to grow and meet the target path.

   --terminationValue <double>
     Minimum difference value to reach before optimizer is terminated.
     (default: 2.0)

   --numberOfIterations <int>
     Maximum number of optimizer iterations (default: 2000)

   --stepLengthFactor <double>
     Optimizer Step Size factor (default: 0.1)

   --stepLengthRelax <double>
     Relaxation factor for the regular step gradient optimizer
     (default: 0.999)

   --extractRadiusFromDistanceMap <boolean>
     Indicates if the tubes radius should be estimated along the path. Extract tube radius along the path using input image.
     To be used when input image is a distance map

   --radiusImage <std::string>
     Original mhd image from where to estimate the tube radius

   --stepRadius <double>
     Step size for radius estimation (default: 0.05)

   --maxRadius <double>
     Maximum Radius for radius estimation.
     (default: 6)

   --startRadius <double>
     Start Radius for radius estimation.
     (default: 1)

   --endPointsFile <string>
     Save the end points of the extracted path.

   --distance <double>
     Distance between the new path and the target path.

   --hardBoundary <bool>
     Indicates wether to stop extracting new path as it approach with in the boundary of the target path.
   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input speed function (mhd file of a real-valued in the
     range [0,1])

   <std::string>
     (required)  Output TRE file of tubes centered along the path.


   Description: Extract Minimal path from an image with path information.

   Author(s): Lucas Gandel (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.

```

#### Algorithm:
[Dan Mueller](http://www.insight-journal.org/browse/publication/213)

----
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by
[Kitware, Inc.](http://www.kitware.com) and licensed under the
[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*
