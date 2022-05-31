TubeTK ComputeTrainingMask Application
=============================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](https://www.kitware.com) and licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).*

#### Overview:

Given an expert vessel mask volume this application defines the vessel and not-vessel volume mask.

Author(s): Sumedha Singla (Kitware), Francois Budin (Kitware)

This work is part of the TubeTK project at Kitware.

#### Algorithm:

1. **Get the input tubes:**
 The first step is to read the input image. The input image is an expert vessel mask,
 where an expert has marked the regions which are vessels.

2. **Compute the centerlines from the masked input**
 All the masked labelled regions in the input are evaluated and a thin centerline is
 generated passing through the center of those masked regions. Binary Thinning Image
 Filter is used.

 * **Threshold and Morphology Image Filter**
 A series of threshold and morphology image filters are applied to the input image.

3. **Write Output image file**
 The final vessel mask is writen as an image file.

#### Commandline Usage:

USAGE:

   ./bin/ComputeTrainingMask  [--returnparameterfile <std::string>]
                              [--processinformationaddress <std::string>]
                              [--xml] [--echo] [--notVesselWidth <double>]
                              [--gap <double>] [--notVesselMask
                              <std::string>] [--] [--version] [-h]
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

   --notVesselWidth <double>
     A double value specifying the width around the region marked as vessel by an expert, that would be regarded as not-vessel expert mask.

   --gap <double>
     A double value specifying the gap between the vessel region marked by an expert and the not-vessel region.

   --notVesselMask <std::string>
     The not-vessel expert mask is saved as an image in the location specified by the string "notVesselMask".

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input image file containing the expert marked vessel.

   <std::string>
     (required)  Output image file containing vessel mask.
