TubeTK VesselTrainMask Application
=============================================

#### Overview:

Given an expert vessel mask volume this application defines the vessel and not-vessel volume mask.

#### Usage:

```
   VesselTrainMask  [--notVesselWidth <double>]
              [--gap <double>]
			  [--notVesselMask]
              [--xml] [--echo]
              [--][--version] [-h]
              <std::string> <std::string>


Where:

   -notVesselWidth <double>
     (required) A double value specifying the width around the vessel
	 region marked by an expert that would be regarded as not vessel expert mask.

   -gap <double>
     (required) A double value specifying the gap between the vessel
	 region marked by expert and the not -vessel region.

   --xml
     Produce xml description of command line arguments (default: 0)

   --echo
     Echo the command line arguments (default: 0)

   -notVesselMask
     The not-vessel mask is saved as an image in the location specified by the
	 string "notVesselMask"

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


   Description: Return image file containing the vessel mask.

   Author(s): Sumedha Singla (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.

```

#### Algorithm:

1. **Get the input tubes:**
 We the input image. The input image is an expert vessel mask,
 where the an expert have marked the regions which are vessels.

2. **Computer the centerlines from the masked input**
 All the masked labelled regions in the input are evaluated and a thin centerline is
 generated passing through the center of those masked regions. Binary Thinning Image
 Filter is used.

 * **Threshold and Morphology Image Filter**
 A series of threshold and morphology image filters are applied to the input image.

3. **Write Output image file**
 the final vessel mask is writen as an image file.

----
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by
[Kitware, Inc.](http://www.kitware.com) and licensed under the
[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*