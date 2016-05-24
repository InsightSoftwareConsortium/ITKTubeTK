TubeTK Tubes to Density Image Application
=========================================

*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](http://www.kitware.com) and licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*

#### Overview:

Given a TRE file containing a set of tubes/vessels, this module generates a
density image highlighting the pixels.


#### Command line usage:

```
USAGE:
		ConvertTubesToDensityImage  [--returnparameterfile <std::string>]
												[--processinformationaddress <std::string>]
												[--xml] [--echo] [--useSquareDistance]
												[--outputSize <std::vector<int>>]
												[--outputSpacing <std::vector<double>>]
												[--inputTemplateImage <std::string>] [--]
												[--version] [-h] <std::string> <std::string>
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

   --useSquareDistance (default: false)
     Use square distance value for the density image.

   --outputSize (optional)
     The image size of the output density, radius and tangent images.
		 User can specify a specific size or can provide a template image.

   --outputSpacing (optional)
     The image spacing of the output density, radius and tangent images.
		 User can specify a specific spacing or can provide a template image.

   --inputTemplateImage
     Template image for determining space and size constraints.

	 --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  File containing the input tubes.

   <std::string>
     (required)  Density output image that is being generated.

   <std::string>
     (required)  Radius output image that is being generated.

   <std::string>
     (required)  Tangent output image that is being generated.

   Description: Computes a density image from tubes

   Author(s): Stephen R. Aylward, Roland Kwitt, Kitware

   Acknowledgements: This work is part of the TubeTK project at Kitware
```
