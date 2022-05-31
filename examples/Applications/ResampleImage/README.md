TubeTK Resample Image Application
=================================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](https://www.kitware.com) and licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).*

#### Overview:

Resample the input image to match another image or to match te required input specifications.

#### Command line usage:

```
USAGE:

   ResampleImage  [--returnparameterfile <std::string>]
                  [--processinformationaddress <std::string>] [--xml]
                  [--echo] [--loadTransform <std::string>] [--interpolator
                  <NearestNeighbor|Linear|BSpline|Sinc>] [--makeHighResIso]
                  [--makeIsotropic] [--resampleFactor
                  <std::vector<double>>] [--index <std::vector<int>>]
                  [--origin <std::vector<double>>] [--spacing
                  <std::vector<double>>] [--matchImage <std::string>] [--]
                  [--version] [-h] <std::string> <std::string>


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

   --loadTransform <std::string>
     Load the transform to be applied.

   --interpolator <NearestNeighbor|Linear|BSpline|Sinc>
     Type of interpolation to perform. (default: Linear)

   --makeHighResIso
     Make spacing isotropic - using smallest voxel size. Overrides other
     matchImage and setSpacing. (default: 0)

   --makeIsotropic
     Make spacing isotropic. Overrides other matchImage and setSpacing.
     (default: 0)

   --resampleFactor <std::vector<double>>
     Factor to increase size. 2,2,2 doubles size in each dimension. Should
     be an N-vector. Overrides matchImage and setSpacing.

   --index <std::vector<int>>
     Index to be used. Should be an N-vector. Overrides matchImage.

   --origin <std::vector<double>>
     Origin to be used. Should be an N-vector. Overrides matchImage.

   --spacing <std::vector<double>>
     Spacing to be used. Should be an N-vector. Overrides matchImage.

   --matchImage <std::string>
     The image from which origin, orientation, and spacing should be taken.

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input volume.

   <std::string>
     (required)  Resampling results.


   Description: Resamples an image to match another image (origin,
   orientation, and spacing), to be isotropic, or to have a specific
   spacing.

   Author(s): Stephen R. Aylward (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.
```
