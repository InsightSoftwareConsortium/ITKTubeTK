TubeTK ConvertTubesToImage Application
=============================================

#### Overview:

Given a TRE file containing a set of tubes/vessels, this module generates a
binary image highlighting the pixels that:

* are anywhere inside a tube when the `-r` flag is specified , or
* are on the centerline of a tube if the `-r` flag is not specified

#### Command line usage:

```
USAGE:

   ConvertTubesToImage  [--returnparameterfile <std::string>]
                        [--processinformationaddress <std::string>] [--xml]
                        [--echo] [-r] [--] [--version] [-h] <std::string>
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

   -r,  --
     Fill-in the radius of the vessels, not just centerlines. (default: 0)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Template image for determining space constraints.

   <std::string>
     (required)  File containing the input tubes.

   <std::string>
     (required)  Binary output image that is being generated.


   Description: Computes a binary image from tubes

   Author(s): Stephen R. Aylward, Roland Kwitt, Kitware

   Acknowledgements: This work is part of the TubeTK project at Kitware
```
