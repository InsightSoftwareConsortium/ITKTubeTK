TubeTK SegmentUsingOtsuThreshold Application
===============================================

### Overview:

Segments an image using the otsu thresholding algorithm and
outputs a segmentation mask.

Optionally, takes a mask as input to exclude regions in the image
from being used to compute an otsu threshold.

#### Command line usage:

```
USAGE:

   SegmentUsingOtsuThreshold  [--returnparameterfile <std::string>]
                              [--processinformationaddress <std::string>]
                              [--xml] [--echo] [--maskValue <double>]
                              [--mask <std::string>] [--] [--version] [-h]
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

   --maskValue <double>
     Only pixels with this mask value will be used to compute the otsu
     threshold (default: 0)

   --mask <std::string>
     Mask volume.

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input volume.

   <std::string>
     (required)  Output volume.


   Description: Run Otsu thresholding to automatically split the
   intensities in an image.

   Author(s): Stephen R. Aylward (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.
```

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](https://www.kitware.com) and licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).*
