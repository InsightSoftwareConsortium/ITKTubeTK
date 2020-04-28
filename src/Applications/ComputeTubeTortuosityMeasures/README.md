TubeTK ComputeTubeTortuosityMeasures Application
=============================================

#### Overview:

Computes the tortuosity measures of all tubes present in a given TRE file and
 writes them to a CSV file.

#### USAGE:

```
   ComputeTubeTortuosityMeasures  [--returnparameterfile <std::string>]
                                  [--processinformationaddress
                                  <std::string>] [--xml] [--echo]
                                  [--histogramMax <double>] [--histogramMin
                                  <double>] [--numberOfHistogramBins <int>]
                                  [--smoothingScale <double>]
                                  [--smoothingMethod
                                  <SMOOTH_TUBE_USING_INDEX_AVERAGE
                                  |SMOOTH_TUBE_USING_INDEX_GAUSSIAN
                                  |SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN>]
                                  [--histogramMetrics] [--curvatureMetrics]
                                  [--oldMetrics] [--basicMetrics] [--]
                                  [--version] [-h] <std::string>
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

   --histogramMax <double>
     Maximum of the range of values the histogram is computed on. (default:
     1)

   --histogramMin <double>
     Minimum of the range of values the histogram is computed on. (default:
     0)

   --numberOfHistogramBins <int>
     Number of bins to use for the curvature histogram (default: 20)

   --smoothingScale <double>
     This parameter has different meanings depending on the smoothing
     method used. If smoothing method is SMOOTH_TUBE_USING_INDEX_AVERAGE
     then smoothing scale is half of the window size. If smoothing method
     is SMOOTH_TUBE_USING_INDEX_GAUSSIAN or
     SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN then smoothing scale is the
     gaussian's standard deviation. (default: 5)

   --smoothingMethod <SMOOTH_TUBE_USING_INDEX_AVERAGE
      |SMOOTH_TUBE_USING_INDEX_GAUSSIAN
      |SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN>
     Specifies the smoothing method to be applied to the vessel before
     computing the tortuosity measures. See SmoothTubeFunctionEnum in
     tubeTubeMath.h for more details. (default:
     SMOOTH_TUBE_USING_INDEX_GAUSSIAN)

   --histogramMetrics
     CurvatureHistogramMetrics (default: 0)

   --curvatureMetrics
     CurvatureScalarMetric, InflectionCount1Metric, InflectionCount2Metric,
     Percentile95Metric, Tau4Metric, TotalCurvatureMetric,
     TotalSquaredCurvatureMetric (default: 0)

   --oldMetrics
     InflectionCountMetric, InflectionPointsMetric, SumOfAnglesMetric.
     (default: 0)

   --basicMetrics
     AverageRadiusMetric, ChordLengthMetric, PathLengthMetric. (default: 0)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input TRE file containing a set of tubes

   <std::string>
     (required)  Output TRE file containing the tortuosity measures


   Description: Computes Tortuosity Measures and writes them to a CSV
   file

   Author(s): Deepak Roy Chittajallu (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.
```

----
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by
[Kitware, Inc.](http://www.kitware.com) and licensed under the
[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*