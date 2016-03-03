/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeImageSegmentationMath_h
#define __itktubeImageSegmentationMath_h

#include "itktubeImageMath.h"

namespace itk
{

namespace tube
{

template< unsigned int VDimension >
class ImageSegmentationMath : public ImageMath<VDimension>
{
public:

  typedef typename ImageMath<VDimension>::PixelType    PixelType;
  typedef typename ImageMath<VDimension>::ImageType    ImageType;

  /** Compute ridgness/vesselness for specified scales. */
  static void EnhanceVessels(
      typename ImageType::Pointer imIn,
      double scaleMin, double scaleMax, double numScales );

  /** Segment using (inclusive) threshold connected components. */
  static void SegmentUsingConnectedThreshold(
      typename ImageType::Pointer & imIn,
      float threshLow, float threshHigh, float labelValue,
      float x, float y, float z );

  /** Run centroid voronoi tessellation on the image. */
  static bool ComputeVoronoiTessellation(
      typename ImageType::Pointer & imIn,
      unsigned int numberOfCentroids,
      unsigned int numberOfIterations,
      unsigned int numberOfSamples,
      const std::string & centroidOutFilePath );

private:
  ImageSegmentationMath();
  ~ImageSegmentationMath();

}; // End class ImageSegmentationMath

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeImageSegmentationMath.hxx"
#endif

#endif // End !defined(__itktubeImageSegmentationMath_h)
