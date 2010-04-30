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
#ifndef __tubePdPfaScoring_h
#define __tubePdPfaScoring_h

#include <map>
#include <vector>
#include <list>
#include <string>
#include <math.h>

#include "itkOrientedImage.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

namespace tube
{

  class BoundingBox;
  class Pixel2D;
  class BoundingBoxChangeList;
  class ScoringOutput;

template< class inputPixelT, class maskPixelT, unsigned int dimensionT >
class PdPfaScoring
{
 public:

  typedef inputPixelT                                            PixelType;
  typedef maskPixelT                                             MaskPixelType;
  typedef itk::Image<PixelType,dimensionT>                       ImageType;
  typedef itk::Image<MaskPixelType,dimensionT>                   MaskType;
  typedef itk::ConnectedComponentImageFilter<MaskType,ImageType> SegmentFilter;
  typedef itk::RelabelComponentImageFilter<ImageType,ImageType>  RelabelFilter;
  typedef itk::ImageRegionConstIterator<ImageType>               ConstIteratorType;
  typedef std::map<Pixel2D,int>                                  ChangeMapType;
  typedef std::multimap<int,Pixel2D>                             LabelMMapType;
  typedef std::list<BoundingBox>                                 BBoxListType;

  void SegmentChanges( ImageType::Pointer changeImage,
                       ChangeMapType& changes,
                       LabelMMapType& changesByLabel );

  void ComputeBoundingBoxes( ImageType::Pointer im,
                             MaskType::Pointer labeledMaskImage,
                             BBoxListType& BBList );

  void ComputeChangeStatistics( MaskType::Pointer trueChangeImage, 
                                MaskType::Pointer foundChangeImage,
                                MaskType::Pointer maskImage,
                                int minPixels,
                                int& totalNumChanges,
                                int& totalNumChangesFound,
                                int& totalNumFalsePositives );

};

// UTILITY CLASSES
class BoundingBox
{
 public:
  BoundingBox() {}
  ~BoundingBox() {}

  int start_pixel_x;
  int end_pixel_x;
  int start_pixel_y;
  int end_pixel_y;

  int numPixels;
};

class Pixel2D
{
 public:
  Pixel2D()  { x = 0; y = 0; }
  ~Pixel2D() {}

  bool operator<(const Pixel2D& rhs) const {
    double norm_lhs = sqrt(x*x + y*y);
    double norm_rhs = sqrt(rhs.x*rhs.x + rhs.y*rhs.y);

    return (norm_lhs < norm_rhs);
  }

  int x;
  int y;
};



// Store the bounding box for all changes
//class BoundingBoxChangeList
class BoundingBoxChangeList
{
 public:
  BoundingBoxChangeList() { }
  ~BoundingBoxChangeList() {}

  int numChanges;

  // A vector specifying the layer on which each
  // of the changes was found
  vector<int>   change_layer;
  vector<int>   change_human_verified;

  vector<double> startGDSCoords_x_um;
  vector<double> endGDSCoords_x_um;
  vector<double> startGDSCoords_y_um;
  vector<double> endGDSCoords_y_um;

  vector<int> startPixelCoords_x;
  vector<int> endPixelCoords_x;
  vector<int> startPixelCoords_y;
  vector<int> endPixelCoords_y;

  vector<double> area_pixels;
  vector<double> area_nm;

};

}  // End namespace tube

#endif
