/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __tubeTreeMathFilters_h
#define __tubeTreeMathFilters_h

#include <itkImageFileReader.h>
#include "itkGroupSpatialObject.h"
#include "itkTubeSpatialObject.h"

namespace tube
{
template< unsigned int VDimension >
class TreeMathFilters
{
public:
  //typedefs
  typedef itk::GroupSpatialObject< VDimension >         TubeGroupType;
  typedef typename TubeGroupType::ChildrenListPointer   TubeListPointerType;
  typedef itk::TubeSpatialObject< VDimension >          TubeType;
  typedef typename TubeType::Pointer                    TubePointerType;
  typedef typename TubeType::TubePointType              TubePointType;
  typedef typename TubeType::PointType                  PositionType;
  typedef itk::IndexValueType                           TubeIdType;
  typedef typename TubeType::TubePointListType          TubePointListType;

  /** Run Fill Gap on the tube-tree. */
  static void FillGap( typename TubeGroupType::Pointer & pTubeGroup,
    char InterpolationMethod );

  static void InterpolatePath(
    typename TubeType::TubePointType * parentNearestPoint,
    typename TubeType::TubePointType * itkNotUsed( childEndPoint ),
    typename TubeType::TubePointListType & newTubePoints,
    char InterpolationMethod );

private:
  TreeMathFilters();
  ~TreeMathFilters();

}; // End class ImageFilters

} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeTreeMathFilters.hxx"
#endif

#endif // End !defined( __tubeTreeMathFilters_h )
