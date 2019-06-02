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
#ifndef __tubeTreeFilters_hxx
#define __tubeTreeFilters_hxx

#include "itkNumericTraits.h"

namespace tube
{

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
TreeFilters< VDimension >::
FillGap( typename TubeGroupType::Pointer & pTubeGroup,
char InterpolationMethod )
{
  char tubeName[] = "Tube";
  TubeListPointerType pTubeList = pTubeGroup->GetChildren(
    pTubeGroup->GetMaximumDepth(), tubeName );

  for( typename TubeGroupType::ChildrenListType::iterator itSourceTubes =
    pTubeList->begin(); itSourceTubes != pTubeList->end(); ++itSourceTubes )
    {
    TubePointerType pCurTube = dynamic_cast< TubeType * >(
      itSourceTubes->GetPointer() );
    TubeIdType curParentTubeId = pCurTube->GetParentId();
    TubePointType* parentNearestPoint = NULL;

    if( pCurTube->GetRoot() == false && curParentTubeId !=
      pTubeGroup->GetId() )
      {
      //find parent target tube
      for( typename TubeGroupType::ChildrenListType::iterator itTubes =
        pTubeList->begin(); itTubes != pTubeList->end(); ++itTubes )
        {
        TubePointerType pTube = dynamic_cast< TubeType * >(
          itTubes->GetPointer() );
        if( pTube->GetId() == curParentTubeId )
          {
          double minDistance = itk::NumericTraits<double>::max();
          int flag =-1;
          for( unsigned int index = 0; index < pTube->GetNumberOfPoints();
            ++index )
            {
            TubePointType* tubePoint = dynamic_cast< TubePointType* >(
              pTube->GetPoint( index ) );
            PositionType tubePointPosition =
              tubePoint->GetPositionInObjectSpace();
            double distance = tubePointPosition.EuclideanDistanceTo(
              pCurTube->GetPoint( 0 )->GetPositionInObjectSpace() );
            if( minDistance > distance )
              {
              minDistance = distance;
              parentNearestPoint = tubePoint;
              flag = 1;
              }
            distance = tubePointPosition.EuclideanDistanceTo(
              pCurTube->GetPoint( pCurTube->GetNumberOfPoints() - 1 )
              ->GetPositionInObjectSpace() );
            if( minDistance > distance )
              {
              minDistance = distance;
              parentNearestPoint = tubePoint;
              flag = 2;
              }
            }

          TubePointListType newTubePoints;
          if( flag == 1 )
            {
            TubePointType* childTubeStartPoint = dynamic_cast<
              TubePointType* >( pCurTube->GetPoint( 0 ) );
            InterpolatePath( parentNearestPoint,
              childTubeStartPoint, newTubePoints, InterpolationMethod );
            TubePointListType targetTubePoints = pCurTube->GetPoints();
            pCurTube->Clear();
            for( unsigned int index = 0; index < newTubePoints.size();
              ++index )
              {
              pCurTube->GetPoints().push_back( newTubePoints[ index ] );
              }
            for( unsigned int i = 0; i < targetTubePoints.size(); ++i )
              {
              pCurTube->GetPoints().push_back( targetTubePoints[ i ] );
              }
            }
          if( flag == 2 )
            {
            TubePointType* childTubeEndPoint =
              dynamic_cast< TubePointType* >
              ( pCurTube->GetPoint( pCurTube->GetNumberOfPoints() - 1 ) );
            InterpolatePath( parentNearestPoint,
              childTubeEndPoint, newTubePoints, InterpolationMethod );
            for( int index = newTubePoints.size() - 1; index >= 0; index-- )
              {
              pCurTube->GetPoints().push_back( newTubePoints[ index ] );
              }
            }
          break;
          }
        }
      }
    }
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
TreeFilters< VDimension >::
InterpolatePath(
  typename TubeType::TubePointType * parentNearestPoint,
  typename TubeType::TubePointType * itkNotUsed( childEndPoint ),
  typename TubeType::TubePointListType & newTubePoints,
  char InterpolationMethod )
{
  if( InterpolationMethod == 'S' )
    {
    newTubePoints.push_back( *parentNearestPoint );
    }
  return;

}
} // End namespace tube

#endif // End !defined( __tubeTreeFilters_hxx )
