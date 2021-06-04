/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __itktubeTubePointWeightsCalculator_hxx
#define __itktubeTubePointWeightsCalculator_hxx

#include "itktubeTubePointWeightsCalculator.h"

namespace itk
{

namespace tube
{

template< unsigned int VDimension,
  class TTubeSpatialObject,
  class TPointWeightFunction,
  class TPointWeights >
TubePointWeightsCalculator< VDimension,
  TTubeSpatialObject,
  TPointWeightFunction,
  TPointWeights >
::TubePointWeightsCalculator( void )
{
}


template< unsigned int VDimension,
  class TTubeSpatialObject,
  class TPointWeightFunction,
  class TPointWeights >
void
TubePointWeightsCalculator< VDimension,
  TTubeSpatialObject,
  TPointWeightFunction,
  TPointWeights >
::Compute( void )
{
  char childName[] = "Tube";
  typename TubeTreeSpatialObjectType::ChildrenListType * tubeList =
    this->m_TubeTreeSpatialObject->GetChildren(
      this->m_TubeTreeSpatialObject->GetMaximumDepth(), childName );

  // Count the tube points.
  SizeValueType tubePoints = 0;
  typedef typename TubeTreeSpatialObjectType::ChildrenListType::iterator
    TubesIteratorType;
  for( TubesIteratorType tubeIterator = tubeList->begin();
       tubeIterator != tubeList->end();
       ++tubeIterator )
    {
    TubeSpatialObjectType * currentTube =
      dynamic_cast< TubeSpatialObjectType * >(
        ( *tubeIterator ).GetPointer() );
    if( currentTube != NULL )
      {
      tubePoints += currentTube->GetNumberOfPoints();
      }
    }
  this->m_PointWeights.SetSize( tubePoints );

  SizeValueType tubeIndex = 0;
  for( TubesIteratorType tubeIterator = tubeList->begin();
       tubeIterator != tubeList->end();
       ++tubeIterator )
    {
    TubeSpatialObjectType * currentTube =
      dynamic_cast< TubeSpatialObjectType * >(
        ( *tubeIterator ).GetPointer() );
    if( currentTube != NULL )
      {
      const typename TubeSpatialObjectType::TubePointListType &
        currentTubePoints = currentTube->GetPoints();
      typedef typename TubeSpatialObjectType::TubePointListType::const_iterator
        TubePointIteratorType;
      for( TubePointIteratorType tubePointIterator =
        currentTubePoints.begin();
        tubePointIterator != currentTubePoints.end();
        ++tubePointIterator )
        {
        this->m_PointWeights[tubeIndex]
          = this->m_PointWeightFunction->Evaluate( *tubePointIterator );
        ++tubeIndex;
        }
      }
    }
  delete tubeList;
}


template< unsigned int VDimension,
  class TTubeSpatialObject,
  class TPointWeightFunction,
  class TPointWeights >
void
TubePointWeightsCalculator< VDimension,
  TTubeSpatialObject,
  TPointWeightFunction,
  TPointWeights >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  if( ! m_TubeTreeSpatialObject.IsNull() )
    {
    os << indent << "TubeTreeSpatialObject: " << m_TubeTreeSpatialObject
      << std::endl;
    }
  else
    {
    os << indent << "TubeTreeSpatialObject: " << "( 0x0 )" << std::endl;
    }
  if( ! m_PointWeightFunction.IsNull() )
    {
    os << indent << "PointWeightFunction: "
       << m_PointWeightFunction << std::endl;
    }
  else
    {
    os << indent << "PointWeightFunction: " << "( 0x0 )" << std::endl;
    }
  os << indent << "PointWeights: " << m_PointWeights << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeTubePointWeightsCalculator_hxx )
