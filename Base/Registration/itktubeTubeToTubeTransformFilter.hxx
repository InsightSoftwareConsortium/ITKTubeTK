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

#ifndef __itktubeTubeToTubeTransformFilter_hxx
#define __itktubeTubeToTubeTransformFilter_hxx

#include "itktubeTubeToTubeTransformFilter.h"

namespace itk
{

namespace tube
{

template< class TTransformType, unsigned int TDimension >
TubeToTubeTransformFilter< TTransformType, TDimension >
::TubeToTubeTransformFilter( void )
{
  m_Output = 0;
  m_OutputIndexToObjectTransform = 0;
  m_Transform = 0;
}

/**
 * Apply the transformation to the tube
 */
template< class TTransformType, unsigned int TDimension >
void
TubeToTubeTransformFilter< TTransformType, TDimension >
::Update( void )
{
  m_Output = GroupType::New();

  // Check if the user set any transform
  if( !m_Transform )
    {
    itkExceptionMacro( << "No transform is set." );
    }

  Point<double, TDimension> inputPoint;
  Point<double, TDimension> inputObjectPoint;
  Point<double, TDimension> worldPoint;
  Point<double, TDimension> transformedWorldPoint;
  Point<double, TDimension> outputObjectPoint;
  Point<double, TDimension> outputPoint;
  CovariantVector< double, TDimension > normal1;
  CovariantVector< double, TDimension > normal2;

  typename TubeType::ChildrenListType::iterator TubeIterator;
  typedef typename TubeType::PointListType      TubePointListType;

  const GroupType * inputGroup = this->GetInput();
  char soTypeName[80];
  strcpy( soTypeName, "VesselTubeSpatialObject" );
  typename TubeType::ChildrenListPointer inputTubeList =
    inputGroup->GetChildren( inputGroup->GetMaximumDepth(), soTypeName );
  for( TubeIterator = inputTubeList->begin();
    TubeIterator != inputTubeList->end();
    TubeIterator++ )
    {
    typename TubeType::Pointer inputTube =
      ((TubeType *)((*TubeIterator).GetPointer()));

    inputTube->ComputeObjectToWorldTransform();
    typename TubeType::TransformType::Pointer
      inputTubeIndexToObjectTransform =
      inputTube->GetIndexToObjectTransform();
    typename TubeType::TransformType::Pointer
      inputTubeIndexToWorldTransform =
      inputTube->GetIndexToWorldTransform();
    typename TubeType::TransformType::Pointer
      inputTubeObjectToWorldTransform =
      inputTube->GetObjectToWorldTransform();

    typename TubeType::Pointer outputTube = TubeType::New();

    outputTube->CopyInformation( inputTube );
    outputTube->Clear();
    outputTube->GetModifiableIndexToObjectTransform()->SetIdentity();
    if( m_OutputIndexToObjectTransform.IsNotNull() )
      {
      outputTube->GetModifiableIndexToObjectTransform()->SetCenter(
        m_OutputIndexToObjectTransform->GetCenter() );
      outputTube->GetModifiableIndexToObjectTransform()->SetMatrix(
        m_OutputIndexToObjectTransform->GetMatrix() );
      outputTube->GetModifiableIndexToObjectTransform()->SetOffset(
        m_OutputIndexToObjectTransform->GetOffset() );
      }
    else
      {
      outputTube->GetModifiableIndexToObjectTransform()->SetCenter(
        inputTubeIndexToObjectTransform->GetCenter() );
      outputTube->GetModifiableIndexToObjectTransform()->SetMatrix(
        inputTubeIndexToObjectTransform->GetMatrix() );
      outputTube->GetModifiableIndexToObjectTransform()->SetOffset(
        inputTubeIndexToObjectTransform->GetOffset() );
      }

    outputTube->ComputeObjectToWorldTransform();

    typename TubeType::TransformType::Pointer
      outputTubeInverseIndexToWorldTransform =
      TubeType::TransformType::New();
    outputTube->GetIndexToWorldTransform()->GetInverse(
      outputTubeInverseIndexToWorldTransform );

    typename TubeType::TransformType::Pointer
      outputTubeInverseObjectToWorldTransform =
      TubeType::TransformType::New();
    outputTube->GetObjectToWorldTransform()->GetInverse(
      outputTubeInverseObjectToWorldTransform );

    TubePointListType tubeList = inputTube->GetPoints();
    typename TubePointListType::const_iterator tubePointIterator =
      tubeList.begin();

    while( tubePointIterator != tubeList.end() )
      {
      inputPoint = (*tubePointIterator).GetPosition();
      inputObjectPoint = inputTubeIndexToObjectTransform
        ->TransformPoint( inputPoint );
      worldPoint = inputTubeIndexToWorldTransform->TransformPoint(
        inputPoint );

      transformedWorldPoint = m_Transform->TransformPoint( worldPoint );

      outputObjectPoint = outputTubeInverseIndexToWorldTransform->
        TransformPoint( transformedWorldPoint );
      outputPoint = outputTubeInverseIndexToWorldTransform->
        TransformPoint( transformedWorldPoint );

      VesselTubeSpatialObjectPoint<TDimension> pnt;
      pnt.SetPosition( outputPoint );

      // get both normals
      typename TubeType::CovariantVectorType n1 = tubePointIterator
        ->GetNormal1();
      typename TubeType::CovariantVectorType n2 = tubePointIterator
        ->GetNormal2();

      // only try transformation of normals if both are non-zero
      if ( !n1.GetVnlVector().is_zero() && !n2.GetVnlVector().is_zero() )
        {
        n1 = inputTubeObjectToWorldTransform->TransformCovariantVector(
          n1, inputObjectPoint );
        n2 = inputTubeObjectToWorldTransform->TransformCovariantVector(
          n2, inputObjectPoint );
        n1 = m_Transform->TransformCovariantVector( n1, worldPoint );
        n2 = m_Transform->TransformCovariantVector( n2, worldPoint );
        n1 = outputTubeInverseObjectToWorldTransform
          ->TransformCovariantVector( n1, transformedWorldPoint );
        n2 = outputTubeInverseObjectToWorldTransform
          ->TransformCovariantVector( n2, transformedWorldPoint );
        n1.Normalize();
        n2.Normalize();
        pnt.SetNormal1( n1 );
        pnt.SetNormal2( n2 );
        }

      for( unsigned int i=0; i<TDimension; ++i )
        {
        inputPoint[i] = tubePointIterator->GetRadius();
        }
      worldPoint = inputTubeIndexToWorldTransform->TransformPoint(
        inputPoint );
      worldPoint = m_Transform->TransformPoint( worldPoint );
      outputPoint = outputTubeInverseIndexToWorldTransform
        ->TransformPoint( worldPoint );
      double radius = 0;
      for( unsigned int i=0; i<TDimension; ++i )
        {
        radius += outputPoint[i] * outputPoint[i];
        }
      radius /= TDimension;
      radius = vcl_sqrt( radius );
      pnt.SetRadius( radius );

      pnt.SetMedialness( (*tubePointIterator).GetMedialness() );
      pnt.SetRidgeness( (*tubePointIterator).GetRidgeness() );
      pnt.SetBranchness( (*tubePointIterator).GetBranchness() );

      outputTube->GetPoints().push_back( pnt );

      ++tubePointIterator;
      }

    m_Output->AddSpatialObject( outputTube );
    }
  delete inputTubeList;
}

template< class TTransformType, unsigned int TDimension >
void
TubeToTubeTransformFilter< TTransformType,TDimension >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Transformation: " << m_Transform << std::endl;
  os << indent << "Transformation: " << m_OutputIndexToObjectTransform
    << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeTubeToTubeTransformFilter_hxx)
