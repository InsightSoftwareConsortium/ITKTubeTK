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
  m_OutputIndexToObjectFrame = 0;
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
      inputIndexToObjectTransform =
      inputTube->GetIndexToObjectTransform();

    typename TubeType::AffineGeometryFrameType::Pointer
      inputIndexToObjectFrame =
      inputTube->GetAffineGeometryFrame();

    typename TubeType::TransformType::Pointer
      inputIndexToWorldTransform =
      inputTube->GetIndexToWorldTransform();

    typename TubeType::TransformType::Pointer
      inputObjectToWorldTransform =
      inputTube->GetObjectToWorldTransform();

    typename TubeType::Pointer outputTube = TubeType::New();

    outputTube->CopyInformation( inputTube );
    outputTube->Clear();
    outputTube->GetModifiableIndexToObjectTransform()->SetIdentity();
    if( m_OutputIndexToObjectFrame.IsNotNull() )
      {
      outputTube->SetAffineGeometryFrame( m_OutputIndexToObjectFrame );
      }
    else
      {
      outputTube->SetAffineGeometryFrame( inputIndexToObjectFrame );
      }

    outputTube->ComputeObjectToWorldTransform();

    typename TubeType::TransformType::Pointer
      outputInverseIndexToWorldTransform =
      TubeType::TransformType::New();
    outputTube->GetIndexToWorldTransform()->GetInverse(
      outputInverseIndexToWorldTransform );

    typename TubeType::TransformType::Pointer
      outputInverseObjectToWorldTransform =
      TubeType::TransformType::New();
    outputTube->GetObjectToWorldTransform()->GetInverse(
      outputInverseObjectToWorldTransform );

    TubePointListType tubeList = inputTube->GetPoints();
    typename TubePointListType::const_iterator tubePointIterator =
      tubeList.begin();

    while( tubePointIterator != tubeList.end() )
      {
      inputPoint = (*tubePointIterator).GetPosition();
      inputObjectPoint = inputIndexToObjectTransform
        ->TransformPoint( inputPoint );
      worldPoint = inputIndexToWorldTransform->TransformPoint(
        inputPoint );

      transformedWorldPoint = m_Transform->TransformPoint( worldPoint );

      outputObjectPoint = outputInverseIndexToWorldTransform->
        TransformPoint( transformedWorldPoint );
      outputPoint = outputInverseIndexToWorldTransform->
        TransformPoint( transformedWorldPoint );

      VesselTubeSpatialObjectPoint<TDimension> pnt;
      pnt.SetPosition( outputPoint );

      // get both normals
      typename TubeType::CovariantVectorType n1 = tubePointIterator
        ->GetNormal1();
      typename TubeType::CovariantVectorType n2 = tubePointIterator
        ->GetNormal2();

      // only try transformation of normals if both are non-zero
      if( !n1.GetVnlVector().is_zero() && !n2.GetVnlVector().is_zero() )
        {
        n1 = inputObjectToWorldTransform->TransformCovariantVector(
          n1, inputObjectPoint );
        n2 = inputObjectToWorldTransform->TransformCovariantVector(
          n2, inputObjectPoint );
        n1 = m_Transform->TransformCovariantVector( n1, worldPoint );
        n2 = m_Transform->TransformCovariantVector( n2, worldPoint );
        n1 = outputInverseObjectToWorldTransform
          ->TransformCovariantVector( n1, transformedWorldPoint );
        n2 = outputInverseObjectToWorldTransform
          ->TransformCovariantVector( n2, transformedWorldPoint );
        n1.Normalize();
        n2.Normalize();
        pnt.SetNormal1( n1 );
        pnt.SetNormal2( n2 );
        }

      typename TubeType::VectorType tang = tubePointIterator->
        GetTangent();
      if( !tang.GetVnlVector().is_zero() )
        {
        tang = inputObjectToWorldTransform->TransformVector(
          tang, inputObjectPoint );
        tang = m_Transform->TransformVector( tang, worldPoint );
        tang = outputInverseObjectToWorldTransform->
          TransformVector( tang, transformedWorldPoint );
        tang.Normalize();
        pnt.SetTangent( tang );
        }

      typename TubeType::VectorType radi;
      for( unsigned int i=0; i<TDimension; ++i )
        {
        radi[i] = tubePointIterator->GetRadius();
        }
      radi = inputIndexToWorldTransform->TransformVector( radi,
        inputPoint );
      radi = m_Transform->TransformVector( radi, worldPoint );
      radi = outputInverseIndexToWorldTransform->
        TransformVector( radi, worldPoint );
      pnt.SetRadius( radi[0] );

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
  os << indent << "Frame: " << m_OutputIndexToObjectFrame << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeTubeToTubeTransformFilter_hxx)
