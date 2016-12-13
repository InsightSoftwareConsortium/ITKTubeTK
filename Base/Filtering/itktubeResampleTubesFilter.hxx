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

#ifndef __itktubeResampleTubesFilter_hxx
#define __itktubeResampleTubesFilter_hxx

#include "itktubeResampleTubesFilter.h"

#include "itkMath.h"

#include "tubeMacro.h"
#include "tubeTubeMath.h"
#include "tubeMatrixMath.h"

#include "itktubeSubSampleTubeTreeSpatialObjectFilter.h"

namespace itk
{
namespace tube
{

//--------------------------------------------------------------------------
template< unsigned int VDimension >
ResampleTubesFilter< VDimension >
::ResampleTubesFilter( void )
{
  m_UseInverseTransform = false;
  m_MatchImage = NULL;
  m_DisplacementField = NULL;
  m_ReadTransformList = NULL;
}

//--------------------------------------------------------------------------
template< unsigned int VDimension >
ResampleTubesFilter< VDimension >
::~ResampleTubesFilter( void )
{
}

//--------------------------------------------------------------------------
template< unsigned int VDimension >
void
ResampleTubesFilter< VDimension >
::SetDisplacementField( DisplacementFieldType * field )
{
  m_DisplacementField = field;
}

//--------------------------------------------------------------------------
template< unsigned int VDimension >
void
ResampleTubesFilter< VDimension >
::SetReadTransformList( BaseTransformListType * tList )
{
  m_ReadTransformList = tList;
}

//--------------------------------------------------------------------------
template< unsigned int VDimension >
void
ResampleTubesFilter< VDimension >
::ReadImageTransform( typename TubeGroupType::TransformType::Pointer &
  outputTransform )
{
  typename ImageType::SpacingType spacing = m_MatchImage->GetSpacing();
  typename ImageType::PointType origin = m_MatchImage->GetOrigin();
  typename ImageType::DirectionType directions =
    m_MatchImage->GetDirection();

  outputTransform = TubeGroupType::TransformType::New();
  outputTransform->SetIdentity();
  itk::Vector< double, VDimension > offset;
  for( unsigned int i = 0; i < VDimension; ++i )
    {
    offset[i] = origin[i];
    }
  outputTransform->SetScale( spacing );
  outputTransform->SetMatrix( directions );
  outputTransform->SetOffset( offset );
}

//--------------------------------------------------------------------------
template< unsigned int VDimension >
typename itk::GroupSpatialObject< VDimension >::Pointer
ResampleTubesFilter< VDimension >
::ApplyDisplacementFieldTransform( typename
  TubeGroupType::TransformType::Pointer outputTransform )
{
  const TubeGroupType * inputTubeGroup = this->GetInput();

  /** Typedefs for Displacement field tranform filter.    */
  typedef itk::tube::TubeToTubeTransformFilter<
    DisplacementFieldTransformType, VDimension >
    DisplacementFieldTransformFilterType;

  // Create new transform
  typename DisplacementFieldTransformType::Pointer transform =
    DisplacementFieldTransformType::New();
  transform->SetDisplacementField( m_DisplacementField );

  // Create the filter and apply
 typename DisplacementFieldTransformFilterType::Pointer filter =
   DisplacementFieldTransformFilterType::New();
  filter->SetInput( inputTubeGroup );
  filter->SetTransform( transform );
  filter->SetOutputIndexToObjectTransform( 
    outputTransform.GetPointer() );
  filter->GraftOutput( this->GetOutput() );
  filter->Update();

  return filter->GetOutput();
}

//--------------------------------------------------------------------------
template< unsigned int VDimension >
typename itk::GroupSpatialObject< VDimension >::Pointer
ResampleTubesFilter< VDimension >
::ApplyInputTransform( typename TubeGroupType::TransformType::Pointer
  outputTransform )
{
  const TubeGroupType * inputTubeGroup = this->GetInput();
  typename TubeGroupType::Pointer tmpTubes;
  tmpTubes = m_InputSpatialObject;

  /** Typedefs for transform read from a file    */
  typedef itk::MatrixOffsetTransformBase< double, VDimension, VDimension >
    MatrixOffsetTransformType;
  typedef itk::tube::TubeToTubeTransformFilter< MatrixOffsetTransformType,
    VDimension >
    MatrixOffsetTransformFilterType;

  BaseTransformListType::const_iterator tListIt;
  for( tListIt = m_ReadTransformList->begin();
    tListIt != m_ReadTransformList->end(); ++tListIt )
    {
    typename MatrixOffsetTransformType::Pointer transform = dynamic_cast<
      MatrixOffsetTransformType * >( ( *tListIt ).GetPointer() );
    typename MatrixOffsetTransformFilterType::Pointer filter =
      MatrixOffsetTransformFilterType::New();
    if( m_UseInverseTransform )
      {
      typename MatrixOffsetTransformType::InverseTransformBaseType::Pointer
        ivT = transform->GetInverseTransform();
      transform = ( MatrixOffsetTransformType * )ivT.GetPointer();
      }

    filter->SetInput( tmpTubes );
    filter->SetTransform( transform );
    filter->SetOutputIndexToObjectTransform( outputTransform );
    filter->GraftOutput( this->GetOutput() );
    filter->Update();
    tmpTubes = filter->GetOutput();
    }
  return tmpTubes;
}

//--------------------------------------------------------------------------
template< unsigned int VDimension >
typename itk::GroupSpatialObject< VDimension >::Pointer
ResampleTubesFilter< VDimension >
::ApplyIdentityAffineTransform( typename
  TubeGroupType::TransformType::Pointer outputTransform )
{
  const TubeGroupType * inputTubeGroup = this->GetInput();

  /** Typedefs for Affine Transform */
  typedef itk::AffineTransform< double, VDimension >
    AffineTransformType;
  typedef itk::tube::TubeToTubeTransformFilter< AffineTransformType,
    VDimension >
    AffineTransformFilterType;

  typename AffineTransformType::Pointer identityAffineTransform =
    AffineTransformType::New();
  identityAffineTransform->SetIdentity();

  typename AffineTransformFilterType::Pointer filter =
    AffineTransformFilterType::New();
  filter->SetInput( inputTubeGroup );
  filter->SetTransform( identityAffineTransform );
  filter->SetOutputIndexToObjectTransform( outputTransform );
  filter->GraftOutput( this->GetOutput() );
  filter->Update();

  return filter->GetOutput();
}

//--------------------------------------------------------------------------
template< unsigned int VDimension >
void
ResampleTubesFilter< VDimension >
::GenerateData( void )
{
  const TubeGroupType * inputTubeGroup = this->GetInput();
  typename TubeGroupType::Pointer tmpTubeGroup = NULL;

  typename TubeGroupType::TransformType::Pointer outputTransform;
  if( m_MatchImage )
    {
    this->ReadImageTransform( outputTransform );
    }
  else
    {
    char soTypeName[80];
    strcpy( soTypeName, "VesselTubeSpatialObject" );
    typename TubeSpatialObjectType::ChildrenListPointer tubeList =
      inputTubeGroup->GetChildren( inputTubeGroup->GetMaximumDepth(),
      soTypeName );
    ( *( tubeList->begin() ) )->ComputeObjectToWorldTransform();
    outputTransform = ( *( tubeList->begin() ) )->
      GetIndexToWorldTransform();
    tubeList->clear();
    delete tubeList;
    }

  if( m_DisplacementField )
    {
    tmpTubeGroup = this->ApplyDisplacementFieldTransform( outputTransform );
    }
  else if( m_ReadTransformList )
    {
    tmpTubeGroup = this->ApplyInputTransform( outputTransform );
    }
  else if( m_MatchImage )
    {
    tmpTubeGroup = this->ApplyIdentityAffineTransform( outputTransform );
    }

  if( m_SamplingFactor != 1 )
    {
    /** Typedefs for Sub samppling filter     */
    typedef itk::tube::SubSampleTubeTreeSpatialObjectFilter< TubeGroupType,
      TubeSpatialObjectType > SubSampleTubeTreeFilterType;

    typename SubSampleTubeTreeFilterType::Pointer subSampleTubeTreeFilter =
      SubSampleTubeTreeFilterType::New();
    if( tmpTubeGroup.IsNotNull() )
      {
      subSampleTubeTreeFilter->SetInput( tmpTubeGroup );
      }
    else
      {
      subSampleTubeTreeFilter->SetInput( inputTubeGroup );
      }

    subSampleTubeTreeFilter->SetSampling( m_SamplingFactor );

    subSampleTubeTreeFilter->GraftOutput( this->GetOutput() );
    try
      {
      subSampleTubeTreeFilter->Update();
      }
    catch( const std::exception &e )
      {
      std::cout << e.what();
      return;
      }
    //this->GraftOutput( subSampleTubeTreeFilter->GetOutput() );
    tmpTubeGroup = subSampleTubeTreeFilter->GetOutput();
    }
  this->GraftOutput( tmpTubeGroup );
  /*
  char soTypeName[80];
  strcpy( soTypeName, "VesselTubeSpatialObject" );
  typename TubeSpatialObjectType::ChildrenListPointer inTubeList =
    inputTubeGroup->GetChildren( inputTubeGroup->GetMaximumDepth(),
      soTypeName );
  typename TubeSpatialObjectType::ChildrenListPointer outTubeList =
    tmpTubeGroup->GetChildren( tmpTubeGroup->GetMaximumDepth(),
      soTypeName );*/
}

//--------------------------------------------------------------------------
template< unsigned int VDimension >
void
ResampleTubesFilter< VDimension >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );
}

} // End namespace tube
} // End namespace itk

#endif // End !defined( __itktubeResampleTubesFilter_hxx )
