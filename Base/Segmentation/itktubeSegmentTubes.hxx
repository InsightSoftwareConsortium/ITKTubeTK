/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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

#ifndef __itktubeSegmentTubes_hxx
#define __itktubeSegmentTubes_hxx

#include "itktubeSegmentTubes.h"

namespace itk
{

namespace tube
{

/**
 * Constructor */
template< class TInputImage >
SegmentTubes<TInputImage>
::SegmentTubes( void )
{
  m_InputImage = NULL;
  m_RadiusInputImage = NULL;
  m_TubeExtractorFilter = TubeExtractorFilterType::New();
  m_SeedMask = NULL;
  m_ScaleMask = NULL;
  m_ExistingTubesMask = NULL;
  m_ExistingTubes = NULL;

  m_UseExistingTubes = false;
  m_Border = 5.0;
  m_TubeGroup = TubeGroupType::New();
}

/**
 * Destructor */
template< class TInputImage >
SegmentTubes<TInputImage>
::~SegmentTubes( void )
{
}

/**
 * Set seed index list */
template< class TInputImage >
void
SegmentTubes<TInputImage>
::SetSeedIndexList( std::vector< ContinuousIndexType > seedI )
{
  double scaleNorm = this->m_InputImage->GetSpacing()[0];
  for( size_t seedNum = 0; seedNum < seedI.size(); ++seedNum )
    {
    m_SeedIndexList.push_back( seedI[seedNum] );
    m_SeedRadiusList.push_back( this->m_Scale / scaleNorm  );
    }
}

/**
 * Set seed index list */
template< class TInputImage >
void
SegmentTubes<TInputImage>
::SetSeedPhysicalCoordinatesList( std::vector< PointType > seedP )
{
  double scaleNorm = this->m_InputImage->GetSpacing()[0];
  for( size_t seedNum = 0; seedNum < seedP.size(); ++seedNum )
    {
    ContinuousIndexType seedIndex;
    bool transformSuccess =
        m_InputImage->TransformPhysicalPointToContinuousIndex
        ( seedP[seedNum], seedIndex );
    if( !transformSuccess )
      {
      std::cerr<<"Could not transform point #"
      << seedNum <<" to seed index." << std::endl;
      continue;
      }
    m_SeedIndexList.push_back( seedIndex );
    m_SeedRadiusList.push_back( this->m_Scale / scaleNorm  );
    }
}

/**
 * Set seed index list */
template< class TInputImage >
void
SegmentTubes<TInputImage>
::SetSeedIndexFromFileList
( std::vector< ContinuousIndexType > seedI, std::vector< ScaleType > seedS)
{
  double scaleNorm = this->m_InputImage->GetSpacing()[0];
  for( size_t seedNum = 0; seedNum < seedI.size(); ++seedNum )
    {
    m_SeedIndexList.push_back( seedI[seedNum] );
    m_SeedRadiusList.push_back( seedS[seedNum] / scaleNorm  );
    }
}

/**
 * Set tube group */
template< class TInputImage >
void
SegmentTubes<TInputImage>
::SetTubeGroup( TubeGroupType * tubes )
{
  m_TubeGroup = tubes;
  m_UseExistingTubes = true;
}

/**
 * Update */
template< class TInputImage >
void
SegmentTubes<TInputImage>
::Update( void )
{
  if( this->m_InputImage )
    {
    this->m_TubeExtractorFilter->SetInputImage( this->m_InputImage );
    }
  else
    {
    throw( "Error: Segment Tubes requires an non-empty input image." );
    }
  if( this->m_RadiusInputImage )
    {
    this->m_TubeExtractorFilter->SetRadiusInputImage( this->m_RadiusInputImage );
    }
  // Set the radius for tube extractor
  double scaleNorm = this->m_InputImage->GetSpacing()[0];
  if( this->m_Scale / scaleNorm < 0.3 )
    {
    throw( "Error: Scale < 0.3 * voxel spacing is unsupported." );
    }
  this->m_TubeExtractorFilter->SetRadius( this->m_Scale / scaleNorm );

  if( this->m_SeedMask )
    {
    itk::ImageRegionConstIteratorWithIndex< TubeMaskImageType > iter(
          this->m_SeedMask, this->m_SeedMask->GetLargestPossibleRegion() );
    itk::ImageRegionConstIterator< ScaleImageType > iterS;
    if( this->m_ScaleMask )
      {
      iterS = itk::ImageRegionConstIterator< ScaleImageType >
        ( this->m_ScaleMask, this->m_ScaleMask->GetLargestPossibleRegion() );
      }

    int count = 0;
    while( !iter.IsAtEnd() )
      {
      if( iter.Get() )
        {
        if( ++count == this->m_SeedMaskStride )
          {
          count = 0;
          m_SeedIndexList.push_back( iter.GetIndex() );
          if( this->m_ScaleMask )
            {
            m_SeedRadiusList.push_back( iterS.Get() / scaleNorm );
            ++iterS;
            }
          else
            {
            m_SeedRadiusList.push_back(  this->m_Scale / scaleNorm );
            }
          }
        }
      ++iter;
      }
    }
  if( this->m_ExistingTubesMask )
    {
    this->m_TubeExtractorFilter->SetTubeMaskImage( this->m_ExistingTubesMask );
    }
  if( m_UseExistingTubes )
    {
    char tubeName[] = "Tube";
    typename TubeGroupType::ChildrenListType * tubeList =
      m_TubeGroup->GetChildren( 9999, tubeName );
    typename TubeGroupType::ChildrenListType::iterator iter = tubeList->begin();
    while( iter != tubeList->end() )
      {
      this->m_TubeExtractorFilter->AddTube( static_cast< TubeType * >(
            iter->GetPointer() ) );
      ++iter;
      }
    }

  if( m_ParameterFile.empty() == false )
    {
    TubeExtractorIOType teReader;
    teReader.SetTubeExtractor( this->m_TubeExtractorFilter );
    teReader.Read(  m_ParameterFile.c_str() );
    }

  this->m_TubeExtractorFilter->SetDebug( false );
  this->m_TubeExtractorFilter->GetRidgeOp()->SetDebug( false );
  this->m_TubeExtractorFilter->GetRadiusOp()->SetDebug( false );

  if( m_Border > 0 )
    {
    typename ImageType::IndexType minIndx = this->m_InputImage->
      GetLargestPossibleRegion().GetIndex();
    typename ImageType::SizeType size = this->m_InputImage->
      GetLargestPossibleRegion().GetSize();
    typename ImageType::IndexType maxIndx = minIndx + size;
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      minIndx[i] += m_Border;
      maxIndx[i] -= m_Border;
      }
    this->m_TubeExtractorFilter->SetExtractBoundMin( minIndx );
    this->m_TubeExtractorFilter->SetExtractBoundMax( maxIndx );
    }

  typename std::vector< ContinuousIndexType >::iterator seedIndexIter =
    this->m_SeedIndexList.begin();
  typename std::vector< ScaleType >::iterator seedRadiusIter =
    this->m_SeedRadiusList.begin();
  unsigned int count = 1;
  bool foundOneTube = false;
  while( seedIndexIter != this->m_SeedIndexList.end() )
    {
    this->m_TubeExtractorFilter->SetRadius( *seedRadiusIter );

    std::cout << "Extracting from index point " << *seedIndexIter
      << " at radius " << *seedRadiusIter << std::endl;
    typename TubeType::Pointer xTube =
      this->m_TubeExtractorFilter->ExtractTube( *seedIndexIter, count, true );
    if( !xTube.IsNull() )
      {
      this->m_TubeExtractorFilter->AddTube( xTube );
      std::cout << "  Extracted " << xTube->GetPoints().size() << " points."
        << std::endl;
      foundOneTube = true;
      }
    else
      {
      std::cout << " Error: Ridge not found for seed # " << count;
      }
    ++seedIndexIter;
    ++seedRadiusIter;
    ++count;
    }
  if ( !foundOneTube )
    {
    std::cout << "No Ridge found at all";
    return;
    }

  //Update tubes transform
  typename TubeTransformType::InputVectorType scaleVector;
  typename TubeTransformType::OffsetType offsetVector;
  typename TubeTransformType::MatrixType directionMatrix;
  typename ImageType::SpacingType spacing = this->m_InputImage->GetSpacing();
  typename ImageType::PointType origin = this->m_InputImage->GetOrigin();
  for (unsigned int i = 0; i < ImageDimension; ++i)
    {
    scaleVector[i] = spacing[i];
    offsetVector[i] = origin[i];
    }

  this->m_TubeExtractorFilter->GetTubeGroup()->GetObjectToParentTransform()->SetScale(
    scaleVector );
  this->m_TubeExtractorFilter->GetTubeGroup()->GetObjectToParentTransform()->SetOffset(
    offsetVector );
  this->m_TubeExtractorFilter->GetTubeGroup()->GetObjectToParentTransform()->SetMatrix(
    this->m_InputImage->GetDirection() );
  this->m_TubeExtractorFilter->GetTubeGroup()->ComputeObjectToWorldTransform();

  std::cout << "Ridge termination code counts:" << std::endl;
  for( unsigned int code = 0; code <
    this->m_TubeExtractorFilter->GetRidgeOp()->GetNumberOfFailureCodes();
    ++code )
    {
    std::cout << "   " << this->m_TubeExtractorFilter->GetRidgeOp()->GetFailureCodeName(
      typename RidgeExtractorFilterType::FailureCodeEnum( code ) ) << " : "
      << this->m_TubeExtractorFilter->GetRidgeOp()->GetFailureCodeCount(
      typename RidgeExtractorFilterType::FailureCodeEnum(
      code ) ) << std::endl;
    }
}

/**
 * Get list of extracted tubes */
template< class TInputImage >
typename SegmentTubes< TInputImage >::TubeGroupType::Pointer
SegmentTubes<TInputImage>
::GetTubeGroup( void )
{
  return this->m_TubeExtractorFilter->GetTubeGroup();
}

/**
 * Get the tube mask image */
template< class TInputImage >
typename SegmentTubes<TInputImage>::TubeMaskImageType::Pointer
SegmentTubes<TInputImage>
::GetTubeMaskImage( void )
{
  return this->m_TubeExtractorFilter->GetTubeMaskImage();
}

/**
 * PrintSelf */
template< class TInputImage >
void SegmentTubes<TInputImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( this->m_InputImage.IsNotNull() )
    {
    os << indent << "Input Image = " << this->m_InputImage << std::endl;
    }
  else
    {
    os << indent << "Input Image = NULL" << std::endl;
    }

  if( this->m_RadiusInputImage.IsNotNull() )
    {
    os << indent << "Radius Input Image = " << this->m_RadiusInputImage
      << std::endl;
    }
  else
    {
    os << indent << "Radius Input Image = NULL" << std::endl;
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeSegmentTubes_hxx)
