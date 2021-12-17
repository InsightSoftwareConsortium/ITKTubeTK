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
#ifndef __itktubeMergeAdjacentImagesFilter_hxx
#define __itktubeMergeAdjacentImagesFilter_hxx

// ITK includes
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageToImageRegistrationHelper.h>
#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include <itkTimeProbesCollectorBase.h>

// TubeTK includes
#include "itkGeneralizedDistanceTransformImageFilter.h"

namespace itk
{

namespace tube
{

template< class TImage >
MergeAdjacentImagesFilter< TImage >
::MergeAdjacentImagesFilter( void )
{
  m_Background = 0;
  m_MaskZero = false;
  m_MaxIterations = 300;
  m_ExpectedOffset = 20;
  m_ExpectedRotation = 0.001;
  m_SamplingRatio = 0.01;
  m_BlendUsingAverage = false;
  m_UseFastBlending = false;

  m_InitialTransform = nullptr;
  m_OutputTransform = nullptr;

  m_Input1 = nullptr;
  m_Input2 = nullptr;
  m_Output = ImageType::New();
}

template< class TImage >
void
MergeAdjacentImagesFilter< TImage >
::SetInput1( const TImage* image )
{
  if( this->m_Input1.GetPointer() != image )
    {
    this->m_Input1 = image;
    }
}

template< class TImage >
void
MergeAdjacentImagesFilter< TImage >
::SetInput2( const TImage* image )
{
  if( this->m_Input2.GetPointer() != image )
    {
    this->m_Input2 = image;
    }
}

template< class TImage >
void
MergeAdjacentImagesFilter< TImage >
::SetPadding( const PaddingType & padding )
{
  m_Padding = padding;
}

template< class TImage >
void
MergeAdjacentImagesFilter< TImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << "Background: " << m_Background << std::endl;
  os << "MaskZero: " << m_MaskZero << std::endl;
  os << "MaxIterations: " << m_MaxIterations << std::endl;
  os << "ExpectedOffset: " << m_ExpectedOffset << std::endl;
  os << "ExpectedRotation: " << m_ExpectedRotation << std::endl;
  os << "SamplingRatio: " << m_SamplingRatio << std::endl;
  os << "BlendUsingAverage: " << m_BlendUsingAverage << std::endl;
  os << "UseFastBlending: " << m_UseFastBlending << std::endl;
}

template< class TImage >
void
MergeAdjacentImagesFilter< TImage >
::LoadInitialTransform( const std::string & filename )
{
  typedef TransformFileReader                    TransformReaderType;
  typedef TransformReaderType::TransformListType TransformListType;

  TransformReaderType::Pointer transformReader = TransformReaderType::New();
  transformReader->SetFileName( filename );
  transformReader->Update();

  const TransformListType *transforms = transformReader->GetTransformList();
  TransformListType::const_iterator transformIt = transforms->begin();
  while( transformIt != transforms->end() )
    {
    if( !strcmp( (*transformIt)->GetNameOfClass(), "AffineTransform") )
      {
      typename TransformType::Pointer affine_read =
        static_cast<TransformType *>( (*transformIt).GetPointer() );
      typename TransformType::ConstPointer affine =
        affine_read.GetPointer();
      SetInitialTransform( affine.GetPointer() );
      break;
      }
    ++transformIt;
    }
}

template< class TImage >
void
MergeAdjacentImagesFilter< TImage >
::SaveOutputTransform( const std::string & filename )
{
  if( m_OutputTransform.IsNotNull() )
    {
    typedef TransformFileWriter TransformWriterType;
  
    TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetFileName( filename );
    transformWriter->SetInput( m_OutputTransform );
    transformWriter->Update();
    }
}

template< class TImage >
void
MergeAdjacentImagesFilter< TImage >
::Update()
{
  // The timeCollector is used to perform basic profiling algorithm components
  // itk::TimeProbesCollectorBase timeCollector;

  // compute min and max coord of input1 with padding
  typename ImageType::IndexType minX1Org;

  minX1Org = m_Input1->GetLargestPossibleRegion().GetIndex();

  if( m_Padding.size() == ImageDimension )
    {
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      minX1Org[i] -= m_Padding[i];
      }
    }
  typename ImageType::SizeType size1 =
    m_Input1->GetLargestPossibleRegion().GetSize();

  typename ImageType::IndexType maxX1Org;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    maxX1Org[i] = minX1Org[i] + size1[i] - 1;
    }
  if( m_Padding.size() == ImageDimension )
    {
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      maxX1Org[i] += 2 * m_Padding[i];
      }
    }

  typename ImageType::IndexType minXOut;
  typename ImageType::IndexType maxXOut;
  typename ImageType::SizeType  sizeOut;

  minXOut = minX1Org;
  maxXOut = maxX1Org;

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    sizeOut[i] = maxXOut[i] - minXOut[i] + 1;
    }

  // compute min-max coord and size of input2 in input1's space with padding
  typename ImageType::IndexType minX2;
  typename ImageType::IndexType minX2Org;

  minX2Org = m_Input2->GetLargestPossibleRegion().GetIndex();
  if( m_Padding.size() == ImageDimension )
    {
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      minX2Org[i] -= m_Padding[i];
      }
    }

  typename ImageType::PointType pointX;
  m_Input2->TransformIndexToPhysicalPoint( minX2Org, pointX );

  if( m_InitialTransform.IsNotNull() )
    {
    pointX = m_InitialTransform->GetInverseTransform()->TransformPoint( pointX );
    }

  m_Input1->TransformPhysicalPointToIndex( pointX, minX2 );

  typename ImageType::SizeType size2 =
    m_Input2->GetLargestPossibleRegion().GetSize();

  typename ImageType::IndexType maxX2;
  typename ImageType::IndexType maxX2Org;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    maxX2Org[i] = minX2Org[i] + size2[i] - 1;
    }
  if( m_Padding.size() == ImageDimension )
    {
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      maxX2Org[i] += 2*m_Padding[i];
      }
    }
  m_Input2->TransformIndexToPhysicalPoint( maxX2Org, pointX );
  if( m_InitialTransform.IsNotNull() )
    {
    pointX = m_InitialTransform->GetInverseTransform()->TransformPoint( pointX );
    }
  m_Input1->TransformPhysicalPointToIndex( pointX, maxX2 );

  // compute min-max coord and size of output
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( minX2[i] < minXOut[i] )
      {
      minXOut[i] = minX2[i];
      }
    if( maxX2[i] < minXOut[i] )
      {
      minXOut[i] = maxX2[i];
      }
    }

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( minX2[i] > maxXOut[i] )
      {
      maxXOut[i] = minX2[i];
      }
    if( maxX2[i] > maxXOut[i] )
      {
      maxXOut[i] = maxX2[i];
      }
    }

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    sizeOut[i] = maxXOut[i] - minXOut[i] + 1;
    }

  // Allocate output image
  // timeCollector.Start( "Allocate output image" );

  typename ImageType::RegionType regionOut;
  regionOut.SetSize( sizeOut );
  regionOut.SetIndex( minXOut );

  m_Output = ImageType::New();
  m_Output->CopyInformation( m_Input1 );
  m_Output->SetRegions( regionOut );
  m_Output->Allocate();
  m_Output->FillBuffer( m_Background );

  itk::ImageRegionConstIteratorWithIndex< ImageType > iter( m_Input1,
    m_Input1->GetLargestPossibleRegion() );
  while( !iter.IsAtEnd() )
    {
    typename ImageType::IndexType indexX = iter.GetIndex();
    double tf = iter.Get();
    if( !m_MaskZero || tf != 0 )
      {
      m_Output->SetPixel( indexX, tf );
      }
    ++iter;
    }

  // timeCollector.Stop( "Allocate output image" );

  // perform registration
  typedef typename itk::ImageToImageRegistrationHelper< ImageType >
    RegFilterType;
  typename RegFilterType::Pointer regOp = RegFilterType::New();
  regOp->SetFixedImage( m_Input1 );
  regOp->SetMovingImage( m_Input2 );
  regOp->SetSampleFromOverlap( true );
  regOp->SetEnableLoadedRegistration( false );
  regOp->SetEnableInitialRegistration( false );
  regOp->SetEnableRigidRegistration( true );
  regOp->SetRigidSamplingRatio( m_SamplingRatio );
  regOp->SetRigidMaxIterations( m_MaxIterations );
  regOp->SetEnableAffineRegistration( false );
  regOp->SetEnableBSplineRegistration( false );
  regOp->SetExpectedOffsetMagnitude( m_ExpectedOffset );
  regOp->SetExpectedRotationMagnitude( m_ExpectedRotation );

  if( m_InitialTransform.IsNotNull() )
    {
    regOp->SetLoadedMatrixTransform( *(m_InitialTransform.GetPointer()) );
    }

  regOp->Initialize();

  if( m_MaxIterations > 0 )
    {
    // timeCollector.Start( "Register images" );

    regOp->SetReportProgress( true );
    regOp->Update();

    // timeCollector.Stop( "Register images" );
    }

  m_OutputTransform = regOp->GetCurrentMatrixTransform();

  // Resample image
  typename ImageType::ConstPointer tmpImage;

  // timeCollector.Start( "Resample Image" );

  regOp->SetFixedImage( m_Output );

  tmpImage = regOp->ResampleImage(
    RegFilterType::OptimizedRegistrationMethodType::LINEAR_INTERPOLATION,
    m_Input2, NULL, NULL, m_Background );

  // timeCollector.Stop( "Resample Image" );
  typedef typename itk::Image<float, ImageDimension> FloatImageType;

  if( m_BlendUsingAverage )
    {
    itk::ImageRegionConstIteratorWithIndex< ImageType > iter2( tmpImage,
      tmpImage->GetLargestPossibleRegion() );

    itk::ImageRegionIteratorWithIndex< ImageType > iterOut( m_Output,
      m_Output->GetLargestPossibleRegion() );

    while( !iter2.IsAtEnd() )
      {
      double iVal = iter2.Get();
      bool image2Point = false;
      if( iVal != m_Background && ( !m_MaskZero || iVal != 0 ) )
        {
        image2Point = true;
        }
      double oVal = iterOut.Get();
      bool imageOutPoint = false;
      if( oVal != m_Background && ( !m_MaskZero || oVal != 0 ) )
        {
        imageOutPoint = true;
        }

      if( imageOutPoint )
        {
        if( image2Point )
          {
          oVal = 0.5 * iVal + 0.5 * oVal;
          }
        }
      else if( image2Point )
        {
        oVal = iVal;
        }

      iterOut.Set( oVal );

      ++iter2;
      ++iterOut;
      }
    }
  else
    {
    // timeCollector.Start( "Resample Image2" );

    typename ImageType::Pointer input2Reg = ImageType::New();
    input2Reg->CopyInformation( tmpImage );
    input2Reg->SetRegions( tmpImage->GetLargestPossibleRegion() );
    input2Reg->Allocate();
    input2Reg->FillBuffer( m_Background );

    typename FloatImageType::Pointer outputMap = FloatImageType::New();
    outputMap->CopyInformation( tmpImage );
    outputMap->SetRegions( tmpImage->GetLargestPossibleRegion() );
    outputMap->Allocate();
    outputMap->FillBuffer( 3 );

    itk::ImageRegionConstIteratorWithIndex< ImageType > iterTmp( tmpImage,
      tmpImage->GetLargestPossibleRegion() );

    itk::ImageRegionIteratorWithIndex< ImageType > iter2( input2Reg,
      input2Reg->GetLargestPossibleRegion() );

    itk::ImageRegionIteratorWithIndex< ImageType > iterOut( m_Output,
      m_Output->GetLargestPossibleRegion() );

    itk::ImageRegionIteratorWithIndex< FloatImageType > iterOutMap( outputMap,
      outputMap->GetLargestPossibleRegion() );

    while( !iter2.IsAtEnd() )
      {
      double iVal = iterTmp.Get();
      iter2.Set( iVal );
      bool image2Point = false;
      if( iVal != m_Background && ( !m_MaskZero || iVal != 0 ) )
        {
        image2Point = true;
        }
      double oVal = iterOut.Get();
      bool imageOutPoint = false;
      if( oVal != m_Background && ( !m_MaskZero || oVal != 0 ) )
        {
        imageOutPoint = true;
        }
      if( image2Point && imageOutPoint )
        {
        iterOutMap.Set( 0 );
        }
      else if( image2Point )
        {
        iterOutMap.Set( 2 );
        }
      else if( imageOutPoint )
        {
        iterOutMap.Set( 1 );
        }
      ++iterTmp;
      ++iter2;
      ++iterOut;
      ++iterOutMap;
      }

    // timeCollector.Stop( "Resample Image2" );

    // timeCollector.Start( "Out Distance Map" );
    typename FloatImageType::Pointer outputDistMap = nullptr;
    typename FloatImageType::Pointer outputVoronoiMap = nullptr;
    if( m_UseFastBlending )
      {
      typedef typename itk::GeneralizedDistanceTransformImageFilter<
        FloatImageType, FloatImageType, FloatImageType >   MapFilterType;
      typename MapFilterType::Pointer mapDistFilter = MapFilterType::New();

      typedef itk::BinaryThresholdImageFilter<FloatImageType, FloatImageType>
        Indicator;
      typename Indicator::Pointer indicator = Indicator::New();

      indicator->SetLowerThreshold( 0 );
      indicator->SetUpperThreshold( 0 );
      indicator->SetOutsideValue( 0 );
      indicator->SetInsideValue(
        mapDistFilter->GetMaximalSquaredDistance() );
      indicator->SetInput( outputMap );
      indicator->Update();

      mapDistFilter->SetInput1( indicator->GetOutput() );
      mapDistFilter->SetInput2( outputMap );
      mapDistFilter->UseImageSpacingOff();
      mapDistFilter->CreateVoronoiMapOn();
      mapDistFilter->Update();
      outputDistMap = mapDistFilter->GetOutput();
      outputVoronoiMap = mapDistFilter->GetVoronoiMap();
      }
    else
      {
      typedef typename itk::DanielssonDistanceMapImageFilter< FloatImageType,
        FloatImageType, FloatImageType >   MapFilterType;
      typename MapFilterType::Pointer mapDistFilter = MapFilterType::New();
      mapDistFilter->SetInput( outputMap );
      mapDistFilter->SetInputIsBinary( false );
      mapDistFilter->SetUseImageSpacing( false );
      mapDistFilter->Update();
      outputDistMap = mapDistFilter->GetDistanceMap();
      outputVoronoiMap = mapDistFilter->GetVoronoiMap();
      }
    // timeCollector.Stop( "Out Distance Map" );

    // timeCollector.Start( "Distance Map Selection" );

    typename FloatImageType::Pointer vorImageMap = FloatImageType::New();
    vorImageMap->CopyInformation( m_Output );
    vorImageMap->SetRegions( m_Output->GetLargestPossibleRegion() );
    vorImageMap->Allocate();
    vorImageMap->FillBuffer( 1 );
    itk::ImageRegionConstIteratorWithIndex< FloatImageType > iterOutVor(
      outputVoronoiMap, outputVoronoiMap->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex< FloatImageType > iterVorMap(
      vorImageMap, vorImageMap->GetLargestPossibleRegion() );
    while( !iterOutVor.IsAtEnd() )
      {
      double tf = iterOutVor.Get();
      if( tf != 1 )
        {
        iterVorMap.Set( 0 );
        }
      ++iterOutVor;
      ++iterVorMap;
      }

    // timeCollector.Stop( "Distance Map Selection" );

    typename FloatImageType::Pointer vorImageDistMap = nullptr;
    if( m_UseFastBlending )
      {
      typedef typename itk::GeneralizedDistanceTransformImageFilter<
        FloatImageType, FloatImageType >   MapFilterType;
      typename MapFilterType::Pointer mapVorFilter = MapFilterType::New();

      // timeCollector.Start( "Voronoi Distance Map" );

      typedef itk::BinaryThresholdImageFilter<FloatImageType, FloatImageType>
        Indicator;
      typename Indicator::Pointer indicator2 = Indicator::New();

      indicator2->SetLowerThreshold( 0 );
      indicator2->SetUpperThreshold( 0 );
      indicator2->SetOutsideValue( 0 );
      indicator2->SetInsideValue(
        mapVorFilter->GetMaximalSquaredDistance() );
      indicator2->SetInput( vorImageMap );
      indicator2->Update();

      mapVorFilter->SetInput1( indicator2->GetOutput() );
      mapVorFilter->SetInput2( vorImageMap );
      mapVorFilter->UseImageSpacingOff();
      mapVorFilter->CreateVoronoiMapOn();
      mapVorFilter->Update();
      vorImageDistMap = mapVorFilter->GetOutput();

      // timeCollector.Stop( "Voronoi Distance Map" );
      }
    else
      {
      // timeCollector.Start( "Voronoi Distance Map" );

      typedef typename itk::SignedDanielssonDistanceMapImageFilter<
        FloatImageType, FloatImageType>  SignedMapFilterType;
      typename SignedMapFilterType::Pointer mapVorFilter =
        SignedMapFilterType::New();

      mapVorFilter->SetInput( vorImageMap );
      mapVorFilter->SetUseImageSpacing( false );
      mapVorFilter->Update();
      vorImageDistMap = mapVorFilter->GetDistanceMap();

      // timeCollector.Stop( "Voronoi Distance Map" );
      }

    // timeCollector.Start( "Voronoi Distance Selection" );

    iter2.GoToBegin();
    iterOut.GoToBegin();
    itk::ImageRegionIteratorWithIndex< FloatImageType > iterOutDistMap(
      outputDistMap, outputDistMap->GetLargestPossibleRegion() );

    itk::ImageRegionIteratorWithIndex< FloatImageType > iterVorDistMap(
      vorImageDistMap, vorImageDistMap->GetLargestPossibleRegion() );

    while( !iter2.IsAtEnd() )
      {
      double iVal = iter2.Get();
      bool image2Point = false;
      if( iVal != m_Background && ( !m_MaskZero || iVal != 0 ) )
        {
        image2Point = true;
        }
      double oVal = iterOut.Get();
      bool imageOutPoint = false;
      if( oVal != m_Background && ( !m_MaskZero || oVal != 0 ) )
        {
        imageOutPoint = true;
        }

      if( imageOutPoint )
        {
        if( image2Point )
          {
          double vDist = iterVorDistMap.Get();
          double oDist = iterOutDistMap.Get();

          if( vDist < 0 && !m_BlendUsingAverage )
            {
            vDist = -vDist;
            double ratio = 0.5 * vDist / ( oDist + vDist ) + 0.5;
            oVal = ratio * oVal + ( 1-ratio )*iVal;
            }
          else if( vDist > 0 && !m_BlendUsingAverage )
            {
            double ratio = 0.5 * vDist / ( oDist + vDist ) + 0.5;
            oVal = ratio * iVal + ( 1 - ratio ) * oVal;
            }
          else
            {
            oVal = 0.5 * iVal + 0.5 * oVal;
            }
          }
        }
      else if( image2Point )
        {
        oVal = iVal;
        }

      iterOut.Set( oVal );

      ++iter2;
      ++iterOut;
      ++iterOutDistMap;
      ++iterVorDistMap;
      }
    // timeCollector.Stop( "Voronoi Distance Selection" );
    }
}

} // End namespace tube

} // End namespace itk

#endif
