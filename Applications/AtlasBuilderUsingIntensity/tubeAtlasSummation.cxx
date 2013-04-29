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


#include "tubeAtlasSummation.h"

#include <cstdio>

namespace tube
{


//-----------------------------------------------------------------------------
AtlasSummation
::AtlasSummation()
{
  m_isStdDeviation = true;
  m_isProcessing = false;

  m_meanBuilder = RobustMeanBuilderType::New();

  m_AdjustResampledImageSize    = false;
  m_AdjustResampledImageOrigin  = false;

  m_OutputSize.Fill(256);
  m_OutputSpacing.Fill(1);
  m_OutputOrigin.Fill(0);
  m_OutSizeSet    = false;
  m_OutSpacingSet = false;
  m_OutOriginSet  = false;

  m_ImageCountThreshold = 1; // Must be > than 0.

  m_image_number = 0;
  m_numOfImages  = 0;  // Variable used for median calculations
  TMedianDefaultPixelValue = itk::NumericTraits<InputPixelType>::max();
  count = 0;
}

//-----------------------------------------------------------------------------
AtlasSummation
::~AtlasSummation()
{}


//-----------------------------------------------------------------------------
void AtlasSummation
::AddImage( InputImageType::Pointer i )
{
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  this->AddImage( i, transform );
}


//-----------------------------------------------------------------------------
void AtlasSummation
::AddImage( InputImageType::Pointer image, TransformType::Pointer t )
{
  InputImageType::Pointer image2 = image;
  if( m_AdjustResampledImageSize || m_AdjustResampledImageOrigin )
    {
    // Form the clipped image and update the base image if necessary
    image = GetClippedImage( image, t );
    SizeType inputSize = image->GetLargestPossibleRegion().GetSize();

    // Did we not start processing
    if( !m_isProcessing )
      {
      Start( image ); // Make sure m_image_nr is reset to 0
      m_meanBuilder->AddImage(image); // Actually add the image
      m_isProcessing = true; // Now, processing has started
      ++m_image_number;
      return;
      }
    else if( UpdateOutputSizeParameter( inputSize ) )
      {
      ::tube::DebugMessage("Updating output image size!");
      m_meanBuilder->UpdateOutputImageSize( inputSize );
      }

    TransformPointer identity = TransformType::New();
    identity->SetIdentity();
    image2 = TransformInputImage( image, identity,
      m_meanBuilder->GetOutputSize(),
      m_meanBuilder->GetOutputSpacing(),
      m_meanBuilder->GetOutputOrigin() );
    }
  else // No adjustments or resampling
    {
    if(!m_isProcessing)
      {
      Start( image );
      m_meanBuilder->AddImage(image);
      m_isProcessing = true;
      ++m_image_number;
      return;
      }
    else
      {
      image2 = TransformInputImage( image, t,
        m_meanBuilder->GetOutputSize(),
        m_meanBuilder->GetOutputSpacing(),
        m_meanBuilder->GetOutputOrigin() );
      }
    }

  m_meanBuilder->AddImage( image2 );
  ++m_image_number;
}


//-----------------------------------------------------------------------------
AtlasSummation::InputImagePointer AtlasSummation
::TransformInputImage( InputImagePointer image,
                       TransformPointer trans,
                       SizeType size,
                       SpacingType spacing,
                       PointType origin )
{
  /*
   * Realign image to the base Image specifications. If no adjustments
   * are done to the image, then add the transform to the filter as well
   */
  typedef itk::ResampleImageFilter< InputImageType, InputImageType >
               ResampleFilterType;

  ResampleFilterType::Pointer transfilter = ResampleFilterType::New();
  transfilter->SetInput( image );
  transfilter->SetOutputSpacing( spacing );
  transfilter->SetSize( size );
  transfilter->SetOutputOrigin( origin  );
  transfilter->SetDefaultPixelValue( DEFAULT_PIXEL_FILL );

  /*
   * For resampling, we need to use the inverse of the transform to get
   * the desired result. Resample needs a FIXED -> MOVING image trans.
   */
  TransformType::Pointer inverse = TransformType::New();
  trans->GetInverse(inverse);
  transfilter->SetTransform(inverse);
  transfilter->Update();

  return transfilter->GetOutput();
}


//-----------------------------------------------------------------------------
void AtlasSummation
::Start( InputImageType::Pointer )
{
  m_image_number = 0;
}


//-----------------------------------------------------------------------------
AtlasSummation::InputImagePointer AtlasSummation
::GetClippedImage( InputImagePointer image, TransformType::Pointer t )
{
  /*
   * The filter does only require output spacing (not size or origin)
   * to run the resampling
   */
  typedef itk::tube::CompleteImageResampleFilter<
    InputImageType, InputImageType, TransformType> ResampleImageFilterType;

  TransformType::Pointer inverse = TransformType::New();
  t->GetInverse(inverse);

  ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();
  resampleFilter->SetInput( image );
  resampleFilter->SetTransform( inverse);
  resampleFilter->SetOutputSpacing( image->GetSpacing() );
  resampleFilter->SetDefaultPixelValue( DEFAULT_PIXEL_FILL );
  resampleFilter->Update();

  // Minimize the image to fit the parameters
  typedef itk::tube::MinimizeImageSizeFilter<InputImageType> MinImageSizeFilterType;

  MinImageSizeFilterType::Pointer minFilter = MinImageSizeFilterType::New();
  minFilter->SetInput( resampleFilter->GetOutput() );
  minFilter->SetThresholdValue( DEFAULT_PIXEL_FILL );
  minFilter->SetDefaultPixelValue( DEFAULT_PIXEL_FILL );

  // Set the clip areas based on the input requests
  if( m_AdjustResampledImageSize )
    {
    minFilter->ClipEndIndicesOn();
    }
  if( m_AdjustResampledImageOrigin )
    {
    minFilter->ClipStartIndicesOn();
    }

  SizeType buffer;
  for( int i = 0; i < TDimensions; i++ )
    {
    buffer[i] = 2;
    }

  minFilter->SetNumberOfBufferPixels( buffer );
  minFilter->Update();
  return minFilter->GetOutput();
}


//-----------------------------------------------------------------------------
bool AtlasSummation
::UpdateOutputSizeParameter( SizeType& inputSize )
{
  SizeType outputSize = m_meanBuilder->GetOutputSize();

  bool sizeChanged = false;
  for( int i = 0; i < TDimensions; i++ )
    {
    if( inputSize[i] > outputSize[i] )
      {
      inputSize[i] = outputSize[i];
      sizeChanged = true;
      }
    }
  return sizeChanged;
}


//-----------------------------------------------------------------------------
void AtlasSummation
::Finalize()
{
  m_meanBuilder->FinalizeOutput();
}


//-----------------------------------------------------------------------------
void AtlasSummation
::WriteImage( MeanImageType::Pointer image, const std::string & file )
{
  typedef itk::ImageFileWriter< MeanImageType > FileWriterType;
  ++count;
  std::stringstream file1;
  file1 << file << count << ".mhd";

  try
    {
    FileWriterType::Pointer writer = FileWriterType::New();
    writer->SetInput( image );
    writer->SetFileName( file1.str().c_str() );
    writer->Update();
    }
  catch( ... )
    {
    ::tube::FmtErrorMessage("Error in writing " + file1.str());
    }
}


//-----------------------------------------------------------------------------
void AtlasSummation
::WriteImage( ProcessImagePointer image, const std::string & file )
{
  typedef itk::ImageFileWriter< ProcessImageType > FileWriterType;

  try
    {
    FileWriterType::Pointer writer = FileWriterType::New();
    writer->SetInput( image );
    writer->SetFileName( file );
    writer->Update();
    }
  catch( ... )
    {
    ::tube::FmtErrorMessage("Error in writing " + file);
    }
}


} // End of namespace tube
