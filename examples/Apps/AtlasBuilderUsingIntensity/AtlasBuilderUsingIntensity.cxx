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

#include "itktubeObjectDocumentToImageFilter.h"
#include "tubeAtlasSummation.h"
#include "tubeMetaObjectDocument.h"

#include "AtlasBuilderUsingIntensityCLP.h"

enum { Dimension = 3 };

typedef itk::Image< float, Dimension > FloatImageType;

void WriteImage( FloatImageType::Pointer i, const std::string & name )
{
  typedef itk::ImageFileWriter< FloatImageType > WriterType;
  WriterType::Pointer writer  = WriterType::New();

  writer->SetInput( i );
  writer->SetFileName( name );
  writer->SetUseCompression( true );
  writer->Update();
}

typedef itk::Image< unsigned short, Dimension > UShortImageType;

void WriteImage( UShortImageType::Pointer i, const std::string & name )
{
  typedef itk::ImageFileWriter< UShortImageType > UShortWriterType;
  UShortWriterType::Pointer writer  = UShortWriterType::New();

  writer->SetInput( i );
  writer->SetFileName( name );
  writer->SetUseCompression( true );
  writer->Update();
}

int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef tube::AtlasSummation                       AtlasSummationType;
  typedef AtlasSummationType::InputPixelType         InputPixelType;
  typedef itk::Image< InputPixelType, Dimension >    InputImageType;
  typedef tube::MetaObjectDocument                   DocumentReaderType;
  typedef itk::tube::ImageDocument                   ImageDocumentType;
  typedef DocumentReaderType::ObjectDocumentListType ImageDocumentListType;
  typedef
    itk::tube::ObjectDocumentToImageFilter< ImageDocumentType, InputImageType >
    DocumentToImageFilter;

  // Reads the ObjectDocuments file
  DocumentReaderType::Pointer reader = new DocumentReaderType();
  reader->SetFileName( imageDocFile.c_str() );

  if( !reader->Read( imageDocFile.c_str() ) )
    {
    tube::ErrorMessage( "Could not read input document!" );

    return EXIT_FAILURE;
    }

  tube::FmtInfoMessage( "Read input document %s", imageDocFile.c_str() );
  ImageDocumentListType imageObjects = reader->GetObjectDocumentList();

  AtlasSummationType * atlasBuilder = new tube::AtlasSummation();

  assert( outputSize.size() == Dimension );

  AtlasSummationType::SizeType size;

  for( unsigned int i = 0; i < outputSpacing.size(); ++i )
    {
    size[i] = outputSpacing[i];
    }

  atlasBuilder->SetOutputSize( size );

  assert( outputSpacing.size() > 0 );

  AtlasSummationType::SpacingType spacing;

  for( unsigned int i = 0; i < outputSpacing.size(); ++i )
    {
    spacing[i] = outputSpacing[i];
    }

  // Configure the atlas builder
  atlasBuilder->SetOutputSpacing( spacing );
  atlasBuilder->SetUseStdDeviation( useStdDeviation );
  atlasBuilder->SetImageCountThreshold( lowerThreshold );
  atlasBuilder->AdjustResampledImageSize( doImageSizeAdjustment );
  atlasBuilder->AdjustResampledImageOrigin( doImageOriginAdjustment );

  ImageDocumentListType::const_iterator it_imgDoc = imageObjects.begin();
  tube::FmtInfoMessage( "Starting image addition..." );

  /* Iteratively add the images to atlas summation method. This is done so that
     only one image must be held in memory at a time. */
  DocumentToImageFilter::Pointer filter = DocumentToImageFilter::New();

  while( it_imgDoc != imageObjects.end() )
    {
    ImageDocumentType::Pointer doc
      = static_cast< ImageDocumentType * >( ( *it_imgDoc ).GetPointer() );
    tube::FmtInfoMessage( "Adding image: %s", doc->GetObjectName().c_str() );

    filter->SetInput( doc );
    filter->SetApplyTransforms( false );
    filter->Update();
    filter->GetComposedTransform();

    if( filter->GetComposedTransformIsIdentity() )
      {
      tube::DebugMessage( "ComposedTransform is IDENTITY" );
      atlasBuilder->AddImage( filter->GetOutput() );
      }
    else
      {
      atlasBuilder->AddImage( filter->GetOutput(),
                              filter->GetComposedTransform().GetPointer() );
      }
    ++it_imgDoc;
    }

  tube::InfoMessage( "Finalizing the images..." );
  atlasBuilder->Finalize();
  tube::InfoMessage( "Done!" );

  // Save output images.
  WriteImage( atlasBuilder->GetMeanImage(), outputMeanAtlas.c_str() );
  WriteImage( atlasBuilder->GetVarianceImage(), outputVarianceAtlas.c_str() );

  if( !outputCountImage.empty() )
    {
    WriteImage( atlasBuilder->GetValidCountImage(), outputCountImage.c_str() );
    }

  delete reader;
  delete atlasBuilder;

  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return DoIt( argc, argv );
}
