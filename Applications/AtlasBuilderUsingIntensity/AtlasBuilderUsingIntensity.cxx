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

#include "itktubeObjectDocumentToImageFilter.h"
#include "tubeAtlasSummation.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeMetaObjectDocument.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>

#include <metaCommand.h>

#include <cstdio>

#include "AtlasBuilderUsingIntensityCLP.h"

using namespace tube;

const unsigned int Dimensions = 3;

/** Processors */
typedef AtlasSummation                                AtlasBuilderType;
typedef itk::AffineTransform< double, Dimensions >    TransformType;

/** Object types */
typedef AtlasBuilderType::InputPixelType              InputPixelType;
typedef itk::Image< InputPixelType, Dimensions >      InputImageType;
typedef AtlasSummation::MeanPixelType                 MeanPixelType;
typedef itk::Image< MeanPixelType, Dimensions >       MeanImageType;
typedef AtlasSummation::VariancePixelType             VariancePixelType;
typedef itk::Image< VariancePixelType, Dimensions >   VarianceImageType;
typedef itk::Image< short, Dimensions >               ShortImageType;
typedef itk::Image< unsigned short, Dimensions >      UShortImageType;
typedef UShortImageType::Pointer                      UShortImagePointer;
typedef itk::Image< float, Dimensions >               FloatImageType;

/** IO */
typedef MetaObjectDocument                            DocumentReaderType;
typedef itk::tube::ImageDocument                      ImageDocumentType;
typedef DocumentReaderType::ObjectListType            ImageDocumentListType;
typedef itk::tube::ObjectDocumentToImageFilter<
  ImageDocumentType, InputImageType >                 DocumentToImageFilter;

/** Function declarations */
void WriteImage( UShortImagePointer image, const std::string & file );
void WriteImage( FloatImageType::Pointer image, const std::string & file );
void SetParameterList( AtlasBuilderType * atlasBuilder , MetaCommand command);
int DoIt( int argc, char *argv[] );

int main( int argc, char* argv[] )
{
  PARSE_ARGS;
  return DoIt( argc, argv );
}

int DoIt( int argc, char *argv[] )
{
  PARSE_ARGS;

  // Reads the ObjectDocuments file
  DocumentReaderType * reader = new DocumentReaderType();
  reader->FileName( imageDocFile.c_str() );
  if( !reader->Read( imageDocFile.c_str() ) )
    {
    tube::ErrorMessage( "Could not read input document!" );
    return EXIT_FAILURE;
    }
  tube::FmtInfoMessage( "Read input document %s", imageDocFile.c_str() );
  ImageDocumentListType * imageObjects = reader->GetObjectList();

  AtlasBuilderType * atlasBuilder = new AtlasBuilderType();

  assert( outputSize.size() == Dimensions );

  AtlasBuilderType::SizeType size;
  for( unsigned int i = 0; i < outputSpacing.size(); ++i )
    {
    size[i] = outputSpacing[i];
    }
  atlasBuilder->SetOutputSize( size );

  assert( outputSpacing.size() > 0 );

  AtlasBuilderType::SpacingType spacing;
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

  ImageDocumentListType::const_iterator it_imgDoc = imageObjects->begin();

  tube::FmtInfoMessage( "Starting image addition..." );

  /*
   * Iteratively add the images to Atlas Summation method
   * (this is done so that only one image must be held in memory at a time )
   */
  DocumentToImageFilter::Pointer filter = DocumentToImageFilter::New();
  while( it_imgDoc != imageObjects->end() )
    {
    ImageDocumentType::Pointer doc =
      static_cast< ImageDocumentType * >( (*it_imgDoc).GetPointer() );

    tube::FmtInfoMessage("Adding image: " + doc->GetObjectName());

    filter->SetInput( doc );
    filter->ApplyTransforms( false );
    filter->Update();

    filter->GetComposedTransform();
    if(filter->ComposedTransformIsIdentity())
      {
      ::tube::DebugMessage("ComposedTransform is IDENTITY");
      atlasBuilder->AddImage( filter->GetOutput() );
      }
    else
      {
      atlasBuilder->AddImage(
        filter->GetOutput(),
        filter->GetComposedTransform().GetPointer() );
      }
    ++it_imgDoc;
    }
  tube::InfoMessage( "Finalizing the images..." );
  atlasBuilder->Finalize();
  tube::InfoMessage( "Done!" );

  // Save output images
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

void WriteImage( UShortImageType::Pointer i, const std::string & name )
{
  typedef itk::ImageFileWriter< UShortImageType > UShortWriterType;
  UShortWriterType::Pointer writer  = UShortWriterType::New();

  writer->SetInput( i );
  writer->SetFileName( name );
  writer->SetUseCompression( true );
  writer->Update();
}

void WriteImage( FloatImageType::Pointer i, const std::string & name )
{
  typedef itk::ImageFileWriter< FloatImageType > WriterType;
  WriterType::Pointer writer  = WriterType::New();

  writer->SetInput( i );
  writer->SetFileName( name );
  writer->SetUseCompression( true );
  writer->Update();
}
