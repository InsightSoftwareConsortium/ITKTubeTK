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
#include "itktubeTubeToTubeTransformFilter.h"

#include <itkImageFileWriter.h>
#include <itkMetaDataObject.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkSpatialObjectWriter.h>
#include <itkEuler3DTransform.h>

/**
 *  This test is a base to generate images and spatial objects for the
 *  registration/metric testing process.
 */

int itktubeSyntheticTubeImageGenerationTest( int argc, char * argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Output_BlurredTubeImage "
              << "Output_VesselTube "
              << "Output_VesselTubeImage "
              << "Input_VesselTubeManuallyModified "
              << "Output_TransformedVesselTubeImage."
              << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<double, 3>                             Image3DType;
  typedef itk::ImageRegionIteratorWithIndex< Image3DType >  Image3DIteratorType;
  typedef itk::TubeSpatialObject<3>                         TubeType;
  typedef itk::TubeSpatialObjectPoint<3>                    TubePointType;
  typedef itk::GroupSpatialObject<3>                        TubeNetType;

  typedef itk::SpatialObjectToImageFilter< TubeNetType, Image3DType >
    SpatialObjectToImageFilterType;
  typedef itk::Euler3DTransform<double>
    TransformType;
  typedef itk::tube::TubeToTubeTransformFilter<TransformType, 3>
    TubeTransformFilterType;

  typedef itk::ImageFileWriter<Image3DType>                 ImageWriterType;
  typedef itk::SpatialObjectWriter<3>                       TubeWriterType;

  Image3DType::SizeType imageSize;
  imageSize[0] = 32;
  imageSize[1] = 32;
  imageSize[2] = 32;

  //------------------------------------------------------------------
  // Generate a simple tube image using Gaussian Filter
  //------------------------------------------------------------------
  std::cout << "Generate a tube blured image..." << std::endl;
  Image3DType::Pointer fixedImage = Image3DType::New();
  fixedImage->SetRegions( imageSize );
  fixedImage->Allocate();
  fixedImage->FillBuffer( 0 );
  fixedImage->Update();

  std::cout << "Start Filling Images ( Square Tube )..." << std::endl;
  Image3DIteratorType fixedIt( fixedImage, fixedImage->GetBufferedRegion() );
  int pixelIndex = 1;
  for( fixedIt.GoToBegin(); !fixedIt.IsAtEnd(); ++fixedIt, ++pixelIndex )
    {
    Image3DType::IndexType index = fixedIt.GetIndex();
    if( ( index[0]>=15 )&&( index[0]<=25 )&&( index[1]>=15 )&&( index[1]<=25 ) )
      {
        fixedIt.Set( 255 - 20 * ( pixelIndex % 5 ) ); // Brighter center
      }
    }

  // Gaussian blur the images to increase the likelihood of vessel
  // spatial object overlapping.
  std::cout << "Apply Gaussian blur..." << std::endl;
  typedef itk::RecursiveGaussianImageFilter<Image3DType, Image3DType>
    GaussianBlurFilterType;
  GaussianBlurFilterType::Pointer blurFilters[3];
  for( int i = 0; i < 3; i++ )
    {
    blurFilters[i] = GaussianBlurFilterType::New();
    blurFilters[i]->SetSigma( 3.0 );
    blurFilters[i]->SetZeroOrder();
    blurFilters[i]->SetDirection( i );
    }
  blurFilters[0]->SetInput( fixedImage );
  blurFilters[1]->SetInput( blurFilters[0]->GetOutput() );
  blurFilters[2]->SetInput( blurFilters[1]->GetOutput() );
  try
    {
    blurFilters[0]->Update();
    blurFilters[1]->Update();
    blurFilters[2]->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // write image
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( argv[1] );
  imageWriter->SetUseCompression( true );
  std::cout << "Write imageFile: " << argv[1] << std::endl;
  imageWriter->SetInput( blurFilters[2]->GetOutput() );
  try
    {
    imageWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  //------------------------------------------------------------------
  // Generate a simple spatial object tube
  //------------------------------------------------------------------
  std::cout << "Create spatial object tube..." << std::endl;
  TubeType::Pointer tube = TubeType::New();

  // Try to add the metaData about the vessel object subtype
  // There is currently some issues on it with ITK.
  // See:
  // http://www.itk.org/Wiki/ITK/Examples/Broken/SimpleOperations/MetaDataDictionary
  // http://public.kitware.com/Bug/view.php?id=12329#bugnotes
  itk::MetaDataDictionary& tubeMetaDictionary = tube->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::string>( tubeMetaDictionary,
                                         "ObjectSubType",
                                         "Vessel" );

  TubePointType point;
  point.SetRadiusInObjectSpace( 2.0 );

  typename TubePointType::PointType pnt;
  for( int i = -550; i < 550; ++i )
    {
    pnt[0] = 15;
    pnt[1] = 15;
    pnt[2] = i/10.0;
    point.SetPositionInObjectSpace( pnt );
    tube->GetPoints().push_back( point );
    }

  TubeNetType::Pointer group = TubeNetType::New();
  group->AddChild( tube );

  std::cout << "Write tubeFile: " << argv[2] << std::endl;
  TubeWriterType::Pointer tubeWriter = TubeWriterType::New();
  tubeWriter->SetFileName( argv[2] );
  tubeWriter->SetInput( group );
  try
    {
    tubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  //------------------------------------------------------------------
  // Write the tube as an image without transformation
  //------------------------------------------------------------------
  std::cout << "Convert the tube into an Image..." << std::endl;
  SpatialObjectToImageFilterType::Pointer imageFilter =
    SpatialObjectToImageFilterType::New();
  imageFilter->SetInput( group );
  imageFilter->SetSize( imageSize );

  Image3DType::PointType origin;
  origin[0] = 0;
  origin[1] = 0;
  origin[2] = 0;
  imageFilter->SetOrigin( origin );
  imageFilter->Update();

  // write image
  ImageWriterType::Pointer imageTubeWriter = ImageWriterType::New();
  imageTubeWriter->SetFileName( argv[3] );
  std::cout << "Write tubeAsImageFile: " << argv[3] << std::endl;
  imageTubeWriter->SetInput( imageFilter->GetOutput() );
  try
    {
    imageTubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  //------------------------------------------------------------------
  // Write the tube as an image with a transformation
  //------------------------------------------------------------------
  std::cout << "Transform and Convert the tube into an Image..." << std::endl;

  // read tube ( spatialObject )
  typedef itk::SpatialObjectReader<3> TubeNetReaderType;
  TubeNetReaderType::Pointer tubeReader = TubeNetReaderType::New();
  std::cout << "Read VesselTube: " << argv[4] << std::endl;
  tubeReader->SetFileName( argv[4] );
  try
    {
    tubeReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  TransformType::Pointer transformTube = TransformType::New();

  TubeType::VectorType translateT;
  translateT[0] = 2.5;
  translateT[1] = 2.5;
  translateT[2] = 2.5;
  transformTube->Translate( translateT );

  TubeType::ScalarType angleX = 0;
  TubeType::ScalarType angleY = 5;
  TubeType::ScalarType angleZ = 0;

  transformTube->Translate( translateT );
  transformTube->SetRotation( angleX, angleY, angleZ );

  itk::Matrix<double,3,3> rotationMatrix;
  itk::Vector<double,3> translation;

  rotationMatrix = transformTube->GetMatrix();
  translation = transformTube->GetTranslation();

  std::cout << rotationMatrix( 0,0 ) << " " << rotationMatrix( 0,1 )
            << " " << rotationMatrix( 0,2 ) << std::endl;
  std::cout << rotationMatrix( 1,0 ) << " " << rotationMatrix( 1,1 )
            << " " << rotationMatrix( 1,2 ) << std::endl;
  std::cout << rotationMatrix( 2,0 ) << " " << rotationMatrix( 2,1 )
            << " " << rotationMatrix( 2,2 ) << std::endl;
  std::cout << translation[0] << " " << translation[1]
            << " " << translation[2] << std::endl;

  // create transform filter
  TubeTransformFilterType::Pointer transformFilter =
    TubeTransformFilterType::New();
  transformFilter->SetInput( tubeReader->GetGroup() );
  transformFilter->SetTransform( transformTube );

  try
    {
    transformFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  SpatialObjectToImageFilterType::Pointer imageFilterTransform =
    SpatialObjectToImageFilterType::New();
  imageFilterTransform->SetInput( transformFilter->GetOutput() );
  imageFilterTransform->SetSize( imageSize );
  imageFilterTransform->SetOrigin( origin );
  imageFilterTransform->Update();

  // write image
  ImageWriterType::Pointer imageTubeWriterT = ImageWriterType::New();
  imageTubeWriterT->SetFileName( argv[5] );
  std::cout << "Write transformedTubeAsImageFile: " << argv[5] << std::endl;
  imageTubeWriterT->SetInput( imageFilterTransform->GetOutput() );
  try
    {
    imageTubeWriterT->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
