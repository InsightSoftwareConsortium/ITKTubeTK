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

#include "tubeMessage.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkOrientImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkTranslationTransform.h>

#include "RegisterUsingImageCentersCLP.h"

template< unsigned int VDimension >
int DoIt( int argc, char * argv[] );

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  ( void )argc;
  ( void )argv;
  return 0;
}

#include "tubeCLIHelperFunctions.h"

// Typedefs independent of image type
typedef itk::ImageIOBase::IOComponentType ScalarPixelType;

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    unsigned int dimension;
    tube::GetImageInformation( 
      inputImageFile.c_str(), componentType, dimension );
    switch( dimension )
      {
      case 3:
        return DoIt< 3 >( argc, argv );
      default:
        tube::ErrorMessage( "Only 3D images supported!" );
        return EXIT_FAILURE;
      }
    }
  catch( const std::exception & exception )
    {
    tube::ErrorMessage( exception.what() );
    return EXIT_FAILURE;
    }
  return EXIT_FAILURE;
}


template< unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  const unsigned int Dimension = VDimension;

  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO( 
      inputImageFile.c_str(), itk::ImageIOFactory::ReadMode );

  imageIO->SetFileName( inputImageFile.c_str() );
  imageIO->ReadImageInformation();

  const ScalarPixelType pixelType = imageIO->GetComponentType();

  switch( pixelType )
    {
    case itk::ImageIOBase::SHORT:
      {
      // Define the image type, the reader type and the nearest-neighbor ( NN )
      // interpolator type that might be used later on.
      typedef itk::Image< short, Dimension >                ImageType;
      typedef itk::ImageFileReader< ImageType >             ReaderType;
      typedef itk::NearestNeighborInterpolateImageFunction<
        ImageType, double >                                 NNInterpolatorType;


      // Read input file
      typename ReaderType::Pointer srcReader = ReaderType::New();
      srcReader->SetFileName( inputImageFile.c_str() );
      srcReader->Update();


      // First, we run an orientation correction filter that ensures that we
      // have a RAI coordinate system.
      typename itk::OrientImageFilter< ImageType, ImageType >::Pointer
        orientationFilter =
        itk::OrientImageFilter<ImageType,ImageType>::New();
      orientationFilter->UseImageDirectionOn();
      orientationFilter->SetDesiredCoordinateOrientation( 
        itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI );
      orientationFilter->SetInput( srcReader->GetOutput() );
      orientationFilter->Update();


      // Get meta-information about the newly oriented image
      typename ImageType::Pointer orientedImage =
        orientationFilter->GetOutput();
      vnl_vector< double > orientedImageOrigin =
        orientedImage->GetOrigin().GetVnlVector();
      typename ImageType::SizeType orientedImageSize =
        orientedImage->GetLargestPossibleRegion().GetSize();


      // Next, we read the reference image. We will use the geometric center
      // of the volume to compute the translation for the input image.
      typename ReaderType::Pointer refReader = ReaderType::New();
      refReader->SetFileName( referenceImageFile.c_str() );
      refReader->Update();
      typename ImageType::Pointer refImage = refReader->GetOutput();
      typename ImageType::SizeType refImageSize =
        refImage->GetLargestPossibleRegion().GetSize();


      // The translation that we use will be a translation of the geometric
      // center of the input image to the geometric center of the reference
      // image.
      typedef itk::TranslationTransform< double, Dimension >
        TranslationTransformType;
      typename TranslationTransformType::Pointer transform =
        TranslationTransformType::New();
      typename TranslationTransformType::OutputVectorType translation;

      for( unsigned int dim=0; dim<Dimension; ++dim )
        {
        if( centerOnZero )
          {
          translation[dim] = orientedImageOrigin[dim];
          }
        else
          {
          translation[dim] = orientedImageOrigin[dim]
            + static_cast<double>( orientedImageSize[dim] )/2.0
            - static_cast<double>( refImageSize[dim] )/2.0;
          }
        }
      transform->Translate( translation );


      // We have to adjust the image origin to reflect the changes induced
      // by the translation transform.
      double outImageOrigin[Dimension];
      for( unsigned int dim=0; dim<Dimension; ++dim )
         {
         outImageOrigin[dim] = orientedImageOrigin[dim]
           - translation[dim];
         }


      // Finally, we run a resampling filter to resample the input image in
      // the new space. We have to take care that, in case of discrete pixel
      // values ( e.g., label map ), we switch from the standard linear
      // interpolator to the NN interpolator.
      typedef itk::ResampleImageFilter<ImageType, ImageType>
        ResampleImageFilterType;
      typename ResampleImageFilterType::Pointer resampleFilter =
        ResampleImageFilterType::New();
      resampleFilter->SetTransform( transform.GetPointer() );
      resampleFilter->SetInput( orientedImage );
      resampleFilter->SetSize( orientedImageSize );
      resampleFilter->SetOutputOrigin( outImageOrigin );
      resampleFilter->SetOutputDirection( orientedImage->GetDirection() );

      if( isLabelMap )
        {
        tube::InfoMessage( "Using NN interpolation!" );
        typename NNInterpolatorType::Pointer nnInterpolator =
          NNInterpolatorType::New();
        resampleFilter->SetInterpolator( nnInterpolator );
        resampleFilter->Update();
        }
      else
        {
        resampleFilter->Update();
        }


      // Write the output image
      typedef itk::ImageFileWriter< ImageType > WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( outputImageFile.c_str() );
      writer->SetInput( resampleFilter->GetOutput() );
      writer->Update();
      break;
      }
    default:
      {
      tube::ErrorMessage( "Only SHORT type supported!" );
      return EXIT_FAILURE;
      }
    }
  return EXIT_SUCCESS;
}
