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

#include <cmath>
#include <iostream>

#include "itkImage.h"
#include "itkImageIOBase.h"
#include "itkImageIOFactory.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkMetaImageIO.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNormalVariateGenerator.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkErodeObjectMorphologyImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkDilateObjectMorphologyImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkNormalizeImageFilter.h"
#include "metaCommand.h"

// tubetk includes
#include "itkTubeRidgeExtractor.h"
#include "itkTubeCVTImageFilter.h"

/** Resamples image a to b if they are different, returns resampled_a */
template< class pixelT, unsigned int dimensionT >
typename itk::Image< pixelT, dimensionT >::Pointer
ResampleImage(
  typename itk::Image< pixelT, dimensionT >::Pointer a,
  typename itk::Image< pixelT, dimensionT >::Pointer b )
{
  typedef itk::Image< pixelT, dimensionT >      ImageType;

  typename ImageType::Pointer output = a;

  for( unsigned int i = 0; i < dimensionT; i++ )
    {
    if( a->GetLargestPossibleRegion().GetSize()[i]
          != b->GetLargestPossibleRegion().GetSize()[i]
        || a->GetLargestPossibleRegion().GetIndex()[i]
            != b->GetLargestPossibleRegion().GetIndex()[i]
        || a->GetSpacing()[i] != b->GetSpacing()[i]
        || a->GetOrigin()[i] != b->GetOrigin()[i]  )
      {
      typedef typename itk::ResampleImageFilter< ImageType,
                ImageType> ResampleFilterType;
      typename ResampleFilterType::Pointer filter =
        ResampleFilterType::New();
      filter->SetInput( a );
      filter->SetSize( b->GetLargestPossibleRegion().GetSize() );
      filter->SetOutputOrigin( b->GetOrigin() );
      filter->SetOutputSpacing( b->GetSpacing() );
      filter->SetDefaultPixelValue( 0 );
      filter->Update();
      output = filter->GetOutput();
      }
    }
  return output;
}

/** Main command */
template< class pixelT, unsigned int dimensionT >
int DoIt( MetaCommand & command )
{
  typedef float                                    PixelType;
  typedef itk::Image< PixelType, dimensionT >      ImageType;
  typedef itk::Image< unsigned char, dimensionT >  ImageTypeUChar;
  typedef itk::Image< unsigned short, dimensionT > ImageTypeUShort;
  typedef itk::Image< short, dimensionT >          ImageTypeShort;

  MetaCommand::OptionVector parsed = command.GetParsedOptions();

  typedef itk::Statistics::NormalVariateGenerator gaussGenType;
  typename gaussGenType::Pointer gaussGen = gaussGenType::New();

  typedef itk::ImageFileReader< ImageType >   VolumeReaderType;

  // Declare a reader
  typename VolumeReaderType::Pointer reader = VolumeReaderType::New();
  reader->SetFileName( command.GetValueAsString( "infile" ).c_str() );
  typename ImageType::Pointer imIn;
  imIn = reader->GetOutput();

  // See if the file can be read - "try" otherwise program will
  //   mysteriously exit on failure in the Object factory
  std::cout << "Reading file ( "
            << command.GetValueAsString( "infile" ).c_str()
            <<" )"<< std::endl;
  try
    {
    reader->Update();
    }
  catch( ... )
    {
    std::cout << "Problems reading file format" << std::endl;
    return EXIT_FAILURE;
    }

  MetaCommand::OptionVector::const_iterator it = parsed.begin();
  while( it != parsed.end() )
    {
    if( ( *it ).name == "Write" )
      {
      std::string outFilename =
        command.GetValueAsString( *it, "filename" );
      std::cout << "Writing output1 ( " << outFilename.c_str()
                << " )" << std::endl;
      typedef itk::ImageFileWriter< ImageType > VolumeWriterType;
      typename VolumeWriterType::Pointer writer = VolumeWriterType::New();
      writer->SetFileName( outFilename.c_str() );
      writer->SetInput( imIn );
      writer->SetUseCompression( true );
      try
        {
        writer->Write();
        }
      catch( itk::ExceptionObject& err )
        {
        std::cout << "WriteImage : " << err << std::endl;
        return EXIT_FAILURE;
        }
      } // end -w

    else if( ( *it ).name == "WriteType" )
      {
      int type = command.GetValueAsInt( "WriteType", "Type" );
      std::string outFilename =
        command.GetValueAsString( *it, "filename" );
      std::cout << "Writing output2 ( " << outFilename.c_str()
                << " )" << std::endl;
      switch( type )
        {
        case 0:
        case 4:
          {
          typedef itk::CastImageFilter< ImageType, ImageTypeUChar>
            CastFilterType;
          typename CastFilterType::Pointer castFilter =
            CastFilterType::New();
          castFilter->SetInput( imIn );

          typedef itk::ImageFileWriter< ImageTypeUChar >
            VolumeWriterType;
          typename VolumeWriterType::Pointer writer =
            VolumeWriterType::New();
          writer->SetFileName( outFilename.c_str() );
          writer->SetInput( castFilter->GetOutput() );
          if( type == 0 )
            {
            writer->SetUseCompression( true );
            }
          writer->Write();
          break;
          }
        case 1:
        case 5:
          {
          typedef itk::CastImageFilter< ImageType, ImageTypeUShort>
            CastFilterType;
          typename CastFilterType::Pointer castFilter =
            CastFilterType::New();
          castFilter->SetInput( imIn );

          typedef itk::ImageFileWriter< ImageTypeUShort >
            VolumeWriterType;
          typename VolumeWriterType::Pointer writer =
            VolumeWriterType::New();
          writer->SetFileName( outFilename.c_str() );
          writer->SetInput( castFilter->GetOutput() );
          if( type == 1 )
            {
            writer->SetUseCompression( true );
            }
          writer->Write();
          break;
          }
        case 2:
        case 6:
          {
          typedef itk::CastImageFilter< ImageType, ImageTypeShort>
            CastFilterType;
          typename CastFilterType::Pointer castFilter =
            CastFilterType::New();
          castFilter->SetInput( imIn );

          typedef itk::ImageFileWriter< ImageTypeShort >
            VolumeWriterType;
          typename VolumeWriterType::Pointer writer =
            VolumeWriterType::New();
          writer->SetFileName( outFilename.c_str() );
          writer->SetInput( castFilter->GetOutput() );
          if( type == 2 )
            {
            writer->SetUseCompression( true );
            }
          writer->Write();
          break;
          }
        case 3:
          {
          typedef itk::CastImageFilter< ImageType, ImageTypeShort>
            CastFilterType;
          typename CastFilterType::Pointer castFilter =
            CastFilterType::New();
          castFilter->SetInput( imIn );

          typedef itk::ImageFileWriter< ImageTypeShort >
            VolumeWriterType;
          typename VolumeWriterType::Pointer writer =
            VolumeWriterType::New();

          itk::MetaImageIO::Pointer metaWriter = itk::MetaImageIO::New();
          writer->SetImageIO( metaWriter );

          writer->SetFileName( outFilename.c_str() );
          writer->SetInput( castFilter->GetOutput() );
          writer->SetUseCompression( false );

          MetaImage * metaImage = metaWriter->GetMetaImagePointer();

          metaImage->ElementSize( 0, imIn->GetSpacing()[0] );
          metaImage->ElementSize( 1, imIn->GetSpacing()[1] );
          metaImage->ElementSize( 2, imIn->GetSpacing()[2] );

          metaImage->AddUserField( "ElementByteOrderMSB",
                                  MET_STRING, strlen( "False" ), "False" );

          writer->Write();
          break;
          }
        case 7:
          {
          typedef itk::ImageFileWriter< ImageType > VolumeWriterType;
          typename VolumeWriterType::Pointer writer = VolumeWriterType::New();
          writer->SetFileName( outFilename.c_str() );
          writer->SetInput( imIn );
          try
            {
            writer->Write();
            }
          catch( itk::ExceptionObject& err )
            {
            std::cout << "WriteImage : " << err << std::endl;
            return EXIT_FAILURE;
            }
          }
        }
      } // end -W

    else if( ( *it ).name == "Intensity" )
      {
      std::cout << "Intensity windowing" << std::endl;
      float valMin = command.GetValueAsFloat( *it, "inValMin" );
      float valMax = command.GetValueAsFloat( *it, "inValMax" );
      float outMin = command.GetValueAsFloat( *it, "outMin" );
      float outMax = command.GetValueAsFloat( *it, "outMax" );
      itk::ImageRegionIterator< ImageType > it2( imIn,
            imIn->GetLargestPossibleRegion() );
      it2.GoToBegin();
      while( !it2.IsAtEnd() )
        {
        double tf = it2.Get();
        tf = ( tf-valMin )/( valMax-valMin );
        if( tf<0 )
          {
          tf = 0;
          }
        if( tf>1 )
          {
          tf = 1;
          }
        tf = ( tf * ( outMax-outMin ) ) + outMin;
        it2.Set( ( PixelType )tf );
        ++it2;
        }
      }

    // IntensityMult
    else if( ( *it ).name == "IntensityMult" )
      {
      std::cout << "Intensity multiplicative bias correct" << std::endl;
      typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
      reader2->SetFileName(
        command.GetValueAsString( *it, "inMeanField" ).c_str() );
      typename ImageType::Pointer imIn2;
      imIn2 = reader2->GetOutput();
      try
        {
        reader2->Update();
        }
      catch( ... )
        {
        std::cout << "Problems reading file format of inFile2."
                  << std::endl;
        return EXIT_FAILURE;
        }
      imIn2 = ResampleImage< PixelType, dimensionT >( imIn2, imIn );
      itk::ImageRegionIterator< ImageType > it2( imIn2,
                     imIn2->GetLargestPossibleRegion() );
      int count = 0;
      double mean = 0;
      it2.GoToBegin();
      while( !it2.IsAtEnd() )
        {
        double tf = it2.Get();
        mean += tf;
        if( tf != 0 )
          {
          ++count;
          }
        ++it2;
        }
      mean /= count;
      itk::ImageRegionIterator< ImageType > it3( imIn,
            imIn->GetLargestPossibleRegion() );
      it3.GoToBegin();
      it2.GoToBegin();
      while( !it3.IsAtEnd() )
        {
        double tf = it3.Get();
        double tf2 = it2.Get();
        if( tf2 != 0 )
          {
          double alpha = mean / tf2;
          tf = tf * alpha;
          it3.Set( ( PixelType )tf );
          }
        ++it3;
        ++it2;
        }
      } // end -I

    // UniformNoise
    else if( ( *it ).name == "UniformNoise" )
      {
      std::cout << "Adding noise" << std::endl;
      float valMin = command.GetValueAsFloat( *it, "inValMin" );
      float valMax = command.GetValueAsFloat( *it, "inValMax" );
      float noiseMean = command.GetValueAsFloat( *it, "noiseMean" );
      float noiseRange = command.GetValueAsFloat( *it, "noiseRange" );
      itk::ImageRegionIterator< ImageType > it2( imIn,
            imIn->GetLargestPossibleRegion() );
      it2.GoToBegin();
      while( !it2.IsAtEnd() )
        {
        double tf = it2.Get();
        if( tf >= valMin && tf <= valMax )
          {
          tf += ( 2.0*( rand()/( double )RAND_MAX )-1 ) * noiseRange
                + noiseMean;
          it2.Set( ( PixelType )tf );
          }
        ++it2;
        }
      } // -N

    // GaussianNoise
    else if( ( *it ).name == "GaussianNoise" )
      {
      std::cout << "Adding noise" << std::endl;
      float valMin = command.GetValueAsFloat( *it, "inValMin" );
      float valMax = command.GetValueAsFloat( *it, "inValMax" );
      float noiseMean = command.GetValueAsFloat( *it, "noiseMean" );
      float noiseStdDev = command.GetValueAsFloat( *it, "noiseStdDev" );
      itk::ImageRegionIterator< ImageType > it2( imIn,
            imIn->GetLargestPossibleRegion() );
      it2.GoToBegin();
      while( !it2.IsAtEnd() )
        {
        double tf = it2.Get();
        if( tf >= valMin && tf <= valMax )
          {
          tf += gaussGen->GetVariate()*noiseStdDev+noiseMean;
          it2.Set( ( PixelType )tf );
          }
        ++it2;
        }
      } // end -n

    // I( x )
    else if( ( *it ).name == "Add" )
      {
      std::cout << "Adding" << std::endl;
      float weight1 = command.GetValueAsFloat( *it, "weight1" );
      float weight2 = command.GetValueAsFloat( *it, "weight2" );
      typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
      reader2->SetFileName( command.GetValueAsString( *it,
          "Infile" ).c_str() );
      typename ImageType::Pointer imIn2;
      imIn2 = reader2->GetOutput();
      try
        {
        reader2->Update();
        }
      catch( ... )
        {
        std::cout << "Problems reading file format of inFile2."
                  << std::endl;
        return EXIT_FAILURE;
        }
      itk::ImageRegionIterator< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType > it2( imIn2,
            imIn2->GetLargestPossibleRegion() );
      it1.GoToBegin();
      it2.GoToBegin();
      while( !it1.IsAtEnd() )
        {
        double tf1 = it1.Get();
        double tf2 = it2.Get();
        double tf = weight1*tf1 + weight2*tf2;
        it1.Set( ( PixelType )tf );
        ++it1;
        ++it2;
        }
      }

    else if( ( *it ).name == "Multiply" )
      {
      std::cout << "Multiplying" << std::endl;
      typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
      reader2->SetFileName( command.GetValueAsString( *it,
          "Infile" ).c_str() );
      typename ImageType::Pointer imIn2;
      imIn2 = reader2->GetOutput();
      try
        {
        reader2->Update();
        }
      catch( ... )
        {
        std::cout << "Problems reading file format of inFile2."
                  << std::endl;
        return EXIT_FAILURE;
        }
      itk::ImageRegionIterator< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType > it2( imIn2,
            imIn2->GetLargestPossibleRegion() );
      it1.GoToBegin();
      it2.GoToBegin();
      while( !it1.IsAtEnd() )
        {
        double tf1 = it1.Get();
        double tf2 = it2.Get();
        double tf = tf1*tf2;
        it1.Set( ( PixelType )tf );
        ++it1;
        ++it2;
        }
      }

    // Normalize
    else if( ( *it ).name == "Normalize" )
      {
      std::cout << "Normalize" << std::endl;
      std::cout << "NOTE: since this filter normalizes the data to lie "
                << "within -1 to 1, integral types will produce an image that "
                << "DOES NOT HAVE a unit variance" << std::endl;
      typedef itk::NormalizeImageFilter< ImageType, ImageType > NormFilterType;
      typename NormFilterType::Pointer normFilter = NormFilterType::New();
      normFilter->SetInput( imIn );
      normFilter->Update();
      imIn = normFilter->GetOutput();
      }

    // I( x )
    else if( ( *it ).name == "Fuse" )
      {
      std::cout << "Fusing" << std::endl;
      float offset2 = command.GetValueAsFloat( *it, "Offset2" );
      typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
      reader2->SetFileName( command.GetValueAsString( *it,
          "Infile2" ).c_str() );
      typename ImageType::Pointer imIn2;
      imIn2 = reader2->GetOutput();
      try
        {
        reader2->Update();
        }
      catch( ... )
        {
        std::cout << "Problems reading file format of inFile2."
                  << std::endl;
        return EXIT_FAILURE;
        }
      itk::ImageRegionIterator< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType > it2( imIn2,
            imIn2->GetLargestPossibleRegion() );
      it1.GoToBegin();
      it2.GoToBegin();
      while( !it1.IsAtEnd() )
        {
        double tf1 = it1.Get();
        double tf2 = it2.Get();
        if( tf2>tf1 )
          {
          double tf = offset2 + tf2;
          it1.Set( ( PixelType )tf );
          }
        ++it1;
        ++it2;
        }
      } // end -a

    // Threshold
    else if( ( *it ).name == "Threshold" )
      {
      std::cout << "Thresholding" << std::endl;

      float threshLow = command.GetValueAsFloat( *it, "threshLow" );
      float threshHigh = command.GetValueAsFloat( *it, "threshHigh" );

      float valTrue = command.GetValueAsFloat( *it, "valTrue" );
      float valFalse = command.GetValueAsFloat( *it, "valFalse" );

      itk::ImageRegionIterator< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      it1.GoToBegin();
      while( !it1.IsAtEnd() )
        {
        double tf1 = it1.Get();
        if( tf1 >= threshLow && tf1 <= threshHigh )
          {
          it1.Set( ( PixelType )valTrue );
          }
        else
          {
          it1.Set( ( PixelType )valFalse );
          }
        ++it1;
        }
      }
    else if( ( *it ).name == "Algorithm" )
      {
      std::cout << "Algorithm" << std::endl;

      float threshLow = command.GetValueAsFloat( *it, "threshLow" );
      float threshHigh = command.GetValueAsFloat( *it, "threshHigh" );

      typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
      reader2->SetFileName( command.GetValueAsString( *it,
          "maskFile" ).c_str() );
      typename ImageType::Pointer imIn2 = reader2->GetOutput();
      try
        {
        reader2->Update();
        }
      catch( ... )
        {
        std::cout << "Problems reading file format of inFile2."
                  << std::endl;
        return EXIT_FAILURE;
        }

      int mode = command.GetValueAsInt( *it, "mode" );

      itk::ImageRegionIterator< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType > it2( imIn2,
            imIn2->GetLargestPossibleRegion() );
      it1.GoToBegin();
      it2.GoToBegin();
      double sum = 0;
      double sumS = 0;
      unsigned int count = 0;
      while( !it1.IsAtEnd() && !it2.IsAtEnd() )
        {
        double maskV = it2.Get();
        if( maskV >= threshLow && maskV <= threshHigh )
          {
          sum += it1.Get();
          sumS += it1.Get() * it1.Get();
          ++count;
          }
        ++it1;
        ++it2;
        }
      double mean = sum/count;
      if( mode == 0 )
        {
        std::cout << "Mean " << mean << std::endl;
        }
      else
        {
        double stdDev = (sumS - (sum*mean))/(count-1);
        std::cout << "StdDev " << stdDev << std::endl;
        }
      }
    else if( ( *it ).name == "Process" )
      {
      std::cout << "Process binary operation" << std::endl;

      typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
      reader2->SetFileName( command.GetValueAsString( *it,
          "file2" ).c_str() );
      typename ImageType::Pointer imIn2 = reader2->GetOutput();
      try
        {
        reader2->Update();
        }
      catch( ... )
        {
        std::cout << "Problems reading file format of inFile2."
                  << std::endl;
        return EXIT_FAILURE;
        }

      int mode = command.GetValueAsInt( *it, "mode" );

      itk::ImageRegionIterator< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType > it2( imIn2,
            imIn2->GetLargestPossibleRegion() );
      it1.GoToBegin();
      it2.GoToBegin();
      if( mode == 0 )
        {
        while( !it1.IsAtEnd() && !it2.IsAtEnd() )
          {
          it1.Set( it1.Get() * it2.Get() );
          ++it1;
          ++it2;
          }
        }
      }
    else if( ( *it ).name == "process" )
      {
      std::cout << "Unary process" << std::endl;

      int mode = command.GetValueAsInt( *it, "mode" );

      itk::ImageRegionIterator< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      it1.GoToBegin();
      if( mode == 0 )
        {
        while( !it1.IsAtEnd() )
          {
          it1.Set( vnl_math_abs( it1.Get() ) );
          ++it1;
          }
        }
      }
    // Masking
    else if( ( *it ).name == "Masking" )
      {
      std::cout << "Masking" << std::endl;

      float threshLow = command.GetValueAsFloat( *it, "threshLow" );
      float threshHigh = command.GetValueAsFloat( *it, "threshHigh" );

      typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
      reader2->SetFileName( command.GetValueAsString( *it,
          "inFile2" ).c_str() );

      float valFalse = command.GetValueAsFloat( *it, "valFalse" );

      typename ImageType::Pointer imIn2 = reader2->GetOutput();
      try
        {
        reader2->Update();
        }
      catch( ... )
        {
        std::cout << "Problems reading file format of inFile2."
                  << std::endl;
        return EXIT_FAILURE;
        }
      imIn2 = ResampleImage< PixelType, dimensionT >( imIn2, imIn );
      itk::ImageRegionIterator< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType > it2( imIn2,
            imIn2->GetLargestPossibleRegion() );
      it1.GoToBegin();
      it2.GoToBegin();
      while( !it1.IsAtEnd() )
        {
        double tf2 = it2.Get();
        if( tf2 >= threshLow && tf2 <= threshHigh )
          {
          //double tf1 = it1.Get();
          //it1.Set( ( PixelType )tf1 );
          }
        else
          {
          it1.Set( ( PixelType )valFalse );
          }
        ++it1;
        ++it2;
        }
      }

    // Morphology
    else if( ( *it ).name == "Morphology" )
      {
      std::cout << "Morphology" << std::endl;

      int mode = command.GetValueAsInt( *it, "mode" );

      float radius = command.GetValueAsFloat( *it, "radius" );

      float foregroundValue = command.GetValueAsFloat( *it,
        "forgroundValue" );
      float backgroundValue = command.GetValueAsFloat( *it,
        "backgroundValue" );

      typedef itk::BinaryBallStructuringElement<PixelType, dimensionT>
        BallType;
      BallType ball;
      ball.SetRadius( 1 );
      ball.CreateStructuringElement();

      typedef itk::ErodeObjectMorphologyImageFilter
                   <ImageType, ImageType, BallType>       ErodeFilterType;
      typedef itk::DilateObjectMorphologyImageFilter
                   <ImageType, ImageType, BallType>       DilateFilterType;
      switch( mode )
        {
        case 0:
          {
          for( int r=0; r<radius; r++ )
            {
            typename ErodeFilterType::Pointer filter =
              ErodeFilterType::New();
            filter->SetBackgroundValue( backgroundValue );
            filter->SetKernel( ball );
            filter->SetObjectValue( foregroundValue );
            filter->SetInput( imIn );
            filter->Update();
            imIn = filter->GetOutput();
            }
          break;
          }
        case 1:
          {
          for( int r=0; r<radius; r++ )
            {
            typename DilateFilterType::Pointer filter =
              DilateFilterType::New();
            filter->SetKernel( ball );
            filter->SetObjectValue( foregroundValue );
            filter->SetInput( imIn );
            filter->Update();
            imIn = filter->GetOutput();
            }
          break;
          }
        }
      }

    else if( ( *it ).name == "overwrite" )
      {
      typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
      reader2->SetFileName( command.GetValueAsString( *it,
          "mask" ).c_str() );
      reader2->Update();
      typename ImageType::Pointer maskIm = reader2->GetOutput();

      float maskKeyVal = command.GetValueAsFloat( *it, "maskKeyVal" );
      float imageKeyVal = command.GetValueAsFloat( *it, "imageKeyVal" );
      float newImageVal = command.GetValueAsFloat( *it, "newImageVal" );

      itk::ImageRegionIterator< ImageType > itIm( imIn,
            imIn->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType > itMask( maskIm,
            maskIm->GetLargestPossibleRegion() );
      while( !itIm.IsAtEnd() )
        {
        if( itMask.Get() == maskKeyVal )
          {
          if( itIm.Get() == imageKeyVal )
            {
            itIm.Set( newImageVal );
            }
          }
        ++itIm;
        ++itMask;
        }

      std::cout << "Overwrite" << std::endl;
      }

    // blur
    else if( ( *it ).name == "blur" )
      {
      std::cout << "Blurring." << std::endl;
      float sigma = command.GetValueAsFloat( *it, "sigma" );

      typename itk::RecursiveGaussianImageFilter< ImageType >::Pointer
        filter;
      typename ImageType::Pointer imTemp;
      for( unsigned int i=0; i<dimensionT; i++ )
        {
        filter = itk::RecursiveGaussianImageFilter< ImageType >::New();
        filter->SetInput( imIn );
        //filter->SetNormalizeAcrossScale( true );
        filter->SetSigma( sigma );

        filter->SetOrder(
                 itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
        filter->SetDirection( i );

        imTemp = filter->GetOutput();
        filter->Update();
        imIn = imTemp;
        }
      }

    // blurOrder
    else if( ( *it ).name == "blurOrder" )
      {
      std::cout << "Blurring." << std::endl;

      float sigma = command.GetValueAsFloat( *it, "sigma" );
      int order = command.GetValueAsInt( *it, "order" );
      int direction = command.GetValueAsInt( *it, "direction" );

      typename itk::RecursiveGaussianImageFilter< ImageType >::Pointer
        filter;
      filter = itk::RecursiveGaussianImageFilter< ImageType >::New();
      filter->SetInput( imIn );
      //filter->SetNormalizeAcrossScale( true );
      filter->SetSigma( sigma );
      filter->SetDirection( direction );
      switch( order )
        {
        case 0:
          filter->SetOrder(
            itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
          break;
        case 1:
          filter->SetOrder(
            itk::RecursiveGaussianImageFilter<ImageType>::FirstOrder );
          break;
        case 2:
          filter->SetOrder(
            itk::RecursiveGaussianImageFilter<ImageType>::SecondOrder );
          break;
        }
      typename ImageType::Pointer imTemp = filter->GetOutput();
      filter->Update();
      imIn = imTemp;
      } // end -B

    // histogram
    else if( ( *it ).name == "histogram" )
      {
      std::cout << "Histogram" << std::endl;

      unsigned int nBins = ( unsigned int )command.GetValueAsInt( *it,
        "nBins" );
      std::string filename =
        command.GetValueAsString( *it, "histOutputFile" );

      itk::ImageRegionIteratorWithIndex< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      it1.GoToBegin();
      double binMin = it1.Get();
      double binMax = it1.Get();
      while( !it1.IsAtEnd() )
        {
        double tf = it1.Get();
        if( tf < binMin )
          {
          binMin = tf;
          }
        else
          {
          if( tf > binMax )
            {
            binMax = tf;
            }
          }
        ++it1;
        }
      std::cout << "  binMin = " << binMin << std::endl;
      std::cout << "  binMax = " << binMax << std::endl;
      it1.GoToBegin();
      itk::Array<double> bin;
      bin.set_size( nBins );
      bin.Fill( 0 );
      while( !it1.IsAtEnd() )
        {
        double tf = it1.Get();
        tf = ( tf-binMin )/( binMax-binMin ) * nBins;
        if( tf>nBins-1 )
          {
          tf = nBins-1;
          }
        else
          {
          if( tf<0 )
            {
            tf = 0;
            }
          }
        bin[( int )tf]++;
        ++it1;
        }
      std::ofstream writeStream;
      writeStream.open( filename.c_str(), std::ios::binary | std::ios::out );
      if( ! writeStream.rdbuf()->is_open() )
        {
        std::cerr << "Cannot write to file : " << filename << std::endl;
        return EXIT_FAILURE;
        }
      for( unsigned int i=0; i<nBins; i++ )
        {
        writeStream << ( i/( double )nBins )*( binMax-binMin )+binMin
                    << " " << bin[i] << std::endl;
        }
      writeStream.close();
      }

    // histogram2
    else if( ( *it ).name == "histogram2" )
      {
      std::cout << "Histogram" << std::endl;

      unsigned int nBins = ( unsigned int )command.GetValueAsInt( *it,
        "nBins" );
      double binMin = command.GetValueAsFloat( *it, "binMin" );
      double binSize = command.GetValueAsFloat( *it, "binSIZE" );
      double binMax = binMin + binSize*nBins;
      std::string filename =
        command.GetValueAsString( *it, "histOutputFile" );

      itk::ImageRegionIteratorWithIndex< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      it1.GoToBegin();
      itk::Array<double> bin;
      bin.set_size( nBins );
      bin.Fill( 0 );
      while( !it1.IsAtEnd() )
        {
        double tf = it1.Get();
        tf = ( tf-binMin )/( binMax-binMin ) * nBins;
        if( tf<nBins && tf>0 )
          {
          bin[( int )tf]++;
          }
        ++it1;
        }
      std::ofstream writeStream;
      writeStream.open( filename.c_str(), std::ios::binary | std::ios::out );
      if( ! writeStream.rdbuf()->is_open() )
        {
        std::cerr << "Cannot write to file : " << filename << std::endl;
        return EXIT_FAILURE;
        }
      for( unsigned int i=0; i<nBins; i++ )
        {
        writeStream << ( i/( double )nBins )*( binMax-binMin )+binMin
                    << " " << bin[i] << std::endl;
        }
      writeStream.close();
      }

    // vessels
    else if( ( *it ).name == "vessels" )
      {
      std::cout << "Vessel Enhancement" << std::endl;

      double scaleMin = command.GetValueAsFloat( *it, "scaleMin" );
      double scaleMax = command.GetValueAsFloat( *it, "scaleMax" );
      double numScales = command.GetValueAsFloat( *it, "numScales" );
      double logScaleStep = (vcl_log(scaleMax) - vcl_log(scaleMin))
        / (numScales-1);

      typedef itk::tube::RidgeExtractor< ImageType > RidgeFuncType;
      typename RidgeFuncType::Pointer imFunc = RidgeFuncType::New();
      imFunc->SetInputImage( imIn );

      typename ImageType::Pointer imIn2 = ImageType::New();
      imIn2->SetRegions( imIn->GetLargestPossibleRegion() );
      imIn2->SetOrigin( imIn->GetOrigin() );
      imIn2->SetSpacing( imIn->GetSpacing() );
      imIn2->Allocate();

      itk::ImageRegionIteratorWithIndex< ImageType > it1( imIn,
            imIn->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType > it2( imIn2,
            imIn2->GetLargestPossibleRegion() );

      double ridgeness = 0;
      double roundness = 0;
      double curvature = 0;
      double scale = scaleMin;
      imFunc->SetScale( scale );
      std::cout << "   Processing scale " << scale << std::endl;
      it1.GoToBegin();
      it2.GoToBegin();
      typename RidgeFuncType::ContinuousIndexType cIndx;
      while( !it1.IsAtEnd() )
        {
        for( unsigned int d=0; d<ImageType::ImageDimension; ++d )
          {
          cIndx[d] = it1.GetIndex()[d];
          }
        ridgeness = imFunc->Ridgeness( cIndx, roundness, curvature  );
        it2.Set( ( PixelType )ridgeness );
        ++it1;
        ++it2;
        }
      for( unsigned int i=1; i<numScales; i++ )
        {
        scale = vcl_exp(vcl_log(scaleMin) + i * logScaleStep);
        imFunc->SetScale( scale );
        std::cout << "   Processing scale " << scale << std::endl;
        it1.GoToBegin();
        it2.GoToBegin();
        while( !it1.IsAtEnd() )
          {
          for( unsigned int d=0; d<ImageType::ImageDimension; ++d )
            {
            cIndx[d] = it1.GetIndex()[d];
            }
          ridgeness = imFunc->Ridgeness( cIndx, roundness, curvature  );
          if( ridgeness > it2.Get() )
            {
            it2.Set( ( PixelType )ridgeness );
            }
          ++it1;
          ++it2;
          }
        }
      it1.GoToBegin();
      it2.GoToBegin();
      while( !it1.IsAtEnd() )
        {
        it1.Set( it2.Get() );
        ++it1;
        ++it2;
        }
      }

    // CorrectionSlice
    else if( ( *it ).name == "CorrectionSlice" )
      {
      std::cout << "Correct intensity slice-by-slice" << std::endl;

      unsigned int numberOfBins =
        ( unsigned int )command.GetValueAsInt( *it, "nBins" );
      unsigned int numberOfMatchPoints =
        ( unsigned int )command.GetValueAsInt( *it, "nMatchPoints" );
      typedef itk::Image<PixelType, 2> ImageType2D;
      typedef itk::HistogramMatchingImageFilter< ImageType2D, ImageType2D >
          HistogramMatchFilterType;
      typename HistogramMatchFilterType::Pointer matchFilter;
      typename ImageType2D::Pointer im2DRef = ImageType2D::New();
      typename ImageType2D::Pointer im2DIn = ImageType2D::New();
      typename ImageType2D::SizeType size2D;
      size2D[0] = imIn->GetLargestPossibleRegion().GetSize()[0];
      size2D[1] = imIn->GetLargestPossibleRegion().GetSize()[1];
      im2DRef->SetRegions( size2D );
      im2DRef->Allocate();
      im2DIn->SetRegions( size2D );
      im2DIn->Allocate();
      itk::ImageRegionIterator< ImageType > it3D( imIn,
            imIn->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType > it3DSliceStart( imIn,
            imIn->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType2D > it2DRef( im2DRef,
            im2DRef->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< ImageType2D > it2DIn( im2DIn,
            im2DIn->GetLargestPossibleRegion() );
      unsigned int z, y, x;
      it3D.GoToBegin();
      unsigned int zMax = 1;
      if( dimensionT == 3 )
        {
        zMax = imIn->GetLargestPossibleRegion().GetSize()[2];
        }
      for( z=0; z<dimensionT && z<zMax; z++ )
        {
        it2DRef.GoToBegin();
        for( y=0; y<imIn->GetLargestPossibleRegion().GetSize()[1]; y++ )
          {
          for( x=0; x<imIn->GetLargestPossibleRegion().GetSize()[0]; x++ )
            {
            it2DRef.Set( it3D.Get() );
            ++it2DRef;
            ++it3D;
            }
          }
        }
      for(; z<zMax; z++ )
        {
        it2DIn.GoToBegin();
        it3DSliceStart = it3D;
        for( y=0; y<imIn->GetLargestPossibleRegion().GetSize()[1]; y++ )
          {
          for( x=0; x<imIn->GetLargestPossibleRegion().GetSize()[0]; x++ )
            {
            it2DIn.Set( it3D.Get() );
            ++it2DIn;
            ++it3D;
            }
          }
        matchFilter = HistogramMatchFilterType::New();
        matchFilter->SetReferenceImage( im2DRef );
        matchFilter->SetInput( im2DIn );
        matchFilter->SetNumberOfHistogramLevels( numberOfBins );
        matchFilter->SetNumberOfMatchPoints( numberOfMatchPoints );
        matchFilter->Update();
        itk::ImageRegionIterator< ImageType2D > it2DOut(
              matchFilter->GetOutput(),
              im2DIn->GetLargestPossibleRegion() );
        it2DRef.GoToBegin();
        it2DOut.GoToBegin();
        it3D = it3DSliceStart;
        for( y=0; y<imIn->GetLargestPossibleRegion().GetSize()[1]; y++ )
          {
          for( x=0; x<imIn->GetLargestPossibleRegion().GetSize()[0]; x++ )
            {
            it2DRef.Set( it2DOut.Get() );
            it3D.Set( it2DOut.Get() );
            ++it2DRef;
            ++it2DOut;
            ++it3D;
            }
          }
        }
      }

    // Correction
    else if( ( *it ).name == "Correction" )
      {
      std::cout << "Correct intensity in the volume" << std::endl;

      unsigned int numberOfBins =
        ( unsigned int )command.GetValueAsInt( *it, "nBins" );
      unsigned int numberOfMatchPoints =
        ( unsigned int )command.GetValueAsInt( *it, "nMatchPoints" );
      typedef itk::HistogramMatchingImageFilter< ImageType, ImageType >
          HistogramMatchFilterType;
      typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
      reader2->SetFileName(
        command.GetValueAsString( *it, "referenceVolume" ).c_str() );
      typename ImageType::Pointer imIn2;
      imIn2 = reader2->GetOutput();
      try
        {
        reader2->Update();
        }
      catch( ... )
        {
        std::cout << "Problems reading file format of inFile2."
                  << std::endl;
        return EXIT_FAILURE;
        }
      typename HistogramMatchFilterType::Pointer matchFilter;
      matchFilter = HistogramMatchFilterType::New();
      matchFilter->SetReferenceImage( imIn2 );
      matchFilter->SetInput( imIn );
      matchFilter->SetNumberOfHistogramLevels( numberOfBins );
      matchFilter->SetNumberOfMatchPoints( numberOfMatchPoints );
      matchFilter->Update();
      imIn = matchFilter->GetOutput();
      }

    // resize
    else if( ( *it ).name == "resize" )
      {
      std::cout << "Resampling." << std::endl;
      double factor = command.GetValueAsFloat( *it, "factor" );

      typename ImageType::Pointer imSub2 = ImageType::New();
      imSub2->SetOrigin( imIn->GetOrigin() );
      typename ImageType::SizeType size;
      typename ImageType::SpacingType spacing;
      if( factor != 0 )
        {
        for( unsigned int i=0; i<dimensionT; i++ )
          {
          size[i] = ( long unsigned int )
                    ( imIn->GetLargestPossibleRegion().GetSize()[i]
                      / factor );
          spacing[i] = imIn->GetSpacing()[i]*factor;
          }
        }
      else
        {
        for( unsigned int i=0; i<dimensionT; i++ )
          {
          spacing[i] = imIn->GetSpacing()[i];
          }

        double meanSpacing = ( spacing[0] + spacing[1] ) / 2;
        if( dimensionT == 3 )
          {
          meanSpacing = ( meanSpacing + spacing[2] ) / 2;
          }
        factor = meanSpacing/spacing[0];
        size[0] = ( long unsigned int )
                  ( imIn->GetLargestPossibleRegion().GetSize()[0]/factor );
        factor = meanSpacing/spacing[1];
        size[1] = ( long unsigned int )
                  ( imIn->GetLargestPossibleRegion().GetSize()[1]/factor );
        factor = meanSpacing/spacing[2];
        spacing[0] = meanSpacing;
        spacing[1] = meanSpacing;
        if( dimensionT == 3 )
          {
          size[2] = ( long unsigned int )
                    ( imIn->GetLargestPossibleRegion().GetSize()[2]
                      / factor );
          spacing[2] = meanSpacing;
          }
        }
      imSub2->SetRegions( size );
      imSub2->SetSpacing( spacing );
      imSub2->Allocate();

      imIn = ResampleImage< PixelType, dimensionT >( imIn, imSub2 );
      }

    // resize2
    else if( ( *it ).name == "resize2" )
      {
      std::cout << "Resampling" << std::endl;
      typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
      reader2->SetFileName( command.GetValueAsString( *it,
          "inFile2" ).c_str() );
      typename ImageType::Pointer imIn2;
      imIn2 = reader2->GetOutput();
      try
        {
        reader2->Update();
        }
      catch( ... )
        {
        std::cout << "Problems reading file format of inFile2."
                  << std::endl;
        return EXIT_FAILURE;
        }
      imIn = ResampleImage< PixelType, dimensionT >( imIn, imIn2 );
      }

    // segment
    else if( ( *it ).name == "segment" )
      {
      std::cout << "Segmenting" << std::endl;

      //int mode = command.GetValueAsInt( *it, "mode" );

      float threshLow = command.GetValueAsFloat( *it, "threshLow" );
      float threshHigh = command.GetValueAsFloat( *it, "threshHigh" );
      float labelValue = command.GetValueAsFloat( *it, "labelValue" );

      float x = command.GetValueAsFloat( *it, "x" );
      float y = command.GetValueAsFloat( *it, "y" );

      typedef itk::ConnectedThresholdImageFilter<ImageType, ImageType>
                 FilterType;
      typename FilterType::Pointer filter = FilterType::New();

      typename ImageType::IndexType seed;
      seed[0] = ( long int )x;
      seed[1] = ( long int )y;

      if( dimensionT == 3 )
        {
        float z = command.GetValueAsFloat( *it, "z" );
        seed[2] = ( long int )z;
        }

      filter->SetInput( imIn );
      filter->SetLower( threshLow );
      filter->SetUpper( threshHigh );
      filter->AddSeed( seed );
      filter->SetReplaceValue( labelValue );
      filter->Update();

      imIn = filter->GetOutput();
      }

    // offset
    else if( ( *it ).name == "offset" )
      {
      double offset[3];
      offset[0] = command.GetValueAsFloat( *it, "offsetX" );
      offset[1] = command.GetValueAsFloat( *it, "offsetY" );
      if( dimensionT == 3 )
        {
        offset[2] = command.GetValueAsFloat( *it, "offsetZ" );
        }
      imIn->SetOrigin( offset );
      }

    // SetRandom
    else if( ( *it ).name == "SetRandom" )
      {
      unsigned int seed = ( unsigned int )command.GetValueAsInt( *it,
        "seedValue" );
      srand( seed );
      gaussGen->Initialize( ( int )seed );
      } // end -S

    // Voronoi
    else if( ( *it ).name == "Voronoi" )
      {
      unsigned int numberOfCentroids =
        ( unsigned int )command.GetValueAsInt( *it, "numCentroids" );
      unsigned int numberOfIterations =
        ( unsigned int )command.GetValueAsInt( *it, "numIters" );
      unsigned int numberOfSamples =
        ( unsigned int )command.GetValueAsInt( *it, "numSamples" );
      std::string filename =
        command.GetValueAsString( *it, "centroidOutFile" );
      typedef itk::tube::CVTImageFilter<ImageType, ImageType> FilterType;
      typename FilterType::Pointer filter = FilterType::New();
      filter->SetInput( imIn );
      filter->SetNumberOfSamples( numberOfSamples );
      filter->SetNumberOfCentroids( numberOfCentroids );
      filter->SetNumberOfIterations( numberOfIterations );
      filter->SetNumberOfSamplesPerBatch( numberOfIterations );
      filter->Update();

      std::ofstream writeStream;
      writeStream.open( filename.c_str(),
        std::ios::binary | std::ios::out );
      if( ! writeStream.rdbuf()->is_open() )
        {
        std::cerr << "Cannot write to file : " << filename << std::endl;
        return EXIT_FAILURE;
        }
      writeStream << numberOfCentroids << std::endl;
      for( unsigned int i=0; i<numberOfCentroids; i++ )
        {
        for( unsigned int j = 0; j<3; j++ )
          {
          writeStream << ( *( filter->GetCentroids() ) )[i][j];
          if( j<2 )
            {
            writeStream << " ";
            }
          }
        writeStream << std::endl;
        }
      writeStream.close();

      imIn = filter->GetOutput();
      typename ImageType::SizeType size =
        imIn->GetLargestPossibleRegion().GetSize();

      filename = filename + ".mat";

      vnl_matrix<int> aMat( numberOfCentroids, numberOfCentroids );
      aMat.fill( 0 );

      itk::Index<dimensionT> indx;
      itk::Index<dimensionT> indx2;
      itk::Index<dimensionT> indx3;
      indx.Fill( 0 );
      bool done = false;
      bool invalid = false;
      bool done2 = false;
      int c, n;
      while( !done )
        {
        c = ( int )( imIn->GetPixel( indx )-1 );
        indx2.Fill( 0 );
        indx2[0] = 1;
        invalid = false;
        done2 = false;
        while( !done2 )
          {
          invalid = false;
          for( unsigned int d=0; d<dimensionT; d++ )
            {
            indx3[d] = indx[d] + indx2[d];
            if( indx3[d] >= ( int )size[d] )
              {
              invalid = true;
              break;
              }
            }
          if( !invalid )
            {
            n = ( int )( imIn->GetPixel( indx3 )-1 );
            if( c != n )
              {
              aMat( c, n ) = 1;
              aMat( n, c ) = 1;
              }
            }
          int i=0;
          indx2[i]++;
          while( !done2 && indx2[i]>=2 )
            {
            indx2[i] = 0;
            i++;
            if( i>2 )
              {
              done2 = true;
              }
            else
              {
              indx2[i]++;
              }
            }
          }
        int i = 0;
        indx[i]++;
        while( !done && indx[i]>=( int )size[i] )
          {
          indx[i] = 0;
          i++;
          if( i>2 )
            {
            done = true;
            }
          else
            {
            if( i == 2 )
              {
              std::cout << "Computing adjacency of slice : " << indx[2]
                        << std::endl;
              }
            indx[i]++;
            }
          }
        }

      writeStream.open( filename.c_str(),
        std::ios::binary | std::ios::out );
      if( ! writeStream.rdbuf()->is_open() )
        {
        std::cerr << "Cannot write to file : " << filename << std::endl;
        return EXIT_FAILURE;
        }
      writeStream << numberOfCentroids << std::endl;
      for( unsigned int i=0; i<numberOfCentroids; i++ )
        {
        for( unsigned int j = 0; j<numberOfCentroids; j++ )
          {
          writeStream << aMat( i, j );
          if( j<numberOfCentroids-1 )
            {
            writeStream << " ";
            }
          }
        writeStream << std::endl;
        }
      writeStream.close();

      } // end -S

    it++;
    }

  return EXIT_SUCCESS;
}

// Description:
// Get the ComponentType and dimension of the image
void GetImageInformation( std::string fileName,
                          itk::ImageIOBase::IOComponentType &componentType,
                          unsigned int & dimension )
{
  // Find out the component type of the image in file
  typedef itk::ImageIOBase::IOComponentType  PixelType;

  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO( fileName.c_str(),
                                        itk::ImageIOFactory::ReadMode );
  if( !imageIO )
    {
    std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
    return;
    }

  // Now that we found the appropriate ImageIO class, ask it to
  // read the meta data from the image file.
  imageIO->SetFileName( fileName.c_str() );
  imageIO->ReadImageInformation();

  componentType = imageIO->GetComponentType();
  dimension = imageIO->GetNumberOfDimensions();
}

int main( int argc, char *argv[] )
{
  MetaCommand command;

  command.SetName( "ImageMath" );
  command.SetVersion( "1.0" );
  command.SetAuthor( "CADDLab @ UNC" );
  command.SetDescription( "Perform several filters on an image" );

  command.SetOption( "Write", "w", false,
    "writes current image to the designated file" );
  command.AddOptionField( "Write", "filename", MetaCommand::STRING, true,
    "", "output filename", MetaCommand::DATA_OUT );

  command.SetOption( "WriteType", "W", false,
    "writes 0=UC 1=US 2=S 3=Old (4-6 uncomp UC,US,S) 7=uncomp F");
  command.AddOptionField( "WriteType", "Type", MetaCommand::INT, true );
  command.AddOptionField( "WriteType", "filename", MetaCommand::STRING,
    true, "", "output filename", MetaCommand::DATA_OUT );

  command.SetOption( "Intensity", "i", false,
    "Intensity window inVal range to outValRange" );
  command.AddOptionField( "Intensity", "inValMin", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Intensity", "inValMax", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Intensity", "outMin", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Intensity", "outMax", MetaCommand::FLOAT,
    true );

  command.SetOption( "IntensityMult", "I", false,
    "Intensity multiplicative correct using inMeanField" );
  command.AddOptionField( "IntensityMult", "inMeanField",
    MetaCommand::INT, true );

  command.SetOption( "GaussianNoise", "n", false,
    "Adds Gaussian noise to all pixels within inVal range" );
  command.AddOptionField( "GaussianNoise", "inValMin", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "GaussianNoise", "inValMax", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "GaussianNoise", "noiseMean",
    MetaCommand::FLOAT, true );
  command.AddOptionField( "GaussianNoise", "noiseStdDev",
    MetaCommand::FLOAT, true );

  command.SetOption( "UniformNoise", "N", false,
    "Adds uniform noise to all pixels within inVal range" );
  command.AddOptionField( "UniformNoise", "inValMin", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "UniformNoise", "inValMax", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "UniformNoise", "noiseMean", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "UniformNoise", "noiseRange",
    MetaCommand::FLOAT, true );

  command.SetOption( "Add", "a", false,
    "I( x ) = weight1*I( x ) + weight2*inFile2( x )" );
  command.AddOptionField( "Add", "weight1", MetaCommand::FLOAT, true );
  command.AddOptionField( "Add", "weight2", MetaCommand::FLOAT, true );
  command.AddOptionField( "Add", "Infile", MetaCommand::STRING, true );

  command.SetOption( "Algorithm", "A", false,
    "Return image value within masked region (mode: 0=mean, 1=stdDev)" );
  command.AddOptionField( "Algorithm", "mode", MetaCommand::INT, true );
  command.AddOptionField( "Algorithm", "threshLow", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Algorithm", "threshHigh", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Algorithm", "maskFile", MetaCommand::STRING,
    true );

  command.SetOption( "blur", "b", false,
    "gaussian blur the image using the given sigma" );
  command.AddOptionField( "blur", "sigma", MetaCommand::FLOAT, true );

  command.SetOption( "blurOrder", "B", false,
    "gaussian blur the image using the given sigma and the order." );
  command.AddOptionField( "blurOrder", "sigma", MetaCommand::FLOAT, true );
  command.AddOptionField( "blurOrder", "order", MetaCommand::INT, true );
  command.AddOptionField( "blurOrder", "direction", MetaCommand::INT,
    true );

  command.SetOption( "CorrectionSlice", "c", false,
    "Correct intensity slice-by-slice using HistogramMatchingFilter" );
  command.AddOptionField( "CorrectionSlice", "nBins", MetaCommand::INT,
    true );
  command.AddOptionField( "CorrectionSlice", "nMatchPoints",
    MetaCommand::INT, true );

  command.SetOption( "Correction", "C", false,
    "Match intensity to another volume using HistogramMatchingFilter" );
  command.AddOptionField( "Correction", "nBins", MetaCommand::INT, true );
  command.AddOptionField( "Correction", "nMatchPoints", MetaCommand::INT,
    true );
  command.AddOptionField(
    "Correction", "referenceVolume", MetaCommand::STRING, true );

  command.SetOption( "Normalize", "d", false,
    "Normalize image by setting mean to zero and standard deviation to one" );

  command.SetOption( "Fuse", "f", false,
    "fuse two images by max, applying offset to second image" );
  command.AddOptionField( "Fuse", "Offset2", MetaCommand::FLOAT, true );
  command.AddOptionField( "Fuse", "Infile2", MetaCommand::STRING, true );

  command.SetOption( "histogram", "l", false,
    "writes the image's histogram to the designated file" );
  command.AddOptionField( "histogram", "nBins", MetaCommand::INT, true );
  command.AddOptionField( "histogram", "histOutputFile",
    MetaCommand::STRING, true );

  command.SetOption( "histogram2", "L", false,
    "writes the image's histogram to the designated file" );
  command.AddOptionField( "histogram2", "nBins", MetaCommand::INT, true );
  command.AddOptionField( "histogram2", "binMin", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "histogram2", "binSIZE", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "histogram2", "histOutputFile",
    MetaCommand::STRING, true );

  command.SetOption( "Masking", "m", false,
    "if tLow<=inFile2(x)<=tHigh then I(x)=I(x) else I(x)=vFalse" );
  command.AddOptionField( "Masking", "threshLow", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Masking", "threshHigh", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Masking", "inFile2", MetaCommand::STRING,
    true );
  command.AddOptionField( "Masking", "valFalse", MetaCommand::FLOAT,
    true );

  command.SetOption( "Morphology", "M", false,
    "Mathematical morphology using a sphere. Mode: 0=erode, 1=dilate." );
  command.AddOptionField( "Morphology", "mode", MetaCommand::INT, true );
  command.AddOptionField( "Morphology", "radius", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Morphology", "forgroundValue",
    MetaCommand::FLOAT, true );
  command.AddOptionField( "Morphology", "backgroundValue",
    MetaCommand::FLOAT, true );

  command.SetOption( "overwrite", "o", false,
    "Replace values within the image, with a mask" );
  command.AddOptionField( "overwrite", "mask", MetaCommand::STRING, true );
  command.AddOptionField( "overwrite", "maskKeyVal", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "overwrite", "imageKeyVal", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "overwrite", "newImageVal", MetaCommand::FLOAT,
    true );

  command.SetOption( "offset", "O", false,
    "Set a new offset for the image" );
  command.AddOptionField( "offset", "offsetX", MetaCommand::FLOAT, true );
  command.AddOptionField( "offset", "offsetY", MetaCommand::FLOAT, true );
  command.AddOptionField( "offset", "offsetZ", MetaCommand::FLOAT, true );

  command.SetOption( "process", "p", false,
    "Process the image using a unary operation (0=abs)" );
  command.AddOptionField( "process", "mode", MetaCommand::INT, true );

  command.SetOption( "Process", "P", false,
    "Process the image using a binary operation (0=multiply)" );
  command.AddOptionField( "Process", "mode", MetaCommand::INT, true );
  command.AddOptionField( "Process", "file2", MetaCommand::STRING, true );

  command.SetOption( "resize", "r", false,
    "Resample to reduce by a factor (factor==0 means make isotropic)" );
  command.AddOptionField( "resize", "factor", MetaCommand::FLOAT, true );

  command.SetOption( "resize2", "R", false, "resample to match inFile2" );
  command.AddOptionField( "resize2", "inFile2", MetaCommand::STRING,
    true );

  command.SetOption( "segment", "s", false,
    "Segment using (inclusive) threshold connected components" );
  command.AddOptionField( "segment", "mode", MetaCommand::INT, true );
  command.AddOptionField( "segment", "threshLow", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "segment", "threshHigh", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "segment", "labelValue", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "segment", "x", MetaCommand::FLOAT, true );
  command.AddOptionField( "segment", "y", MetaCommand::FLOAT, true );
  command.AddOptionField( "segment", "z", MetaCommand::FLOAT, true );

  command.SetOption( "SetRandom", "S", false,
    "Sets the random number seed - to repeat experiments" );
  command.AddOptionField( "SetRandom", "seedValue", MetaCommand::FLOAT,
    true );

  command.SetOption( "Threshold", "t", false,
    "if tLow<=I(x)<=tHigh then I(x)=vTrue else I(x)=vFalse" );
  command.AddOptionField( "Threshold", "threshLow", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Threshold", "threshHigh", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Threshold", "valTrue", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Threshold", "valFalse", MetaCommand::FLOAT,
    true );

  command.SetOption( "Multiply","u",false,
    "I( x ) = I( x ) * inFile2( x )" );
  command.AddOptionField( "Multiply", "Infile", MetaCommand::STRING,
    true );

  command.SetOption( "vessels", "z", false,
    "Compute ridgness/vesselness for specified scales" );
  command.AddOptionField( "vessels", "scaleMin", MetaCommand::INT, true );
  command.AddOptionField( "vessels", "scaleMax", MetaCommand::INT, true );
  command.AddOptionField( "vessels", "numScales", MetaCommand::INT, true );

  command.SetOption( "Voronoi", "Z", false,
    "Run centroid voronoi tessellation on the image" );
  command.AddOptionField( "Voronoi", "numCentroids",
    MetaCommand::FLOAT, true );
  command.AddOptionField( "Voronoi", "numIters",
    MetaCommand::FLOAT, true );
  command.AddOptionField( "Voronoi", "numSamples",
    MetaCommand::FLOAT, true );
  command.AddOptionField( "Voronoi", "centroidOutFile",
    MetaCommand::STRING, true );

  command.AddField( "infile", "infile filename",
    MetaCommand::STRING, MetaCommand::DATA_IN );

  // Parsing
  if( !command.Parse( argc, argv ) )
    {
    return EXIT_FAILURE;
    }

  itk::ImageIOBase::IOComponentType componentType;
  unsigned int dimension;

  try
    {
    GetImageInformation( command.GetValueAsString( "infile" ),
                         componentType, dimension );
    if( dimension == 2 )
      {
      switch( componentType )
        {
        case itk::ImageIOBase::UCHAR:
          return DoIt<unsigned char, 2>( command );
        case itk::ImageIOBase::CHAR:
          return DoIt<char, 2>( command );
        case itk::ImageIOBase::USHORT:
          return DoIt<unsigned short, 2>( command );
        case itk::ImageIOBase::SHORT:
          return DoIt<short, 2>( command );
        case itk::ImageIOBase::UINT:
          return DoIt<unsigned int, 2>( command );
        case itk::ImageIOBase::INT:
          return DoIt<int, 2>( command );
        case itk::ImageIOBase::ULONG:
          return DoIt<unsigned long, 2>( command );
        case itk::ImageIOBase::LONG:
          return DoIt<long, 2>( command );
        case itk::ImageIOBase::FLOAT:
          return DoIt<float, 2>( command );
        case itk::ImageIOBase::DOUBLE:
          return DoIt<double, 2>( command );
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
          std::cout << "unknown component type" << std::endl;
          return EXIT_FAILURE;
        }
      }
    else if( dimension == 3 )
      {
      switch( componentType )
        {
        case itk::ImageIOBase::UCHAR:
          return DoIt<unsigned char, 3>( command );
        case itk::ImageIOBase::CHAR:
          return DoIt<char, 3>( command );
        case itk::ImageIOBase::USHORT:
          return DoIt<unsigned short, 3>( command );
        case itk::ImageIOBase::SHORT:
          return DoIt<short, 3>( command );
        case itk::ImageIOBase::UINT:
          return DoIt<unsigned int, 3>( command );
        case itk::ImageIOBase::INT:
          return DoIt<int, 3>( command );
        case itk::ImageIOBase::ULONG:
          return DoIt<unsigned long, 3>( command );
        case itk::ImageIOBase::LONG:
          return DoIt<long, 3>( command );
        case itk::ImageIOBase::FLOAT:
          return DoIt<float, 3>( command );
        case itk::ImageIOBase::DOUBLE:
          return DoIt<double, 3>( command );
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
          std::cout << "unknown component type" << std::endl;
          return EXIT_FAILURE;
        }
      }
    }
  catch( itk::ExceptionObject &excep )
    {
    std::cerr << argv[0] << ": itk exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

}
