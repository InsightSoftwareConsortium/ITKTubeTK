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

#include "tubeImageFilters.h"

#include <itkImageFileWriter.h>

#include <metaCommand.h>
#include "ImageMathCLP.h"

/** Main command */
template< class TPixel, unsigned int VDimension >
int DoIt( MetaCommand & command )
{
  typedef float                                    PixelType;
  typedef itk::Image< PixelType, VDimension >      ImageType;
  typedef itk::Image< unsigned char, VDimension >  ImageTypeUChar;
  typedef itk::Image< unsigned short, VDimension > ImageTypeUShort;
  typedef itk::Image< short, VDimension >          ImageTypeShort;

  MetaCommand::OptionVector parsed = command.GetParsedOptions();


  int CurrentSeed = 42;

  typedef itk::ImageFileReader< ImageType >       VolumeReaderType;

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

          for( unsigned int d=0; d<ImageType::ImageDimension; ++d )
            {
            metaImage->ElementSize( d, imIn->GetSpacing()[d] );
            }

          metaImage->AddUserField( "ElementByteOrderMSB",
            MET_STRING, std::strlen( "False" ), "False" );

          writer->Write();
          break;
          }
        case 7:
          {
          typedef itk::ImageFileWriter< ImageType > VolumeWriterType;
          typename VolumeWriterType::Pointer writer =
            VolumeWriterType::New();
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
      tube::ImageFilters< VDimension >::ApplyIntensityWindowing( 
        imIn,
        command.GetValueAsFloat( *it, "inValMin" ),
        command.GetValueAsFloat( *it, "inValMax" ),
        command.GetValueAsFloat( *it, "outMin" ),
        command.GetValueAsFloat( *it, "outMax" ) );
      }

    // IntensityMult
    else if( ( *it ).name == "IntensityMult" )
      {
      std::cout << "Intensity multiplicative bias correct" << std::endl;
      bool success = tube::ImageFilters< VDimension >::
        ApplyIntensityMultiplicativeWithBiasCorrection( imIn,
        command.GetValueAsString( *it, "inMeanField" ) );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      } // end -I

    // UniformNoise
    else if( ( *it ).name == "UniformNoise" )
      {
      std::cout << "Adding noise" << std::endl;
      tube::ImageFilters< VDimension >::AddUniformNoise( imIn,
        command.GetValueAsFloat( *it, "inValMin" ),
        command.GetValueAsFloat( *it, "inValMax" ),
        command.GetValueAsFloat( *it, "noiseMean" ),
        command.GetValueAsFloat( *it, "noiseRange" ), CurrentSeed );
      } // -N

    // GaussianNoise
    else if( ( *it ).name == "GaussianNoise" )
      {
      std::cout << "Adding noise" << std::endl;
      tube::ImageFilters< VDimension >::AddGaussianNoise( imIn,
        command.GetValueAsFloat( *it, "inValMin" ),
        command.GetValueAsFloat( *it, "inValMax" ),
        command.GetValueAsFloat( *it, "noiseMean" ),
        command.GetValueAsFloat( *it, "noiseStdDev" ), CurrentSeed );
      } // end -n

    // I( x )
    else if( ( *it ).name == "Add" )
      {
      std::cout << "Adding" << std::endl;
      bool success = tube::ImageFilters< VDimension >::AddImages( imIn,
        command.GetValueAsString( *it, "Infile" ),
        command.GetValueAsFloat( *it, "weight1" ),
        command.GetValueAsFloat( *it, "weight2" ) );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      }

    else if( ( *it ).name == "Multiply" )
      {
      std::cout << "Multiplying" << std::endl;
      bool success = tube::ImageFilters< VDimension >::MultiplyImages( 
        imIn, command.GetValueAsString( *it, "Infile" ) );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      }

    // Mirror pad
    else if( ( *it ).name == "MirrorPad" )
      {
      tube::ImageFilters< VDimension >::MirrorAndPadImage( imIn,
        command.GetValueAsInt( *it, "numPadVoxels" ) );
      }

    // Normalize
    else if( ( *it ).name == "Normalize" )
      {
      std::cout << "Normalize" << std::endl;
      std::cout << "NOTE: since this filter normalizes the data to lie "
        << "within -1 to 1, integral types will produce an image that "
        << "DOES NOT HAVE a unit variance" << std::endl;

      tube::ImageFilters< VDimension >::template NormalizeImage<
        TPixel >( imIn, command.GetValueAsInt( *it, "type" ) );
      }

    // I( x )
    else if( ( *it ).name == "Fuse" )
      {
      std::cout << "Fusing" << std::endl;
      bool success = tube::ImageFilters< VDimension >::FuseImages( 
        imIn, command.GetValueAsString( *it, "Infile2" ),
        command.GetValueAsFloat( *it, "Offset2" ) );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      } // end -a

    // Median
    else if( ( *it ).name == "Median" )
      {
      std::cout << "Median filtering" << std::endl;
      bool success = tube::ImageFilters< VDimension >::MedianImage( 
        imIn, command.GetValueAsInt( *it, "Size" ) );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      } // end -g

    // Threshold
    else if( ( *it ).name == "Threshold" )
      {
      std::cout << "Thresholding" << std::endl;
      tube::ImageFilters< VDimension >::ThresholdImage( imIn,
        command.GetValueAsFloat( *it, "threshLow" ),
        command.GetValueAsFloat( *it, "threshHigh" ),
        command.GetValueAsFloat( *it, "valTrue" ),
        command.GetValueAsFloat( *it, "valFalse" ) );
      }
    else if( ( *it ).name == "Algorithm" )
      {
      std::cout << "Algorithm" << std::endl;
      bool success = false;
      int mode = command.GetValueAsInt( *it, "mode" );
      double value = tube::ImageFilters< VDimension >::
        ComputeImageStdDevOrMeanWithinRangeUsingMask( imIn,
        command.GetValueAsString( *it, "maskFile" ),
        command.GetValueAsFloat( *it, "threshLow" ),
        command.GetValueAsFloat( *it, "threshHigh" ),
        mode, success );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      std::cout << ( mode == 0 ? "Mean " : "StdDev " ) << value
        << std::endl;
      }
    else if( ( *it ).name == "Process" )
      {
      std::cout << "Process binary operation" << std::endl;
      int mode = command.GetValueAsInt( *it, "mode" );
      if( mode == 0 )
        {
        bool success =
          tube::ImageFilters< VDimension >::MultiplyImages( imIn,
          command.GetValueAsString( *it, "file2" ) );
        if( !success )
          {
          return EXIT_FAILURE;
          }
        }
      }
    else if( ( *it ).name == "process" )
      {
      std::cout << "Unary process" << std::endl;
      int mode = command.GetValueAsInt( *it, "mode" );
      if( mode == 0 )
        {
        tube::ImageFilters< VDimension >::AbsoluteImage( imIn );
        }
      }
    // Masking
    else if( ( *it ).name == "Masking" )
      {
      std::cout << "Masking" << std::endl;
      bool success = tube::ImageFilters< VDimension >::
        MaskImageWithValueIfNotWithinSecondImageRange( imIn,
        command.GetValueAsString( *it, "inFile2" ),
        command.GetValueAsFloat( *it, "threshLow" ),
        command.GetValueAsFloat( *it, "threshHigh" ),
        command.GetValueAsFloat( *it, "valFalse" ) );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      }

    // Morphology
    else if( ( *it ).name == "Morphology" )
      {
      std::cout << "Morphology" << std::endl;

      tube::ImageFilters< VDimension >::MorphImage( imIn,
        command.GetValueAsInt( *it, "mode" ),
        command.GetValueAsFloat( *it, "radius" ),
        command.GetValueAsFloat( *it, "forgroundValue" ),
        command.GetValueAsFloat( *it, "backgroundValue" ) );
      }

    else if( ( *it ).name == "overwrite" )
      {
      tube::ImageFilters< VDimension >::OverwriteImage( imIn,
        command.GetValueAsString( *it, "mask" ),
        command.GetValueAsFloat( *it, "maskKeyVal" ),
        command.GetValueAsFloat( *it, "imageKeyVal" ),
        command.GetValueAsFloat( *it, "newImageVal" ) );
      std::cout << "Overwrite" << std::endl;
      }

    // blur
    else if( ( *it ).name == "blur" )
      {
      std::cout << "Blurring." << std::endl;
      tube::ImageFilters< VDimension >::BlurImage( imIn,
        command.GetValueAsFloat( *it, "sigma" ) );
      }

    // blurOrder
    else if( ( *it ).name == "blurOrder" )
      {
      std::cout << "Blurring." << std::endl;
      tube::ImageFilters< VDimension >::BlurOrderImage( imIn,
        command.GetValueAsFloat( *it, "sigma" ),
        command.GetValueAsInt( *it, "order" ),
        command.GetValueAsInt( *it, "direction" ) );
      } // end -B

    // histogram
    else if( ( *it ).name == "histogram" )
      {
      std::cout << "Histogram" << std::endl;
      bool success = tube::ImageFilters< VDimension >::
        ComputeImageHistogram( imIn,
        static_cast<unsigned int>( command.GetValueAsInt( *it, "nBins" ) ),
        command.GetValueAsString( *it, "histOutputFile" ) );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      }

    // histogram2
    else if( ( *it ).name == "histogram2" )
      {
      std::cout << "Histogram" << std::endl;
      bool success = tube::ImageFilters< VDimension >::
        ComputeImageHistogram2( imIn,
        static_cast<unsigned int>( command.GetValueAsInt( *it, "nBins" ) ),
        command.GetValueAsFloat( *it, "binMin" ),
        command.GetValueAsFloat( *it, "binSIZE" ),
        command.GetValueAsString( *it, "histOutputFile" ) );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      }

    // vessels
    else if( ( *it ).name == "vessels" )
      {
      std::cout << "Vessel Enhancement" << std::endl;
      tube::ImageFilters< VDimension >::EnhanceVessels( imIn,
        command.GetValueAsFloat( *it, "scaleMin" ),
        command.GetValueAsFloat( *it, "scaleMax" ),
        command.GetValueAsFloat( *it, "numScales" ) );
      }

    // CorrectionSlice
    else if( ( *it ).name == "CorrectionSlice" )
      {
      std::cout << "Correct intensity slice-by-slice" << std::endl;
      tube::ImageFilters< VDimension >::
        CorrectIntensitySliceBySliceUsingHistogramMatching( imIn,
        static_cast<unsigned int>( command.GetValueAsInt( *it, "nBins" ) ),
        static_cast<unsigned int>( command.GetValueAsInt( *it,
        "nMatchPoints" ) ) );
      }

    // Correction
    else if( ( *it ).name == "Correction" )
      {
      std::cout << "Correct intensity in the volume" << std::endl;
      bool success = tube::ImageFilters< VDimension >::
        CorrectIntensityUsingHistogramMatching( imIn,
        static_cast<unsigned int>( command.GetValueAsInt( *it, "nBins" ) ),
        static_cast<unsigned int>( command.GetValueAsInt( *it,
        "nMatchPoints" ) ), command.GetValueAsString( *it,
        "referenceVolume" ) );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      }

    // resize
    else if( ( *it ).name == "resize" )
      {
      std::cout << "Resampling." << std::endl;
      tube::ImageFilters< VDimension >::Resize( imIn,
        command.GetValueAsFloat( *it, "factor" ) );
      }

    // resize2
    else if( ( *it ).name == "resize2" )
      {
      std::cout << "Resampling" << std::endl;
      bool success = tube::ImageFilters< VDimension >::Resize( imIn,
        command.GetValueAsString( *it, "inFile2" ) );
      if( !success )
        {
        return EXIT_FAILURE;
        }
      }

    // segment
    else if( ( *it ).name == "segment" )
      {
      std::cout << "Segmenting" << std::endl;
      tube::ImageFilters< VDimension >::
        SegmentUsingConnectedThreshold( imIn,
        command.GetValueAsFloat( *it, "threshLow" ),
        command.GetValueAsFloat( *it, "threshHigh" ),
        command.GetValueAsFloat( *it, "labelValue" ),
        command.GetValueAsFloat( *it, "x" ),
        command.GetValueAsFloat( *it, "y" ),
        VDimension == 3 ? command.GetValueAsFloat( *it, "z" ) : 0.0 );
      }

    // offset
    else if( ( *it ).name == "offset" )
      {
      double offset[VDimension];
      offset[0] = command.GetValueAsFloat( *it, "offsetX" );
      offset[1] = command.GetValueAsFloat( *it, "offsetY" );
      if( VDimension == 3 )
        {
        offset[VDimension-1] = command.GetValueAsFloat( *it, "offsetZ" );
        }
      imIn->SetOrigin( offset );
      }

    // SetRandom
    else if( ( *it ).name == "SetRandom" )
      {
      CurrentSeed = ( unsigned int )command.GetValueAsInt( *it,
        "seedValue" );
      } // end -S

    // Voronoi
    else if( ( *it ).name == "Voronoi" )
      {
      tube::ImageFilters< VDimension >::
        ComputeVoronoiTessellation( imIn,
        static_cast<unsigned int>( command.GetValueAsInt( *it,
        "numCentroids" ) ),
        static_cast<unsigned int>( command.GetValueAsInt( *it,
        "numIters" ) ),
        static_cast<unsigned int>( command.GetValueAsInt( *it,
        "numSamples" ) ),
        command.GetValueAsString( *it,
        "centroidOutFile" ) );
      if( !imIn )
        {
        return EXIT_FAILURE;
        }
      } // end -S

    // ExtractSlice
    else if( ( *it ).name == "ExtractSlice" )
      {
      tube::ImageFilters< VDimension >::ExtractSlice( imIn,
        static_cast<unsigned int>( command.GetValueAsInt( *it,
          "dimension" ) ),
        static_cast<unsigned int>( command.GetValueAsInt( *it,
          "slice" ) ) );
      } // end -e

    ++it;
    }

  return EXIT_SUCCESS;
}

// Description:
// Get the ComponentType and dimension of the image
void GetImageInformation( std::string fileName,
  itk::ImageIOBase::IOComponentType &componentType,
  unsigned int & dimension )
{
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO( fileName.c_str(),
                                        itk::ImageIOFactory::ReadMode );
  if( !imageIO )
    {
    std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
    dimension = 0;
    return;
    }

  // Now that we found the appropriate ImageIO class, ask it to
  // read the meta data from the image file.
  imageIO->SetFileName( fileName.c_str() );
  imageIO->ReadImageInformation();

  componentType = imageIO->GetComponentType();
  dimension = imageIO->GetNumberOfDimensions();
}

int main( int argc, char * argv[] )
{
  //PARSE_ARGS;
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
    "write UChar,UShort,Short,Old Meta,U-UChar,U-UShort,U-Short,U-Float )"
    );
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
    "Return image value within masked region ( mode: 0=mean, 1=stdDev )" );
  command.AddOptionField( "Algorithm", "mode", MetaCommand::INT, true );
  command.AddOptionField( "Algorithm", "threshLow", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Algorithm", "threshHigh", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Algorithm", "maskFile", MetaCommand::STRING,
    true );

  command.SetOption( "blur", "b", false,
    "Gaussian blur the image using the given sigma" );
  command.AddOptionField( "blur", "sigma", MetaCommand::FLOAT, true );

  command.SetOption( "blurOrder", "B", false,
    "Gaussian blur the image using the given sigma and the order." );
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
    "Normalize: 0 = mean/std; 1 = FWHM ; 2 = FWHM mean ( shift ) only" );
  command.AddOptionField( "Normalize", "type", MetaCommand::INT, true );

  command.SetOption( "Fuse", "f", false,
    "fuse two images by max, applying offset to second image" );
  command.AddOptionField( "Fuse", "Offset2", MetaCommand::FLOAT, true );
  command.AddOptionField( "Fuse", "Infile2", MetaCommand::STRING, true );

  command.SetOption( "Median", "g", false,
    "Apply a median filter to the image" );
  command.AddOptionField( "Median", "Size", MetaCommand::INT, true );

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
    "if inFile( x ) in [tLow, tHigh] then I( x )=I( x ) else I( x )=vFalse" );
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
    "Process the image using a unary operation ( 0=abs )" );
  command.AddOptionField( "process", "mode", MetaCommand::INT, true );

  command.SetOption( "Process", "P", false,
    "Process the image using a binary operation ( 0=multiply )" );
  command.AddOptionField( "Process", "mode", MetaCommand::INT, true );
  command.AddOptionField( "Process", "file2", MetaCommand::STRING, true );

  command.SetOption( "resize", "r", false,
    "Resample to reduce by a factor ( factor==0 means make isotropic )" );
  command.AddOptionField( "resize", "factor", MetaCommand::FLOAT, true );

  command.SetOption( "resize2", "R", false, "resample to match inFile2" );
  command.AddOptionField( "resize2", "inFile2", MetaCommand::STRING,
    true );

  command.SetOption( "segment", "s", false,
    "Segment using ( inclusive ) threshold connected components" );
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
    "if I( x ) in [tLow,tHigh] then I( x )=vTrue else I( x )=vFalse" );
  command.AddOptionField( "Threshold", "threshLow", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Threshold", "threshHigh", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Threshold", "valTrue", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Threshold", "valFalse", MetaCommand::FLOAT,
    true );

  command.SetOption( "Multiply", "u", false,
    "I( x ) = I( x ) * inFile2( x )" );
  command.AddOptionField( "Multiply", "Infile", MetaCommand::STRING,
    true );

  command.SetOption( "MirrorPad", "x", false,
    "Use mirroring to pad an image" );
  command.AddOptionField( "MirrorPad", "numPadVoxels", MetaCommand::INT,
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

  command.SetOption( "ExtractSlice", "e", false,
    "Extract a single slice from the image" );
  command.AddOptionField( "ExtractSlice", "dimension",
    MetaCommand::INT, true );
  command.AddOptionField( "ExtractSlice", "slice",
    MetaCommand::INT, true );
  // Parsing
  if( !command.Parse( argc, argv ) )
    {
    return EXIT_FAILURE;
    }

  itk::ImageIOBase::IOComponentType componentType;


  try
    {
    unsigned int dimension;
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
    else
      {
      std::cout << "only images of dimension 2 or 3 allowed !" << std::endl;
      return EXIT_FAILURE;
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
