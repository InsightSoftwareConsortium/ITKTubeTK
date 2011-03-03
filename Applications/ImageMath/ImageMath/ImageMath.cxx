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

// Meta Command include
#include "metaCommand.h"

// SimpleITK includes
#include "sitkImageFileReader.h"
#include "sitkImageFileWriter.h"
#include "sitkBinaryThresholdImageFilter.h"
#include "sitkIntensityWindowingImageFilter.h"
#include "sitkMultiplyImageFilter.h"
#include "sitkConnectedThresholdImageFilter.h"
#include "sitkExtractImageFilter.h"

// sitkIM includes
#include "sitkIMMorphologyFilter.h"
#include "sitkIMAddFilter.h"
#include "sitkIMBlurFilter.h"
#include "sitkIMBlurOrderFilter.h"
#include "sitkIMUnaryProcessFilter.h"
#include "sitkIMMaskingFilter.h"
#include "sitkIMResampleFilter.h"
#include "sitkIMIntensityMultFilter.h"
#include "sitkIMUniformNoiseFilter.h"
#include "sitkIMGaussianNoiseFilter.h"
#include "sitkIMFuseFilter.h"
#include "sitkIMAlgorithmReporter.h"
#include "sitkIMBinaryProcessFilter.h"
#include "sitkIMHistogramReporter.h"
#include "sitkIMHistogram2Reporter.h"
#include "sitkIMVesselsFilter.h"
#include "sitkIMCorrectionSliceFilter.h"
#include "sitkIMCorrectionFilter.h"
#include "sitkIMResizeFilter.h"
#include "sitkIMVoronoiReporter.h"
#include "sitkIMTypeImageWriter.h"


/******************************************************************************
 * Main body function - Dispatch the proper sequence of commands
 */
int DoIt( MetaCommand & command )
{
  // Parse the incoming command line
  MetaCommand::OptionVector parsed = command.GetParsedOptions();

  // Global seed for the gaussian generator
  int globalRandSeed = sitkIM::GaussianNoiseFilter::NO_SEED;


  // Read the file
  std::cout << "Reading file ( "
            << command.GetValueAsString( "infile" ).c_str()
            <<" )"<< std::endl;
  itk::simple::ImageFileReader reader;
  reader.SetFileName( command.GetValueAsString( "infile" ) );
  std::auto_ptr<itk::simple::Image> im( reader.Execute() );

  // Cast to float type
  itk::simple::CastImageFilter castFilter;
  castFilter.SetOutputPixelType( itk::simple::sitkFloat32 );
  im.reset( castFilter.Execute( im.get() ) );

  //////
  // Run the commands in sequence
  //////
  MetaCommand::OptionVector::const_iterator it = parsed.begin();
  while( it != parsed.end() )
    {
    if( ( *it ).name == "Write" )
      {
      std::string outFilename =
        command.GetValueAsString( *it, "filename" );
      std::cout << "Writing output1 ( " << outFilename.c_str()
                << " )" << std::endl;

      // Save the file
      itk::simple::ImageFileWriter writer;
      writer.SetFileName( outFilename );
      writer.Execute( im.get() );

      } // end -w

    else if( ( *it ).name == "WriteType" )
      {
      unsigned int type = (unsigned int)command.GetValueAsInt( "WriteType", "Type" );
      std::string outFilename =
        command.GetValueAsString( *it, "filename" );
      std::cout << "Writing output2 ( " << outFilename.c_str()
                << " )" << std::endl;

      // Use the sitIM type writer
      sitkIM::TypeImageWriter writer;
      writer.Execute(im.get(), type, outFilename);
      
      } // end -W

    else if( ( *it ).name == "Intensity" )
      {
      std::cout << "Intensity windowing" << std::endl;
      float valMin = command.GetValueAsFloat( *it, "inValMin" );
      float valMax = command.GetValueAsFloat( *it, "inValMax" );
      float outMin = command.GetValueAsFloat( *it, "outMin" );
      float outMax = command.GetValueAsFloat( *it, "outMax" );

      // Do the filtering
      itk::simple::IntensityWindowingImageFilter filter;
      im.reset( filter.Execute(im.get(), valMin, valMax, outMin, outMax) );

      }

    // IntensityMult
    else if( ( *it ).name == "IntensityMult" )
      {
      std::cout << "Intensity multiplicative bias correct" << std::endl;

      // Open the second image
      std::string inFile2 = command.GetValueAsString( *it, "inMeanField" );
      itk::simple::ImageFileReader reader;
      reader.SetFileName( inFile2 );
      itk::simple::Image* im2 = reader.Execute();

      // Do the filtering
      sitkIM::IntensityMultFilter filter;
      im.reset( filter.Execute(im.get(), im2) );

      // Clean up
      delete im2;

      } // end -I

    // UniformNoise
    else if( ( *it ).name == "UniformNoise" )
      {
      std::cout << "Adding noise" << std::endl;
      float valMin = command.GetValueAsFloat( *it, "inValMin" );
      float valMax = command.GetValueAsFloat( *it, "inValMax" );
      float noiseMean = command.GetValueAsFloat( *it, "noiseMean" );
      float noiseRange = command.GetValueAsFloat( *it, "noiseRange" );

      // Do the filtering
      sitkIM::UniformNoiseFilter filter;
      im.reset( filter.Execute( im.get(), valMin, valMax, noiseMean, noiseRange, globalRandSeed ) );
      
      } // -N

    // GaussianNoise
    else if( ( *it ).name == "GaussianNoise" )
      {
      std::cout << "Adding noise" << std::endl;
      float valMin = command.GetValueAsFloat( *it, "inValMin" );
      float valMax = command.GetValueAsFloat( *it, "inValMax" );
      float noiseMean = command.GetValueAsFloat( *it, "noiseMean" );
      float noiseStdDev = command.GetValueAsFloat( *it, "noiseStdDev" );

      // Do the filtering
      sitkIM::GaussianNoiseFilter filter;
      im.reset( filter.Execute( im.get(), valMin, valMax, noiseMean, noiseStdDev, globalRandSeed ) );

      } // end -n

    // Add
    else if( ( *it ).name == "Add" )
      {
      std::cout << "Adding" << std::endl;
      float weight1 = command.GetValueAsFloat( *it, "weight1" );
      float weight2 = command.GetValueAsFloat( *it, "weight2" );

      // Read the second image
      itk::simple::ImageFileReader reader;
      reader.SetFileName( command.GetValueAsString( *it, "Infile" ) );
      itk::simple::Image* im2 = reader.Execute();

      // Call the library filter
      sitkIM::AddFilter filter;
      im.reset( filter.Execute( im.get(), im2, weight1, weight2 ) );

      // Clean up
      delete im2;

      }

    // Multiply
    else if( ( *it ).name == "Multiply" )
      {
      std::cout << "Multiplying" << std::endl;

      // Read the second image
      itk::simple::ImageFileReader reader;
      reader.SetFileName( command.GetValueAsString( *it, "Infile" ) );
      itk::simple::Image* im2 = reader.Execute();

      // Call the library filter
      itk::simple::MultiplyImageFilter filter;
      im.reset( filter.Execute( im.get(), im2 ) );

      // Clean up
      delete im2;
      }

    // Fuse
    else if( ( *it ).name == "Fuse" )
      {
      std::cout << "Fusing" << std::endl;
      float offset = command.GetValueAsFloat( *it, "Offset2" );

      // Read the second image
      itk::simple::ImageFileReader reader;
      reader.SetFileName( command.GetValueAsString( *it, "Infile2" ) );
      itk::simple::Image* im2 = reader.Execute();

      // Call the library filter
      sitkIM::FuseFilter filter;
      im.reset( filter.Execute( im.get(), im2, offset ) );

      // Clean up
      delete im2;

      } // end -a

    // Threshold
    else if( ( *it ).name == "Threshold" )
      {
      std::cout << "Thresholding" << std::endl;

      float threshLow = command.GetValueAsFloat( *it, "threshLow" );
      float threshHigh = command.GetValueAsFloat( *it, "threshHigh" );

      float valTrue = command.GetValueAsFloat( *it, "valTrue" );
      float valFalse = command.GetValueAsFloat( *it, "valFalse" );

      // Filter it
      itk::simple::BinaryThresholdImageFilter filter;
      im.reset( filter.Execute( im.get(), threshLow, threshHigh, valTrue, valFalse ) );

      }

    // Algorithm
    else if( ( *it ).name == "Algorithm" )
      {
      std::cout << "Algorithm" << std::endl;

      float threshLow = command.GetValueAsFloat( *it, "threshLow" );
      float threshHigh = command.GetValueAsFloat( *it, "threshHigh" );

      // Read the second image
      itk::simple::ImageFileReader reader;
      reader.SetFileName( command.GetValueAsString( *it, "maskFile" ) );
      itk::simple::Image* im2 = reader.Execute();

      // Get the mode
      int mode = command.GetValueAsInt( *it, "mode" );

      // Call the reporter
      sitkIM::AlgorithmReporter reporter;
      double out = reporter.Execute( im.get(), im2, threshLow, threshHigh,
                    static_cast<sitkIM::AlgorithmReporter::ModeType>( mode ) );

      if( mode == 0 )
        {
        std::cout << "Mean " << out << std::endl;
        }
      else
        {
        std::cout << "StdDev " << out << std::endl;
        }

      // Clean up
      delete im2;
      }

    // Binary Process
    else if( ( *it ).name == "Process" )
      {
      std::cout << "Process binary operation" << std::endl;

      // Read the second image
      itk::simple::ImageFileReader reader;
      reader.SetFileName( command.GetValueAsString( *it, "file2" ) );
      itk::simple::Image* im2 = reader.Execute();

      // Get the mode
      int mode = command.GetValueAsInt( *it, "mode" );

      // Call the filter
      sitkIM::BinaryProcessFilter filter;
      im.reset( filter.Execute( im.get(), im2, mode ) );

      // Clean up
      delete im2;
      }

    // Unary Process
    else if( ( *it ).name == "process" )
      {
      std::cout << "Unary process" << std::endl;

      int mode = command.GetValueAsInt( *it, "mode" );

      // Filter it
      sitkIM::UnaryProcessFilter filter;
      im.reset( filter.Execute( im.get(), mode ) );
      }

    // Masking
    else if( ( *it ).name == "Masking" )
      {
      std::cout << "Masking" << std::endl;

      float threshLow = command.GetValueAsFloat( *it, "threshLow" );
      float threshHigh = command.GetValueAsFloat( *it, "threshHigh" );

      std::string inFile2 = command.GetValueAsString( *it, "inFile2" );

      float valFalse = command.GetValueAsFloat( *it, "valFalse" );

      // Open the second image
      itk::simple::ImageFileReader reader;
      reader.SetFileName( inFile2 );
      itk::simple::Image* im2 = reader.Execute();

      // Do the masking
      sitkIM::MaskingFilter filter;
      im.reset( filter.Execute( im.get(), im2, threshLow, threshHigh, valFalse ) );

      // Clean up
      delete im2;
      }

    // Morphology
    else if( ( *it ).name == "Morphology" )
      {
      std::cout << "Morphology" << std::endl;

      // Get the parameters
      int mode = command.GetValueAsInt( *it, "mode" );
      float radius = command.GetValueAsFloat( *it, "radius" );
      float foregroundValue = command.GetValueAsFloat( *it,
        "forgroundValue" );
      float backgroundValue = command.GetValueAsFloat( *it,
        "backgroundValue" );

      // Set up the sitkIM library filter
      sitkIM::MorphologyFilter filter;
      im.reset( filter.Execute( im.get(), static_cast<sitkIM::MorphologyFilter::OperationType>(mode),
                         foregroundValue, backgroundValue, radius ) );

      }

    // blur
    else if( ( *it ).name == "blur" )
      {
      std::cout << "Blurring." << std::endl;
      float sigma = command.GetValueAsFloat( *it, "sigma" );

      // Run library filter
      sitkIM::BlurFilter filter;
      im.reset( filter.Execute( im.get(), sigma ) );
      }

    // blurOrder
    else if( ( *it ).name == "blurOrder" )
      {
      std::cout << "Blurring." << std::endl;

      float sigma = command.GetValueAsFloat( *it, "sigma" );
      int order = command.GetValueAsInt( *it, "order" );
      int direction = command.GetValueAsInt( *it, "direction" );

      // Set up the filter
      sitkIM::BlurOrderFilter filter;
      im.reset( filter.Execute( im.get(), sigma, order, direction ) );
      } // end -B

    // histogram
    else if( ( *it ).name == "histogram" )
      {
      std::cout << "Histogram" << std::endl;

      unsigned int nBins = ( unsigned int )command.GetValueAsInt( *it,
        "nBins" );
      std::string filename =
        command.GetValueAsString( *it, "histOutputFile" );

      // Set up the reporter (ignore the output)
      sitkIM::HistogramReporter reporter;
      reporter.Execute( im.get(), nBins, filename );

      }

    // histogram2
    else if( ( *it ).name == "histogram2" )
      {
      std::cout << "Histogram" << std::endl;

      unsigned int nBins = ( unsigned int )command.GetValueAsInt( *it,
        "nBins" );
      double binMin = command.GetValueAsFloat( *it, "binMin" );
      double binSize = command.GetValueAsFloat( *it, "binSIZE" );
      std::string filename =
        command.GetValueAsString( *it, "histOutputFile" );

      // Set up the reporter (ignore the output)
      sitkIM::Histogram2Reporter reporter;
      reporter.Execute( im.get(), nBins, binMin, binSize, filename );
      }

    // vessels
    else if( ( *it ).name == "vessels" )
      {
      std::cout << "Vessel Enhancement" << std::endl;

      double scaleMin = command.GetValueAsFloat( *it, "scaleMin" );
      double scaleMax = command.GetValueAsFloat( *it, "scaleMax" );
      unsigned int numScales = (unsigned int)command.GetValueAsInt( *it, "numScales" );

      // Set up the filter
      sitkIM::VesselsFilter filter;
      im.reset( filter.Execute( im.get(), scaleMin, scaleMax, numScales ) );
      }

    // CorrectionSlice
    else if( ( *it ).name == "CorrectionSlice" )
      {
      std::cout << "Correct intensity slice-by-slice" << std::endl;

      unsigned int numberOfBins =
        ( unsigned int )command.GetValueAsInt( *it, "nBins" );
      unsigned int numberOfMatchPoints =
        ( unsigned int )command.GetValueAsInt( *it, "nMatchPoints" );

      // Set up the filter
      sitkIM::CorrectionSliceFilter filter;
      im.reset( filter.Execute( im.get(), numberOfBins, numberOfMatchPoints ) );
      }

    // Correction
    else if( ( *it ).name == "Correction" )
      {
      std::cout << "Correct intensity in the volume" << std::endl;

      unsigned int numberOfBins =
        ( unsigned int )command.GetValueAsInt( *it, "nBins" );
      unsigned int numberOfMatchPoints =
        ( unsigned int )command.GetValueAsInt( *it, "nMatchPoints" );
      std::string inFile2 = command.GetValueAsString( *it, "referenceVolume" );

      // Open the second image
      itk::simple::ImageFileReader reader;
      reader.SetFileName( inFile2 );
      itk::simple::Image* im2 = reader.Execute();

      // Do the masking
      sitkIM::CorrectionFilter filter;
      im.reset( filter.Execute( im.get(), im2, numberOfBins, numberOfMatchPoints ) );

      // Clean up
      delete im2;
      }

    // resize
    else if( ( *it ).name == "resize" )
      {
      std::cout << "Resampling." << std::endl;
      double factor = command.GetValueAsFloat( *it, "factor" );

      // Set up the filter
      sitkIM::ResizeFilter filter;
      im.reset( filter.Execute( im.get(), factor ) );
      }

    // resize2
    else if( ( *it ).name == "resize2" )
      {
      std::cout << "Resampling" << std::endl;
      std::string inFile2 = command.GetValueAsString( *it, "inFile2" );

      // Open the second image
      itk::simple::ImageFileReader reader;
      reader.SetFileName( inFile2 );
      itk::simple::Image* im2 = reader.Execute();

      // Set up the filter
      sitkIM::ResampleFilter filter;
      im.reset( filter.Execute( im.get(), im2 ) );

      // Clean up
      delete im2;
      }

    // segment
    else if( ( *it ).name == "segment" )
      {
      std::cout << "Segmenting" << std::endl;

      // IGNORING MODE????
      //int mode = command.GetValueAsInt( *it, "mode" );

      float threshLow = command.GetValueAsFloat( *it, "threshLow" );
      float threshHigh = command.GetValueAsFloat( *it, "threshHigh" );
      float labelValue = command.GetValueAsFloat( *it, "labelValue" );

      float x = command.GetValueAsFloat( *it, "x" );
      float y = command.GetValueAsFloat( *it, "y" );
      float z = command.GetValueAsFloat( *it, "z" );

      // Set up the filter
      itk::simple::ConnectedThresholdImageFilter filter;
      filter.SetLower( threshLow ).SetUpper( threshHigh ).SetReplaceValue( labelValue );
      std::vector<unsigned int> seed;
      seed.push_back(x);
      seed.push_back(y);
      if (im->GetDimension() == 3)
        {
        seed.push_back(z);
        }
      filter.SetSeed( seed );
      im.reset( filter.Execute( im.get() ) );
      }

    // offset
    else if( ( *it ).name == "offset" )
      {
      std::vector<double> offset;
      offset.push_back( command.GetValueAsFloat( *it, "offsetX" ) );
      offset.push_back( command.GetValueAsFloat( *it, "offsetY" ) );
      if( im->GetDimension() == 3 )
        {
        offset.push_back( command.GetValueAsFloat( *it, "offsetZ" ) );
        }

      // Update offset
      im->SetOrigin( offset );
      }

    // SetRandom
    else if( ( *it ).name == "SetRandom" )
      {
      int seed = command.GetValueAsInt( *it, "seedValue" );
      globalRandSeed = seed;
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

      // Set up the reporter
      sitkIM::VoronoiReporter reporter;
      im.reset( reporter.Execute( im.get(), numberOfCentroids, numberOfIterations,
                           numberOfSamples, filename ) );
      } // end -Z

    // ExtractSlice
    else if( ( *it ).name == "ExtractSlice" )
      {
      unsigned int dimension =
        ( unsigned int )command.GetValueAsInt( *it, "dimension" );
      unsigned int slice =
        ( unsigned int )command.GetValueAsInt( *it, "slice" );

      // Set up the slice extractor
      itk::simple::ExtractImageFilter filter;
      im.reset( filter.Execute(im.get(), slice, dimension) );
      } // end -e

    it++;
    }

  return EXIT_SUCCESS;
}


/******************************************************************************
 * Main entry point - Set up the command line parsing
 */
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

  command.SetOption( "Fuse", "f", false,
    "fuse two images by max, applying offset to second image" );
  command.AddOptionField( "Fuse", "Offset2", MetaCommand::FLOAT, true );
  command.AddOptionField( "Fuse", "Infile2", MetaCommand::STRING, true );

  command.SetOption( "Add", "a", false,
    "I( x ) = weight1*I( x ) + weight2*inFile2( x )" );
  command.AddOptionField( "Add", "weight1", MetaCommand::FLOAT, true );
  command.AddOptionField( "Add", "weight2", MetaCommand::FLOAT, true );
  command.AddOptionField( "Add", "Infile", MetaCommand::STRING, true );

  command.SetOption( "Multiply","u",false,
    "I( x ) = I( x ) * inFile2( x )" );
  command.AddOptionField( "Multiply", "Infile", MetaCommand::STRING,
    true );

  command.SetOption( "Algorithm", "A", false,
    "Return image value within masked region (mode: 0=mean, 1=stdDev)" );
  command.AddOptionField( "Algorithm", "mode", MetaCommand::INT, true );
  command.AddOptionField( "Algorithm", "threshLow", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Algorithm", "threshHigh", MetaCommand::FLOAT,
    true );
  command.AddOptionField( "Algorithm", "maskFile", MetaCommand::STRING,
    true );

  command.SetOption( "process", "p", false,
    "Process the image using a unary operation (0=abs)" );
  command.AddOptionField( "process", "mode", MetaCommand::INT, true );

  command.SetOption( "Process", "P", false,
    "Process the image using a binary operation (0=multiply)" );
  command.AddOptionField( "Process", "mode", MetaCommand::INT, true );
  command.AddOptionField( "Process", "file2", MetaCommand::STRING, true );

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

  command.SetOption( "blur", "b", false,
    "gaussian blur the image using the given sigma" );
  command.AddOptionField( "blur", "sigma", MetaCommand::FLOAT, true );

  command.SetOption( "blurOrder", "B", false,
    "gaussian blur the image using the given sigma and the order." );
  command.AddOptionField( "blurOrder", "sigma", MetaCommand::FLOAT, true );
  command.AddOptionField( "blurOrder", "order", MetaCommand::INT, true );
  command.AddOptionField( "blurOrder", "direction", MetaCommand::INT,
    true );

  command.SetOption( "vessels", "z", false,
    "Compute ridgness/vesselness for specified scales" );
  command.AddOptionField( "vessels", "scaleMin", MetaCommand::INT, true );
  command.AddOptionField( "vessels", "scaleMax", MetaCommand::INT, true );
  command.AddOptionField( "vessels", "numScales", MetaCommand::INT, true );

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

  command.SetOption( "offset", "O", false,
    "Set a new offset for the image" );
  command.AddOptionField( "offset", "offsetX", MetaCommand::FLOAT, true );
  command.AddOptionField( "offset", "offsetY", MetaCommand::FLOAT, true );
  command.AddOptionField( "offset", "offsetZ", MetaCommand::FLOAT, true );

  command.SetOption( "SetRandom", "S", false,
    "Sets the random number seed - to repeat experiments" );
  command.AddOptionField( "SetRandom", "seedValue", MetaCommand::FLOAT,
    true );

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

  command.SetOption( "ExtractSlice", "e", false,
    "Extract a single slice from the image" );
  command.AddOptionField( "ExtractSlice", "dimension",
    MetaCommand::INT, true );
  command.AddOptionField( "ExtractSlice", "slice",
    MetaCommand::INT, true );

  command.AddField( "infile", "infile filename",
    MetaCommand::STRING, MetaCommand::DATA_IN );

  // Parsing
  if( !command.Parse( argc, argv ) )
    {
    return EXIT_FAILURE;
    }

  //
  // Run the main body function
  //
  return DoIt( command );

}
