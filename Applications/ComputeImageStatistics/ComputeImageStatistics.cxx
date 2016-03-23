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

#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <vnl/vnl_math.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "ComputeImageStatisticsCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"


template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "tubeMaskToStats",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef itk::Image< TPixel,  VDimension >        MaskType;
  typedef itk::Image< unsigned int,  VDimension >  ConnCompType;
  typedef itk::Image< float,  VDimension >         VolumeType;
  typedef itk::ImageFileReader< VolumeType >       VolumeReaderType;
  typedef itk::ImageFileReader< MaskType >         MaskReaderType;

  //
  //
  //
  timeCollector.Start("Load mask");
  typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
  maskReader->SetFileName( inputMask.c_str() );
  try
    {
    maskReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading mask: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load mask");
  double progress = 0.1;
  progressReporter.Report( progress );

  typename MaskType::Pointer curMask = maskReader->GetOutput();

  //
  //
  //
  timeCollector.Start("Load volume");
  typename VolumeReaderType::Pointer volumeReader = VolumeReaderType::New();
  volumeReader->SetFileName( inputVolume.c_str() );
  try
    {
    volumeReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load volume");
  progress = 0.2;
  progressReporter.Report( progress );

  typename VolumeType::Pointer curVolume = volumeReader->GetOutput();

  //
  //
  //
  typename ConnCompType::Pointer curConnComp = ConnCompType::New();
  curConnComp->CopyInformation( curVolume );
  curConnComp->SetRegions( curVolume->GetLargestPossibleRegion() );
  curConnComp->Allocate();

  timeCollector.Start("Connected Components");

  typedef std::map< TPixel, unsigned int > MapType;

  MapType maskMap;

  itk::ImageRegionIterator< MaskType > maskIter( curMask,
    curMask->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ConnCompType > connCompIter( curConnComp,
    curConnComp->GetLargestPossibleRegion() );
  typename MapType::iterator mapIter;
  unsigned int numberOfComponents = 0;
  unsigned int id = 0;
  TPixel lastMaskV = maskIter.Get() + 1;
  TPixel maskV = 0;
  while( !connCompIter.IsAtEnd() )
    {
    maskV = maskIter.Get();

    if( maskV != lastMaskV )
      {
      lastMaskV = maskV;
      mapIter = maskMap.find( maskV );
      if( mapIter != maskMap.end() )
        {
        id = mapIter->second;
        }
      else
        {
        id = numberOfComponents;
        maskMap[ maskV ] = id;
        ++numberOfComponents;
        }
      }
    ++maskIter;
    ++connCompIter;
    }

  itk::Array< double > compMean( numberOfComponents );
  compMean.Fill( 0 );
  itk::Array< double > compMin( numberOfComponents );
  compMin.Fill( 0 );
  itk::Array< double > compMax( numberOfComponents );
  compMax.Fill( 0 );
  itk::Array< double > compStdDev( numberOfComponents );
  compStdDev.Fill( 0 );
  itk::Array< double > compCount( numberOfComponents );
  compCount.Fill( 0 );
  itk::Array< TPixel > compValue( numberOfComponents );
  compValue.Fill( 0 );

  unsigned int maxNumBins = 2000;
  vnl_matrix< unsigned int > compHisto( numberOfComponents, maxNumBins );
  compHisto.fill( 0 );


  maskIter.GoToBegin();
  connCompIter.GoToBegin();
  itk::ImageRegionIterator< VolumeType > volumeIter( curVolume,
    curVolume->GetLargestPossibleRegion() );
  id = 0;
  lastMaskV = maskIter.Get() + 1;
  maskV = 0;
  while( !connCompIter.IsAtEnd() )
    {
    maskV = maskIter.Get();

    if( maskV != lastMaskV )
      {
      lastMaskV = maskV;
      mapIter = maskMap.find( maskV );
      if( mapIter != maskMap.end() )
        {
        id = mapIter->second;
        }
      else
        {
        id = numberOfComponents;
        maskMap[ maskV ] = id;
        compValue[ id ] = maskV;
        ++numberOfComponents;
        }
      }

    double volumeV = volumeIter.Get();
    compMean[ id ] += volumeV;
    compStdDev[ id ] += volumeV * volumeV;
    if( compCount[ id ] == 0 )
      {
      compMin[ id ] = volumeV;
      compMax[ id ] = volumeV;
      }
    else
      {
      if( volumeV < compMin[ id ] )
        {
        compMin[ id ] = volumeV;
        }
      else
        {
        if( volumeV > compMax[ id ] )
          {
          compMax[ id ] = volumeV;
          }
        }
      }
    ++compCount[ id ];

    connCompIter.Set( id );

    ++maskIter;
    ++volumeIter;
    ++connCompIter;
    }

  maskIter.GoToBegin();
  volumeIter.GoToBegin();
  while( !volumeIter.IsAtEnd() )
    {
    maskV = maskIter.Get();
    mapIter = maskMap.find( maskV );
    id = mapIter->second;
    unsigned int bin = static_cast< unsigned int >( ( ( volumeIter.Get() -
      compMin[id] ) /
      ( compMax[id] - compMin[id] ) ) * maxNumBins + 0.5 );
    ++compHisto[id][bin];
    ++maskIter;
    ++volumeIter;
    }

  unsigned int numberOfQuantiles = quantiles.size();
  vnl_matrix< double > quantileValue( numberOfComponents,
    numberOfQuantiles );
  for( unsigned int quantileNumber=0; quantileNumber<numberOfQuantiles;
    ++quantileNumber )
    {
    double quant = quantiles[ quantileNumber ];
    for( unsigned int comp=0; comp<numberOfComponents; ++comp )
      {
      unsigned int targetCount = static_cast< unsigned int >( quant *
        compCount[ comp ] );
      unsigned int bin = 0;
      unsigned int binCount = 0;
      while( binCount + compHisto[ comp ][ bin ] < targetCount )
        {
        binCount += compHisto[ comp ][ bin ];
        ++bin;
        }
      double binPortion = static_cast< double >( targetCount - binCount )
        / compHisto[ comp ][ bin ];
      double binSize = ( compMax[ comp ] - compMin[ comp ] ) / maxNumBins;
      quantileValue[ comp ][ quantileNumber ] = ( bin + binPortion )
        * binSize + compMin[ comp ];
      }
    }

  std::cout << "Number of components = " << numberOfComponents << std::endl;
  std::ofstream writeStream;
  if( ! csvStatisticsFile.empty() )
    {
    writeStream.open( csvStatisticsFile.c_str(),
      std::ios::binary | std::ios::out );
    if( !writeStream.rdbuf()->is_open() )
      {
      std::cerr << "Cannot write to file " << csvStatisticsFile << std::endl;
      return EXIT_FAILURE;
      }
    }

  std::cout << "id, Value, Count, Mean, StdDev, Min, Max";
  for( unsigned int q=0; q<numberOfQuantiles; ++q )
    {
    std::cout << ", " << quantiles[ q ];
    }
  std::cout << std::endl;
  if( ! csvStatisticsFile.empty() )
    {
    writeStream << "id, Value, Count, Mean, StdDev, Min, Max";
    for( unsigned int q=0; q<numberOfQuantiles; ++q )
      {
      writeStream << ", " << quantiles[ q ];
      }
    writeStream << std::endl;
    }
  for( unsigned int i=0; i<numberOfComponents; ++i )
    {
    std::cout << i << ", ";
    std::cout << static_cast< double >( compValue[ i ] ) << ", ";
    std::cout << compCount[ i ] << ", ";
    if( ! csvStatisticsFile.empty() )
      {
      writeStream << i << ", ";
      writeStream << static_cast< double >( compValue[ i ] ) << ", ";
      writeStream << compCount[ i ] << ", ";
      }
    if( compCount[ i ] > 0 )
      {
      compStdDev[ i ] = std::sqrt( ( compStdDev[ i ]
        - ( ( compMean[ i ] * compMean[ i ] ) / compCount[ i ] ) )
        / ( compCount[ i ] - 1 ) );
      compMean[ i ] /= compCount[ i ];
      }
    std::cout << compMean[ i ] << ", ";
    std::cout << compStdDev[ i ] << ", ";
    std::cout << compMin[ i ] << ", ";
    std::cout << compMax[ i ];
    for( unsigned int q=0; q<numberOfQuantiles; ++q )
      {
      std::cout << ", " << quantileValue[ i ][ q ];
      }
    std::cout << std::endl;
    if( ! csvStatisticsFile.empty() )
      {
      writeStream << compMean[ i ] << ", ";
      writeStream << compStdDev[ i ] << ", ";
      writeStream << compMin[ i ] << ", ";
      writeStream << compMax[ i ];
      for( unsigned int q=0; q<numberOfQuantiles; ++q )
        {
        writeStream << ", " << quantileValue[ i ][ q ];
        }
      writeStream << std::endl;
      }
    }
  if( ! csvStatisticsFile.empty() )
    {
    writeStream.close();
    }

  connCompIter.GoToBegin();
  volumeIter.GoToBegin();
  while( !connCompIter.IsAtEnd() )
    {
    id = connCompIter.Get();

    volumeIter.Set( compMean[ id ] );

    ++connCompIter;
    ++volumeIter;
    }
  timeCollector.Stop("Connected Components");

  typedef itk::ImageFileWriter< VolumeType  >   ImageWriterType;

  timeCollector.Start("Save data");
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputMeanImage.c_str() );
  writer->SetInput( curVolume );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume: Exception caught: "
      + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Save data");
  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputMask, argc, argv );
}
