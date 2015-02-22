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

#include "ComputeImageStatisticsUsingMaskCLP.h"

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

  unsigned int maxNumComponents = curVolume->GetLargestPossibleRegion().
    GetNumberOfPixels();
  itk::Array< double > compMean(maxNumComponents);
  compMean.Fill( 0 );
  itk::Array< double > compMin(maxNumComponents);
  compMin.Fill( 0 );
  itk::Array< double > compMax(maxNumComponents);
  compMax.Fill( 0 );
  itk::Array< double > compStdDev(maxNumComponents);
  compStdDev.Fill( 0 );
  itk::Array< double > compCount(maxNumComponents);
  compCount.Fill( 0 );

  itk::ImageRegionIterator< MaskType > maskIter( curMask,
    curMask->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< VolumeType > volumeIter( curVolume,
    curVolume->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ConnCompType > connCompIter( curConnComp,
    curConnComp->GetLargestPossibleRegion() );
  typename MapType::iterator mapIter;
  unsigned int numberOfComponents = 0;
  unsigned int id = 0;
  TPixel maskV = 0;
  double volumeV = 0;
  while( !connCompIter.IsAtEnd() )
    {
    maskV = maskIter.Get();

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

    volumeV = volumeIter.Get();
    compMean[id] += volumeV;
    compStdDev[id] += volumeV * volumeV;
    if( compCount[id] == 0 )
      {
      compMin[id] = volumeV;
      compMax[id] = volumeV;
      }
    else
      {
      if( volumeV < compMin[id] )
        {
        compMin[id] = volumeV;
        }
      else
        {
        if( volumeV > compMax[id] )
          {
          compMax[id] = volumeV;
          }
        }
      }
    ++compCount[id];

    connCompIter.Set( id );

    ++maskIter;
    ++volumeIter;
    ++connCompIter;
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

  std::cout << "id, Value, Count, Mean, StdDev, Min, Max" << std::endl;
  if( ! csvStatisticsFile.empty() )
    {
    writeStream << "id, Value, Count, Mean, StdDev, Min, Max" << std::endl;
    }
  for( unsigned int i=0; i<numberOfComponents; i++ )
    {
    std::cout << i << ", ";
    std::cout << maskMap[i] << ", ";
    std::cout << compCount[i] << ", ";
    if( ! csvStatisticsFile.empty() )
      {
      writeStream << i << ", ";
      writeStream << maskMap[i] << ", ";
      writeStream << compCount[i] << ", ";
      }
    if( compCount[i] > 0 )
      {
      compStdDev[i] = vcl_sqrt( ( compStdDev[i]
        - ( ( compMean[i] * compMean[i] ) / compCount[i] ) )
        / ( compCount[i] - 1 ) );
      compMean[i] /= compCount[i];
      }
    std::cout << compMean[i] << ", ";
    std::cout << compStdDev[i] << ", ";
    std::cout << compMin[i] << ", ";
    std::cout << compMax[i] << std::endl;
    if( ! csvStatisticsFile.empty() )
      {
      writeStream << compMean[i] << ", ";
      writeStream << compStdDev[i] << ", ";
      writeStream << compMin[i] << ", ";
      writeStream << compMax[i] << std::endl;
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

    volumeIter.Set( compMean[id] );

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
