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

#include "itkTimeProbesCollectorBase.h"

#include "itkImageFileReader.h"

#include "itkSpatialObjectReader.h"
#include "itkGroupSpatialObject.h"
#include "itkTubeSpatialObject.h"

#include "ComputeTubeFlyThroughImageCLP.h"

template< class TPixel, unsigned int TDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "<ModuleName>CLP.h"
//   and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;
  
  // Load Input Image
  typedef TPixel                                        PixelType;
  typedef itk::Image< PixelType, VDimension >           ImageType;
  typedef itk::ImageFileReader< ImageType >             ImageReaderType;
  
  timeCollector.Start("Loading Input Image");
  
  typename ImageReaderType::Pointer imageReader = ReaderType::New();
  
  try
  {
    imageReader->SetFileName( inputImageFile.c_str() );    
    imageReader->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    tube::ErrorMessage( "Error loading input image: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
  }
  
  typename ImageType::Pointer inputImage = imageReader.GetOutput();
  
  timeCollector.Stop("Loading Input Image");
  
  // Load TRE File
  typedef itk::SpatialObjectReader< VDimension >      	TubesReaderType;  
  typedef itk::GroupSpatialObject< Dimension >        	TubeGroupType;

  timeCollector.Start("Loading Input TRE File");
  
  typename TubesReaderType::Pointer tubesReader = TubesReaderType::New();
  
  try
  {
    tubesReader->SetFileName( file );
    tubesReader->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    tube::ErrorMessage( "Error loading TRE File: "
			+ std::string(err.GetDescription()));
    timeCollector.Report();
    return EXIT_FAILURE;
  }

  TubeGroupType::Pointer tubeGroup = tubesReader->GetGroup(); 
  
  timeCollector.Stop("Loading Input TRE File");

  // Find the user-specified tube
  typedef typename TubeGroupType::ChildrenListType	TubeListType;
  typedef itk::TubeSpatialObject< VDimension >   	TubeType;
  typename TubeType::Pointer inputTube = ITK_NULLPTR;
  
  timeCollector.Start("Finding the user specified tube");
  
  TubeListType * tubeList = tubeGroup->GetChildren();
  typename TubeListType::const_iterator itTubes = tubeList->begin();
  
  while( itTubes != tubeList->end() )
  {
      if( (*itTubes)->GetId() == inputTubeId )
      {	
	inputTube = dynamic_cast<TubeType *>( (*itTubes).GetPointer() );
	break;
      }    
      
      itTubes++;
  }
  
  delete tubeList;
  
  if( inputTube == ITK_NULLPTR )
  {
      tube::ErrorMessage("Unable to find the specified tube");
      timeCollector.Report();
      return EXIT_FAILURE;
  }
  
  inputTube->ComputeTangentsAndNormals();
  
  timeCollector.Stop("Finding the user specified tube");  
  
  // Get list of tube points
  typedef typename TubeType::PointListType		TubePointListType;     
  
  TubePointListType tubePointList = inputTube->GetPoints();
  if( tubePointList.size() <= 0 )
  {
      tube::ErrorMessage("The specified tube does not contain any points");
      timeCollector.Report();
      return EXIT_FAILURE;
  }
    
  // Determine maximum radius among all tube points 
  typename TubePointListType::const_iteraror itPts = tubePointList.begin();
  double maxTubeRadius = itPts->GetRadius();
  
  while( itPts != tubePointList.end() )
  {
      if( itPts->GetRadius() > maxTubeRadius )
      {
	  maxTubeRadius = itPts->GetRadius();	
      }    
  }  
  
  // Determine minimum spacing accross all dimensions
  typename ImageType::SpacingType inputSpacing = inputImage->GetSpacing();  
  double minInputSpacing = inputSpacing[0];
  
  for(unsigned int i = 1; i < VDimension; i++)
  {
    if( inputSpacing[i] < minInputSpacing )
    {
      minInputSpacing = inputSpacing[i];
    }  
  }
    
  // Setup the output image and allocate space for it
  typename ImageType::Pointer outputImage = ImageType::New();  

    // set start index
    typename ImageType::IndexType startIndex;
    startIndex.Fill(0);

    // set size
    typename ImageType::SizeType size;  
    for(unsigned int i = 0; i < VDimension-1; i++)
    {
	size[i] = 2 * (typename ImageType::SizeType)
		      (maxTubeRadius / minInputSpacing) + 1;
    }
    size[VDimension-1] = tubePointList.size();	

    // set regions
    typename ImageType::RegionType region;  
    region.SetIndex( startIndex );  
    region.SetSize( size );    
    outputImage.SetRegions( region );  
    
    // set spacing
    typename ImageType::SpacingType outputSpacing;
     
    for(unsigned int i = 0; i < VDimension; i++)
    {
      outputSpacing[i] = minInputSpacing;
    }
    outputImage->SetSpacing(outputSpacing);
      
    // Allocate and initialize 
    outputImage->Allocate();
    outputImage->FillBuffer( 0 );
  
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  try
  {
    PARSE_ARGS;
  }
  catch( const std::exception & err )
  {
    tube::ErrorMessage( err.what() );
    return EXIT_FAILURE;
  }
  
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputImageFile, argc, argv );
  
}
