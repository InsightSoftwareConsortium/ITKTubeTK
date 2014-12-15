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
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMinimumMaximumImageFilter.h"

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

  // Ensure that the input image dimension is valid
  // We only support 2D and 3D Images due to the 
  // limitation of itkTubeSpatialObject
  if(VDimension != 2 || VDimension != 3)
  {
    tube::ErrorMessage( "Error: Only 2D and 3D Images are currently supported." );
    return EXIT_FAILURE;
  }
  
  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;
  
  // Load Input Image
  typedef TPixel                                        PixelType;
  typedef itk::Image< PixelType, VDimension >           ImageType;
  typedef itk::ImageFileReader< ImageType >             ImageReaderType;
  
  timeCollector.Start("Loading Input Image");
  
  typename ImageReaderType::Pointer pImageReader = ImageReaderType::New();
  
  try
  {
    pImageReader->SetFileName( inputImageFile.c_str() );    
    pImageReader->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    tube::ErrorMessage( "Error loading input image: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
  }
  
  typename ImageType::Pointer pInputImage = pImageReader->GetOutput();
  
  timeCollector.Stop("Loading Input Image");
  
  // Load TRE File
  typedef itk::SpatialObjectReader< VDimension >      	TubesReaderType;  
  typedef itk::GroupSpatialObject< VDimension >        	TubeGroupType;

  timeCollector.Start("Loading Input TRE File");
  
  typename TubesReaderType::Pointer pTubeFileReader = TubesReaderType::New();
  
  try
  {
    pTubeFileReader->SetFileName( inputTREFile.c_str() );
    pTubeFileReader->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    tube::ErrorMessage( "Error loading TRE File: "
			+ std::string(err.GetDescription()));
    timeCollector.Report();
    return EXIT_FAILURE;
  }

  typename TubeGroupType::Pointer pTubeGroup = pTubeFileReader->GetGroup(); 
  
  timeCollector.Stop("Loading Input TRE File");

  // Find the user-specified tube
  typedef typename TubeGroupType::ChildrenListType	TubeListType;
  typedef itk::TubeSpatialObject< VDimension >   	TubeType;
  typename TubeType::Pointer pInputTube;
  bool blnTubeFound = false;
  
  timeCollector.Start("Finding the user specified tube");
  
  TubeListType * pTubeList = pTubeGroup->GetChildren();
  typename TubeListType::const_iterator itTubes = pTubeList->begin();
  
  while( itTubes != pTubeList->end() )
  {
      if( (*itTubes)->GetId() == inputTubeId )
      {	
	pInputTube = dynamic_cast<TubeType *>( (*itTubes).GetPointer() );
	blnTubeFound = true;
	break;
      }    
      
      itTubes++;
  }
  
  delete pTubeList;  
  
  if( !blnTubeFound )
  {
      tube::ErrorMessage("Unable to find the specified tube");
      timeCollector.Report();
      return EXIT_FAILURE;
  }
  
  pInputTube->RemoveDuplicatePoints();
  pInputTube->ComputeTangentAndNormals();
  
  timeCollector.Stop("Finding the user specified tube");  
  
  // Get list of tube points
  typedef typename TubeType::PointListType		TubePointListType;     
  
  TubePointListType tubePointList = pInputTube->GetPoints();
  if( tubePointList.size() <= 0 )
  {
      tube::ErrorMessage("The specified tube does not contain any points");
      timeCollector.Report();
      return EXIT_FAILURE;
  }
    
  // Determine maximum radius among all tube points 
  typename TubePointListType::const_iterator itPts = tubePointList.begin();
  double maxTubeRadius = itPts->GetRadius();
  
  while( itPts != tubePointList.end() )
  {
      if( itPts->GetRadius() > maxTubeRadius )
      {
	  maxTubeRadius = itPts->GetRadius();	
      }    
      
      itPts++;
  }  
  
  // Determine the mean distance between consecutive tube points
  double meanTubePointDist = 0; 
  
  for(unsigned int pid = 1; pid < tubePointList.size(); pid++)
  {
    // compute distance between current and previous tube point
    double curDist = 0;
    typename TubeType::PointType p1 = tubePointList[pid-1].GetPosition();
    typename TubeType::PointType p2 = tubePointList[pid].GetPosition();     
    
    for(unsigned int i = 0; i < VDimension; i++)
    {
      curDist += (p2[i] - p1[i]) * (p2[i] - p1[i]);
    }
    
    // update mean
    meanTubePointDist += curDist;
  }

  meanTubePointDist /= tubePointList.size();
  
  // Determine minimum spacing of the input image
  typename ImageType::SpacingType inputSpacing = pInputImage->GetSpacing();  
  double minInputSpacing = inputSpacing[0];
  
  for(unsigned int i = 1; i < VDimension; i++)
  {
    if( inputSpacing[i] < minInputSpacing )
    {
      minInputSpacing = inputSpacing[i];
    }  
  }
    
  // Setup the output image and allocate space for it
  typename ImageType::Pointer pOutputImage = ImageType::New();  

    // set spacing
    // For last dimension its set to mean consecutive point distance  
    // For other dimensions its set to the minimum input spacing 
    typename ImageType::SpacingType outputSpacing;
    
    for(unsigned int i = 0; i < VDimension-1; i++)
    {
      outputSpacing[i] = minInputSpacing;
    }
    
    outputSpacing[VDimension-1] = meanTubePointDist;
    
    pOutputImage->SetSpacing(outputSpacing);
      
    // set start index
    typename ImageType::IndexType startIndex;
    startIndex.Fill(0);

    // set size
    typename ImageType::SizeType size;  
    for(unsigned int i = 0; i < VDimension-1; i++)
    {
	size[i] = 2 * (typename ImageType::SizeValueType)
		      (0.5 + maxTubeRadius / outputSpacing[i]) + 1;
    }
    size[VDimension-1] = tubePointList.size();	

    // set regions
    typename ImageType::RegionType region;  
    region.SetIndex( startIndex );  
    region.SetSize( size );    
    pOutputImage->SetRegions( region );  
    
    // Allocate and initialize 
    pOutputImage->Allocate();
    pOutputImage->FillBuffer( 0 );
  
  // For each tube point, extract normal plane image
  // and fill into corresponding slice in the output image
  typedef typename TubeType::TubePointType TubePointType;
  typedef typename TubePointType::CovariantVectorType TubeNormalType;  
  
  typedef itk::ImageRegionIteratorWithIndex< ImageType > ImageIteratorType;  
  
  typedef itk::LinearInterpolateImageFunction< ImageType, double > InterpolatorType;

  typedef itk::MinimumMaximumImageFilter< ImageType > MinMaxImageFilterType;

  timeCollector.Start("Generating output image");

  unsigned int ptInd = 0;

  typename MinMaxImageFilterType::Pointer minmaxFilter = MinMaxImageFilterType::New();
  minmaxFilter->SetInput( pInputImage );
  minmaxFilter->Update();
  TPixel outsideVal = minmaxFilter->GetMinimum();
  
  typename InterpolatorType::Pointer pInterpolator = InterpolatorType::New(); 
  pInterpolator->SetInputImage( pInputImage );
  
  for( itPts = tubePointList.begin(); 
	itPts != tubePointList.end(); 
	itPts++, ptInd++ 
      )
  { 
    // Get position, radius and frenet-serret basis of current tube point
    // in the world coordinate system
    typename TubeType::PointType curTubePosition = (*itPts).GetPosition();
    double curTubeRadius = (*itPts).GetRadius();
    TubeNormalType curTubeNormal1 = (*itPts).GetNormal1();
    TubeNormalType curTubeNormal2 = (*itPts).GetNormal2();
    
    // Define slice region in the output image
    typename ImageType::RegionType sliceRegion;
    
      typename ImageType::IndexType sliceStartIndex;
      sliceStartIndex.Fill(0);
      sliceStartIndex[VDimension-1] = ptInd;
      
      typename ImageType::SizeType sliceSize;
      sliceSize = pOutputImage->GetLargestPossibleRegion().GetSize();
      sliceSize[VDimension-1] = 0;
    
    sliceRegion.SetIndex( sliceStartIndex );
    sliceRegion.SetSize( sliceSize );
    
    // Iterate through corresponding slice of output image and fill each pixel
    ImageIteratorType itOutSlice(pOutputImage, sliceRegion);
    
    for( itOutSlice.GoToBegin(); !itOutSlice.IsAtEnd(); ++itOutSlice )
    {
      // get index of the current output pixel
      typename ImageType::IndexType curOutIndex = itOutSlice.GetIndex();
      
      // compute corresponding position in the input image
      typename ImageType::PointType curInputPoint = curTubePosition;
      double distToCenter = 0;
      
      if(VDimension == 2) 
      {
	double stepN1 = curOutIndex[0] - 0.5 * sliceSize[0] * outputSpacing[0];
	
	for(unsigned int i = 0; i < VDimension; i++)
	{
	    curInputPoint[i] += stepN1 * curTubeNormal1[i]; 	  
	}
	
	distToCenter = stepN1;
      }
      else if (VDimension == 3)
      {
	double stepN1 = curOutIndex[0] - 0.5 * sliceSize[0] * outputSpacing[0];
	double stepN2 = curOutIndex[1] - 0.5 * sliceSize[1] * outputSpacing[1];
	
	for(unsigned int i = 0; i < VDimension; i++)
	{
	    curInputPoint[i] += stepN1 * curTubeNormal1[i]; 	  
	    curInputPoint[i] += stepN2 * curTubeNormal2[i]; 	  
	}	
	
	distToCenter = std::sqrt( stepN1 * stepN1 + stepN2 * stepN2  );
      }
      
      // Get the intensity value from the input image using interpolation
      if( distToCenter <= curTubeRadius &&
	  pInterpolator->IsInsideBuffer(curInputPoint) )
      {
	itOutSlice.Set( pInterpolator->Evaluate(curInputPoint) );
      }
      else
      {
	itOutSlice.Set( outsideVal );
      }
    }
  }  
    
  timeCollector.Stop("Generating output image");
  
  // Write output image
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;
  
  timeCollector.Start("Writing output image");
  
  typename ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  
  try
  {
    imageWriter->SetFileName( outputImageFile.c_str() );    
    imageWriter->SetInput( pOutputImage );
    imageWriter->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    tube::ErrorMessage( "Error writing output image: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
  }
  
  timeCollector.Stop("Writing output image");
  
  // All done
  timeCollector.Report();
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
