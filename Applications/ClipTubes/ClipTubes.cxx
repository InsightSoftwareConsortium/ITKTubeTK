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
#include <iostream>
#include <sstream>

#include "itkTimeProbesCollectorBase.h"
#include "tubeMessage.h"

#include "tubeMacro.h"
#include "metaScene.h"

#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkGroupSpatialObject.h"

#include "ClipTubesCLP.h"

template< unsigned int VDimension >
bool isInside (itk::Point<double,VDimension> pointPos, double tubeRadius,
		  std::vector<double> boxPos, std::vector<double> boxSize)
{
  // Return a boolean indicating if any slice of the tube is included in the box
  // A slice is considered as a point and an associated radius
  bool hasXInside=false,hasYInside=false,hasZInside=false;
  
  if( (pointPos[0]+tubeRadius)>=boxPos[0] 
      && (pointPos[0]-tubeRadius)<=(boxPos[0]+boxSize[0]) )
    {
    hasXInside=true;
    }
  if( pointPos[1]<=boxPos[1] 
      && pointPos[1]>=(boxPos[1]-boxSize[1]) )
    {
    hasYInside=true;
    }
  switch( VDimension )
    {
    case 2:
      hasZInside=true;
      break;

    case 3:
      if( (pointPos[2]+tubeRadius)>=boxPos[2] 
	  && (pointPos[2]-tubeRadius)<=(boxPos[2]+boxSize[2]) )
	{
	hasZInside=true;
	}
      break;

    default:
      tubeErrorMacro(<< "Error: Only 2D and 3D data is currently supported.");
      return EXIT_FAILURE;
    }
 
  return (hasXInside && hasYInside && hasZInside);
}


template< unsigned int VDimension >
int DoIt (int argc, char * argv[])
{
  PARSE_ARGS;

  // Ensure that the input image dimension is valid
  // We only support 2D and 3D Images due to the
  // limitation of itkTubeSpatialObject
  if( VDimension != 2 && VDimension != 3 )
    {
    tube::ErrorMessage("Error: Only 2D and 3D data is currently supported.");
    return EXIT_FAILURE;
    }
    
  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;
  
  // Load TRE File
  tubeStandardOutputMacro( << "\n>> Loading TRE File" );

  typedef itk::SpatialObjectReader< VDimension > TubesReaderType;
  typedef itk::GroupSpatialObject< VDimension >  TubeGroupType;
//   typedef itk::TubeSpatialObject< VDimension >  TubeType;    //WARNING : dynamic_cast on "typedef itk::TubeSpatialObject< VDimension > TubeType" (before GetPoints() call) causes SEGFAULT 
  typedef itk::VesselTubeSpatialObject< VDimension >  TubeType; //          so TubeType corresponds to VesselTubeType to prevent issues
  typedef itk::VesselTubeSpatialObjectPoint< VDimension >  TubePointType;// so TubePointType corresponds to VesselTubePointType WARNING

  timeCollector.Start( "Loading Input TRE File" );
  
  typename TubesReaderType::Pointer tubeFileReader = TubesReaderType::New();
  
  try
    {
     tubeFileReader->SetFileName( inputTREFile.c_str() );
     tubeFileReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
     tube::ErrorMessage( "Error loading TRE File: "
                          + std::string( err.GetDescription() ) );
     timeCollector.Report();
     return EXIT_FAILURE;
    }
  
  typename TubeGroupType::Pointer pSourceTubeGroup = tubeFileReader->GetGroup();
  typename TubeGroupType::ChildrenListPointer pSourceTubeList =
  pSourceTubeGroup->GetChildren();
  
  timeCollector.Stop( "Loading Input TRE File" );
  
  // Compute clipping   
  tubeStandardOutputMacro( << "\n>> Finding Tubes for Clipping" );
  
  timeCollector.Start( "Selecting Tubes" );
  
  typename TubeGroupType::Pointer pTargetTubeGroup = TubeGroupType::New();//Target Group to save desired tubes
  
  int targetTubeId=0;
  
  for( typename TubeGroupType::ChildrenListType::iterator
       tubeList_it = pSourceTubeList->begin();
       tubeList_it != pSourceTubeList->end(); ++tubeList_it)
    { 
      //**** Source Tube **** :
      typename TubeType::Pointer pCurSourceTube = 
      dynamic_cast< TubeType* >( tubeList_it->GetPointer() ); 
      //dynamic_cast verification
      if(!pCurSourceTube)
	return EXIT_FAILURE;
      
      typename TubeType::PointListType TargetPointList;//Point List for TargetTube
      
      //Get points in current source tube
      typename TubeType::PointListType pointList = pCurSourceTube->GetPoints(); 
      for( typename TubeType::PointListType::const_iterator
	   pointList_it = pointList.begin();
	   pointList_it != pointList.end(); ++pointList_it)
	{
	  TubePointType curSourcePoint = *pointList_it;
	  typename TubePointType::PointType curSourcePos = 
	  curSourcePoint.GetPosition();
	  //Save point in target tube if it belongs to the box  
	  if(isInside(curSourcePos,curSourcePoint.GetRadius(),boxCorner,boxSize))
	  { 
	    if(ClipTubes)
	    {
	      TargetPointList.push_back(curSourcePoint);
	    }
	    else
	    {
	       pCurSourceTube->SetId(targetTubeId);
	       ++targetTubeId;
	       pTargetTubeGroup->AddSpatialObject(pCurSourceTube);
	       break;
	    }  
	  }  	 
	  else
	  {
	    if(TargetPointList.size()>0)
	    {
	      //**** Target Tube **** :
	      typename TubeType::Pointer pTargetTube = TubeType::New();//Target Tube to add to Group
	      pTargetTube->SetId(targetTubeId);
	      ++targetTubeId;
	      //Save clipped tube
	      pTargetTube->SetPoints(TargetPointList);
	      pTargetTubeGroup->AddSpatialObject(pTargetTube);
	      
	      TargetPointList.clear();
	    }
	  }
	} 
     
    }
    
  timeCollector.Stop( "Selecting Tubes" );
  
  // Write output TRE file
  tubeStandardOutputMacro(
    << "\n>> Writing TRE file" );

  timeCollector.Start( "Writing output TRE file" );
  
  typedef itk::SpatialObjectWriter< VDimension > TubeWriterType;
  typename TubeWriterType::Pointer tubeWriter = TubeWriterType::New();

  try
    {
    tubeWriter->SetFileName( outputTREFile.c_str() );
    tubeWriter->SetInput(pTargetTubeGroup);
    tubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing TRE file: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Writing output TRE file" );
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
  
  if(boxCorner.empty() || boxSize.empty())
  {
    tube::ErrorMessage(
      "Error: longflags --boxCorner and --boxSize are both required");
    return EXIT_FAILURE;
  }
  
  MetaScene *mScene = new MetaScene;
  mScene->Read( inputTREFile.c_str() );
  
  if( mScene->GetObjectList()->empty() )
    {
    tubeWarningMacro( << "Input TRE file has no spatial objects" );
    return EXIT_SUCCESS;
    }

  switch( mScene->GetObjectList()->front()->NDims() )
    {
    case 2:
      return DoIt<2>( argc, argv );
      break;

    case 3:
      return DoIt<3>( argc, argv );
      break;

    default:
      tubeErrorMacro(
	<< "Error: Only 2D and 3D data is currently supported.");
      return EXIT_FAILURE;
    }
}