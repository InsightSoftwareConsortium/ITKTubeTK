
#include <iostream>
#include <sstream>


#include "tubeMessage.h"

#include "tubeMacro.h"
#include "metaScene.h"

#include "itkSpatialObjectReader.h"
#include "itkGroupSpatialObject.h"
// #include "itkTubeSpatialObject.h"               // WARNING include?
// #include "itkTubeSpatialObjectPoint.h"
// #include "itkVesselTubeSpatialObject.h"
// #include "itkVesselTubeSpatialObjectPoint.h"

#include "ClipTubesCLP.h"

template< unsigned int VDimension >
bool isInside (itk::Point<double,VDimension> pointPos, double tubeRadius, std::vector<double> boxPos, std::vector<double> boxSize)
{
  // Return a boolean indicating if a slice of the tube is included in the box
  // A slice of the tube is considered as a point and an associated radius
  bool hasXInside=false,hasYInside=false,hasZInside=false;
  
  if((pointPos[0]+tubeRadius)>=boxPos[0] && (pointPos[0]-tubeRadius)<=(boxPos[0]+boxSize[0]))
    hasXInside=true;
  if(pointPos[1]<=boxPos[1] && pointPos[1]>=(boxPos[1]-boxSize[1]))
    hasYInside=true;
  if((pointPos[2]+tubeRadius)>=boxPos[2] && (pointPos[2]-tubeRadius)<=(boxPos[2]+boxSize[2]))
    hasZInside=true;
  
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
  
  // Load TRE File
  tubeStandardOutputMacro( << "\n>> Loading TRE File" );

  typedef itk::SpatialObjectReader< VDimension > TubesReaderType;
  typedef itk::GroupSpatialObject< VDimension >  TubeGroupType;
//   typedef itk::TubeSpatialObject< VDimension >  TubeType;    //WARNING : dynamic_cast on "typedef itk::TubeSpatialObject< VDimension > TubeType" (before GetPoints() call) causes SEGFAULT 
  typedef itk::VesselTubeSpatialObject< VDimension >  TubeType; //          so TubeType corresponds to VesselTubeType to prevent issues
  typedef itk::VesselTubeSpatialObjectPoint< VDimension >  TubePointType;// so TubePointType corresponds to VesselTubePointType WARNING

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
     return EXIT_FAILURE;
    }
    
  // Get source tubes   
  typename TubeGroupType::Pointer pTubeGroup = tubeFileReader->GetGroup();
  typename TubeGroupType::ChildrenListPointer pTubeList = pTubeGroup->GetChildren();

  for( typename TubeGroupType::ChildrenListType::iterator
       tubeList_it = pTubeList->begin();
       tubeList_it != pTubeList->end(); ++tubeList_it)
    {
      typename TubeType::Pointer pCurSourceTube = dynamic_cast< TubeType* >( tubeList_it->GetPointer() ); 
      //dynamic_cast verification
      if(!pCurSourceTube)
	return EXIT_FAILURE;
      //Get points in current source tube
      typename TubeType::PointListType pointList = pCurSourceTube->GetPoints(); 
      for( typename TubeType::PointListType::const_iterator
	   pointList_it = pointList.begin();
	   pointList_it != pointList.end(); ++pointList_it)
	{
	  TubePointType curSourcePoint = *pointList_it;
	  typename TubePointType::PointType curSourcePos = curSourcePoint.GetPosition();
	  if(isInside(curSourcePos,curSourcePoint.GetRadius(),boxCorner,boxSize))
	  {
	    std::cout<<pCurSourceTube->GetId()<<std::endl;
	    break;
	  }  
	 
	}
    }
  
  
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
      tubeErrorMacro(<< "Error: Only 2D and 3D data is currently supported.");
      return EXIT_FAILURE;
    }
}