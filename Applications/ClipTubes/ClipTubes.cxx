
#include <iostream>
#include <sstream>


#include "tubeMessage.h"

#include "tubeMacro.h"
#include "metaScene.h"

#include "itkSpatialObjectReader.h"
#include "itkGroupSpatialObject.h"
// #include "itkImageFileReader.h"
// #include "itkTubeSpatialObject.h"
// #include "itkTubeSpatialObjectPoint.h"
// #include "itkVesselTubeSpatialObject.h"
// #include "itkVesselTubeSpatialObjectPoint.h"



#include "ClipTubesCLP.h"

template< unsigned int VDimension >
int DoIt (int argc, char * argv[])
{
  PARSE_ARGS;

  // Load TRE File
  tubeStandardOutputMacro( << "\n>> Loading TRE File" );

  typedef itk::SpatialObjectReader< VDimension > TubesReaderType;//ToDo: Template DoIt<VDimension>
  typedef itk::GroupSpatialObject< VDimension >  TubeGroupType;

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

   
   
  std::cout << "Hello TubeTK World !" << std::endl;
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