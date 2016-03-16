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
#include "tubeMacro.h"
#include "metaScene.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkPolyLineParametricPath.h"
#include "itkSingleImageCostFunction.h"
#include "itkGradientDescentOptimizer.h"
#include "itkPathIterator.h"
#include "itkGroupSpatialObject.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itktubeRadiusExtractor2.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"

#include "itkSpeedFunctionToPathFilter.h"
#include "itkSpeedFunctionPathInformation.h"
#include "itkIterateNeighborhoodOptimizer.h"

#include <sstream>

#include "SegmentTubeUsingMinimalPathCLP.h"

template< class TPixel, unsigned int DimensionT >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h"
//   and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int DimensionT >
bool
IsPointTooNear( typename itk::GroupSpatialObject< DimensionT >::Pointer &sourceTubeGroup,
              itk::Point< double, DimensionT > outsidePoint,
              itk::Point< double, DimensionT > &nearestPoint,
              double thresholdDistance)
{
  typedef itk::GroupSpatialObject< DimensionT >           TubeGroupType;
  typedef itk::VesselTubeSpatialObject< DimensionT >      TubeType;
  typedef itk::VesselTubeSpatialObjectPoint< DimensionT > TubePointType;

  double minDistance = itk::NumericTraits<double>::max();
  double nearestPointRadius = 0.0;

  typename TubeGroupType::ChildrenListPointer sourceTubeList =
    sourceTubeGroup->GetChildren();
  for( typename TubeGroupType::ChildrenListType::iterator
    tubeList_it = sourceTubeList->begin();
    tubeList_it != sourceTubeList->end(); ++tubeList_it)
    {
    //**** Source Tube **** :
    typename TubeType::Pointer pCurSourceTube =
      dynamic_cast< TubeType* >( tubeList_it->GetPointer() );
    //dynamic_cast verification
    if(!pCurSourceTube)
      {
      return EXIT_FAILURE;
      }
    pCurSourceTube->ComputeObjectToWorldTransform();
    //Get points in current source tube
    typename TubeType::PointListType pointList =
      pCurSourceTube->GetPoints();
    //Get Index to World Transformation
    typename TubeType::TransformType * pTubeIndexPhysTransform =
      pCurSourceTube->GetIndexToWorldTransform();
    for( typename TubeType::PointListType::const_iterator
      pointList_it = pointList.begin();
      pointList_it != pointList.end(); ++pointList_it )
      {
      TubePointType curSourcePoint = *pointList_it;
      //Transform parameters in physical space
      typename TubePointType::PointType curSourcePos =
        pTubeIndexPhysTransform->TransformPoint(
          curSourcePoint.GetPosition() );
      double distance = curSourcePos.SquaredEuclideanDistanceTo( outsidePoint );
      if( minDistance > distance )
        {
        minDistance = distance;
        for( unsigned int i = 0; i<DimensionT; i++ )
          {
          nearestPoint[i] = curSourcePos[i];
          }
        nearestPointRadius = curSourcePoint.GetRadius();
        }
      }
    }
  if( thresholdDistance > 0 )
    {
    if( minDistance < thresholdDistance*thresholdDistance)
      {
      return true;
      }
    else
      {
      return false;
      }
    }
  if( minDistance < nearestPointRadius*nearestPointRadius)
    {
    return true;
    }
  else
    {
    return false;
    }
}

template< unsigned int DimensionT >
void WriteFCSVFile( std::string filename,
  std::vector< itk::Point< double, DimensionT > > &pointList )
{
  std::fstream of;
  of.open(filename.c_str(), std::fstream::out);
  if( !of.is_open() )
    {
    std::cout << "WriteData: unable to open file " << filename.c_str()
      << " for writing";
    return;
    }
  // put down a header
  of << "# Markups fiducial file version = " << "4.5 \n";
  of << "# CoordinateSystem = " << 0 << "\n";
  // label the columns
  // id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID
  // orientation is a quaternion, angle and axis
  // associatedNodeID and description and label can be empty strings
  // id,x,y,z,ow,ox,oy,oz,vis,sel,lock,,,
  // label can have spaces, everything up to next comma is used, no quotes
  // necessary, same with the description
  of << "# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,"
    "associatedNodeID" << "\n";

  typedef itk::Point< double, DimensionT >                PointType;

  for( typename std::vector< PointType >::iterator it = pointList.begin();
    it != pointList.end(); ++it )
    {
    PointType curPoint = *it;
    of << ","; //As we don't know the id.
    of << -1 * curPoint[0]  << "," << -1 * curPoint[1] << ",";
    if( DimensionT == 3 )
      {
      of << curPoint[2] << ",";
      }
    of << "0,0,0,0,"; // As we don't know the orientation
    of << "1,"; //visibility is true
    of << "1,"; //selection is true
    of << ",";
    of << ",,";
    of << "\n";
    }
  of.close();
}

template< class TPixel, unsigned int DimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;
  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter( "Extract Minimal Path",
    CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                                    PixelType;
  typedef itk::Image< PixelType, DimensionT >       ImageType;
  typedef itk::ImageFileReader< ImageType >         ReaderType;

  typedef itk::SpatialObjectReader< DimensionT >          TubesReaderType;
  typedef itk::GroupSpatialObject< DimensionT >           TubeGroupType;
  typedef itk::VesselTubeSpatialObject< DimensionT >      TubeType;
  typedef itk::VesselTubeSpatialObjectPoint< DimensionT > TubePointType;

  timeCollector.Start( "Load data" );

  //Read input Image
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( InputImage.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::stringstream out;
    out << "ExceptionObject caught !" << std::endl;
    out << err << std::endl;
    tube::ErrorMessage( out.str() );
    timeCollector.Stop( "Load data" );
    return EXIT_FAILURE;
    }

  typename ImageType::Pointer speed = reader->GetOutput();
  speed->DisconnectPipeline();

  //Read radius extraction Image
  if( !RadiusImage.empty() )
    {
    reader->SetFileName( RadiusImage.c_str() );
    try
      {
      reader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::stringstream out;
      out << "ExceptionObject caught !" << std::endl;
      out << err << std::endl;
      tube::ErrorMessage( out.str() );
      timeCollector.Stop( "Load data" );
      return EXIT_FAILURE;
      }
    }
  typename ImageType::Pointer radiusExtractorInput = reader->GetOutput();

  //Read input centerline
  typename TubesReaderType::Pointer tubeFileReader = TubesReaderType::New();
  if( !InputPathFile.empty() )
    {
    try
      {
      tubeFileReader->SetFileName( InputPathFile.c_str() );
      tubeFileReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Error loading TRE File: "
        + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    }

  timeCollector.Stop( "Load data" );
  progressReporter.Report( 0.1 );

  timeCollector.Start( "Set parameters" );
  typedef itk::PolyLineParametricPath<DimensionT>               PathType;
  typedef itk::SpeedFunctionToPathFilter< ImageType, PathType > PathFilterType;
  typedef itk::Point< double, DimensionT >                      PointType;

  typedef itk::LinearInterpolateImageFunction< ImageType, double >
    InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

  typedef itk::SingleImageCostFunction< ImageType > CostFunctionType;
  typename CostFunctionType::Pointer costFunction = CostFunctionType::New();
  costFunction->SetInterpolator( interpolator );

  //Get Input image information
  typedef typename TubeType::TransformType TransformType;
  typename TransformType::InputVectorType scaleVector;
  typename TransformType::OffsetType offsetVector;
  typename ImageType::SpacingType spacing = speed->GetSpacing();
  typename ImageType::PointType origin = speed->GetOrigin();
  double tubeSpacing[DimensionT];

  for ( unsigned int k = 0; k < DimensionT; ++k )
    {
    scaleVector[k] = spacing[k];
    offsetVector[k] = origin[k];
    tubeSpacing[k] = spacing[k];
    }

  // Create path information
  typedef itk::SpeedFunctionPathInformation< PointType > PathInformationType;
  typename PathInformationType::Pointer pathInfo = PathInformationType::New();

  if( Path.size() >= 2 )
    {
    for( unsigned int k = 0; k<Path.size(); k++ )
      {
      PointType path;
      for( unsigned int i = 0; i<DimensionT; i++ )
        {
        path[i]=Path[k][i];
        }
      if( k == 0 )
        {
        pathInfo->SetStartPoint( path );
        }
      else if( k >= Path.size() - 1 )
        {
        pathInfo->SetEndPoint( path );
        }
      else
        {
        pathInfo->AddWayPoint( path );
        }
      }
    }
  else if ( Path.size() == 1)
    {
    PointType startPositionPoint;
    PointType path;
    for( unsigned int i = 0; i<DimensionT; i++ )
      {
      path[i] = Path[0][i];
      startPositionPoint[i] = Path[0][i];
      }
    pathInfo->SetStartPoint( path );
    typename TubeGroupType::Pointer sourceTubeGroup =
      tubeFileReader->GetGroup();
    PointType pointPath;
    IsPointTooNear< TPixel, DimensionT >
      ( sourceTubeGroup, startPositionPoint, pointPath, -1 );
    pathInfo->SetEndPoint( pointPath );
    }
  else
    {
    tubeErrorMacro(
      << "Error: Path should contain at least a Start and an End Point" );
    timeCollector.Stop( "Set parameters" );
    return EXIT_FAILURE;
    }

  // Create path filter
  typename PathFilterType::Pointer pathFilter = PathFilterType::New();
  pathFilter->SetInput( speed );
  pathFilter->SetCostFunction( costFunction );
  pathFilter->SetTerminationValue( TerminationValue );
  pathFilter->AddPathInformation( pathInfo );

  // Set Optimizer
  if( Optimizer == "Iterate_Neighborhood" )
    {
    // Create IterateNeighborhoodOptimizer
    typedef itk::IterateNeighborhoodOptimizer OptimizerType;
    typename OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->MinimizeOn();
    optimizer->FullyConnectedOn();
    typename OptimizerType::NeighborhoodSizeType size( DimensionT );
    for( unsigned int i = 0; i < DimensionT; i++ )
      {
      size[i] = speed->GetSpacing()[i] * StepLengthFactor;
      }
    optimizer->SetNeighborhoodSize( size );
    pathFilter->SetOptimizer( optimizer );
    }
  else if( Optimizer == "Gradient_Descent" )
    {
    // Create GradientDescentOptimizer
    typedef itk::GradientDescentOptimizer OptimizerType;
    typename OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetNumberOfIterations( NumberOfIterations );
    pathFilter->SetOptimizer( optimizer );
    }
  else if( Optimizer == "Regular_Step_Gradient_Descent" )
    {
    // Compute the minimum spacing
    double minspacing = spacing[0];
    for( unsigned int dim = 0; dim < DimensionT; dim++ )
      {
      if ( spacing[dim] < minspacing )
        {
        minspacing = spacing[dim];
        }
      }
    // Create RegularStepGradientDescentOptimizer
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typename OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetNumberOfIterations( NumberOfIterations );
    optimizer->SetMaximumStepLength( 1.0 * StepLengthFactor*minspacing );
    optimizer->SetMinimumStepLength( 0.5 * StepLengthFactor*minspacing );
    optimizer->SetRelaxationFactor( StepLengthRelax );
    pathFilter->SetOptimizer( optimizer );
    }
  else
    {
    tubeErrorMacro(
      << "Error: Optimizer not known" );
    timeCollector.Stop( "Set parameters" );
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Set parameters" );
  progressReporter.Report( 0.2 );

  timeCollector.Start( "Extract minimal path" );

  try
    {
    pathFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::stringstream out;
    out << "ExceptionObject caught !" << std::endl;
    out << err << std::endl;
    tube::ErrorMessage( out.str() );
    timeCollector.Stop( "Extract minimal path" );
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Extract minimal path" );
  progressReporter.Report( 0.4 );

  // Create output TRE file
  typename TubeGroupType::Pointer pTubeGroup = TubeGroupType::New();

  // Update tubes transform
  pTubeGroup->GetObjectToParentTransform()->SetScale(
    scaleVector );
  pTubeGroup->GetObjectToParentTransform()->SetOffset(
    offsetVector );
  pTubeGroup->GetObjectToParentTransform()->SetMatrix(
    speed->GetDirection() );
  pTubeGroup->ComputeObjectToWorldTransform();

  timeCollector.Start( "Rasterizing path" );

  for ( unsigned int i = 0; i<pathFilter->GetNumberOfOutputs(); i++ )
    {
    // Get the path
    typename PathType::Pointer path = pathFilter->GetOutput( i );
    // Check path is valid
    if ( path->GetVertexList()->Size() == 0 )
      {
      std::cout << "WARNING: Path " << ( i + 1 )
        << " contains no points!" << std::endl;
      continue;
      }

    // Output centerline in TRE file
    typename TubeType::PointListType tubePointList;
    typename PathType::VertexListType * vertexList = path->GetVertexList();
    double cost = 0;
    for( unsigned int k = 0; k < vertexList->Size(); k++ )
      {
      PointType pathPoint;
      speed->TransformContinuousIndexToPhysicalPoint(
        vertexList->GetElement(k), pathPoint );
      typename ImageType::IndexType imageIndex;
      if ( speed->TransformPhysicalPointToIndex( pathPoint, imageIndex ) )
        {
        cost = cost + speed->GetPixel( imageIndex );
        }
      if( !InputPathFile.empty() && HardBoundary )
        {
        typename TubeGroupType::Pointer sourceTubeGroup =
        tubeFileReader->GetGroup();
        PointType nearPoint;
        bool isNear = IsPointTooNear< TPixel, DimensionT >
          ( sourceTubeGroup, pathPoint, nearPoint, Distance );
        if( isNear )
          {
          continue;
          }
        else
          {
          HardBoundary = false;
          }
        }
      for( unsigned int d = 0; d < DimensionT; d++ )
        {
        pathPoint[d]=(pathPoint[d]-origin[d])/spacing[d];
        }
      TubePointType tubePoint;
      tubePoint.SetPosition( pathPoint );
      tubePoint.SetID( k );
      if( ExtractRadiusUsingInputImage )
        {
        tubePoint.SetRadius( speed->GetPixel( imageIndex ) );
        }
      tubePointList.push_back( tubePoint );
      }
    typename TubeType::Pointer pTube = TubeType::New();
    pTube->SetPoints( tubePointList );
    pTube->ComputeTangentAndNormals();
    pTube->SetSpacing( tubeSpacing );
    pTube->SetId( i );

    // Extract Radius
    if( !ExtractRadiusUsingInputImage && !RadiusImage.empty() )
      {
      typedef itk::tube::RadiusExtractor2<ImageType> RadiusExtractorType;
      typename RadiusExtractorType::Pointer radiusExtractor
        = RadiusExtractorType::New();
      radiusExtractor->SetInputImage( radiusExtractorInput );
      radiusExtractor->SetRadiusStart( StartRadius );
      radiusExtractor->SetRadiusMin( 0.2 );
      radiusExtractor->SetRadiusMax( MaxRadius );
      radiusExtractor->SetRadiusStep( StepRadius );
      radiusExtractor->SetRadiusTolerance( 0.025 );
      radiusExtractor->SetDebug( false );
      radiusExtractor->ExtractRadii( pTube );
      }

    pTubeGroup->AddSpatialObject( pTube );
    pTubeGroup->ComputeObjectToWorldTransform();
    std::cout << "Cost associated with centerline: " << cost
        << std::endl;
    }

  timeCollector.Stop( "Rasterizing path" );
  progressReporter.Report( 0.9 );

  timeCollector.Start( "Write output data" );

  // Write output TRE file
  typedef itk::SpatialObjectWriter< DimensionT > TubeWriterType;
  typename TubeWriterType::Pointer tubeWriter = TubeWriterType::New();
  try
    {
    tubeWriter->SetFileName( OutputTREFile.c_str() );
    tubeWriter->SetInput( pTubeGroup );
    tubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing TRE file: "
      + std::string( err.GetDescription() ) );
    timeCollector.Stop( "Write output data" );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  if( !PathEndPoints.empty() )
    {
    std::vector< PointType > endPointList;
    typename TubeGroupType::ChildrenListPointer pChildrenList =
      pTubeGroup->GetChildren();
    for( typename TubeGroupType::ChildrenListType::iterator
    tubeList_it = pChildrenList->begin();
    tubeList_it != pChildrenList->end(); ++tubeList_it )
      {
      //**** Source Tube **** :
      typename TubeType::Pointer pCurSourceTube =
        dynamic_cast< TubeType* >( tubeList_it->GetPointer() );
      //dynamic_cast verification
      if( !pCurSourceTube )
        {
        return EXIT_FAILURE;
        }
      pCurSourceTube->ComputeObjectToWorldTransform();
      typename TubeType::PointListType pointList =
        pCurSourceTube->GetPoints();
      TubePointType startSourcePoint = pointList.front();
      TubePointType endSourcePoint = pointList.back();
      typename TubePointType::PointType startSourcePosIndexSpace =
        pCurSourceTube->GetIndexToWorldTransform()->TransformPoint(
        startSourcePoint.GetPosition() );
      typename TubePointType::PointType endSourcePosIndexSpace =
        pCurSourceTube->GetIndexToWorldTransform()->TransformPoint(
        endSourcePoint.GetPosition() );
      PointType startPoint;
      for( unsigned int i = 0; i < DimensionT; i++ )
        {
        startPoint[i] = startSourcePosIndexSpace[i];
        }
      PointType endPoint;
      for( unsigned int i = 0; i < DimensionT; i++ )
        {
        endPoint[i] = endSourcePosIndexSpace[i];
        }
      endPointList.push_back( startPoint );
      endPointList.push_back( endPoint );
      }
    WriteFCSVFile( PathEndPoints, endPointList );
    }

  timeCollector.Stop( "Write output data" );
  progressReporter.Report( 1.0 );

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputImage.
  return tube::ParseArgsAndCallDoIt( InputImage, argc, argv );

}
