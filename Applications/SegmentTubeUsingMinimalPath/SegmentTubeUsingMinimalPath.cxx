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
              itk::Point< double, DimensionT > &nearestPoint )
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
  if( minDistance < nearestPointRadius*nearestPointRadius)
    {
    return true;
    }
  else
    {
    return false;
    }
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
  if( !TargetTubeFileName.empty() )
    {
    try
      {
      tubeFileReader->SetFileName( TargetTubeFileName.c_str() );
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
  PointType startPathPoint;
  if( StartPoint.size() == 1 )
    {
    for( unsigned int i = 0; i < DimensionT; i++ )
      {
      startPathPoint[i]=StartPoint[0][i];
      }
    pathInfo->SetStartPoint( startPathPoint );
    }
  else
    {
    tubeErrorMacro(
      << "Error: Path must contain at only one Start Point" );
    timeCollector.Stop( "Set parameters" );
    return EXIT_FAILURE;
    }
  if( IntermediatePoints.size() >= 1 )
  {
  for( unsigned int k = 0; k < IntermediatePoints.size(); k++ )
    {
    PointType pathPoint;
    for( unsigned int i = 0; i < DimensionT; i++ )
      {
       pathPoint[i]=IntermediatePoints[k][i];
      }
    pathInfo->AddWayPoint( pathPoint );
    }
  }
  if( EndPoint.size() == 1 )
    {
    PointType pathPoint;
    for( unsigned int i = 0; i < DimensionT; i++ )
      {
      pathPoint[i]=EndPoint[0][i];
      }
    pathInfo->SetEndPoint( pathPoint );
    }
  else if ( !TargetTubeFileName.empty() && EndPoint.size() == 0 )
    {
    typename TubeGroupType::Pointer sourceTubeGroup =
      tubeFileReader->GetGroup();
    PointType pointPath;
    IsPointTooNear< TPixel, DimensionT >
      ( sourceTubeGroup, startPathPoint, pointPath );
    pathInfo->SetEndPoint( pointPath );
    }
  else
    {
    tubeErrorMacro(
      << "Error: Atmost one End/Target Point or Target Tube should be provided. " );
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
      if( !TargetTubeFileName.empty() &&
        ConnectionOption == "Connect To Target Tube Surface" )
        {
        typename TubeGroupType::Pointer sourceTubeGroup =
        tubeFileReader->GetGroup();
        PointType nearPoint;
        bool isNear = IsPointTooNear< TPixel, DimensionT >
          ( sourceTubeGroup, pathPoint, nearPoint );
        if( isNear )
          {
          continue;
          }
        }
      for( unsigned int d = 0; d < DimensionT; d++ )
        {
        pathPoint[d]=(pathPoint[d]-origin[d])/spacing[d];
        }
      TubePointType tubePoint;
      tubePoint.SetPosition( pathPoint );
      tubePoint.SetID( k );
      tubePointList.push_back( tubePoint );
      }
    typename TubeType::Pointer pTube = TubeType::New();
    pTube->SetPoints( tubePointList );
    pTube->ComputeTangentAndNormals();
    pTube->SetSpacing( tubeSpacing );
    pTube->SetId( i );

    // Extract Radius
    if( !RadiusImage.empty() )
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
