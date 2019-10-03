/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __itktubeSegmentTubeUsingMinimalPathFilter_hxx
#define __itktubeSegmentTubeUsingMinimalPathFilter_hxx

#include "itktubeSegmentTubeUsingMinimalPathFilter.h"

// MinimalPathExtraction Imports
#include "itkSpeedFunctionToPathFilter.h"
#include "itkIterateNeighborhoodOptimizer.h"
#include "itkSingleImageCostFunction.h"

namespace itk
{
namespace tube
{
template< unsigned int Dimension, class TInputPixel >
SegmentTubeUsingMinimalPathFilter< Dimension, TInputPixel >
::SegmentTubeUsingMinimalPathFilter( void )
{
  m_SpeedImage = NULL;
  m_RadiusImage = NULL;
  m_TargetTubeGroup = NULL;
  m_ConnectToTargetTubeSurface = false;
  m_OptimizationMethod = "Regular_Step_Gradient_Descent";
  m_OptimizerTerminationValue = 2.0;
  m_OptimizerNumberOfIterations = 2000;
  m_OptimizerStepLengthFactor = 0.1;
  m_OptimizerStepLengthRelax = 0.999;
  m_StartRadius = 1;
  m_MaxRadius = 6;
  m_StepSizeForRadiusEstimation = 0.5;
  m_CostAssociatedWithExtractedTube = 0.0;
  m_Output = NULL;
}

template< unsigned int Dimension, class TInputPixel >
void
SegmentTubeUsingMinimalPathFilter< Dimension, TInputPixel >
::SetIntermediatePoints( std::vector< PointType > pathPoints )
{
  m_IntermediatePoints = pathPoints;
}

template< unsigned int Dimension, class TInputPixel >
void
SegmentTubeUsingMinimalPathFilter< Dimension, TInputPixel >
::Update( void )
{
  typedef itk::PolyLineParametricPath< Dimension > PathType;
  typedef itk::SpeedFunctionToPathFilter
    < InputImageType, PathType > PathFilterType;
  typedef itk::LinearInterpolateImageFunction< InputImageType, double >
    InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

  typedef itk::SingleImageCostFunction< InputImageType > CostFunctionType;
  typename CostFunctionType::Pointer costFunction = CostFunctionType::New();
  costFunction->SetInterpolator( interpolator );

  //Get Input image information
  typedef typename TubeType::TransformType TransformType;
  typename TransformType::InputVectorType scaleVector;
  typename TransformType::OffsetType offsetVector;
  typename InputImageType::SpacingType spacing = m_SpeedImage->GetSpacing();
  typename InputImageType::PointType origin = m_SpeedImage->GetOrigin();

  for( unsigned int k = 0; k < Dimension; ++k )
    {
    scaleVector[k] = spacing[k];
    offsetVector[k] = origin[k];
    }
  // Create path information
  typedef itk::SpeedFunctionPathInformation< PointType > PathInformationType;
  typename PathInformationType::Pointer pathInfo = PathInformationType::New();
  pathInfo->SetStartPoint( m_StartPoint );
  for( unsigned int i = 0; i < m_IntermediatePoints.size(); i++ )
    {
    pathInfo->AddWayPoint( m_IntermediatePoints[i] );
    }
  if( m_TargetTubeGroup )
    {
    PointType pointPath;
    this->IsPointTooNear( m_TargetTubeGroup, m_StartPoint, pointPath );
    pathInfo->SetEndPoint( pointPath );
    }
  else
    {
    pathInfo->SetEndPoint( m_EndPoint );
    }

  // Create path filter
  typename PathFilterType::Pointer pathFilter = PathFilterType::New();
  pathFilter->SetInput( m_SpeedImage );
  pathFilter->SetCostFunction( costFunction );
  pathFilter->SetTerminationValue( m_OptimizerTerminationValue );
  pathFilter->AddPathInformation( pathInfo );

  // Set Optimizer
  if( m_OptimizationMethod == "Iterate_Neighborhood" )
    {
    // Create IterateNeighborhoodOptimizer
    typedef itk::IterateNeighborhoodOptimizer OptimizerType;
    typename OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->MinimizeOn();
    optimizer->FullyConnectedOn();
    typename OptimizerType::NeighborhoodSizeType size( Dimension );
    for( unsigned int i = 0; i < Dimension; i++ )
      {
      size[i] = m_SpeedImage->GetSpacing()[i] * m_OptimizerStepLengthFactor;
      }
    optimizer->SetNeighborhoodSize( size );
    pathFilter->SetOptimizer( optimizer );
    }
  else if( m_OptimizationMethod == "Gradient_Descent" )
    {
    // Create GradientDescentOptimizer
    typedef itk::GradientDescentOptimizer OptimizerType;
    typename OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetNumberOfIterations( m_OptimizerNumberOfIterations );
    pathFilter->SetOptimizer( optimizer );
    }
  else if( m_OptimizationMethod == "Regular_Step_Gradient_Descent" )
    {
    // Compute the minimum spacing
    double minspacing = spacing[0];
    for( unsigned int dim = 0; dim < Dimension; dim++ )
      {
      if( spacing[dim] < minspacing )
        {
        minspacing = spacing[dim];
        }
      }
    // Create RegularStepGradientDescentOptimizer
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typename OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetNumberOfIterations( m_OptimizerNumberOfIterations );
    optimizer->SetMaximumStepLength
      ( 1.0 * m_OptimizerStepLengthFactor * minspacing );
    optimizer->SetMinimumStepLength
      ( 0.5 * m_OptimizerStepLengthFactor * minspacing );
    optimizer->SetRelaxationFactor( m_OptimizerStepLengthRelax );
    pathFilter->SetOptimizer( optimizer );
    }
  else
    {
    return;
    }

  try
    {
    pathFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::stringstream out;
    out << "ExceptionObject caught !" << std::endl;
    out << err << std::endl;
    return;
    }

  // Create output TRE file
  m_Output =InputSpatialObjectType::New();

  // Update tubes transform
  m_Output->GetObjectToParentTransform()->Scale(
    scaleVector );
  m_Output->GetObjectToParentTransform()->SetOffset(
    offsetVector );
  m_Output->GetObjectToParentTransform()->SetMatrix(
    m_SpeedImage->GetDirection() );
  m_Output->Update();
  m_CostAssociatedWithExtractedTube = 0.0;
  for( unsigned int i = 0; i < pathFilter->GetNumberOfOutputs(); i++ )
    {
    // Get the path
    typename PathType::Pointer path = pathFilter->GetOutput( i );
    // Check path is valid
    if( path->GetVertexList()->Size() == 0 )
      {
      std::cout << "WARNING: Path " << ( i + 1 )
        << " contains no points!" << std::endl;
      continue;
      }

    // Output centerline in TRE file
    typename TubeType::TubePointListType tubePointList;
    typename PathType::VertexListType * vertexList = path->GetVertexList();
    for( unsigned int k = 0; k < vertexList->Size(); k++ )
      {
      PointType pathPoint;
      m_SpeedImage->TransformContinuousIndexToPhysicalPoint(
        vertexList->GetElement( k ), pathPoint );
      typename InputImageType::IndexType imageIndex;
      if( m_SpeedImage->TransformPhysicalPointToIndex
        ( pathPoint, imageIndex ) )
        {
        m_CostAssociatedWithExtractedTube +=
          m_SpeedImage->GetPixel( imageIndex );
        }
      if( m_ConnectToTargetTubeSurface )
        {
        PointType nearPoint;
        bool isNear = this->IsPointTooNear
          ( m_TargetTubeGroup, pathPoint, nearPoint );
        if( isNear )
          {
          continue;
          }
        }
      TubePointType tubePoint;
      tubePoint.SetPositionInObjectSpace( vertexList->GetElement( k ) );
      tubePoint.SetId( k );
      tubePointList.push_back( tubePoint );
      }
    typename TubeType::Pointer pTube = TubeType::New();
    pTube->SetPoints( tubePointList );
    pTube->ComputeTangentAndNormals();
    pTube->SetId( i );

    // Extract Radius
    if( m_RadiusImage )
      {
      typedef itk::tube::RadiusExtractor2< InputImageType >
        RadiusExtractorType;
      typename RadiusExtractorType::Pointer radiusExtractor
        = RadiusExtractorType::New();
      radiusExtractor->SetInputImage( m_RadiusImage );
      radiusExtractor->SetRadiusStart( m_StartRadius );
      radiusExtractor->SetRadiusMin( 0.2 );
      radiusExtractor->SetRadiusMax( m_MaxRadius );
      radiusExtractor->SetRadiusStep( m_StepSizeForRadiusEstimation );
      radiusExtractor->SetRadiusTolerance( 0.025 );
      radiusExtractor->SetDebug( false );
      radiusExtractor->ExtractRadii( pTube );
      }

    m_Output->AddChild( pTube );
    m_Output->Update();
    }
}

template< unsigned int Dimension, class TInputPixel >
bool
SegmentTubeUsingMinimalPathFilter< Dimension, TInputPixel >
::IsPointTooNear( const InputSpatialObjectType * sourceTubeGroup,
              PointType outsidePoint,
              PointType &nearestPoint )
{

  double minDistance = itk::NumericTraits<double>::max();
  double nearestPointRadius = 0.0;

  typename InputSpatialObjectType::ChildrenListPointer sourceTubeList =
    sourceTubeGroup->GetChildren();
  for( typename InputSpatialObjectType::ChildrenListType::iterator
    tubeList_it = sourceTubeList->begin();
    tubeList_it != sourceTubeList->end(); ++tubeList_it )
    {
    //**** Source Tube **** :
    typename TubeType::Pointer pCurSourceTube =
      dynamic_cast< TubeType* >( tubeList_it->GetPointer() );
    //dynamic_cast verification
    if( !pCurSourceTube )
      {
      return EXIT_FAILURE;
      }
    pCurSourceTube->Update();
    //Get points in current source tube
    typename TubeType::TubePointListType pointList =
      pCurSourceTube->GetPoints();
    //Get Index to World Transformation
    typename TubeType::TransformType * pTubeObjectPhysTransform =
      pCurSourceTube->GetObjectToWorldTransform();
    for( typename TubeType::TubePointListType::const_iterator
      pointList_it = pointList.begin();
      pointList_it != pointList.end(); ++pointList_it )
      {
      TubePointType curSourcePoint = *pointList_it;
      //Transform parameters in physical space
      typename TubePointType::PointType curSourcePos =
        pTubeObjectPhysTransform->TransformPoint(
          curSourcePoint.GetPositionInObjectSpace() );
      double distance =
        curSourcePos.SquaredEuclideanDistanceTo( outsidePoint );
      if( minDistance > distance )
        {
        minDistance = distance;
        for( unsigned int i = 0; i < Dimension; i++ )
          {
          nearestPoint[i] = curSourcePos[i];
          }
        nearestPointRadius = curSourcePoint.GetRadiusInObjectSpace();
        }
      }
    }
  if( minDistance < nearestPointRadius*nearestPointRadius )
    {
    return true;
    }
  else
    {
    return false;
    }
}

template< unsigned int Dimension, class TInputPixel >
void
SegmentTubeUsingMinimalPathFilter< Dimension, TInputPixel >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

}

} // end namespace tube
} // end namespace itk

#endif // End !defined( __itktubeSegmentTubeUsingMinimalPathFilter_hxx )
