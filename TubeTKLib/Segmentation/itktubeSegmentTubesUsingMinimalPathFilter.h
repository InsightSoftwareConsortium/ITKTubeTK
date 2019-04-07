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
#ifndef __itktubeSegmentTubesUsingMinimalPathFilter_h
#define __itktubeSegmentTubesUsingMinimalPathFilter_h

#include "itkImage.h"

//ITK imports
#include "itkGroupSpatialObject.h"
#include "itkTubeSpatialObjectPoint.h"
#include "itkTubeSpatialObject.h"
#include "itkPathIterator.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkPolyLineParametricPath.h"
#include "itkGradientDescentOptimizer.h"
#include "itkNumericTraits.h"
#include "itkObject.h"
//TubeTK imports
#include "itktubeRadiusExtractor2.h"

namespace itk
{
namespace tube
{
/** \class SegmentTubesUsingMinimalPathFilter
 * \brief Segmenta a tube between the start and
 * end point using an input speed image.
 * This filter uses itk::minimumPathExtraction filter to per the minimum
 * path between the end points.
 *
 */

template< unsigned int Dimension, class TInputPixel >
class SegmentTubesUsingMinimalPathFilter: public Object
{
public:
  /** Standard class typedefs. */
  typedef SegmentTubesUsingMinimalPathFilter Self;
  typedef Object                             Superclass;
  typedef SmartPointer< Self >               Pointer;
  typedef SmartPointer< const Self >         ConstPointer;

  typedef itk::Image< TInputPixel, Dimension > InputImageType;
  typedef itk::GroupSpatialObject< Dimension > InputSpatialObjectType;

  typedef typename InputSpatialObjectType::Pointer TubeGroupPointer;

  typedef itk::Point< double, Dimension >           PointType;
  typedef itk::TubeSpatialObjectPoint< Dimension >  TubePointType;
  typedef itk::TubeSpatialObject< Dimension >       TubeType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro
    ( SegmentTubesUsingMinimalPathFilter, Object );

  itkSetMacro( SpeedImage, typename InputImageType::Pointer );
  itkGetMacro( SpeedImage, typename InputImageType::Pointer );
  itkSetMacro( RadiusImage, typename InputImageType::Pointer );
  itkGetMacro( RadiusImage, typename InputImageType::Pointer );
  itkSetMacro( StartPoint, PointType );
  itkGetMacro( StartPoint, PointType );
  itkSetMacro( EndPoint, PointType );
  itkGetMacro( EndPoint, PointType );
  itkSetMacro( ConnectToTargetTubeSurface, bool );
  itkGetMacro( ConnectToTargetTubeSurface, bool );
  itkSetMacro( OptimizationMethod, std::string );
  itkGetMacro( OptimizationMethod, std::string );
  itkSetMacro( OptimizerTerminationValue, double );
  itkGetMacro( OptimizerTerminationValue, double );
  itkSetMacro( OptimizerNumberOfIterations, int );
  itkGetMacro( OptimizerNumberOfIterations, int );
  itkSetMacro( OptimizerStepLengthFactor, double );
  itkGetMacro( OptimizerStepLengthFactor, double );
  itkSetMacro( OptimizerStepLengthRelax, double );
  itkGetMacro( OptimizerStepLengthRelax, double );
  itkSetMacro( StartRadius, double );
  itkGetMacro( StartRadius, double );
  itkSetMacro( MaxRadius, double );
  itkGetMacro( MaxRadius, double );
  itkSetMacro( StepSizeForRadiusEstimation, double );
  itkGetMacro( StepSizeForRadiusEstimation, double );
  itkGetMacro( CostAssociatedWithExtractedTube, double );
  itkSetMacro( CostAssociatedWithExtractedTube, double );
  /** Sets the input tubes */
  itkSetMacro( TargetTubeGroup, TubeGroupPointer );
  itkGetMacro( TargetTubeGroup, TubeGroupPointer );
  itkGetMacro( Output, TubeGroupPointer );

  void SetIntermediatePoints( std::vector< PointType > );
  void Update( void );
protected:
  SegmentTubesUsingMinimalPathFilter( void );
  ~SegmentTubesUsingMinimalPathFilter() {}

  void PrintSelf( std::ostream & os, Indent indent ) const;
  bool IsPointTooNear( const InputSpatialObjectType * sourceTubeGroup,
              PointType outsidePoint,
              PointType &nearestPoint );
private:
  SegmentTubesUsingMinimalPathFilter( const Self & );
  void operator=( const Self & );

  typename InputImageType::Pointer  m_SpeedImage;
  typename InputImageType::Pointer  m_RadiusImage;
  PointType                         m_StartPoint;
  PointType                         m_EndPoint;
  std::vector< PointType >          m_IntermediatePoints;
  TubeGroupPointer                  m_TargetTubeGroup;
  bool                              m_ConnectToTargetTubeSurface;
  std::string                       m_OptimizationMethod;
  double                            m_OptimizerTerminationValue;
  int                               m_OptimizerNumberOfIterations;
  double                            m_OptimizerStepLengthFactor;
  double                            m_OptimizerStepLengthRelax;
  double                            m_StartRadius;
  double                            m_MaxRadius;
  double                            m_StepSizeForRadiusEstimation;
  double                            m_CostAssociatedWithExtractedTube;
  TubeGroupPointer                  m_Output;

}; //End class SegmentTubesUsingMinimalPathFilter
} // End namespace tube
} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSegmentTubesUsingMinimalPathFilter.hxx"
#endif

#endif // End !defined( __itktubeSegmentTubesUsingMinimalPathFilter_h )
