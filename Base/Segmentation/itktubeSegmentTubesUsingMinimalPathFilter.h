/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
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
#include "itkLinearInterpolateImageFunction.h"
#include "itkPolyLineParametricPath.h"
#include "itkGradientDescentOptimizer.h"
#include "itkNumericTraits.h"

#include "itktubeSpatialObjectToSpatialObjectFilter.h"

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

template< class TInputSpatialObject, class TInputImage >
class SegmentTubesUsingMinimalPathFilter:
  public SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TInputSpatialObject >
{
public:
  /** Standard class typedefs. */
  typedef SegmentTubesUsingMinimalPathFilter              Self;
  typedef SpatialObjectToSpatialObjectFilter
    < TInputSpatialObject, TInputSpatialObject >          Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  typedef TInputSpatialObject                      InputSpatialObjectType;
  typedef TInputImage                              InputImageType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SegmentTubesUsingMinimalPathFilter, SpatialObjectToSpatialObjectFilter);

  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  typedef itk::Point< double, ImageDimension >                PointType;
  typedef itk::VesselTubeSpatialObjectPoint< ImageDimension > TubePointType;
  typedef itk::VesselTubeSpatialObject< ImageDimension >      TubeType;

  itkSetMacro( SpeedImage, typename InputImageType::Pointer );
  itkGetMacro( SpeedImage, typename InputImageType::Pointer );
  itkSetMacro( RadiusImage, typename InputImageType::Pointer );
  itkGetMacro( RadiusImage, typename InputImageType::Pointer );
  itkSetMacro( StartPoint, PointType );
  itkGetMacro( StartPoint, PointType );
  itkSetMacro( EndPoint, PointType );
  itkGetMacro( EndPoint, PointType );
  itkSetMacro( ExtractEndPointFromTargetTube, bool );
  itkGetMacro( ExtractEndPointFromTargetTube, bool );
  itkSetMacro( ConnectToTargetTubeSurface, bool );
  itkGetMacro( ConnectToTargetTubeSurface, bool );
  itkSetMacro( ExtractRadius, bool );
  itkGetMacro( ExtractRadius, bool );
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

  void SetIntermediatePoints( std::vector< PointType > );
protected:
  SegmentTubesUsingMinimalPathFilter( void );
  ~SegmentTubesUsingMinimalPathFilter() {}

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  virtual void GenerateData( void );
  bool IsPointTooNear( const InputSpatialObjectType * sourceTubeGroup,
              PointType outsidePoint,
              PointType &nearestPoint );
private:

  SegmentTubesUsingMinimalPathFilter(const Self &);
  void operator=(const Self &);

  typename InputImageType::Pointer  m_SpeedImage;
  typename InputImageType::Pointer  m_RadiusImage;
  PointType                         m_StartPoint;
  PointType                         m_EndPoint;
  std::vector< PointType >          m_IntermediatePoints;
  bool                              m_ExtractEndPointFromTargetTube;
  bool                              m_ConnectToTargetTubeSurface;
  bool                              m_ExtractRadius;
  std::string                       m_OptimizationMethod;
  double                            m_OptimizerTerminationValue;
  int                               m_OptimizerNumberOfIterations;
  double                            m_OptimizerStepLengthFactor;
  double                            m_OptimizerStepLengthRelax;
  double                            m_StartRadius;
  double                            m_MaxRadius;
  double                            m_StepSizeForRadiusEstimation;
  double                            m_CostAssociatedWithExtractedTube;

}; //End class SegmentTubesUsingMinimalPathFilter
} // End namespace tube
} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSegmentTubesUsingMinimalPathFilter.hxx"
#endif

#endif // End !defined( __itktubeSegmentTubesUsingMinimalPathFilter_h )
