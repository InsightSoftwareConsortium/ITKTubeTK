/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __tubeSegmentTubesUsingMinimalPath_h
#define __tubeSegmentTubesUsingMinimalPath_h

// ITK includes
#include <itkGroupSpatialObject.h>
#include <itkMacro.h>
#include <itkProcessObject.h>


// TubeTK includes
#include "itktubeSegmentTubesUsingMinimalPathFilter.h"
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class SegmentTubesUsingMinimalPath
 *
 *  \ingroup TubeTK
 */

template< unsigned int Dimension, class TInputPixel >
class SegmentTubesUsingMinimalPath:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef SegmentTubesUsingMinimalPath               Self;
  typedef itk::ProcessObject                         Superclass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef itk::tube::SegmentTubesUsingMinimalPathFilter
    < Dimension, TInputPixel >                       FilterType;

  typedef typename FilterType::InputImageType         InputImageType;
  typedef typename InputImageType::Pointer            InputImagePointer;
  typedef typename FilterType::InputSpatialObjectType TubeGroupType;
  typedef typename TubeGroupType::Pointer             TubeGroupPointer;
  typedef typename FilterType::PointType              PointType;


  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( SegmentTubesUsingMinimalPath, Object );

  /* Set target tubes */
  tubeWrapSetMacro( TargetTubeGroup, TubeGroupPointer, Filter );

  /** Set speed Image */
  tubeWrapSetMacro( SpeedImage, InputImagePointer, Filter );
  tubeWrapGetMacro( SpeedImage, InputImagePointer, Filter );

  /* Set radius image */
  tubeWrapSetMacro( RadiusImage, InputImagePointer, Filter );
  tubeWrapGetMacro( RadiusImage, InputImagePointer, Filter );

  /* Set start point for the path */
  tubeWrapSetMacro( StartPoint, PointType, Filter );

  /* Set end point for the path */
  tubeWrapSetMacro( EndPoint, PointType, Filter );

  /* Set if the extract path connects to the surface of the target tube */
  tubeWrapSetMacro( ConnectToTargetTubeSurface, bool, Filter );

  /* Set Optimization method parameters. */
  tubeWrapSetMacro( OptimizationMethod, std::string, Filter );
  tubeWrapSetMacro( OptimizerTerminationValue, double, Filter );
  tubeWrapSetMacro( OptimizerNumberOfIterations, int, Filter );
  tubeWrapSetMacro( OptimizerStepLengthFactor, double, Filter );
  tubeWrapSetMacro( OptimizerStepLengthRelax, double, Filter );

  /*Set radius extraction parameters. */
  tubeWrapSetMacro( StartRadius, double, Filter );
  tubeWrapSetMacro( MaxRadius, double, Filter );
  tubeWrapSetMacro( StepSizeForRadiusEstimation, double, Filter );
  tubeWrapGetMacro( CostAssociatedWithExtractedTube, double, Filter );
  /* Get the extracted minimum path tube */
  tubeWrapGetMacro( Output, TubeGroupPointer, Filter );

  void SetIntermediatePoints( std::vector< PointType > );
  /* Runs tubes to image conversion */
  tubeWrapUpdateMacro( Filter );

protected:
  SegmentTubesUsingMinimalPath( void );
  ~SegmentTubesUsingMinimalPath() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkSegmentTubesUsingMinimalPathFilter parameters **/
  SegmentTubesUsingMinimalPath( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename FilterType::Pointer m_Filter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeSegmentTubesUsingMinimalPath.hxx"
#endif

#endif // End !defined( __tubeSegmentTubesUsingMinimalPath_h )
