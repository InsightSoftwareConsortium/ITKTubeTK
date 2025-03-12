/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef tubeSegmentTubeUsingMinimalPath_h
#define tubeSegmentTubeUsingMinimalPath_h

// ITK includes
#include <itkGroupSpatialObject.h>
#include <itkMacro.h>
#include <itkProcessObject.h>


// TubeTK includes
#include "itktubeSegmentTubeUsingMinimalPathFilter.h"
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class SegmentTubeUsingMinimalPath
 *
 *  \ingroup TubeTK
 */

template <unsigned int Dimension, class TInputPixel>
class SegmentTubeUsingMinimalPath : public itk::ProcessObject
{
public:
  /** Standard class type alias. */
  using Self = SegmentTubeUsingMinimalPath;
  using Superclass = itk::ProcessObject;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using FilterType = itk::tube::SegmentTubeUsingMinimalPathFilter<Dimension, TInputPixel>;

  using InputImageType = typename FilterType::InputImageType;
  using InputImagePointer = typename InputImageType::Pointer;
  using TubeGroupType = typename FilterType::InputSpatialObjectType;
  using TubeGroupPointer = typename TubeGroupType::Pointer;
  using PointType = typename FilterType::PointType;


  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(SegmentTubeUsingMinimalPath);

  /* Set target tubes */
  tubeWrapSetMacro(TargetTubeGroup, TubeGroupPointer, Filter);

  /** Set speed Image */
  tubeWrapSetMacro(SpeedImage, InputImagePointer, Filter);
  tubeWrapGetMacro(SpeedImage, InputImagePointer, Filter);

  /* Set radius image */
  tubeWrapSetMacro(RadiusImage, InputImagePointer, Filter);
  tubeWrapGetMacro(RadiusImage, InputImagePointer, Filter);

  /* Set start point for the path */
  tubeWrapSetMacro(StartPoint, PointType, Filter);

  /* Set end point for the path */
  tubeWrapSetMacro(EndPoint, PointType, Filter);

  /* Set if the extract path connects to the surface of the target tube */
  tubeWrapSetMacro(ConnectToTargetTubeSurface, bool, Filter);

  /* Set Optimization method parameters. */
  tubeWrapSetMacro(OptimizationMethod, std::string, Filter);
  tubeWrapSetMacro(OptimizerTerminationValue, double, Filter);
  tubeWrapSetMacro(OptimizerNumberOfIterations, int, Filter);
  tubeWrapSetMacro(OptimizerStepLengthFactor, double, Filter);
  tubeWrapSetMacro(OptimizerStepLengthRelax, double, Filter);

  /*Set radius extraction parameters. */
  tubeWrapSetMacro(StartRadius, double, Filter);
  tubeWrapSetMacro(MaxRadius, double, Filter);
  tubeWrapGetMacro(CostAssociatedWithExtractedTube, double, Filter);
  /* Get the extracted minimum path tube */
  tubeWrapGetMacro(Output, TubeGroupPointer, Filter);

  void SetIntermediatePoints(std::vector<PointType>);
  /* Runs tubes to image conversion */
  tubeWrapUpdateMacro(Filter);

protected:
  SegmentTubeUsingMinimalPath(void);
  ~SegmentTubeUsingMinimalPath() {}
  void
  PrintSelf(std::ostream & os, itk::Indent indent) const override;

private:
  /** itkSegmentTubeUsingMinimalPathFilter parameters **/
  SegmentTubeUsingMinimalPath(const Self &);
  void
  operator=(const Self &);

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void
  SetInput(const DataObjectIdentifierType &, itk::DataObject *) override {};

  typename FilterType::Pointer m_Filter;
};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#  include "tubeSegmentTubeUsingMinimalPath.hxx"
#endif

#endif // End !defined( __tubeSegmentTubeUsingMinimalPath_h )
