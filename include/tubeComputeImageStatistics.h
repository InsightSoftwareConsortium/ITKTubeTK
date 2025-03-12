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
#ifndef tubeComputeImageStatistics_h
#define tubeComputeImageStatistics_h

// ITK includes
#include <itkProcessObject.h>
#include <itkImage.h>
#include <itkImageToImageFilter.h>

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "itktubeComputeImageStatistics.h"

namespace tube
{

/** \class ComputeImageStatistics
 * \brief Computes image statistics
 *
 * \ingroup TubeTK
 */

template <class TPixel, unsigned int VDimension>
class ComputeImageStatistics : public itk::ProcessObject
{
public:
  /** Standard class type alias. */
  using Self = ComputeImageStatistics;
  using Superclass = itk::ProcessObject;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using FilterType = itk::tube::ComputeImageStatistics<TPixel, VDimension>;
  using MaskType = typename FilterType::MaskType;
  using VolumeType = typename FilterType::VolumeType;


  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkTypeMacro(ComputeImageStatistics, ProcessObject);

  /** Get statistics components */
  tubeWrapGetMacro(CompMean, std::vector<double>, Filter);
  tubeWrapGetMacro(CompMin, std::vector<double>, Filter);
  tubeWrapGetMacro(CompMax, std::vector<double>, Filter);
  tubeWrapGetMacro(CompStdDev, std::vector<double>, Filter);
  tubeWrapGetMacro(CompCount, std::vector<double>, Filter);
  tubeWrapGetMacro(CompValue, std::vector<TPixel>, Filter);
  tubeWrapGetMacro(NumberOfComponents, unsigned int, Filter);

  /** Set/Get quantiles */
  tubeWrapSetMacro(Quantiles, std::vector<float>, Filter);
  tubeWrapGetMacro(Quantiles, std::vector<float>, Filter);

  /** Set/Get input mask */
  tubeWrapSetObjectMacro(InputMask, MaskType, Filter);
  tubeWrapGetObjectMacro(InputMask, MaskType, Filter);

  /** Compute image statistics */
  tubeWrapUpdateMacro(Filter);

  /** Write statistics to a CSV formatted file */
  void
  WriteCSVStatistics(std::string csvStatisticsFile) const;

  /** Set/Get input image */
  tubeWrapSetConstObjectMacro(Input, VolumeType, Filter);
  tubeWrapGetConstObjectMacro(Input, VolumeType, Filter);

  /** Get output image */
  tubeWrapGetObjectMacro(Output, VolumeType, Filter);

protected:
  ComputeImageStatistics(void);
  ~ComputeImageStatistics() {}

  void
  PrintSelf(std::ostream & os, itk::Indent indent) const;

private:
  ComputeImageStatistics(const Self &);

  void
  operator=(const Self &);

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void
  SetInput(const DataObjectIdentifierType &, itk::DataObject *) {};

  typename FilterType::Pointer m_Filter;


}; // End class ComputeImageStatistics

} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#  include "tubeComputeImageStatistics.hxx"
#endif

#endif // End !defined( __tubeComputeImageStatistics_h )
