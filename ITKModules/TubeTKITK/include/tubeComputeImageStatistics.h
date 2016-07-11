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
#ifndef __tubeComputeImageStatistics_h
#define __tubeComputeImageStatistics_h

// ITK includes
#include <itkImage.h>
#include <itkImageToImageFilter.h>

// TubeTK includes
#include "itktubeComputeImageStatistics.h"
#include "tubeWrappingMacros.h"

namespace tube
{

/** \class ComputeImageStatistics
 * \brief Computes image statistics
 *
 * \ingroup TubeTKITK
 */

template< class TPixel, unsigned int VDimension >
class ComputeImageStatistics
  : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef ComputeImageStatistics           Self;
  typedef itk::Object                      Superclass;
  typedef itk::SmartPointer< Self >        Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  typedef itk::tube::ComputeImageStatistics< TPixel, VDimension > FilterType;
  typedef typename FilterType::MaskType                           MaskType;
  typedef typename FilterType::VolumeType                         VolumeType;


  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ComputeImageStatistics, Object );

  /** Get statistics components */
  tubeWrapGetMacro( CompMean, itk::Array< double >, Filter );
  tubeWrapGetMacro( CompMin, itk::Array< double >, Filter );
  tubeWrapGetMacro( CompMax, itk::Array< double >, Filter );
  tubeWrapGetMacro( CompStdDev, itk::Array< double >, Filter );
  tubeWrapGetMacro( CompCount, itk::Array< double >, Filter );
  tubeWrapGetMacro( CompValue, itk::Array< TPixel >, Filter );
  tubeWrapGetMacro( NumberOfComponents, unsigned int, Filter );

  /** Set/Get quantiles */
  tubeWrapSetMacro( Quantiles, std::vector<float>, Filter );
  tubeWrapGetMacro( Quantiles, std::vector<float>, Filter );

  /** Set/Get input mask */
  tubeWrapSetObjectMacro( InputMask, MaskType, Filter );
  tubeWrapGetObjectMacro( InputMask, MaskType, Filter );

  /** Compute image statistics */
  tubeWrapUpdateMacro(Filter);

  /** Write statistics to a CSV formatted file */
  void WriteCSVStatistics( std::string csvStatisticsFile ) const;

  /** Set/Get input image */
  tubeWrapSetConstObjectMacro( Input, VolumeType, Filter );
  tubeWrapGetConstObjectMacro( Input, VolumeType, Filter );

  /** Get output image */
  tubeWrapGetObjectMacro(Output, VolumeType, Filter);

protected:
  ComputeImageStatistics( void );
  ~ComputeImageStatistics() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  ComputeImageStatistics(const Self &);
  void operator=(const Self &);

  typename FilterType::Pointer m_Filter;


}; // End class ComputeImageStatistics

} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeComputeImageStatistics.hxx"
#endif

#endif // End !defined(__tubeComputeImageStatistics_h)
