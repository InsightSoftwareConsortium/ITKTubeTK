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
#ifndef __tubeComputeImageSimilarityMetrics_h
#define __tubeComputeImageSimilarityMetrics_h

// ITK includes
#include <itkProcessObject.h>

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "itktubeComputeImageSimilarityMetrics.h"

namespace tube
{
/** \class ComputeImageSimilarityMetrics
 *  \brief Computes similarity between two images using correlation or
 * mutual information
 *
 *  \ingroup TubeTKITK
 */

template< class TInputImage >
class ComputeImageSimilarityMetrics:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ComputeImageSimilarityMetrics              Self;
  typedef itk::ProcessObject                         Superclass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef itk::tube::ComputeImageSimilarityMetrics<
    TInputImage >                                    FilterType;

  typedef typename FilterType::ImageType             ImageType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ComputeImageSimilarityMetrics, Object );

  /** Set/Get use of correlation or mutual information to compute similarity */
  tubeWrapSetMacro( UseCorrelation, bool, Filter );
  tubeWrapGetMacro( UseCorrelation, bool, Filter );

  /** Set/Get portion of voxels used to compute image similarity */
  tubeWrapSetMacro( SamplingRate, double, Filter );
  tubeWrapGetMacro( SamplingRate, double, Filter );

  /** Set/Get input image 1 */
  tubeWrapSetConstObjectMacro( Input1, ImageType, Filter );
  tubeWrapGetConstObjectMacro( Input1, ImageType, Filter );

  /** Set/Get input image 2 */
  tubeWrapSetConstObjectMacro( Input2, ImageType, Filter );
  tubeWrapGetConstObjectMacro( Input2, ImageType, Filter );

  /** Compute image similarity */
  tubeWrapUpdateMacro( Filter );

  /** Get image similarity */
  tubeWrapGetMacro( Output, double, Filter );

protected:
  ComputeImageSimilarityMetrics( void );
  ~ComputeImageSimilarityMetrics() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkComputeImageSimilarityMetricsFilter parameters **/
  ComputeImageSimilarityMetrics( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename FilterType::Pointer m_Filter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeComputeImageSimilarityMetrics.hxx"
#endif

#endif // End !defined( __tubeComputeImageSimilarityMetrics_h )
