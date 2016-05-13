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

#ifndef __itktubeComputeImageSimilarityMetrics_h
#define __itktubeComputeImageSimilarityMetrics_h

// ITK includes
#include <itkObject.h>
#include <itkObjectFactory.h>

namespace itk
{

namespace tube
{

/** \class ComputeImageSimilarityMetrics
 * \brief Computes similarity between two images using correlation or
 * mutual information
 */

template< class TInputImage >
class ComputeImageSimilarityMetrics
  : public Object
{
public:

  /** Standard class typedefs. */
  typedef ComputeImageSimilarityMetrics               Self;
  typedef Object                                      SuperClass;
  typedef SmartPointer< Self >                        Pointer;
  typedef SmartPointer< const Self >                  ConstPointer;

  /** custom typedefs */
  typedef TInputImage                                 ImageType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ComputeImageSimilarityMetrics, Object );

  /** Set/Get use of correlation or mutual information to compute similarity */
  itkSetMacro( UseCorrelation, bool );
  itkGetMacro( UseCorrelation, bool );

  /** Set/Get sampling rate */
  itkSetMacro( SamplingRate, double );
  itkGetMacro( SamplingRate, double );

  /** Set/Get input image 1 */
  itkSetConstObjectMacro( Input1, ImageType );
  itkGetConstObjectMacro( Input1, ImageType );

  /** Set/Get input image 2 */
  itkSetConstObjectMacro( Input2, ImageType );
  itkGetConstObjectMacro( Input2, ImageType );

  /** Compute image similarity */
  void Update( void );

  /** Get image similarity */
  itkGetMacro(Output, double);

protected:

  ComputeImageSimilarityMetrics( void );
  ~ComputeImageSimilarityMetrics( void ) {};

  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  typename ImageType::ConstPointer            m_Input1;
  typename ImageType::ConstPointer            m_Input2;
  bool                                        m_UseCorrelation;
  double                                      m_SamplingRate;
  double                                      m_Output;

}; // End class ComputeImageSimilarityMetrics

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeComputeImageSimilarityMetrics.hxx"
#endif

#endif // End !defined(__itktubeComputeImageSimilarityMetrics_h)
