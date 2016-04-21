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
#ifndef __tubeComputeBinaryImageSimilarityMetrics_h
#define __tubeComputeBinaryImageSimilarityMetrics_h

#include "itktubeLabelOverlapMeasuresImageFilter.h"
#include "itkObject.h"
#include "itkNumericTraits.h"

namespace tube
{
/** \class ComputeBinaryImageSimilarityMetrics
 *
 *  \ingroup TubeTKITK
 */

template< typename TInputImage >
class ComputeBinaryImageSimilarityMetrics:
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef ComputeBinaryImageSimilarityMetrics     Self;
  typedef itk::SmartPointer< Self >               Pointer;
  typedef itk::SmartPointer< const Self >         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ComputeBinaryImageSimilarityMetrics, Object);

  typedef TInputImage                                      ImageType;
  typedef typename TInputImage::PixelType                  LabelType;
  /** Type to use form computations. */
 // typedef typename NumericTraits< LabelType >::RealType    RealType;

  /** Set the source image. */
  void SetSourceImage( const ImageType * image );

  /** Set the target image. */
  void SetTargetImage( const ImageType * image );

  /** measures over all labels */
  float GetTotalOverlap( void ) const;
  float GetUnionOverlap( void ) const;
  float GetMeanOverlap( void ) const;
  float GetSimilarity( void ) const;
  float GetFalseNegativeError( void ) const;
  float GetFalsePositiveError( void ) const;

  void Update();

protected:
  ComputeBinaryImageSimilarityMetrics( void );
  ~ComputeBinaryImageSimilarityMetrics() {}
  void PrintSelf(std::ostream & os, itk::Indent indent) const;

private:
  /** itktubeLabelOverlapMeasuresImageFilter parameters **/
  ComputeBinaryImageSimilarityMetrics(const Self &);
  void operator=(const Self &);

  typedef itk::tube::LabelOverlapMeasuresImageFilter< ImageType >
    MetricFilterType;
  typename MetricFilterType::Pointer m_MetricFilter;

};
} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeComputeBinaryImageSimilarityMetrics.hxx"
#endif

#endif // End !defined( __tubeComputeBinaryImageSimilarityMetrics_h )
