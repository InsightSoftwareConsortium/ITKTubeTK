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
//#include "itkNumericTraits.h"
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class ComputeBinaryImageSimilarityMetrics
 *
 *  \ingroup TubeTKITK
 */

template< class TInputImage >
class ComputeBinaryImageSimilarityMetrics:
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef ComputeBinaryImageSimilarityMetrics     Self;
  typedef itk::SmartPointer< Self >               Pointer;
  typedef itk::SmartPointer< const Self >         ConstPointer;
  typedef TInputImage                             InputImageType;
  typedef typename TInputImage::PixelType         LabelType;

  typedef itk::tube::LabelOverlapMeasuresImageFilter< InputImageType >
    FilterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ComputeBinaryImageSimilarityMetrics, Object);

  /** Set the source image. */
  tubeWrapSetConstObjectMacro( SourceImage, InputImageType, Filter );

  /** Set the target image. */
  tubeWrapSetConstObjectMacro( TargetImage, InputImageType, Filter );

  /** measures over all labels */
  tubeWrapGetMacro( TotalOverlap, float, Filter );
  tubeWrapGetMacro( UnionOverlap, float, Filter );
  tubeWrapGetMacro( MeanOverlap, float, Filter );
  tubeWrapGetMacro( Similarity, float, Filter );
  tubeWrapGetMacro( FalseNegativeError, float, Filter );
  tubeWrapGetMacro( FalsePositiveError, float, Filter );

  /** Compute image similarity */
  tubeWrapUpdateMacro(Filter);

protected:
  ComputeBinaryImageSimilarityMetrics( void );
  ~ComputeBinaryImageSimilarityMetrics() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itktubeLabelOverlapMeasuresImageFilter parameters **/
  ComputeBinaryImageSimilarityMetrics( const Self & );
  void operator=( const Self & );

  typename FilterType::Pointer m_Filter;

};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeComputeBinaryImageSimilarityMetrics.hxx"
#endif

#endif // End !defined( __tubeComputeBinaryImageSimilarityMetrics_h )
