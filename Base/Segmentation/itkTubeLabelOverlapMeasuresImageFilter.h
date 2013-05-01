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
#ifndef __itkTubeLabelOverlapMeasuresImageFilter_h
#define __itkTubeLabelOverlapMeasuresImageFilter_h

#include "itkInPlaceImageFilter.h"
#include "itkFastMutexLock.h"
#include "itkNumericTraits.h"

#include "itksys/hash_map.hxx"

namespace itk {

namespace tube {

/** \class LabelOverlapMeasuresImageFilter
 * \brief Computes overlap measures between the set same set of labels of
 * pixels of two images.  Background is assumed to be 0.
 *
 * \sa LabelOverlapMeasuresImageFilter
 *
 * \ingroup MultiThreaded
 */
template<class TLabelImage>
class ITK_EXPORT LabelOverlapMeasuresImageFilter :
    public InPlaceImageFilter< TLabelImage >
{
public:
  /** Standard Self typedef */
  typedef LabelOverlapMeasuresImageFilter                Self;
  typedef ImageToImageFilter<TLabelImage,TLabelImage>    Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Runtime information support. */
  itkTypeMacro( LabelOverlapMeasuresImageFilter, ImageToImageFilter );

  /** Image related typedefs. */
  typedef TLabelImage                                   LabelImageType;
  typedef typename TLabelImage::Pointer                 LabelImagePointer;
  typedef typename TLabelImage::ConstPointer            LabelImageConstPointer;

  typedef typename TLabelImage::RegionType              RegionType;
  typedef typename TLabelImage::SizeType                SizeType;
  typedef typename TLabelImage::IndexType               IndexType;

  typedef typename TLabelImage::PixelType               LabelType;

  /** Type to use form computations. */
  typedef typename NumericTraits<LabelType>::RealType RealType;

  /** \class LabelLabelOverlapMeasuress
   * \brief Metrics stored per label */
  class LabelSetMeasures
    {
    public:
    // default constructor
    LabelSetMeasures( void )
      {
      m_Source = 0;
      m_Target = 0;
      m_Union = 0;
      m_Intersection = 0;
      m_SourceComplement = 0;
      m_TargetComplement = 0;
      }

    // added for completeness
    LabelSetMeasures& operator=( const LabelSetMeasures& l )
      {
      m_Source = l.m_Source;
      m_Target = l.m_Target;
      m_Union = l.m_Union;
      m_Intersection = l.m_Intersection;
      m_SourceComplement = l.m_SourceComplement;
      m_TargetComplement = l.m_TargetComplement;
      return *this;
      }

    unsigned long m_Source;
    unsigned long m_Target;
    unsigned long m_Union;
    unsigned long m_Intersection;
    unsigned long m_SourceComplement;
    unsigned long m_TargetComplement;
    };

  /** Type of the map used to store data per label */
  typedef itksys::hash_map<LabelType, LabelSetMeasures> MapType;
  typedef typename MapType::iterator                    MapIterator;
  typedef typename MapType::const_iterator              MapConstIterator;

  /** Image related typedefs. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TLabelImage::ImageDimension );

  /** Set the source image. */
  void SetSourceImage( const LabelImageType * image )
    { this->SetNthInput( 0, const_cast<LabelImageType *>( image ) ); }

  /** Set the target image. */
  void SetTargetImage( const LabelImageType * image )
    { this->SetNthInput( 1, const_cast<LabelImageType *>( image ) ); }

  /** Get the source image. */
  const LabelImageType * GetSourceImage( void )
    { return this->GetInput( 0 ); }

  /** Get the target image. */
  const LabelImageType * GetTargetImage( void )
    { return this->GetInput( 1 ); }

  /** Get the label set measures */
  MapType GetLabelSetMeasures( void )
    { return this->m_LabelSetMeasures; }

  /**
   * tric overlap measures
   */
  /** measures over all labels */
  RealType GetTotalOverlap( void );
  RealType GetUnionOverlap( void );
  RealType GetMeanOverlap( void );
  RealType GetVolumeSimilarity( void );
  RealType GetFalseNegativeError( void );
  RealType GetFalsePositiveError( void );
  /** measures over individual labels */
  RealType GetTargetOverlap( LabelType );
  RealType GetUnionOverlap( LabelType );
  RealType GetMeanOverlap( LabelType );
  RealType GetVolumeSimilarity( LabelType );
  RealType GetFalseNegativeError( LabelType );
  RealType GetFalsePositiveError( LabelType );
  /** alternative names */
  RealType GetJaccardCoefficient( void )
    { return this->GetUnionOverlap(); }
  RealType GetJaccardCoefficient( LabelType label )
    { return this->GetUnionOverlap( label ); }
  RealType GetDiceCoefficient( void )
    { return this->GetMeanOverlap(); }
  RealType GetDiceCoefficient( LabelType label )
    { return this->GetMeanOverlap( label ); }


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( Input1HasNumericTraitsCheck,
    ( Concept::HasNumericTraits<LabelType> ) );
  /** End concept checking */
#endif

protected:
  LabelOverlapMeasuresImageFilter( void );
  ~LabelOverlapMeasuresImageFilter( void ) {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void BeforeThreadedGenerateData( void );

  void AfterThreadedGenerateData( void );

  /** Multi-thread version GenerateData. */
  void ThreadedGenerateData( const RegionType&, ThreadIdType );

  // Override since the filter needs all the data for the algorithm
  void GenerateInputRequestedRegion( void );

  // Override since the filter produces all of its output
  void EnlargeOutputRequestedRegion( DataObject *data );

private:
  //purposely not implemented
  LabelOverlapMeasuresImageFilter( const Self& );
  //purposely not implemented
  void operator=( const Self& );

  std::vector<MapType>                          m_LabelSetMeasuresPerThread;
  MapType                                       m_LabelSetMeasures;

  SimpleFastMutexLock                           m_Mutex;

}; // end of class

} // end namespace tube

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeLabelOverlapMeasuresImageFilter.txx"
#endif

#endif
