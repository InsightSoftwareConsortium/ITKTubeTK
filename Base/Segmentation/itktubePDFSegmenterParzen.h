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

#ifndef __itktubePDFSegmenterParzen_h
#define __itktubePDFSegmenterParzen_h

#include "itktubePDFSegmenterBase.h"

#include <itkImage.h>
#include <itkListSample.h>

#include <vector>

namespace itk
{

namespace tube
{

#define PARZEN_MAX_NUMBER_OF_FEATURES 4

template< class TImage, class TLabelMap >
class PDFSegmenterParzen : public PDFSegmenterBase< TImage, TLabelMap >
{
public:

  typedef PDFSegmenterParzen                         Self;
  typedef PDFSegmenterBase< TImage, TLabelMap >      Superclass;
  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;

  itkTypeMacro( PDFSegmenterParzen, PDFSegmenterBase );

  itkNewMacro( Self );

  //
  // Template Args Typedefs
  //
  typedef TImage                               InputImageType;
  typedef TLabelMap                            LabelMapType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    TImage::ImageDimension );

  //
  // Superclass Typedefs
  //
  typedef typename Superclass::FeatureVectorGeneratorType
                                                 FeatureVectorGeneratorType;
  typedef typename Superclass::FeatureValueType  FeatureValueType;
  typedef typename Superclass::FeatureVectorType FeatureVectorType;
  typedef typename Superclass::FeatureImageType  FeatureImageType;

  typedef typename Superclass::LabelMapPixelType LabelMapPixelType;

  typedef typename Superclass::ObjectIdType      ObjectIdType;
  typedef typename Superclass::ObjectIdListType  ObjectIdListType;

  typedef typename Superclass::ProbabilityPixelType
                                                 ProbabilityPixelType;

  typedef typename Superclass::ProbabilityVectorType
                                                 ProbabilityVectorType;
  typedef typename Superclass::ProbabilityImageType
                                                 ProbabilityImageType;

  typedef typename Superclass::VectorDoubleType  VectorDoubleType;
  typedef typename Superclass::VectorIntType     VectorIntType;
  typedef typename Superclass::VectorUIntType    VectorUIntType;

  //
  // Custom Typedefs
  //
  typedef float                                HistogramPixelType;

  typedef Image< HistogramPixelType, PARZEN_MAX_NUMBER_OF_FEATURES >
    HistogramImageType;

  typedef HistogramPixelType                   PDFPixelType;
  typedef HistogramImageType                   PDFImageType;

  typedef Image< LabelMapPixelType, PARZEN_MAX_NUMBER_OF_FEATURES >
    LabeledFeatureSpaceType;

  //
  // Methods
  //
  itkSetMacro( HistogramSmoothingStandardDeviation, double );
  itkGetMacro( HistogramSmoothingStandardDeviation, double );
  itkSetMacro( OutlierRejectPortion, double );
  itkGetMacro( OutlierRejectPortion, double );

  typename PDFImageType::Pointer GetClassPDFImage(
    unsigned int classNum ) const;

  void SetClassPDFImage( unsigned int classNum,
    typename PDFImageType::Pointer classPDF );

  const VectorUIntType & GetNumberOfBinsPerFeature( void ) const;
  void             SetNumberOfBinsPerFeature( const VectorUIntType & nBin );
  const VectorDoubleType & GetBinMin( void ) const;
  void             SetBinMin( const VectorDoubleType & binMin );
  const VectorDoubleType & GetBinSize( void ) const;
  void             SetBinSize( const VectorDoubleType & binMin );

  /** Given one PDF per class, generate a labelmap of feature space */
  void GenerateLabeledFeatureSpace( void );

  void SetLabeledFeatureSpace( typename LabeledFeatureSpaceType::Pointer
    labeledFeatureSpace );

  typename LabeledFeatureSpaceType::Pointer GetLabeledFeatureSpace( void )
    const;

  virtual void Update( void );

  //
  // Must overwrite
  //
  virtual ProbabilityVectorType GetProbabilityVector( const
    FeatureVectorType & fv ) const;

protected:

  PDFSegmenterParzen( void );
  virtual ~PDFSegmenterParzen( void );

  virtual void GeneratePDFs( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  PDFSegmenterParzen( const Self & );          // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  // Superclass typedefs
  typedef std::vector< typename ProbabilityImageType::Pointer >
    ProbabilityImageVectorType;

  typedef std::vector< ProbabilityPixelType > ListVectorType;
  typedef std::vector< ListVectorType >       ListSampleType;
  typedef std::vector< ListSampleType >       ClassListSampleType;

  // Custom typedefs
  typedef std::vector< typename HistogramImageType::Pointer >
    ClassHistogramImageType;

  ClassHistogramImageType         m_InClassHistogram;
  VectorDoubleType                m_HistogramBinMin;
  VectorDoubleType                m_HistogramBinSize;
  VectorUIntType                  m_HistogramNumberOfBin;

  double                          m_OutlierRejectPortion;

  double                          m_HistogramSmoothingStandardDeviation;

  typename LabeledFeatureSpaceType::Pointer m_LabeledFeatureSpace;

}; // End class PDFSegmenterParzen

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePDFSegmenterParzen.hxx"
#endif

#endif // End !defined( __itktubePDFSegmenterParzen_h )
