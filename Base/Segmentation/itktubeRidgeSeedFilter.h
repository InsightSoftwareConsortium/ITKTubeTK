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

#ifndef __itktubeRidgeSeedFilter_h
#define __itktubeRidgeSeedFilter_h

#include "itktubeBasisFeatureVectorGenerator.h"
#include "itktubePDFSegmenter.h"
#include "itktubeRidgeFFTFeatureVectorGenerator.h"

#include <itkImage.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <vector>

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap >
class RidgeSeedFilter : public Object
{
public:

  typedef RidgeSeedFilter                            Self;
  typedef Object                                     Superclass;
  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;

  itkTypeMacro( RidgeSeedFilter, ImageToImageFilter );

  itkNewMacro( Self );

  typedef TImage                                  ImageType;
  typedef TImage                                  InputImageType;
  typedef Image< float, TImage::ImageDimension >  OutputImageType;

  typedef TLabelMap                               LabelMapType;
  typedef typename LabelMapType::PixelType        LabelMapPixelType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    TImage::ImageDimension );

  typedef RidgeFFTFeatureVectorGenerator< TImage >
    RidgeFeatureGeneratorType;

  typedef typename RidgeFeatureGeneratorType::FeatureValueType
    FeatureValueType;
  typedef typename RidgeFeatureGeneratorType::FeatureVectorType
    FeatureVectorType;

  typedef typename RidgeFeatureGeneratorType::IndexType
    IndexType;
  typedef typename RidgeFeatureGeneratorType::RidgeScalesType
    RidgeScalesType;
  typedef typename RidgeFeatureGeneratorType::ValueListType
    WhitenMeansType;
  typedef typename RidgeFeatureGeneratorType::ValueListType
    WhitenStdDevsType;

  typedef BasisFeatureVectorGenerator< TImage, LabelMapType >
    SeedFeatureGeneratorType;

  typedef typename SeedFeatureGeneratorType::ObjectIdType ObjectIdType;
  typedef typename SeedFeatureGeneratorType::VectorType   VectorType;
  typedef typename SeedFeatureGeneratorType::MatrixType   MatrixType;
  typedef typename SeedFeatureGeneratorType::FeatureImageType
    BasisImageType;

  typedef PDFSegmenter< OutputImageType, 4, LabelMapType >
    PDFSegmenterType;
  typedef typename  PDFSegmenterType::ProbabilityPixelType
    ProbabilityPixelType;
  typedef typename  PDFSegmenterType::ProbabilityImageType
    ProbabilityImageType;

  void SetInput( typename ImageType::Pointer img );
  void AddInput( typename ImageType::Pointer img );

  void SetLabelMap( typename LabelMapType::Pointer img );

  typename SeedFeatureGeneratorType::Pointer
    GetSeedFeatureGenerator( void );
  typename RidgeFeatureGeneratorType::Pointer
    GetRidgeFeatureGenerator( void );

  typename PDFSegmenterType::Pointer GetPDFSegmenter( void );

  void            SetScales( const RidgeScalesType & Scales );
  RidgeScalesType GetScales( void ) const;

  // Basis
  void   SetWhitenMeans( const WhitenMeansType & means );
  void   SetWhitenStdDevs( const WhitenStdDevsType & stdDevs );
  const WhitenMeansType &   GetWhitenMeans( void ) const;
  const WhitenStdDevsType & GetWhitenStdDevs( void ) const;

  unsigned int GetNumberOfBasis( void ) const;

  double       GetBasisValue( unsigned int basisNum ) const;
  VectorType   GetBasisVector( unsigned int basisNum ) const;
  MatrixType   GetBasisMatrix( void ) const;
  VectorType   GetBasisValues( void ) const;

  typename BasisImageType::Pointer GetBasisImage( unsigned int num = 0 )
    const;

  void   SetBasisValue( unsigned int basisNum, double value );
  void   SetBasisVector( unsigned int basisNum, const VectorType & vec );
  void   SetBasisMatrix( const MatrixType & mat );
  void   SetBasisValues( const VectorType & values );

  // PDFSegmenter
  typename ProbabilityImageType::Pointer
    GetClassProbabilityForInput( unsigned int objectNum ) const;

  typename ProbabilityImageType::Pointer
    GetClassProbabilityDifferenceForInput( unsigned int objectNum ) const;

  // Ridge, Basis, and PDFSegmenter
  itkSetMacro( RidgeId, ObjectIdType );
  itkGetMacro( RidgeId, ObjectIdType );
  itkSetMacro( BackgroundId, ObjectIdType );
  itkGetMacro( BackgroundId, ObjectIdType );
  itkSetMacro( UnknownId, ObjectIdType );
  itkGetMacro( UnknownId, ObjectIdType );

  itkSetMacro( SeedTolerance, double );
  itkGetMacro( SeedTolerance, double );

  itkSetMacro( Skeletonize, bool );
  itkGetMacro( Skeletonize, bool );

  itkSetMacro( UseIntensityOnly, bool );
  itkGetMacro( UseIntensityOnly, bool );

  itkSetMacro( TrainClassifier, bool );
  itkGetMacro( TrainClassifier, bool );

  // Local
  void   Update();
  void   ClassifyImages();

  typename LabelMapType::Pointer GetOutput( void );

  typename OutputImageType::Pointer GetOutputSeedScales( void );

protected:

  RidgeSeedFilter( void );
  virtual ~RidgeSeedFilter( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  RidgeSeedFilter( const Self & );    // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  typename RidgeFeatureGeneratorType::Pointer     m_RidgeFeatureGenerator;
  typename SeedFeatureGeneratorType::Pointer      m_SeedFeatureGenerator;
  typename PDFSegmenterType::Pointer              m_PDFSegmenter;

  ObjectIdType   m_RidgeId;
  ObjectIdType   m_BackgroundId;
  ObjectIdType   m_UnknownId;

  double         m_SeedTolerance;

  bool           m_Skeletonize;

  bool           m_UseIntensityOnly;

  bool           m_TrainClassifier;

  typename LabelMapType::Pointer m_LabelMap;

}; // End class RidgeSeedFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRidgeSeedFilter.hxx"
#endif

#endif // End !defined(__itktubeRidgeSeedFilter_h)
