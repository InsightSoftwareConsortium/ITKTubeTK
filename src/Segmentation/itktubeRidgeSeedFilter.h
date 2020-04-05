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

#ifndef __itktubeRidgeSeedFilter_h
#define __itktubeRidgeSeedFilter_h

#include "itktubeBasisFeatureVectorGenerator.h"
#include "itktubePDFSegmenterBase.h"
#include "itktubePDFSegmenterParzen.h"
#include "itktubeRidgeFFTFeatureVectorGenerator.h"

#include <itkImage.h>
#include <itkImageToImageFilter.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <vector>

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap >
class RidgeSeedFilter : public ImageToImageFilter< TImage, TLabelMap >
{
public:

  typedef RidgeSeedFilter                            Self;
  typedef ImageToImageFilter< TImage, TLabelMap >    Superclass;
  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;

  itkTypeMacro( RidgeSeedFilter, ImageToImageFilter );

  itkNewMacro( Self );

  typedef TImage                                  InputImageType;
  typedef Image< float, TImage::ImageDimension >  OutputImageType;

  typedef TLabelMap                               LabelMapType;
  typedef typename LabelMapType::PixelType        LabelMapPixelType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    TImage::ImageDimension );

  typedef RidgeFFTFeatureVectorGenerator< InputImageType >
    RidgeFeatureGeneratorType;

  typedef typename RidgeFeatureGeneratorType::FeatureValueType
    FeatureValueType;
  typedef typename RidgeFeatureGeneratorType::FeatureVectorType
    FeatureVectorType;
  typedef typename RidgeFeatureGeneratorType::FeatureImageType
    FeatureImageType;

  typedef typename RidgeFeatureGeneratorType::IndexType
    IndexType;
  typedef typename RidgeFeatureGeneratorType::RidgeScalesType
    RidgeScalesType;
  typedef typename RidgeFeatureGeneratorType::ValueListType
    WhitenMeansType;
  typedef typename RidgeFeatureGeneratorType::ValueListType
    WhitenStdDevsType;

  typedef BasisFeatureVectorGenerator< InputImageType, LabelMapType >
    SeedFeatureGeneratorType;

  typedef typename SeedFeatureGeneratorType::ObjectIdType ObjectIdType;
  typedef typename SeedFeatureGeneratorType::VectorType   VectorType;
  typedef typename SeedFeatureGeneratorType::MatrixType   MatrixType;

  typedef PDFSegmenterBase< InputImageType, LabelMapType >
    PDFSegmenterType;
  typedef PDFSegmenterParzen< InputImageType, LabelMapType >
    PDFSegmenterParzenType;
  typedef typename  PDFSegmenterType::ProbabilityPixelType
    ProbabilityPixelType;
  typedef typename  PDFSegmenterType::ProbabilityImageType
    ProbabilityImageType;

  virtual void SetInput( const InputImageType * img ) override;
  virtual void SetInput( unsigned int id, const InputImageType * img ) override;

  using Superclass::AddInput;
  virtual void AddInput( const InputImageType * img );

  void SetLabelMap( LabelMapType * img );
  typename LabelMapType::Pointer GetLabelMap( void )
  { return this->GetOutput(); };

  typename SeedFeatureGeneratorType::Pointer
    GetSeedFeatureGenerator( void );
  typename RidgeFeatureGeneratorType::Pointer
    GetRidgeFeatureGenerator( void );

  typename PDFSegmenterType::Pointer GetPDFSegmenter( void );
  void SetPDFSegmenter( PDFSegmenterType * pdfSegmenter );

  void            SetScales( const RidgeScalesType & Scales );
  RidgeScalesType GetScales( void ) const;

  // Basis
  void         SetInputWhitenMeans( const WhitenMeansType & means );
  const WhitenMeansType &   GetInputWhitenMeans( void ) const;

  void         SetInputWhitenStdDevs( const WhitenStdDevsType & stdDevs );
  const WhitenStdDevsType & GetInputWhitenStdDevs( void ) const;

  void         SetOutputWhitenMeans( const WhitenMeansType & means );
  const WhitenMeansType &   GetOutputWhitenMeans( void ) const;

  void         SetOutputWhitenStdDevs( const WhitenStdDevsType & stdDevs );
  const WhitenStdDevsType & GetOutputWhitenStdDevs( void ) const;

  unsigned int GetNumberOfBasis( void ) const;

  void         SetBasisValue( unsigned int basisNum, double value );
  double       GetBasisValue( unsigned int basisNum ) const;
  void         SetBasisVector( unsigned int basisNum, const VectorType & vec );
  VectorType   GetBasisVector( unsigned int basisNum ) const;
  void         SetBasisMatrix( const MatrixType & mat );
  MatrixType   GetBasisMatrix( void ) const;
  void         SetBasisValues( const VectorType & values );
  VectorType   GetBasisValues( void ) const;

  typename FeatureImageType::Pointer GetBasisImage( unsigned int num = 0 )
    const;


  // PDFSegmenter
  typename ProbabilityImageType::Pointer
    GetClassProbabilityImage( unsigned int objectNum ) const;

  typename ProbabilityImageType::Pointer
    GetClassLikelihoodRatioImage( unsigned int objectNum );

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
  void   Update() override;
  void   ClassifyImages();

  typename LabelMapType::Pointer GetOutput( void );

  typename OutputImageType::Pointer GetOutputSeedScales( void );

protected:

  RidgeSeedFilter( void );
  virtual ~RidgeSeedFilter( void );

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  RidgeSeedFilter( const Self & );    // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const typename Superclass::DataObjectIdentifierType &,
    itk::DataObject * ) override
    {};

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

  mutable typename LabelMapType::Pointer m_LabelMap;

  std::vector< typename ProbabilityImageType::Pointer > m_RatioImageVector;

}; // End class RidgeSeedFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRidgeSeedFilter.hxx"
#endif

#endif // End !defined( __itktubeRidgeSeedFilter_h )
