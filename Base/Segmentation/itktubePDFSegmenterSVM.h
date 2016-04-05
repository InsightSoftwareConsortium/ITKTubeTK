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

#ifndef __itktubePDFSegmenterSVM_h
#define __itktubePDFSegmenterSVM_h

#include "itktubePDFSegmenterBase.h"

#include "svm.h"

#include <itkImage.h>
#include <itkListSample.h>

#include <vector>

namespace itk
{

namespace tube
{

#define MAX_NUMBER_OF_FEATURES 5

template< class TImage, class TLabelMap >
class PDFSegmenterSVM : public PDFSegmenterBase< TImage, TLabelMap >
{
public:

  typedef PDFSegmenterSVM                            Self;
  typedef PDFSegmenterBase< TImage, TLabelMap >      Superclass;
  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;

  itkTypeMacro( PDFSegmenterSVM, PDFSegmenterBase );

  itkNewMacro( Self );

  //
  // Superclass Typedefs
  //
  typedef TImage                               InputImageType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    TImage::ImageDimension );

  typedef FeatureVectorGenerator< InputImageType >
                                               FeatureVectorGeneratorType;
  typedef typename FeatureVectorGeneratorType::FeatureValueType
                                               FeatureValueType;
  typedef typename FeatureVectorGeneratorType::FeatureVectorType
                                               FeatureVectorType;
  typedef typename FeatureVectorGeneratorType::FeatureImageType
                                               FeatureImageType;

  typedef TLabelMap                            LabelMapType;
  typedef typename LabelMapType::PixelType     LabelMapPixelType;

  typedef int                                  ObjectIdType;
  typedef std::vector< ObjectIdType >          ObjectIdListType;

  typedef float                                ProbabilityPixelType;
  typedef Image< ProbabilityPixelType, TImage::ImageDimension >
                                               ProbabilityImageType;

  typedef std::vector< double >                VectorDoubleType;
  typedef std::vector< int >                   VectorIntType;
  typedef std::vector< unsigned int >          VectorUIntType;

  //
  // Custom Typedefs
  //

  //
  // Methods
  //
  itkGetMacro( TrainingDataStride, unsigned int );
  itkSetMacro( TrainingDataStride, unsigned int );

  svm_model * GetModel( void );
  void SetModel( svm_model * model );

  svm_parameter * GetParameter( void );
  void SetParameter( svm_parameter * parameter );

  //
  // Must overwrite
  //
  virtual ProbabilityPixelType GetClassProbability( unsigned int
    classNum, const FeatureVectorType & fv ) const;

protected:

  PDFSegmenterSVM( void );
  virtual ~PDFSegmenterSVM( void );

  //
  // Must overwrite
  //
  virtual void GeneratePDFs( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  PDFSegmenterSVM( const Self & );       // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  // Superclass typedefs
  typedef std::vector< typename ProbabilityImageType::Pointer >
                                              ProbabilityImageVectorType;
  typedef std::vector< ProbabilityPixelType > ListVectorType;
  typedef std::vector< ListVectorType >       ListSampleType;
  typedef std::vector< ListSampleType >       ClassListSampleType;

  // Custom typedefs
  svm_model          * m_Model;
  svm_parameter        m_Parameter;

  unsigned int         m_TrainingDataStride;

}; // End class PDFSegmenterSVM

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePDFSegmenterSVM.hxx"
#endif

#endif // End !defined(__itktubePDFSegmenterSVM_h)
