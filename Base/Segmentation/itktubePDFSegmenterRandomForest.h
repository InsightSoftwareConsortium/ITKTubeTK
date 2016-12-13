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

#ifndef __itktubePDFSegmenterRandomForest_h
#define __itktubePDFSegmenterRandomForest_h

#include "itktubePDFSegmenterBase.h"

#include "andres/marray.hxx"
#include "andres/ml/decision-trees.hxx"

#include <itkImage.h>
#include <itkListSample.h>

#include <vector>

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap >
class PDFSegmenterRandomForest
: public PDFSegmenterBase< TImage, TLabelMap >
{
public:

  typedef PDFSegmenterRandomForest                   Self;
  typedef PDFSegmenterBase< TImage, TLabelMap >      Superclass;
  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;

  itkTypeMacro( PDFSegmenterRandomForest, PDFSegmenterBase );

  itkNewMacro( Self );

  //
  // Template Args Typesdefs
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
  typedef andres::ml::DecisionForest< FeatureValueType, unsigned int,
    ProbabilityPixelType >                       DecisionForestType;

  //
  // Methods
  //
  itkGetMacro( TrainingDataStride, unsigned int );
  itkSetMacro( TrainingDataStride, unsigned int );

  DecisionForestType & GetModel( void );
  void SetModel( DecisionForestType & model );

  itkGetMacro( NumberOfDecisionTrees, unsigned int );
  itkSetMacro( NumberOfDecisionTrees, unsigned int );

  //
  // Must overwrite
  //
  virtual ProbabilityVectorType GetProbabilityVector( const
    FeatureVectorType & fv ) const;

protected:

  PDFSegmenterRandomForest( void );
  virtual ~PDFSegmenterRandomForest( void );

  //
  // Must overwrite
  //
  virtual void GeneratePDFs( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  PDFSegmenterRandomForest( const Self & ); // Purposely not implemented
  void operator = ( const Self & );         // Purposely not implemented

  // Superclass typedefs
  typedef std::vector< ProbabilityPixelType > ListVectorType;
  typedef std::vector< ListVectorType >       ListSampleType;
  typedef std::vector< ListSampleType >       ClassListSampleType;

  // Custom typedefs
  DecisionForestType            m_Model;

  unsigned int                  m_TrainingDataStride;

  unsigned int                  m_NumberOfDecisionTrees;

}; // End class PDFSegmenterRandomForest

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePDFSegmenterRandomForest.hxx"
#endif

#endif // End !defined( __itktubePDFSegmenterRandomForest_h )
