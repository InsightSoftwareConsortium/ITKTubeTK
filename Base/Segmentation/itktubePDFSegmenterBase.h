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

#ifndef __itktubePDFSegmenterBase_h
#define __itktubePDFSegmenterBase_h

#include "itktubeFeatureVectorGenerator.h"

#include <itkImage.h>
#include <itkListSample.h>

#include <vector>

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap >
class PDFSegmenterBase : public Object
{
public:

  typedef PDFSegmenterBase                     Self;
  typedef Object                               Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( PDFSegmenterBase, Object );

  itkNewMacro( Self );

  //
  // Template Args Typedefs
  //
  typedef TImage                               InputImageType;

  typedef TLabelMap                            LabelMapType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    TImage::ImageDimension );

  //
  // Base Typedefs
  //
  typedef FeatureVectorGenerator< InputImageType >
                                               FeatureVectorGeneratorType;
  typedef typename FeatureVectorGeneratorType::FeatureValueType
                                               FeatureValueType;
  typedef typename FeatureVectorGeneratorType::FeatureVectorType
                                               FeatureVectorType;
  typedef typename FeatureVectorGeneratorType::FeatureImageType
                                               FeatureImageType;

  typedef typename LabelMapType::PixelType     LabelMapPixelType;

  typedef int                                  ObjectIdType;
  typedef std::vector< ObjectIdType >          ObjectIdListType;

  typedef float                                ProbabilityPixelType;
  typedef std::vector< ProbabilityPixelType >  ProbabilityVectorType;
  typedef Image< ProbabilityPixelType, TImage::ImageDimension >
                                               ProbabilityImageType;

  typedef std::vector< double >                VectorDoubleType;
  typedef std::vector< int >                   VectorIntType;
  typedef std::vector< unsigned int >          VectorUIntType;

  //
  // Methods
  //
  void SetFeatureVectorGenerator( typename
    FeatureVectorGeneratorType::Pointer fvg );

  void ClearObjectIds( void );
  void SetObjectId( ObjectIdType objectId );
  void AddObjectId( ObjectIdType objectId );
  void SetObjectId( const ObjectIdListType _objectId );
  const VectorIntType & GetObjectId( void ) const;

  unsigned int GetNumberOfClasses( void ) const;
  unsigned int GetNumberOfObjectIds( void ) const;

  unsigned int GetNumberOfFeatures( void ) const;

  unsigned int GetObjectNumberFromId( ObjectIdType id ) const;

  void   SetObjectPDFWeight( unsigned int num, double weight );
  void   SetObjectPDFWeight( const VectorDoubleType & weight );
  const VectorDoubleType & GetObjectPDFWeight( void ) const;

  itkSetMacro( VoidId, ObjectIdType );
  itkGetMacro( VoidId, ObjectIdType );

  itkSetObjectMacro( LabelMap, LabelMapType );
  itkGetObjectMacro( LabelMap, LabelMapType );

  itkSetMacro( ErodeRadius, int );
  itkGetMacro( ErodeRadius, int );
  itkSetMacro( DilateFirst, bool );
  itkGetMacro( DilateFirst, bool );
  itkSetMacro( HoleFillIterations, int );
  itkGetMacro( HoleFillIterations, int );
  itkSetMacro( ProbabilityImageSmoothingStandardDeviation, double );
  itkGetMacro( ProbabilityImageSmoothingStandardDeviation, double );

  /** Copy the input object mask to the output mask, overwritting the
   *   classification assigned to those voxels. Default is false. */
  itkSetMacro( ReclassifyObjectLabels, bool );
  itkGetMacro( ReclassifyObjectLabels, bool );

  /** Copy the input not-object mask to the output mask, overwritting the
   *   classification assigned to those voxels. Default is false. */
  itkSetMacro( ReclassifyNotObjectLabels, bool );
  itkGetMacro( ReclassifyNotObjectLabels, bool );

  /** All object, void, and notObject pixels are force to being classified
   * as object or notObject. Default is false. */
  itkSetMacro( ForceClassification, bool );
  itkGetMacro( ForceClassification, bool );

  /** Reduce sample size per class to match that of the class with the
   * sample size. Default is true. */
  itkSetMacro( BalanceClassSampleSize, bool );
  itkGetMacro( BalanceClassSampleSize, bool );

  void SetProgressProcessInformation( void * processInfo, double fraction,
    double start );

  virtual void Update( void );

  virtual void ClassifyImages( void );

  // Overwrite for speedup
  virtual typename ProbabilityImageType::Pointer GetClassProbabilityImage(
    unsigned int classNum ) const;

  //
  // Must overwrite
  //
  virtual ProbabilityVectorType GetProbabilityVector( const
    FeatureVectorType & fv ) const;



protected:

  PDFSegmenterBase( void );
  virtual ~PDFSegmenterBase( void );

  void BalanceClassSampleSize( void );

  virtual void GenerateSample( void );

  //
  // Must overwrite
  //
  virtual void GeneratePDFs( void );

  virtual void ApplyPDFs( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

  typedef std::vector< typename ProbabilityImageType::Pointer >
                                              ProbabilityImageVectorType;
  typedef std::vector< ProbabilityPixelType > ListVectorType;
  typedef std::vector< ListVectorType >       ListSampleType;
  typedef std::vector< ListSampleType >       ClassListSampleType;

  ClassListSampleType                           m_InClassList;
  ListSampleType                                m_OutClassList;

  typename FeatureVectorGeneratorType::Pointer  m_FeatureVectorGenerator;

  typename LabelMapType::Pointer                m_LabelMap;

  ProbabilityImageVectorType                    m_ProbabilityImageVector;

  bool                m_SampleUpToDate;
  bool                m_PDFsUpToDate;
  bool                m_ClassProbabilityImagesUpToDate;

  ObjectIdListType    m_ObjectIdList;
  ObjectIdType        m_VoidId;

  void              * m_ProgressProcessInfo;
  double              m_ProgressFraction;
  double              m_ProgressStart;

private:

  PDFSegmenterBase( const Self & );      // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  VectorDoubleType    m_PDFWeightList;

  int                 m_ErodeRadius;
  bool                m_DilateFirst;
  int                 m_HoleFillIterations;
  double              m_ProbabilityImageSmoothingStandardDeviation;
  bool                m_ReclassifyObjectLabels;
  bool                m_ReclassifyNotObjectLabels;
  bool                m_ForceClassification;

  bool                m_BalanceClassSampleSize;

}; // End class PDFSegmenterBase

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePDFSegmenterBase.hxx"
#endif

#endif // End !defined( __itktubePDFSegmenterBase_h )
