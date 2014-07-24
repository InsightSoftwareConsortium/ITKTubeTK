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

#ifndef __itktubePDFSegmenter_h
#define __itktubePDFSegmenter_h

#include <itkImage.h>
#include <itkListSample.h>

#include <vector>

namespace itk
{

namespace tube
{

#define MAX_NUMBER_OF_FEATURES 5

template< class TImage, unsigned int N, class TLabelMap >
class PDFSegmenter : public Object
{
public:

  typedef PDFSegmenter                         Self;
  typedef Object                               Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( PDFSegmenter, Object );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  typedef TImage                               ImageType;
  typedef typename ImageType::PixelType        PixelType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    TImage::ImageDimension );

  typedef TLabelMap                            LabelMapType;
  typedef typename LabelMapType::PixelType     LabelMapPixelType;

  typedef int                                  ObjectIdType;
  typedef std::vector< ObjectIdType >          ObjectIdListType;

  typedef float                                ProbabilityPixelType;
  typedef Image< ProbabilityPixelType, TImage::ImageDimension >
                                               ProbabilityImageType;

  typedef float                                HistogramPixelType;
  typedef Image< HistogramPixelType, N >       HistogramImageType;

  typedef HistogramPixelType                   PDFPixelType;
  typedef HistogramImageType                   PDFImageType;

  typedef Image< LabelMapPixelType, N >        LabeledFeatureSpaceType;

  typedef std::vector< double >                VectorDoubleType;
  typedef std::vector< int >                   VectorIntType;
  typedef std::vector< unsigned int >          VectorUIntType;

  //
  // Methods
  //
  void SetInput( unsigned int featureNumber,
    typename ImageType::Pointer vol );

  void ClearObjectIds( void );
  void SetObjectId( ObjectIdType objectId );
  void AddObjectId( ObjectIdType objectId );
  void SetObjectId( const ObjectIdListType _objectId );
  const VectorIntType & GetObjectId( void ) const;

  unsigned int GetNumberOfClasses( void ) const;
  unsigned int GetNumberOfObjectIds( void ) const;

  unsigned int GetNumberOfFeatures( void ) const
    { return N; };

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
  itkSetMacro( HistogramSmoothingStandardDeviation, double );
  itkGetMacro( HistogramSmoothingStandardDeviation, double );
  itkSetMacro( OutlierRejectPortion, double );
  itkGetMacro( OutlierRejectPortion, double );
  itkSetMacro( Draft, bool );
  itkGetMacro( Draft, bool );

  typename ProbabilityImageType::Pointer
    GetClassProbabilityForInput( unsigned int classNum ) const;

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

  void SetProgressProcessInformation( void * processInfo, double fraction,
    double start );

  /** Given one PDF per class, generate a labelmap of feature space */
  void GenerateLabeledFeatureSpace( void );

  void SetLabeledFeatureSpace( typename LabeledFeatureSpaceType::Pointer
    labeledFeatureSpace );

  typename LabeledFeatureSpaceType::Pointer GetLabeledFeatureSpace( void )
    const;

  void Update( void );
  void ClassifyImages( void );

protected:

  PDFSegmenter( void );
  virtual ~PDFSegmenter( void );

  void GenerateSample( void );
  void GeneratePDFs( void );
  void ApplyPDFs( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  PDFSegmenter( const Self & );          // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  typedef std::vector< typename ProbabilityImageType::Pointer >
    ProbabilityImageVectorType;

  typedef std::vector< typename HistogramImageType::Pointer >
    ClassHistogramImageType;

  typedef itk::Vector< ProbabilityPixelType, N+ImageDimension >      ListVectorType;
  typedef itk::Statistics::ListSample< ListVectorType >   ListSampleType;
  typedef std::vector< typename ListSampleType::Pointer > ClassListSampleType;

  bool                                     m_SampleUpToDate;
  bool                                     m_PDFsUpToDate;
  bool                                     m_ImagesUpToDate;

  ClassListSampleType                      m_InClassList;
  typename ListSampleType::Pointer         m_OutClassList;

  ClassHistogramImageType                  m_InClassHistogram;
  VectorDoubleType                         m_HistogramBinMin;
  VectorDoubleType                         m_HistogramBinSize;
  VectorUIntType                           m_HistogramNumberOfBin;

  //  Data
  std::vector< typename ImageType::Pointer > m_InputImageList;

  typename LabelMapType::Pointer  m_LabelMap;

  ObjectIdListType                m_ObjectIdList;
  ObjectIdType                    m_VoidId;

  VectorDoubleType                m_PDFWeightList;

  int                             m_ErodeRadius;
  bool                            m_DilateFirst;
  int                             m_HoleFillIterations;
  double                          m_ProbabilityImageSmoothingStandardDeviation;
  double                          m_HistogramSmoothingStandardDeviation;
  double                          m_OutlierRejectPortion;
  bool                            m_Draft;
  bool                            m_ReclassifyObjectLabels;
  bool                            m_ReclassifyNotObjectLabels;
  bool                            m_ForceClassification;

  ProbabilityImageVectorType      m_ProbabilityImageVector;

  typename LabeledFeatureSpaceType::Pointer    m_LabeledFeatureSpace;

  void                          * m_ProgressProcessInfo;
  double                          m_ProgressFraction;
  double                          m_ProgressStart;

}; // End class PDFSegmenter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePDFSegmenter.hxx"
#endif

#endif // End !defined(__itktubePDFSegmenter_h)
