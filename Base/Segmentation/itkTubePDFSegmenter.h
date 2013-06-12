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

#ifndef __itkTubePDFSegmenter_h
#define __itkTubePDFSegmenter_h

#include <itkImage.h>
#include <itkListSample.h>

#include <vector>

namespace itk
{

namespace tube
{

#define MAX_NUMBER_OF_FEATURES 4

template< class ImageT, unsigned int N, class LabelmapT >
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
  typedef ImageT                               ImageType;
  typedef typename ImageType::PixelType        PixelType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef LabelmapT                            MaskImageType;
  typedef typename MaskImageType::PixelType    MaskPixelType;

  typedef int                                  ObjectIdType;
  typedef std::vector< ObjectIdType >          ObjectIdListType;

  typedef float                                ProbabilityPixelType;
  typedef Image< ProbabilityPixelType,
    ImageT::ImageDimension >
                                               ProbabilityImageType;
  typedef std::vector< ProbabilityPixelType >  ProbabilityListType;

  typedef float                                HistogramPixelType;
  typedef Image< HistogramPixelType, N >       HistogramImageType;

  typedef HistogramPixelType                   PDFPixelType;
  typedef HistogramImageType                   PDFImageType;

  typedef Vector< double, N >                  VectorDoubleNType;

  //
  // Methods
  //
  void SetInputVolume( unsigned int featureNumber,
    typename ImageType::Pointer vol );

  void ClearObjectIds( void );
  void SetObjectId( ObjectIdType objectId );
  void AddObjectId( ObjectIdType objectId );

  ObjectIdType GetObjectId( unsigned int num = 0 ) const;
  unsigned int GetObjectNumberFromId( ObjectIdType id ) const;

  void   SetObjectPDFWeight( unsigned int num, double weight );
  double GetObjectPDFWeight( unsigned int num ) const;

  itkSetMacro( VoidId, ObjectIdType );
  itkGetMacro( VoidId, ObjectIdType );

  itkSetObjectMacro( Labelmap, MaskImageType );
  itkGetObjectMacro( Labelmap, MaskImageType );

  itkSetMacro( ErodeRadius, int );
  itkSetMacro( HoleFillIterations, int );
  itkSetMacro( ProbabilityImageSmoothingStandardDeviation, double );
  itkSetMacro( HistogramSmoothingStandardDeviation, double );
  itkSetMacro( OutlierRejectPortion, double );
  itkSetMacro( Draft, bool );

  int GetNumberOfClasses( void );

  const typename ProbabilityImageType::Pointer GetClassProbabilityVolume(
    unsigned int classNum ) const;

  const typename PDFImageType::Pointer GetClassPDFImage(
    unsigned int classNum );
  void SetClassPDFImage( unsigned int classNum,
    typename PDFImageType::Pointer classPDF );

  double GetPDFBinMin( unsigned int featureNum );
  void   SetPDFBinMin( unsigned int featureNum, double val );
  double GetPDFBinScale( unsigned int featureNum );
  void   SetPDFBinScale( unsigned int featureNum, double val );

  /** Copy the input object mask to the output mask, overwritting the
   *   classification assigned to those voxels. Default is false. */
  itkSetMacro( ReclassifyObjectMask, bool );
  itkGetMacro( ReclassifyObjectMask, bool );

  /** Copy the input not-object mask to the output mask, overwritting the
   *   classification assigned to those voxels. Default is false. */
  itkSetMacro( ReclassifyNotObjectMask, bool );
  itkGetMacro( ReclassifyNotObjectMask, bool );

  /** All object, void, and notObject pixels are force to being classified
   * as object or notObject. Default is false. */
  itkSetMacro( ForceClassification, bool );
  itkGetMacro( ForceClassification, bool );

  void SetProgressProcessInformation( void * processInfo, double fraction,
    double start );

  void Update( void );
  void ClassifyImages( void );
  void GenerateLabeledFeatureSpace( void );

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

  typedef itk::Vector< PixelType, N+ImageDimension >
    ListVectorType;
  typedef itk::Statistics::ListSample< ListVectorType >
    ListSampleType;
  typedef std::vector< typename ListSampleType::Pointer >
    ClassListSampleType;

  bool                                     m_SampleUpToDate;
  bool                                     m_PDFsUpToDate;
  bool                                     m_ImagesUpToDate;

  ClassListSampleType                      m_InClassList;
  typename ListSampleType::Pointer         m_OutClassList;

  ClassHistogramImageType                  m_InClassHistogram;
  typename HistogramImageType::Pointer     m_OutHistogram;
  VectorDoubleNType                        m_HistogramBinMin;
  VectorDoubleNType                        m_HistogramBinMax;
  VectorDoubleNType                        m_HistogramBinScale;
  unsigned int                             m_HistogramNumBinsND;
  unsigned int                             m_HistogramNumBins1D;

  //  Data
  std::vector< typename ImageType::Pointer > m_InputVolumeList;

  typename MaskImageType::Pointer m_Labelmap;

  ObjectIdListType                m_ObjectIdList;
  ObjectIdType                    m_VoidId;

  ProbabilityListType             m_PDFWeightList;

  int                             m_ErodeRadius;
  int                             m_HoleFillIterations;
  double                          m_ProbabilityImageSmoothingStandardDeviation;
  double                          m_HistogramSmoothingStandardDeviation;
  double                          m_OutlierRejectPortion;
  bool                            m_Draft;
  bool                            m_ReclassifyObjectMask;
  bool                            m_ReclassifyNotObjectMask;
  bool                            m_ForceClassification;

  ProbabilityImageVectorType      m_ProbabilityImageVector;

  void                          * m_ProgressProcessInfo;
  double                          m_ProgressFraction;
  double                          m_ProgressStart;

}; // End class PDFSegmenter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubePDFSegmenter.txx"
#endif

#endif // End !defined(__itkTubePDFSegmenter_h)
