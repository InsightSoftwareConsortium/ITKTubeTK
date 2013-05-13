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

#include <vector>

#include "itkImage.h"
#include "itkListSample.h"

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

  typedef Image<float, N>                      HistogramImageType;
  typedef HistogramImageType                   PDFImageType;
  typedef Vector<double, N>                    ListDoubleType;

  //
  // Methods
  //
  void SetInputVolume( unsigned int featureNumber,
    typename ImageType::Pointer vol );

  void ClearObjectIds( void )
    {
    m_ObjectIdList.clear();
    }

  void AddObjectId( ObjectIdType objectId )
    {
    m_ObjectIdList.push_back( objectId );
    }

  ObjectIdType GetObjectId( int num = 0 )
    {
    return m_ObjectIdList[ num ];
    }

  itkSetMacro( VoidId, ObjectIdType );
  itkGetMacro( VoidId, ObjectIdType );

  itkSetObjectMacro( Labelmap, MaskImageType );
  itkGetObjectMacro( Labelmap, MaskImageType );

  itkSetMacro( ErodeRadius, int );
  itkSetMacro( HoleFillIterations, int );
  itkSetMacro( FprWeight, double );
  itkSetMacro( ProbabilitySmoothingStandardDeviation, double );
  itkSetMacro( Draft, bool );

  int GetNumberOfClasses( void );

  int GetClassId( unsigned int classNum );

  const typename ProbabilityImageType::Pointer GetClassProbabilityVolume(
    unsigned int classNum );

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
  typename ListSampleType::Pointer         m_OutList;

  ClassHistogramImageType                  m_InClassHisto;
  typename HistogramImageType::Pointer     m_OutHisto;
  ListDoubleType                           m_HistoBinMin;
  ListDoubleType                           m_HistoBinMax;
  ListDoubleType                           m_HistoBinScale;
  unsigned int                             m_HistoNumBinsND;
  unsigned int                             m_HistoNumBins1D;

  //  Data
  std::vector< typename ImageType::Pointer > m_InputVolumeList;

  typename MaskImageType::Pointer m_Labelmap;

  ObjectIdListType                m_ObjectIdList;
  ObjectIdType                    m_VoidId;

  int                             m_ErodeRadius;
  int                             m_HoleFillIterations;
  double                          m_ProbabilitySmoothingStandardDeviation;
  double                          m_FprWeight;
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
