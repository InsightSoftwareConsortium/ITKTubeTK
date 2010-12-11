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
#include "itkOrientedImage.h"

namespace itk
{

namespace tube
{

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
  typedef std::vector< int >                   ObjectIdListType;

  typedef float                                ProbabilityPixelType;
  typedef itk::OrientedImage< ProbabilityPixelType,
    ::itk::GetImageDimension< ImageT >::ImageDimension >
                                               ProbabilityImageType;


  //
  // Methods
  //
  itkSetObjectMacro( InputVolume1, ImageType );
  itkSetObjectMacro( InputVolume2, ImageType );
  itkSetObjectMacro( InputVolume3, ImageType );

  void SetObjectId( ObjectIdType objectId )
    {
    m_ObjectIdList.clear();
    m_ObjectIdList.push_back( objectId );
    };

  void AddObjectId( ObjectIdType objectId )
    {
    m_ObjectIdList.push_back( objectId );
    };

  ObjectIdType GetObjectId( int num = 0 )
    {
    return m_ObjectIdList[ num ];
    };

  itkSetMacro( VoidId, ObjectIdType );
  itkGetMacro( VoidId, ObjectIdType );

  itkSetObjectMacro( Labelmap, MaskImageType );
  itkGetObjectMacro( Labelmap, MaskImageType );

  itkSetMacro( UseTexture, bool );
  itkSetMacro( ErodeRadius, int );
  itkSetMacro( HoleFillIterations, int );
  itkSetMacro( FprWeight, double );
  itkSetMacro( ProbabilitySmoothingStandardDeviation, double );
  itkSetMacro( Draft, bool );

  const typename ProbabilityImageType::Pointer * GetProbabilityImage(
    unsigned int classNum );

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

protected:

  typedef std::vector< typename ProbabilityImageType::Pointer >
    ProbabilityImageVectorType;

  PDFSegmenter( void );
  virtual ~PDFSegmenter( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

  typename ImageType::Pointer GenerateTextureImage(
    const ImageType * im );

private:

  PDFSegmenter( const Self & );          // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  //  Data
  typename ImageType::Pointer     m_InputVolume1;
  typename ImageType::Pointer     m_InputVolume2;
  typename ImageType::Pointer     m_InputVolume3;

  typename MaskImageType::Pointer m_Labelmap;

  ObjectIdListType                m_ObjectIdList;
  ObjectIdType                    m_VoidId;

  bool                            m_UseTexture;
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

};

}

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubePDFSegmenter.txx"
#endif

#endif
