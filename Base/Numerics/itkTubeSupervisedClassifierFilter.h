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
#ifndef __itkTubeSupervisedClassifierFilter_h
#define __itkTubeSupervisedClassifierFilter_h

#include <vector>

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "itkImage.h"

#include "itkTubeImageToImageFilter.h"

#include "itkTubeFeatureVectorGenerator.h"

namespace itk
{

namespace tube
{

template< class ImageT, class LabelmapT >
class SupervisedClassifierFilter : public ImageToImageFilter< ImageT,
  Image< ImageT::ImageDimension, LabelmapT::PixelType > >
{
public:

  typedef SupervisedClassifierFilter           Self;
  typedef ImageToImageFilter<  ImageT, Image< ImageT::ImageDimension,
   unsigned short > >                          Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkNewMacro( Self );

  itkTypeMacro( SupervisedClassifierFilter, ImageToImageFilter );

  //
  // Custom Typedefs
  //
  typedef ImageT                                           ImageType;

  typedef Image< ImageT::ImageDimension, unsigned short >  OutputImageType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef LabelmapT                                  MaskImageType;

  typedef FeatureVectorGenerator::FeatureVectorType  FeatureVectorType;

  typedef LabelmapT::PixelType                       ObjectIdType;

  typedef std::vector< ObjectIdType >                ObjectIdListType;

  //
  // Methods
  //
  void SetFeatureVectorGenerator( FeatureVectorGenerator::Pointer fGen );
  FeatureVectorGenerator::Pointer GetFeatureVectorGenerator( void );

  void             SetObjectId( ObjectIdType objectId );
  void             AddObjectId( ObjectIdType objectId );
  ObjectIdType     GetObjectId( unsigned int num = 0 );
  unsigned int     GetNumberOfObjectIds( void );

  itkSetObjectMacro( Labelmap, MaskImageType );
  itkGetObjectMacro( Labelmap, MaskImageType );

  void SetProgressProcessInformation( void * processInfo, double fraction,
    double start );

protected:

  SupervisedClassifierGenerator( void );
  virtual ~SupervisedClassifierGenerator( void );

  typename MaskImageType::Pointer m_Labelmap;

  ObjectIdListType                m_ObjectIdList;

  FeatureVectorGenerator::Pointer m_FeatureVectorGenerator;

  // ITK Filter Stuff
  void PrintSelf( std::ostream & os, Indent indent ) const;

  void ThreadedGenerateData( const OutputImageRegionType & threadRegion,
    ThreadIdType threadId );

  void BeforeThreadedGenerateData( void );
  void AfterThreadedGenerateData( void );

private:

  // Purposely not implemented
  SupervisedClassifierGenerator( const Self & );
  void operator = ( const Self & );      // Purposely not implemented

  //  Data
  void                          * m_ProgressProcessInfo;
  double                          m_ProgressFraction;
  double                          m_ProgressStart;

};

}

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeSupervisedClassifierGenerator.txx"
#endif

#endif
