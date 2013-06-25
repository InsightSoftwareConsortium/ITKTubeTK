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
#include "itktubeRidgeFeatureVectorGenerator.h"

#include <itkImage.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <vector>

namespace itk
{

namespace tube
{

template< class ImageT, class LabelmapT >
class ITK_EXPORT RidgeSeedFilter :
  public Object
{
public:

  typedef RidgeSeedFilter                            Self;
  typedef Object                                     Superclass;
  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;

  itkTypeMacro( RidgeSeedFilter, ImageToImageFilter );

  itkNewMacro( Self );

  typedef ImageT                                  ImageType;
  typedef ImageT                                  InputImageType;
  typedef Image< float, ImageT::ImageDimension >  OutputImageType;

  typedef LabelmapT                               LabelmapType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef RidgeFeatureVectorGenerator< ImageT >   RidgeFeatureGeneratorType;

  typedef BasisFeatureVectorGenerator< ImageT, LabelmapType >
                                                  SeedFeatureGeneratorType;

  typedef typename RidgeFeatureGeneratorType::FeatureValueType
                                                  FeatureValueType;
  typedef typename RidgeFeatureGeneratorType::FeatureVectorType
                                                  FeatureVectorType;

  typedef typename RidgeFeatureGeneratorType::IndexType
                                                  IndexType;

  typedef typename RidgeFeatureGeneratorType::RidgeScalesType
                                                  RidgeScalesType;

  typedef typename SeedFeatureGeneratorType::ObjectIdType
                                                  ObjectIdType;

  typedef PDFSegmenter< OutputImageType, 3, LabelmapType >
                                                  PDFSegmenterType;

  void SetInput( typename ImageType::Pointer img );
  typename ImageType::Pointer GetInput( void );

  void SetLabelmap( typename LabelmapType::Pointer img );
  typename LabelmapType::Pointer GetLabelmap( void );

  typename SeedFeatureGeneratorType::Pointer
    GetSeedFeatureGenerator( void );

  typename RidgeFeatureGeneratorType::Pointer
    GetRidgeFeatureGenerator( void );

  typename PDFSegmenterType::Pointer GetPDFSegmenter( void );

  // Ridge
  void  SetIntensityRange( float intensityMin, float intensityMax );
  void  SetIntensityMin( float intensityMin );
  float GetIntensityMin( void );
  void  SetIntensityMax( float intensityMax );
  float GetIntensityMax( void );

  void  SetScales( const RidgeScalesType & Scales );
  RidgeScalesType GetScales( void );

  // Ridge and PDFSegmenter
  void         SetObjectId( ObjectIdType objectId );
  void         AddObjectId( ObjectIdType objectId );
  ObjectIdType GetObjectId( unsigned int num = 0 ) const;
  unsigned int GetNumberOfObjectIds( void ) const;

  // Local
  void Update();
  void ClassifyImages();

  typename LabelmapType::Pointer GetOutput( void );

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

}; // End class RidgeSeedFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRidgeSeedFilter.hxx"
#endif

#endif // End !defined(__itktubeRidgeSeedFilter_h)
