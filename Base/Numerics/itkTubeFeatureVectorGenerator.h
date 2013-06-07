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
#ifndef __itkTubeFeatureVectorGenerator_h
#define __itkTubeFeatureVectorGenerator_h

#include <vector>

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "itkImage.h"

namespace itk
{

namespace tube
{

template< class ImageT >
class FeatureVectorGenerator : public ProcessObject
{
public:

  typedef FeatureVectorGenerator               Self;
  typedef Object                               Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( FeatureVectorGenerator, ProcessObject );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  typedef ImageT                                        ImageType;
  typedef std::vector< typename ImageType::Pointer >    ImageListType;

  typedef typename ImageT::IndexType                    IndexType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef typename ImageT::PixelType           FeatureValueType;
  typedef vnl_vector< FeatureValueType >       FeatureVectorType;
  typedef Image< FeatureValueType, ImageT::ImageDimension >
                                               FeatureImageType;

  typedef double                               ValueType;
  typedef std::vector< ValueType >             ValueListType;

  //
  // Methods
  //
  void SetInputImage( typename ImageType::Pointer img );
  void AddInputImage( typename ImageType::Pointer img );

  virtual typename ImageType::Pointer GetInputImage(
    unsigned int num ) const;

  unsigned int GetNumberOfInputImages() const;

  void UpdateWhitenFeatureImageStats( void );

  void SetWhitenMeans( const ValueListType & means );
  const ValueListType & GetWhitenMeans( void ) const;

  void SetWhitenStdDevs( const ValueListType & stdDevs );
  const ValueListType & GetWhitenStdDevs( void ) const;

  void SetWhitenFeatureImageMean( unsigned int num, double mean );
  double GetWhitenFeatureImageMean( unsigned int num ) const;

  void SetWhitenFeatureImageStdDev( unsigned int num, double stdDev );
  double GetWhitenFeatureImageStdDev( unsigned int num ) const;

  virtual unsigned int GetNumberOfFeatures( void ) const;

  virtual FeatureVectorType GetFeatureVector(
    const IndexType & indx ) const;

  virtual FeatureValueType GetFeatureVectorValue(
    const IndexType & indx, unsigned int fNum ) const;

  virtual typename FeatureImageType::Pointer GetFeatureImage(
    unsigned int num ) const;

protected:

  FeatureVectorGenerator( void );
  virtual ~FeatureVectorGenerator( void );

  ImageListType                   m_InputImageList;

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  // Purposely not implemented
  FeatureVectorGenerator( const Self & );
  void operator = ( const Self & );      // Purposely not implemented

  //  Data
  ValueListType                   m_WhitenFeatureImageMean;
  ValueListType                   m_WhitenFeatureImageStdDev;

};

}

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeFeatureVectorGenerator.txx"
#endif

#endif
