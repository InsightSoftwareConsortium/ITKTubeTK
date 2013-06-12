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

#ifndef __itkTubeRidgeFeatureVectorGenerator_h
#define __itkTubeRidgeFeatureVectorGenerator_h

#include "itkTubeFeatureVectorGenerator.h"

#include <itkImage.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <vector>

namespace itk
{

namespace tube
{

template< class ImageT >
class RidgeFeatureVectorGenerator :
  public FeatureVectorGenerator<
    Image< typename ImageT::PixelType, ImageT::ImageDimension > >
{
public:

  typedef RidgeFeatureVectorGenerator          Self;
  typedef FeatureVectorGenerator< ImageT >     Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( RidgeFeatureVectorGenerator, FeatureVectorGenerator );

  itkNewMacro( Self );

  //
  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef typename Superclass::FeatureValueType      FeatureValueType;

  typedef typename Superclass::FeatureVectorType     FeatureVectorType;

  typedef typename Superclass::ImageType             ImageType;

  typedef typename Superclass::IndexType             IndexType;

  typedef std::vector< double >                      RidgeScalesType;

  //
  virtual unsigned int GetNumberOfFeatures( void ) const;

  void SetIntensityRange( float intensityMin, float intensityMax );
  float GetIntensityMin( void ) const;
  float GetIntensityMax( void ) const;

  void SetIntensityRangeByPercentile( float percentile,
    bool findBrightPoints=true );

  void SetScales( const RidgeScalesType & scales );
  const RidgeScalesType & GetScales( void ) const;

  virtual FeatureVectorType GetFeatureVector(
    const IndexType & indx ) const;

  virtual FeatureValueType GetFeatureVectorValue( const IndexType & indx,
    unsigned int fNum ) const;

protected:

  RidgeFeatureVectorGenerator( void );
  virtual ~RidgeFeatureVectorGenerator( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  // Purposely not implemented
  RidgeFeatureVectorGenerator( const Self & );
  void operator = ( const Self & );      // Purposely not implemented

  double                             m_IntensityMin;
  double                             m_IntensityMax;

  RidgeScalesType                    m_Scales;

}; // class

}  // tube namespace

}  // itk namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeRidgeFeatureVectorGenerator.txx"
#endif

#endif  // __itkTubeRidgeFeatureVectorGenerator_h
