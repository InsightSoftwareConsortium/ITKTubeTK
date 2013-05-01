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
#ifndef __itkTubeRidgeSeedGenerator_h
#define __itkTubeRidgeSeedGenerator_h

#include <vector>

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "itkImage.h"

#include "itkTubeLDAGenerator.h"

namespace itk
{

namespace tube
{

template< class ImageT, class LabelmapT >
class RidgeSeedGenerator :
  public LDAGenerator< itk::Image< float, ImageT::ImageDimension >,
  LabelmapT >
{
public:

  typedef RidgeSeedGenerator                   Self;
  typedef LDAGenerator< ImageT, LabelmapT >    Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( RidgeSeedGenerator, LDAGenerator );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  typedef ImageT                                     RidgeImageType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef typename Superclass::MaskImageType         MaskImageType;

  typedef typename Superclass::FeatureType           FeatureType;
  typedef typename Superclass::FeatureVectorType     FeatureVectorType;
  typedef typename Superclass::ObjectIdType          ObjectIdType;
  typedef typename Superclass::LDAValuesType         LDAValuesType;
  typedef typename Superclass::LDAVectorType         LDAVectorType;
  typedef typename Superclass::LDAMatrixType         LDAMatrixType;
  typedef typename Superclass::LDAImageType          LDAImageType;

  typedef std::vector< double >                      RidgeScalesType;

  //
  // Methods
  //
  void SetRidgeImage( typename RidgeImageType::Pointer img );
  typename RidgeImageType::Pointer GetRidgeImage( void );

  virtual unsigned int GetNumberOfFeatures( void );

  void SetIntensityRange( float intensityMin, float intensityMax );
  float GetIntensityMin( void );
  float GetIntensityMax( void );

  void SetIntensityRangeByPercentile( float percentile,
    bool findBrightPoints=true );

  void SetScales( const RidgeScalesType & scales );

  RidgeScalesType & GetScales( void );

  void Update( void );

  void UpdateLDAImages( void );

protected:

  RidgeSeedGenerator( void );
  virtual ~RidgeSeedGenerator( void );

  typedef ContinuousIndex< double, ImageDimension > ContinuousIndexType;

  //
  // Methods from LDAGenerator
  //
  virtual void GenerateLDA( void );

  //
  // Methods
  //
  void GenerateFeatureImages( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  RidgeSeedGenerator( const Self & );    // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  typename RidgeImageType::Pointer   m_RidgeImage;

  double                             m_IntensityMin;
  double                             m_IntensityMax;

  RidgeScalesType                    m_Scales;
};

}

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeRidgeSeedGenerator.txx"
#endif

#endif
