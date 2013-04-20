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
class ITK_EXPORT RidgeSeedGenerator :
  public ImageToImageFilter< ImageT,
           Image< float, ImageT::ImageDimension > >
{
public:

  typedef RidgeSeedGenerator                   Self;
  typedef ImageToImageFilter< ImageT, LabelmapT >    Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( RidgeSeedGenerator, ImageToImageFilter );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  typedef ImageT                                     InputImageType;
  typedef Image< float, Image::ImageDimension >      OutputImageType;

  typedef typename LabelmapT                         MaskImageType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef LDAGenerator< OutputImageType, LabelmapT >    LDAGeneratorType;

  typedef typename LDAGeneratorType::FeatureType        FeatureType;
  typedef typename LDAGeneratorType::FeatureVectorType  FeatureVectorType;
  typedef typename LDAGeneratorType::ObjectIdType       ObjectIdType;

  typedef typename LDAGeneratorType::LDAValuesType      LDAValuesType;
  typedef typename LDAGeneratorType::LDAVectorType      LDAVectorType;
  typedef typename LDAGeneratorType::LDAMatrixType      LDAMatrixType;
  typedef typename LDAGeneratorType::LDAImageType       LDAImageType;

  typedef std::vector< double >                         RidgeScalesType;

  //
  // Methods
  //
  void SetInputImage( typename RidgeImageType::Pointer img );
  typename RidgeImageType::Pointer GetRidgeImage( void );

  virtual unsigned int GetNumberOfFeatures( void );

  void SetIntensityRange( float intensityMin, float intensityMax );
  itkSetMacro( IntensityMin, float );
  itkGetMacro( IntensityMin, float );
  itkSetMacro( IntensityMax, float );
  itkGetMacro( IntensityMax, float );

  void SetIntensityRangeByPercentile( float percentile,
    bool findBrightPoints=true );

  itkSetMacro( Scales, RidgeScalesType );
  itkGetMacro( Scales, RidgeScalesType );

  void UpdateLDAImages();

protected:

  RidgeSeedGenerator( void );
  virtual ~RidgeSeedGenerator( void );

  typedef ContinuousIndex< double, ImageDimension > ContinuousIndexType;

  void GenerateInputRequestedRegion();

  void GenerateData();

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  RidgeSeedGenerator( const Self & );    // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  LDAGenerator::Pointer              m_LDAGenerator;

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
