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
#ifndef __itkTubeTubeSeedGenerator_h
#define __itkTubeTubeSeedGenerator_h

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
class TubeSeedGenerator : public LDAGenerator< ImageT, LabelmapT >
{
public:

  typedef TubeSeedGenerator                    Self;
  typedef LDAGenerator< ImageT, LabelmapT >    Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( TubeSeedGenerator, LDAGenerator );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  typedef ImageT                                        ImageType;
  typedef std::vector< typename ImageType::Pointer >    ImageListType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef LabelmapT                            MaskImageType;

  typedef double                               FeatureType;
  typedef vnl_vector< FeatureType >            FeatureVectorType;

  typedef int                                  ObjectIdType;
  typedef std::vector< ObjectIdType >          ObjectIdListType;

  typedef vnl_vector< double >                 ObjectMeanType;
  typedef std::vector< ObjectMeanType >        ObjectMeanListType;

  typedef vnl_matrix< double >                 ObjectCovarianceType;
  typedef std::vector< ObjectCovarianceType >  ObjectCovarianceListType;

  typedef vnl_vector< double >                 LDAValuesType;
  typedef vnl_vector< double >                 LDAVectorType;
  typedef vnl_matrix< double >                 LDAMatrixType;

  typedef std::vector< double >                NJetScalesType;

  typedef itk::Image< float, ImageDimension >           LDAImageType;
  typedef std::vector< typename LDAImageType::Pointer > LDAImageListType;

  //
  // Methods
  //
  unsigned int GetNumberOfFeatures( void );

  void SetIntensityRange( float intensityMin, float intensityMax );
  float GetIntensityMin( void );
  float GetIntensityMax( void );

  void SetIntensityThresholdByPercentile( float percentile,
    bool findBrightPoints=true );

  void SetScales( const NJetScalesType & scales );

  NJetScalesType & GetScales( void );

  const typename LDAImageType::Pointer & GetNJetFeatureImage(
    unsigned int num );

  void Update();

  void UpdateLDAImages();

protected:

  TubeSeedGenerator( void );
  virtual ~TubeSeedGenerator( void );

  typedef ContinuousIndex< double, ImageDimension > ContinuousIndexType;

  void GenerateNJetFeatureImages( void );

  virtual LDAValuesType GetFeatureVector( const ContinuousIndexType & indx );

  virtual void GenerateLDA( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  TubeSeedGenerator( const Self & );          // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  NJetScalesType m_Scales;

  LDAImageListType  m_NJetFeatureImageList;

  FeatureVectorType m_NJetFeatureVector;
};

}

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeTubeSeedGenerator.txx"
#endif

#endif
