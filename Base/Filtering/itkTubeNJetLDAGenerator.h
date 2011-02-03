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
#ifndef __itkTubeNJetLDAGenerator_h
#define __itkTubeNJetLDAGenerator_h

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
class NJetLDAGenerator : public LDAGenerator< ImageT, LabelmapT >
{
public:

  typedef NJetLDAGenerator                     Self;
  typedef LDAGenerator< ImageT, LabelmapT >    Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( NJetLDAGenerator, LDAGenerator );

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

  void SetZeroScales( const NJetScalesType & scales );
  void SetFirstScales( const NJetScalesType & scales );
  void SetSecondScales( const NJetScalesType & scales );
  void SetRidgeScales( const NJetScalesType & scales );

  NJetScalesType & GetZeroScales( void );
  NJetScalesType & GetFirstScales( void );
  NJetScalesType & GetSecondScales( void );
  NJetScalesType & GetRidgeScales( void );

  void SetForceIntensityConsistency( bool _forceIntensity );
  bool GetForceIntensityConsistency( void );
  void SetForceOrientationInsensitivity( bool _forceOrientation );
  bool GetForceOrientationInsensitivity( void );

  const typename LDAImageType::Pointer & GetNJetFeatureImage(
    unsigned int num );

  const typename LDAImageType::Pointer & GetNJetKernelImage(
    unsigned int num );

  void Update();

  void UpdateLDAImages();

protected:

  NJetLDAGenerator( void );
  virtual ~NJetLDAGenerator( void );

  typedef ContinuousIndex< double, ImageDimension > ContinuousIndexType;

  void GenerateNJetFeatureImages( void );

  virtual LDAValuesType GetFeatureVector( const ContinuousIndexType & indx );

  virtual void GenerateLDA( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  NJetLDAGenerator( const Self & );          // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  NJetScalesType m_ZeroScales;
  NJetScalesType m_FirstScales;
  NJetScalesType m_SecondScales;
  NJetScalesType m_RidgeScales;

  bool m_ForceIntensityConsistency;
  bool m_ForceOrientationInsensitivity;

  typename LDAImageType::Pointer m_NJetKernelImage;

  LDAImageListType m_NJetFeatureImageList;

  FeatureVectorType m_NJetFeatureVector;
};

}

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeNJetLDAGenerator.txx"
#endif

#endif
