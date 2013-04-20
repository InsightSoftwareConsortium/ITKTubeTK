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
#ifndef __itkTubeBasisFeatureVectorGenerator_h
#define __itkTubeBasisFeatureVectorGenerator_h

#include <vector>

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "itkImage.h"

#include "itkTubeFeatureVectorGenerator.h"

namespace itk
{

namespace tube
{

template< class ImageT >
class BasisFeatureVectorGenerator
: public FeatureVectorGenerator< ImageT >
{
public:

  typedef BasisFeatureVectorGenerator          Self;
  typedef FeatureVectorGenerator< ImageT >     Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( BasisFeatureVectorGenerator, FeatureVectorGenerator );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  typedef Superclass::ImageType                ImageType;

  typedef LabelmapT                            MaskImageType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef Superclass::FeatureType              FeatureType;
  typedef Superclass::FeatureVectorType        FeatureVectorType;


  typedef LabelmapT::PixelType                 ObjectIdType;

  typedef std::vector< ObjectIdType >          ObjectIdListType;

  typedef vnl_vector< double >                 BasisValuesType;
  typedef vnl_vector< double >                 BasisVectorType;
  typedef vnl_matrix< double >                 BasisMatrixType;

  //
  // Methods
  //
  void         SetInputFeatureVectorGenerator(
                 FeatureVectorGeneratorType::Pointer & fGen );

  void         SetObjectId( ObjectIdType objectId );
  void         AddObjectId( ObjectIdType objectId );
  ObjectIdType GetObjectId( unsigned int num = 0 ) const;
  unsigned int GetNumberOfObjectIds( void ) const;

  ObjectMeanType       * GetObjectMean( ObjectIdType objectId ) const;
  ObjectCovarianceType * GetObjectCovariance( ObjectIdType objectId ) const;

  ObjectMeanType       * GetGlobalMean( void ) const;
  ObjectCovarianceType * GetGlobalCovariance( void ) const;

  unsigned int    GetNumberOfBasis( void ) const;

  BasisVectorType GetBasisVector( unsigned int basisNum ) const;
  double          GetBasisValue( unsigned int basisNum ) const;
  BasisMatrixType GetBasisMatrix( void ) const;
  BasisValuesType GetBasisValues( void ) const;

  void            SetBasisVector( unsigned int basisNum,
                    const BasisVectorType & vec );
  void            SetBasisValue( unsigned int basisNum, double value );
  void            SetBasisMatrix( const BasisMatrixType & mat );
  void            SetBasisValues( const BasisValuesType & values );

  itkSetMacro( PerformLDA, bool );
  itkGetMacro( PerformLDA, bool );
  itkSetMacro( PerformPCA, bool );
  itkGetMacro( PerformPCA, bool );

  virtual void GenerateBasis( void );

  void SetNumberOfBasisToUseAsFeatures( int numBasisUsed );
  int  GetNumberOfFeatures( void ) const;

  virtual void GetFeatureVector( const IndexType & indx ) const;

  virtual void GetFeatureVectorValue( const IndexType & indx,
    unsigned int fNum ) const;

protected:

  BasisFeatureVectorGenerator( void );
  virtual ~BasisFeatureVectorGenerator( void );

  virtual void GenerateStatistics( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  // Purposely not implemented
  BasisFeatureVectorGenerator( const Self & );
  void operator = ( const Self & );      // Purposely not implemented

  //  Data
  FeatureVectorGeneratorType::Pointer
                                  m_InputFeatureVectorGenerator;

  bool                            m_PerformLDA;
  bool                            m_PerformPCA;

  unsigned int                    m_NumberOfBasis;
  unsigned int                    m_NumberOfBasisToUseAsFeatures;

  BasisMatrixType                 m_BasisMatrix;
  BasisValuesType                 m_BasisValues;
};

}

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeBasisFeatureVectorGenerator.txx"
#endif

#endif
